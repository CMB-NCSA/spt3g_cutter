import argparse
import logging
import pandas
import time
from pyaml_env import parse_config
import multiprocessing as mp
import spt3g_cutter.fitsfinder as fitsfinder
import spt3g_cutter.cutterlib as cutterlib


def cmdline():

    # Make a proto-parse use to read in the default yaml configuration
    # file, Turn off help, so we print all options in response to -h
    conf_parser = argparse.ArgumentParser(add_help=False)
    conf_parser.add_argument("-c", "--configfile", help="SPT3G config file")
    args, remaining_argv = conf_parser.parse_known_args()
    # If we have -c or --config, then we proceed to read it
    if args.configfile:
        conf_defaults = parse_config(args.configfile)
    else:
        conf_defaults = {}

    # 2. This is the main parser
    parser = argparse.ArgumentParser(description="SPT3G thumbnail cutter tool",
                                     # Inherit options from config_parser
                                     parents=[conf_parser])
    # Location (RA,DEC, XSIZE, YSIZE)
    parser.add_argument("--inputList",
                        help="Input CSV file with positions (RA,DEC) and optional (XSIZE,YSIZE) in arcmins")
    parser.add_argument("--xsize", type=float, action="store", default=None,
                        help="Length of x-side in arcmins of image [default = 10]")
    parser.add_argument("--ysize", type=float, action="store", default=None,
                        help="Length of y-side of in arcmins image [default = 10]")
    # File location
    parser.add_argument("--outdir", type=str, action='store', default=None,
                        help="Location for output files")
    parser.add_argument("--prefix", type=str, action='store', default='spt3g',
                        help="Prefix for thumbnail filenames [default='spt3g']")
    # DB options
    parser.add_argument("--dbname", type=str, action='store', default=None,
                        help="Database (file) to connect")
    parser.add_argument("--tablename", type=str, action='store', default='file_info_v1',
                        help="Name of tablw with file informatiom")
    parser.add_argument("--bands", nargs="*", default=['90GHz', '150GHz', '220GHz'],
                        help="The bands to select from: 90GHz, 150GHz and 220GHz")
    parser.add_argument("--filetypes", nargs="*", default=['raw', 'filtered'],
                        help="The filetype to select: 'raw/filtered'")
    parser.add_argument("--date_start", type=str, action='store', default=None,
                        help="The START date to search for files formatted [YYYY-MM-DD]")
    parser.add_argument("--date_end", type=str, action='store', default=None,
                        help="The END date to search for files formatted [YYYY-MM-DD]")
    parser.add_argument("--yearly", type=str, action='store', default=None,
                        help="The yearly tag or tags to use [i.e. yearly_winter_2020]")

    # Logging options (loglevel/log_format/log_format_date)
    default_log_format = '[%(asctime)s.%(msecs)03d][%(levelname)s][%(name)s][%(funcName)s] %(message)s'
    default_log_format_date = '%Y-%m-%d %H:%M:%S'
    parser.add_argument("--loglevel", action="store", default='INFO', type=str.upper,
                        choices=["DEBUG", "INFO", "WARNING", "ERROR", "CRITICAL"],
                        help="Logging Level [DEBUG/INFO/WARNING/ERROR/CRITICAL]")
    parser.add_argument("--log_format", action="store", type=str, default=default_log_format,
                        help="Format for logging")
    parser.add_argument("--log_format_date", action="store", type=str, default=default_log_format_date,
                        help="Format for date section of logging")

    # Use multiprocessing
    parser.add_argument("--np", action="store", default=1, type=int,
                        help="Run using multi-process, 0=automatic, 1=single-process [default]")

    # args = parser.parse_args()
    # Set the defaults of argparse using the values in the yaml config file
    parser.set_defaults(**conf_defaults)
    args = parser.parse_args(args=remaining_argv)
    args.loglevel = getattr(logging, args.loglevel)

    print(args)

    # Make sure that both date_start/end are defined or both are None
    if args.date_start is None and args.date_end is None:
        pass
    elif isinstance(args.date_start, str) and isinstance(args.date_end, str):
        pass
    else:
        raise ValueError('Both --date_start and --date_end must be defined')

    return args


def run(args):

    # Create logger
    cutterlib.create_logger(level=args.loglevel,
                            log_format=args.log_format,
                            log_format_date=args.log_format_date)
    logger = logging.getLogger(__name__)

    # Read in CSV file with pandas
    logger.info(f"Reading: {args.inputList}")
    df = pandas.read_csv(args.inputList)

    # Get the arrays with ra, dec, xsize, ysize
    (xsize, ysize) = fitsfinder.check_xysize(df, xsize=args.xsize, ysize=args.ysize)
    args.ra = df.RA.values
    args.dec = df.DEC.values

    # Connect, get query and run query
    dbhandle = fitsfinder.connect_db(args.dbname)
    query = fitsfinder.get_query(args.tablename,
                                 bands=args.bands,
                                 filetypes=args.filetypes,
                                 date_start=args.date_start,
                                 date_end=args.date_end,
                                 yearly=args.yearly)
    logger.info(f"Running query: {query}")
    rec = fitsfinder.query2rec(query, dbhandle)

    cutout_names = {}

    # Get the number of processors to use
    NP = cutterlib.get_NP(args.np)
    if NP > 1:
        p = mp.Pool(processes=NP)
        logger.info(f"Will use {NP} processors to convert and ingest")
        manager = mp.Manager()
        return_dict = manager.dict()
    else:
        return_dict = None

    # Loop over all files
    args.files = rec['FILE'].tolist()
    Nfiles = len(args.files)
    logger.info(f"Found {Nfiles} files")
    k = 1
    t0 = time.time()
    for file in args.files:
        counter = f"{k}/{Nfiles} files"
        ar = (file, args.ra, args.dec, return_dict)
        kw = {'xsize': xsize, 'ysize': ysize, 'units': 'arcmin',
              'prefix': args.prefix, 'outdir': args.outdir, 'counter': counter}
        if NP > 1:
            p.apply_async(cutterlib.fitscutter, args=ar, kwds=kw)
        else:
            res = cutterlib.fitscutter(*ar, **kw)
            cutout_names.update(res)
        k += 1

    if NP > 1:
        p.close()
        p.join()
        # Update with returned dictionary
        cutout_names = dict(return_dict)

    # Store the dict with all of the cutout names
    args.cutout_names = cutout_names

    args = cutterlib.capture_job_metadata(args)

    # Write the manifest yaml file
    cutterlib.write_manifest(args)

    logger.info(f"Grand Total time: {cutterlib.elapsed_time(t0)}")
