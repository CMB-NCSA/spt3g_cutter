import argparse
import logging
import pandas
import time
import sys
from pyaml_env import parse_config
import multiprocessing as mp
import spt3g_cutter
import spt3g_cutter.fitsfinder as fitsfinder
import spt3g_cutter.cutterlib as cutterlib
import os
import psutil
import copy


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
    parser.add_argument("--tablename", type=str, action='store', default='file_info_v2',
                        help="Name of tablw with file informatiom")
    parser.add_argument("--bands", nargs="*", default=['90GHz', '150GHz', '220GHz'],
                        help="The bands to select from: 90GHz, 150GHz and 220GHz")
    parser.add_argument("--filetypes", nargs="*", default=['passthrough', 'filtered'],
                        help="The filetype to select: 'passthrough/filtered'")
    parser.add_argument("--date_start", type=str, action='store', default=None,
                        help="The START date to search for files formatted [YYYY-MM-DD]")
    parser.add_argument("--date_end", type=str, action='store', default=None,
                        help="The END date to search for files formatted [YYYY-MM-DD]")
    parser.add_argument("--yearly", nargs="*", default=None,
                        help="The yearly tag or tags to use [i.e. yearly_winter_2020]")
    parser.add_argument("--get_lightcurve", action='store_true', default=False,
                        help="Extract light curve at pixel position for each object")
    parser.add_argument("--version", action="version", version=f"spt3g_cutter: {spt3g_cutter.__version__}",
                        help="Print version and exit")
    parser.add_argument("--get_uniform_coverage", action='store_true', default=False,
                        help="Get only objects within the uniform coverage")
    parser.add_argument("--nofits", action='store_true', default=False,
                        help="Do not create fits files for stamps")
    # Read options
    parser.add_argument("--stage", action='store_true', default=False,
                        help="Stage input files before operanting on them.")
    parser.add_argument("--stage_path", action='store', default=None,
                        help="Path for indirect write.")

    # Logging options (loglevel/log_format/log_format_date)
    if 'LOG_LEVEL' in os.environ:
        default_log_level = os.environ['LOG_LEVEL']
    else:
        default_log_level = 'INFO'
    default_log_format = '[%(asctime)s.%(msecs)03d][%(levelname)s][%(name)s][%(funcName)s] %(message)s'
    default_log_format_date = '%Y-%m-%d %H:%M:%S'
    parser.add_argument("--loglevel", action="store", default=default_log_level, type=str.upper,
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

    # In case we pass and yearly tags
    if args.yearly is not None:
        if len(args.yearly) == 0:
            args.yearly = None

    # Make sure that both date_start/end are defined or both are None
    if args.date_start is None and args.date_end is None:
        pass
    elif isinstance(args.date_start, str) and isinstance(args.date_end, str):
        pass
    else:
        raise ValueError('Both --date_start and --date_end must be defined')

    if args.stage_path is None and args.stage is True:
        if 'SPT3G_INGEST_STAGE_PATH' in os.environ:
            args.stage_path = os.environ['SPT3G_INGEST_STAGE_PATH']
        else:
            args.stage_path = '/tmp'

    # Define the prefix for staging input files
    if args.stage is True:
        args.stage_prefix = os.path.join(args.stage_path, 'spt3g_cutter-stage-')
    else:
        args.stage_prefix = None
    return args


def run(args):

    # Create logger
    cutterlib.create_logger(level=args.loglevel,
                            log_format=args.log_format,
                            log_format_date=args.log_format_date)
    logger = logging.getLogger(__name__)

    logger.info(f"Received command call:\n{' '.join(sys.argv[0:-1])}")
    logger.info(f"Running spt3g_cutter:{spt3g_cutter.__version__}")
    logger.info(f"Running with args: \n{args}")

    # Read in CSV file with pandas
    logger.info(f"Reading: {args.inputList}")
    df = pandas.read_csv(args.inputList)

    # Get the arrays with ra, dec, xsize, ysize
    (xsize, ysize) = fitsfinder.check_xysize(df, xsize=args.xsize, ysize=args.ysize)
    args.ra = df.RA.values.tolist()
    args.dec = df.DEC.values.tolist()
    if 'OBJID' in df.keys():
        args.objID = df.OBJID.values.tolist()
    else:
        args.objID = None

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
    rejected_names = {}
    lightcurve = {}

    # Get the number of processors to use
    NP = cutterlib.get_NP(args.np)
    if NP > 1:
        p = mp.Pool(processes=NP)
        logger.info(f"Will use {NP} processors for process")
        manager = mp.Manager()
        cutout_dict = manager.dict()
        rejected_dict = manager.dict()
        lightcurve_dict = manager.dict()
        results = []
    else:
        cutout_dict = None
        rejected_dict = None
        lightcurve_dict = None

    # Create the outdir if it does not exists
    if not os.path.exists(args.outdir):
        os.makedirs(args.outdir, mode=0o755, exist_ok=True)

    # Loop over all files
    args.files = rec['FILE'].tolist()
    Nfiles = len(args.files)
    logger.info(f"Found {Nfiles} files")
    k = 1
    t0 = time.time()
    for file in args.files:
        counter = f"{k}/{Nfiles} files"

        # Make a copy of objID if not None:
        if args.objID is None:
            objID = None
        else:
            objID = copy.deepcopy(args.objID)

        ar = (file, args.ra, args.dec, cutout_dict, rejected_dict, lightcurve_dict)
        kw = {'xsize': xsize, 'ysize': ysize, 'units': 'arcmin', 'objID': objID,
              'prefix': args.prefix, 'outdir': args.outdir, 'counter': counter,
              'get_lightcurve': args.get_lightcurve,
              'get_uniform_coverage': args.get_uniform_coverage,
              'nofits': args.nofits,
              'stage': args.stage,
              'stage_prefix': args.stage_prefix}

        if NP > 1:
            # Get result to catch exceptions later, after close()
            s = p.apply_async(cutterlib.fitscutter, args=ar, kwds=kw)
            results.append(s)
        else:
            names, pos, lc = cutterlib.fitscutter(*ar, **kw)
            cutout_names.update(names)
            rejected_names.update(pos)
            lightcurve.update(lc)
        k += 1

    if NP > 1:
        p.close()
        # Check for exceptions
        for r in results:
            r.get()
        p.join()

        # Update with returned dictionary, we need to make them real
        # dictionaries, instead DictProxy objects returned from multiprocessing
        logger.info("Updating returned dictionaries")
        cutout_names = cutout_dict.copy()
        rejected_names = rejected_dict.copy()
        lightcurve = lightcurve_dict.copy()
        p.terminate()
        del p

    # Time it took to just cut
    logger.info(f"Total cutting time: {cutterlib.elapsed_time(t0)}")

    # Store the dict with all of the cutout names and rejects
    args.cutout_names = cutout_names
    args.rejected_names = rejected_names

    args = cutterlib.capture_job_metadata(args)

    # Report total memory usage
    logger.info(f"Memory: {psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3} Gb")
    process = psutil.Process(os.getpid())
    logger.info(f"Memory percent: {process.memory_percent()} %")

    # # Clean up
    if NP > 1:
        logger.info("Deleting variables -- probably futile")
        del manager
        del cutout_names
        del cutout_dict
        del rejected_names
        del rejected_dict
        del lightcurve_dict

    if args.get_lightcurve:

        # Get the observations dictionary
        args.obs_dict = cutterlib.get_obs_dictionary(lightcurve)
        logger.info(f"Size of lightcurve: {sys.getsizeof(lightcurve)/1024/1024}")
        logger.info(f"Size of args.obs_dict: {sys.getsizeof(args.obs_dict)/1024/1024}")

        # Create new pool
        # NP = 1  # Remove this line once we have a machine with more memory
        NP = len(args.bands)
        logger.info(f"Creating pool with: {NP} processes for repack lightcurve")
        if NP > 1:
            p = mp.Pool(processes=NP)
            results = []

        for BAND in args.bands:
            for FILETYPE in args.filetypes:
                logger.info(f"Memory: {psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3} Gb")
                ar = (lightcurve, BAND, FILETYPE, args)
                logger.info(f"Memory percent: {process.memory_percent()} %")
                if NP > 1:
                    s = p.apply_async(cutterlib.repack_lightcurve_band_filetype, args=ar)
                    results.append(s)
                else:
                    cutterlib.repack_lightcurve_band_filetype(*ar)

        if NP > 1:
            p.close()
            # Check for exceptions
            for r in results:
                r.get()
            p.join()
            p.terminate()
            del p

    # Write the manifest file
    cutterlib.write_manifest(args)
    logger.info(f"Grand Total time: {cutterlib.elapsed_time(t0)}")
