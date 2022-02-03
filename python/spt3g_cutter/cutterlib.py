#!/usr/bin/env python

"""
A set of simple proto-function to make postage stamps using fitsio
Based on desthumbs code circa 2015
"""

import fitsio
import os
import sys
import spt3g_cutter
import spt3g_cutter.astrometry as astrometry
import time
import numpy
import copy
from collections import OrderedDict
from astropy.wcs import WCS
from astropy.utils.exceptions import AstropyWarning
import logging
from logging.handlers import RotatingFileHandler
import warnings
import multiprocessing
import yaml
import datetime
import subprocess

# To avoid header warning from astropy
warnings.filterwarnings('ignore', category=AstropyWarning, append=True)

# Logger
LOGGER = logging.getLogger(__name__)

# Naming template
PREFIX = 'SPT3G'
OBJ_ID = "{prefix}J{ra}{dec}"
FITS_OUTNAME = "{outdir}/{objID}_{filter}_{obsid}_{filetype_ext}.{ext}"
LOG_OUTNAME = "{outdir}/{objID}.{ext}"
BASE_OUTNAME = "{objID}"
BASEDIR_OUTNAME = "{outdir}/{objID}"
FILETYPE_EXT = {'passthrough': 'psth', 'filtered': 'fltd'}


def configure_logger(logger, logfile=None, level=logging.NOTSET, log_format=None, log_format_date=None):
    """
    Configure an existing logger
    """
    # Define formats
    if log_format:
        FORMAT = log_format
    else:
        FORMAT = '[%(asctime)s.%(msecs)03d][%(levelname)s][%(name)s][%(funcName)s] %(message)s'
    if log_format_date:
        FORMAT_DATE = log_format_date
    else:
        FORMAT_DATE = '%Y-%m-%d %H:%M:%S'
    formatter = logging.Formatter(FORMAT, FORMAT_DATE)

    # Need to set the root logging level as setting the level for each of the
    # handlers won't be recognized unless the root level is set at the desired
    # appropriate logging level. For example, if we set the root logger to
    # INFO, and all handlers to DEBUG, we won't receive DEBUG messages on
    # handlers.
    logger.setLevel(level)

    handlers = []
    # Set the logfile handle if required
    if logfile:
        fh = RotatingFileHandler(logfile, maxBytes=2000000, backupCount=10)
        fh.setFormatter(formatter)
        fh.setLevel(level)
        handlers.append(fh)
        logger.addHandler(fh)

    # Set the screen handle
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)
    sh.setLevel(level)
    handlers.append(sh)
    logger.addHandler(sh)

    return


def create_logger(logfile=None, level=logging.NOTSET, log_format=None, log_format_date=None):
    """
    Simple logger that uses configure_logger()
    """
    logger = logging.getLogger(__name__)
    configure_logger(logger, logfile=logfile, level=level,
                     log_format=log_format, log_format_date=log_format_date)
    logging.basicConfig(handlers=logger.handlers, level=level)
    logger.propagate = False
    return logger


def elapsed_time(t1, verb=False):
    """
    Returns the time between t1 and the current time now
    I can can also print the formatted elapsed time.
    ----------
    t1: float
        The initial time (in seconds)
    verb: bool, optional
        Optionally print the formatted elapsed time
    returns
    -------
    stime: float
        The elapsed time in seconds since t1
    """
    t2 = time.time()
    stime = "%dm %2.2fs" % (int((t2-t1)/60.), (t2-t1) - 60*int((t2-t1)/60.))
    if verb:
        print("Elapsed time: {}".format(stime))
    return stime


def get_fits_hdu_extensions_byfilename(filename):
    """
    Return the HDU extension for coadds (old-school) based on the extension
    name. Check if dealing with .fz or .fits files
    """
    if os.path.basename(os.path.splitext(filename)[-1]) == '.fz':
        sci_hdu = 1
        wgt_hdu = 2
    elif os.path.basename(os.path.splitext(filename)[-1]) == '.fits':
        sci_hdu = 0
        wgt_hdu = 1
    else:
        raise NameError("ERROR: No .fz or .fits files found")
    return sci_hdu, wgt_hdu


def update_wcs_matrix(header, x0, y0, naxis1, naxis2, ra, dec, proj='ZEA'):
    """
    Update the wcs header object with the right CRPIX[1, 2] CRVAL[1, 2] for a
    given subsection

    Parameters:
    header: fits style header
        The header to work with
    x0, y0: float
        The new center of the image
    naxis1, naxis2: int
        The number of pixels on each axis.

    Returns:
        fits style header with the new center.
    """

    # We need to make a deep copy/otherwise if fails
    h = copy.deepcopy(header)
    # Get the astropy.wcs object
    wcs = WCS(h)
    # Recompute CRVAL1/2 on the new center x0,y0
    CRVAL1, CRVAL2 = wcs.wcs_pix2world(x0, y0, 1)
    # Recast numpy objects as floats
    CRVAL1 = float(CRVAL1)
    CRVAL2 = float(CRVAL2)
    # Asign CRPIX1/2 on the new image
    CRPIX1 = int(naxis1/2.0)
    CRPIX2 = int(naxis2/2.0)
    # Update the values
    h['CRVAL1'] = CRVAL1
    h['CRVAL2'] = CRVAL2
    h['CRPIX1'] = CRPIX1
    h['CRPIX2'] = CRPIX2

    if proj == 'TAN':
        h['CTYPE1'] = 'RA---TAN'
        h['CTYPE2'] = 'DEC--TAN'
        # Delete some key that are not needed
        dkeys = ['PROJ', 'LONPOLE', 'LATPOLE', 'POLAR', 'ALPHA0', 'DELTA0', 'X0', 'Y0']
        for k in dkeys:
            h.delete(k)
    elif proj == 'ZEA':
        h['LATPOLE'] = CRVAL2
    else:
        raise NameError(f"Projection: {proj} not implemented")

    # New record for RA/DEC center
    recs = [{'name': 'racut', 'value': ra, 'comment': 'RA of cutout'},
            {'name': 'deccut', 'value': dec, 'comment': 'DEC of cutout'}]
    for rec in recs:
        h.add_record(rec)

    return h


def check_inputs(ra, dec, xsize, ysize, objID=None):

    """ Check and fix inputs for cutout"""
    # Make sure that RA,DEC are the same type
    if type(ra) != type(dec):
        raise TypeError('RA and DEC need to be the same type()')

    if objID is not None and type(objID) != type(ra):
        raise TypeError('objID needs to be the same type() as RA,DEC')

    if type(xsize) != type(ysize):
        raise TypeError('XSIZE and YSIZE need to be the same type()')
    # Make them iterable and proper length
    if hasattr(ra, '__iter__') is False and hasattr(dec, '__iter__') is False:
        ra = [ra]
        dec = [dec]
    if objID is not None and hasattr(objID, '__iter__') is False:
        objID = [objID]
    if hasattr(xsize, '__iter__') is False and hasattr(ysize, '__iter__') is False:
        xsize = [xsize]*len(ra)
        ysize = [ysize]*len(ra)
    # Make sure they are all of the same length
    if len(ra) != len(dec):
        raise TypeError('RA and DEC need to be the same length')
    if objID is not None and len(objID) != len(ra):
        raise TypeError('objID needs to be the same length as RA, DEC')
    if len(xsize) != len(ysize):
        raise TypeError('XSIZE and YSIZE need to be the same length')
    if (len(ra) != len(xsize)) or (len(ra) != len(ysize)):
        raise TypeError('RA, DEC and XSIZE and YSIZE need to be the same length')
    # Now make sure that all objID are unique
    if objID is not None and len(set(objID)) != len(objID):
        raise TypeError('Elements in objID are not unique')
    # If objID is None, return a list of None of the same lenght as ra,dec
    if objID is None:
        objID = [objID]*len(ra)

    return ra, dec, xsize, ysize, objID


def get_thumbFitsName(ra, dec, filter, obsid, filetype_ext,
                      objID=None, prefix=PREFIX, ext='fits', outdir=os.getcwd()):
    """ Common function to set the Fits thumbnail name """
    ra = astrometry.dec2deg(ra/15., sep="", plussign=False)
    dec = astrometry.dec2deg(dec, sep="", plussign=True)
    if objID is None:
        objID = OBJ_ID.format(ra=ra, dec=dec, prefix=prefix)
    # Locals need to be captured at the end
    kw = locals()
    outname = FITS_OUTNAME.format(**kw)
    return outname


def get_thumbBaseDirName(ra, dec, objID=None, prefix=PREFIX, outdir=os.getcwd()):
    """ Common function to set the Fits thumbnail name """
    ra = astrometry.dec2deg(ra/15., sep="", plussign=False)
    dec = astrometry.dec2deg(dec, sep="", plussign=True)
    if objID is None:
        objID = OBJ_ID.format(ra=ra, dec=dec, prefix=prefix)
    # Locals need to be captured at the end
    kw = locals()
    basedir = BASEDIR_OUTNAME.format(**kw)
    return basedir


def get_thumbLogName(ra, dec, objID=None, prefix=PREFIX, ext='log', outdir=os.getcwd()):
    """ Common function to set the Fits thumbnail name """
    ra = astrometry.dec2deg(ra/15., sep="", plussign=False)
    dec = astrometry.dec2deg(dec, sep="", plussign=True)
    if objID is None:
        objID = OBJ_ID.format(ra=ra, dec=dec, prefix=prefix)
    # Locals need to be captured at the end
    kw = locals()
    outname = LOG_OUTNAME.format(**kw)
    return outname


def get_thumbBaseName(ra, dec, objID=None, prefix=PREFIX):
    """ Common function to set the Fits thumbnail name """
    ra = astrometry.dec2deg(ra/15., sep="", plussign=False)
    dec = astrometry.dec2deg(dec, sep="", plussign=True)
    if objID is None:
        objID = OBJ_ID.format(ra=ra, dec=dec, prefix=prefix)
    # Locals need to be captured at the end
    kw = locals()
    outname = BASE_OUTNAME.format(**kw)
    return outname


def get_headers_hdus(filename):

    header = OrderedDict()
    hdu = OrderedDict()

    is_compressed = False
    with fitsio.FITS(filename) as fits:
        # Case 1 -- for well-defined fitsfiles with EXTNAME
        for k in range(len(fits)):
            h = fits[k].read_header()
            # Is compressed
            if h.get('ZIMAGE'):
                is_compressed = True
            # Make sure that we can get the EXTNAME
            if not h.get('EXTNAME'):
                continue
            extname = h['EXTNAME'].strip()
            if extname == 'COMPRESSED_IMAGE':
                is_compressed = True
                continue
            header[extname] = h
            hdu[extname] = k

        # Case 2 -- files without EXTNAME
        if len(header) < 1:
            LOGGER.debug("Getting EXTNAME by compression")
            if is_compressed:
                sci_hdu = 1
                wgt_hdu = 2
            else:
                sci_hdu = 0
                wgt_hdu = 1
            # Assign headers and hdus
            header['SCI'] = fits[sci_hdu].read_header()
            hdu['SCI'] = sci_hdu
            try:
                header['WGT'] = fits[wgt_hdu].read_header()
                hdu['WGT'] = wgt_hdu
            except IOError:
                LOGGER.warning(f"No WGT HDU for: {filename}")

    return header, hdu


def get_NP(MP):

    """ Get the number of processors in the machine
    if MP == 0, use all available processor
    """
    # For it to be a integer
    MP = int(MP)
    if MP == 0:
        NP = multiprocessing.cpu_count()
    elif isinstance(MP, int):
        NP = MP
    else:
        raise ValueError('MP is wrong type: %s, integer type' % MP)
    return NP


def fitscutter(filename, ra, dec, cutout_names, rejected_positions, lightcurve,
               objID=None, xsize=1.0, ysize=1.0, units='arcmin', get_lightcurve=False,
               prefix=PREFIX, outdir=None, clobber=True, logger=None, counter=''):
    """
    Makes cutouts around ra, dec for a give xsize and ysize
    ra,dec can be scalars or lists/arrays

    Input Parameters
    ----------------
    filename: string
        The fitsfile (fits or fits.fz) to cut from
    ra, dec: float or array of floats
        The position in decimal degrees where we can to cut

    Input/Output Parameters
    -----------------------
    These are dictionary that need to be inputs to be passed back when
    running using multiprocessing

    cutout_names: dictionary
        Dict story the names of the cutout files
    rejected_positions: dictionary
        Dict of rejected positions with centers outside the FITS files
    lightcurve: dictionary
        Dict with lightcurve information

    Optional inputs
    ---------------
    xsize, ysize: float or array of floats
        The x,y size of the stamps in arcmins
    objID: string or list of string
        The list of object ID for each ra, dec pair
    units: string
        The units for the xsise, ysize
    prefix: string
        The prefix to prepend to the output names (i.e.: spt3g)
    outdir: string
        The path for the output directory
    clobber: Bool
        Overwrite file if they exist
    logger: logging object
        Optional logging object
    counter: string
        Optional counter to pass on for tracking flow

    Returns:
        cutout_names, rejected_positions
    """

    # global timer for function
    t1 = time.time()
    if not logger:
        logger = LOGGER

    if not outdir:
        outdir = os.getcwd()

    # Check and fix inputs
    ra, dec, xsize, ysize, objID = check_inputs(ra, dec, xsize, ysize, objID)

    logger.info(f"Working on FITS file: {filename} -- {counter}")
    logger.info(f"Will cut: {len(ra)} stamps")

    # Check for the units
    if units == 'arcsec':
        scale = 1
    elif units == 'arcmin':
        scale = 60
    elif units == 'degree':
        scale = 3600
    else:
        raise Exception("ERROR: must define units as arcses/arcmin/degree only")

    # Get header/extensions/hdu
    t0 = time.time()
    header, hdunum = get_headers_hdus(filename)
    logger.debug(f"Done Getting header, hdus: {elapsed_time(t0)}")
    extnames = header.keys()  # Gets SCI and WGT
    logger.debug(f"Found EXTNAMES:{extnames}")

    # Get the pixel-scale of the input image
    pixelscale = astrometry.get_pixelscale(header['SCI'], units='arcsec')

    # Read in wcs with astropy
    wcs = WCS(header['SCI'])

    # Get the dimensions of the parent image
    if 'EXTNAME' in header['SCI'] and header['SCI']['EXTNAME'].strip() == 'COMPRESSED_IMAGE':
        NAXIS1 = header['SCI']['ZNAXIS1']
        NAXIS2 = header['SCI']['ZNAXIS2']
    else:
        NAXIS1 = header['SCI']['NAXIS1']
        NAXIS2 = header['SCI']['NAXIS2']

    # Extract the band/filter from the header
    if 'BAND' in header['SCI']:
        band = header['SCI']['BAND'].strip()
    elif 'FILTER' in header['SCI']:
        band = header['SCI']['FILTER'].strip()
    else:
        raise Exception("ERROR: Cannot provide suitable BAND/FILTER from SCI header")

    # Extract OBSID from the header
    if 'OBS-ID' in header['SCI']:
        obsid = str(header['SCI']['OBS-ID']).strip()
    else:
        raise Exception("ERROR: Cannot provide suitable OBS-ID from SCI header")

    # Extract FILETYPE from the header
    if 'FILETYPE' in header['SCI']:
        filetype = str(header['SCI']['FILETYPE']).strip()
    else:
        # Try to get it from the filename
        raise Exception("ERROR: Cannot provide suitable FILETYPE from SCI header")

    if 'DATE-BEG' in header['SCI']:
        date_beg = str(header['SCI']['DATE-BEG']).strip()
    else:
        raise Exception("ERROR: Cannot provide suitable DATE_BEG from SCI header")

    # The extension to use for FILETYPE
    filetype_ext = FILETYPE_EXT[filetype]

    # Intitialize the FITS object
    t0 = time.time()
    ifits = fitsio.FITS(filename, 'r')
    logger.debug(f"Done loading fitsio.FITS({filename}): {elapsed_time(t0)}")

    if cutout_names is None:
        cutout_names = {}
    if rejected_positions is None:
        rejected_positions = {}
    if lightcurve is None:
        lightcurve = {}

    # Local lists/dicts
    outnames = []
    rejected = []
    lc_local = {}

    # Define the ID for the lightcurve information from this filename
    if get_lightcurve:
        lcID = f'{obsid}_{filetype}'
        lc_local['DATE-BEG'] = date_beg
        lc_local['BAND'] = band
        lc_local['objID'] = objID
        lc_local['FILETYPE'] = filetype

    ######################################
    # Loop over ra/dec and xsize,ysize
    for k in range(len(ra)):

        # The basename for the (ra,dec)
        if objID[k] is None:
            objID[k] = get_thumbBaseName(ra[k], dec[k], prefix=prefix)

        # Define the geometry of the thumbnail
        x0, y0 = wcs.wcs_world2pix(ra[k], dec[k], 1)
        x0 = round(float(x0))
        y0 = round(float(y0))
        dx = int(0.5*xsize[k]*scale/pixelscale)
        dy = int(0.5*ysize[k]*scale/pixelscale)
        naxis1 = 2*dx  # +1
        naxis2 = 2*dy  # +1
        y1 = y0-dy
        y2 = y0+dy
        x1 = x0-dx
        x2 = x0+dx

        # Make sure the (x0,y0) is contained within the image
        if x0 < 0 or y0 < 0 or x0 > NAXIS1 or y0 > NAXIS2:
            LOGGER.warning(f"(RA,DEC):{ra[k]},{dec[k]} outside {filename}")
            LOGGER.warning(f"(x0,y0):{x0},{y0} > {NAXIS1},{NAXIS2}")
            rejected.append(f"{ra[k]}, {dec[k]}, {objID[k]}")
            continue

        # Make sure we are not going beyond the limits
        # if negative set it to zero
        if y1 < 0:
            y1 = 0
        if y2 > NAXIS2:
            y2 = NAXIS2
        if x1 < 0:
            x1 = 0
        if x2 > NAXIS1:
            x2 = NAXIS1

        im_section = OrderedDict()
        h_section = OrderedDict()
        LOGGER.debug(f"Found naxis1,naxis2: {naxis1},{naxis2}")
        LOGGER.debug(f"Found x1,x2: {x1},{x2}")
        LOGGER.debug(f"Found y1,y2: {y1},{y2}")

        for EXTNAME in extnames:
            # The hdunum for that extname
            HDUNUM = hdunum[EXTNAME]
            # Append data from (x0, y0) pixel from EXTNAME
            if get_lightcurve:
                data_extname = float(ifits[HDUNUM][int(y0), int(x0)][0][0])
                lc_local.setdefault(f'flux_{EXTNAME}', []).append(data_extname)

            # Create a canvas
            im_section[EXTNAME] = numpy.zeros((naxis1, naxis2))
            # Read in the image section we want for SCI/WGT
            im_section[EXTNAME] = ifits[HDUNUM][int(y1):int(y2), int(x1):int(x2)]
            # Correct NAXIS1 and NAXIS2
            naxis1 = numpy.shape(im_section[EXTNAME])[1]
            naxis2 = numpy.shape(im_section[EXTNAME])[0]
            # Update the WCS in the headers and make a copy
            h_section[EXTNAME] = update_wcs_matrix(header[EXTNAME], x0, y0, naxis1, naxis2, ra[k], dec[k])
            # Add the objID to the header of the thumbnail
            rec = {'name': 'OBJECT', 'value': objID[k], 'comment': 'Name of the objID'}
            h_section[EXTNAME].add_record(rec)

        # Get the basedir
        basedir = get_thumbBaseDirName(ra[k], dec[k], objID=objID[k], prefix=prefix, outdir=outdir)
        if not os.path.exists(basedir):
            os.makedirs(basedir, mode=0o755, exist_ok=True)

        # Construct the name of the Thumbmail using BAND/FILTER/prefix/etc
        outname = get_thumbFitsName(ra[k], dec[k], band, obsid, filetype_ext,
                                    objID=objID[k], prefix=prefix, outdir=basedir)
        # Save the outnames without the output directory
        outnames.append(outname.replace(f"{outdir}/", ''))
        # Write out the file
        t0 = time.time()
        ofits = fitsio.FITS(outname, 'rw', clobber=clobber)
        for EXTNAME in extnames:
            ofits.write(im_section[EXTNAME], extname=EXTNAME, header=h_section[EXTNAME])
        ofits.close()
        logger.debug(f"Done writing {outname}: {elapsed_time(t0)}")

    ifits.close()
    logger.info(f"Done {filename} in {elapsed_time(t1)} -- {counter}")

    # Assing internal lists/dict to managed dictionalks
    cutout_names[filename] = outnames
    if get_lightcurve:
        lightcurve[lcID] = lc_local

    if len(rejected) > 0:
        rejected_positions[filename] = rejected
        logger.info(f"{len(rejected)} positions for {filename}")

    return cutout_names, rejected_positions, lightcurve


def get_id_names(ra, dec, prefix):
    "Get the ID names associated with every position"
    names = []
    for k in range(len(ra)):
        names.append(get_thumbBaseName(ra[k], dec[k], prefix=prefix))
    return names


def get_size_on_disk(outdir):
    size = subprocess.check_output(['du', '-sh', outdir]).split()[0].decode('ascii')
    return size


def get_job_info(args):
    " Get the JOB_ID and JOB_OUTPUT_DIR from the environment "
    JOB_ID = None
    JOB_OUTPUT_DIR = None
    if 'JOB_ID' in os.environ:
        JOB_ID = os.environ['JOB_ID']
    if 'JOB_OUTPUT_DIR' in os.environ:
        JOB_OUTPUT_DIR = os.environ['JOB_OUTPUT_DIR']
    return JOB_ID, JOB_OUTPUT_DIR


def get_positions_idnames(args):
    "Get the id names for all positions"
    positions_idnames = []
    for k in range(len(args.ra)):
        positions_idnames.append(f"{args.ra[k]}, {args.dec[k]}, {args.id_names[k]}")
    return positions_idnames


def capture_job_metadata(args):
    """ Get more information abot this job for the manifest"""

    # Get the ID names for each ra,dec pair and store them
    if args.objID is None:
        args.id_names = get_id_names(args.ra, args.dec, args.prefix)
    else:
        args.id_names = args.objID

    # Get the positions and id_names
    args.input_positions = get_positions_idnames(args)

    # Make a list of all of the cutout cutout_names
    cutout_files = []
    for file in args.cutout_names.keys():
        cutout_files.extend(args.cutout_names[file])
    args.cutout_files = cutout_files

    # Get the size on disk
    args.size_on_disk = get_size_on_disk(args.outdir)
    args.files_on_disk = len(args.cutout_files)

    # Get the job information from k8s
    (args.JOB_ID, args.JOB_OUTPUT_DIR) = get_job_info(args)
    return args


def repack_lightcurve(lightcurve, args):
    "Repack the lightcurve dictionary keyed by objID"

    LC = {}
    for objID in args.id_names:
        dates = {}
        flux_SCI = {}
        flux_WGT = {}
        # Loop over the observations (OBS-ID + filetype)
        for obs in lightcurve:

            # Check if we have weight flux:
            got_WGT = False
            if 'flux_WGT' in lightcurve[obs]:
                got_WGT = True

            # Get filetype and band
            FILETYPE = lightcurve[obs]['FILETYPE']
            BAND = lightcurve[obs]['BAND']
            DATE_BEG = lightcurve[obs]['DATE-BEG']

            # Initialize dictionary for per band
            if BAND not in dates:
                LOGGER.debug(f"Initializing dates/flux for {BAND}")
                dates[BAND] = {}
                flux_SCI[BAND] = {}
                if got_WGT:
                    flux_WGT[BAND] = {}

            # Initialize dictionary for nested dict for filetype
            for band in dates.keys():
                if FILETYPE not in dates[band]:
                    LOGGER.debug(f"Initializing dates/flux for {band}/{FILETYPE}")
                    dates[band][FILETYPE] = []
                    flux_SCI[band][FILETYPE] = []
                    if got_WGT:
                        flux_WGT[band][FILETYPE] = []

            # Get the date and store in list for [band][filter]
            dates[BAND][FILETYPE].append(DATE_BEG)
            # Store data in list keyed to filetype

            # Get the index for objID
            idx = lightcurve[obs]['objID'].index(objID)
            flux_sci = lightcurve[obs]['flux_SCI'][idx]
            flux_SCI[BAND][FILETYPE].append(flux_sci)
            if 'flux_WGT' in lightcurve[obs]:
                flux_wgt = lightcurve[obs]['flux_WGT'][idx]
                flux_WGT[BAND][FILETYPE].append(flux_wgt)

        # Put everything into a main dictionary
        LC[objID] = {}
        LC[objID]['dates'] = dates
        LC[objID]['flux_SCI'] = flux_SCI
        LC[objID]['flux_WGT'] = flux_WGT

    return LC


def write_lightcurve(args):

    d = datetime.datetime.today()
    date = d.isoformat('T', 'seconds')
    comment = f"# Lightcurve file created by: spt3g_cutter-{spt3g_cutter.__version__} on {date}\n"

    yaml_file = os.path.join(args.outdir, 'lightcurve.yaml')
    with open(yaml_file, 'w') as lightcurve_file:
        lightcurve_file.write(comment)
        yaml.dump(args.lc, lightcurve_file, sort_keys=False, default_flow_style=False)
    LOGGER.info(f"Wrote lightcurve file to: {yaml_file}")


def write_manifest(args):

    """Write YAML file with files created and input options"""

    ordered = ['bands', 'date_start', 'date_end', 'tablename', 'dbname', 'np', 'outdir',
               'inputList', 'yearly', 'files', 'id_names', 'size_on_disk',
               'JOB_ID', 'JOB_OUTPUT_DIR', 'files_on_disk', 'cutout_files',
               'rejected_positions', 'input_positions']
    manifest = {}

    d = args.__dict__
    for key in ordered:
        manifest[key] = d[key]

    d = datetime.datetime.today()
    date = d.isoformat('T', 'seconds')
    comment = f"# Manifest file created by: spt3g_cutter-{spt3g_cutter.__version__} on {date}\n"

    yaml_file = os.path.join(args.outdir, 'manifest.yaml')
    with open(yaml_file, 'w') as manifest_file:
        manifest_file.write(comment)
        yaml.dump(manifest, manifest_file, sort_keys=False, default_flow_style=False)
    LOGGER.info(f"Wrote manifest file to: {yaml_file}")


if __name__ == "__main__":

    # Example of inputs:
    # ra,dec can be list or scalars
    filename = 'yearly_150GHz_winter_2020_tonly.fits.fz'
    ra = [358.3406816, 355.8253677]
    dec = [-58.9660379, -57.1623017]

    xsize = [10]*len(ra)
    ysize = [10]*len(ra)

    # Create logger
    create_logger()
    logger = logging.getLogger(__name__)

    t0 = time.time()
    fitscutter(filename, ra, dec, xsize=xsize, ysize=ysize, units='arcmin',
               prefix=PREFIX)
    logger.info(f"Done: {elapsed_time(t0)}")
