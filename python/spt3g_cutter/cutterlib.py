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
import spt3g_cutter.fields as fields
import time
import numpy
import copy
from collections import OrderedDict
from astropy.wcs import WCS
# from astropy.utils.exceptions import AstropyWarning
import logging
from logging.handlers import RotatingFileHandler
# import warnings
import multiprocessing
import yaml
import datetime
import subprocess
import numpy as np
import pandas
import dateutil
from tempfile import mkdtemp
import errno
import shutil
import psutil
from astropy.io import fits
from astropy.time import Time
import json
from multiprocessing.managers import DictProxy

# To avoid header warning from astropy
# warnings.filterwarnings('ignore', category=AstropyWarning, append=True)

# Logger
LOGGER = logging.getLogger(__name__)

# Naming template
PREFIX = 'SPT3G'
OBJ_ID = "{prefix}J{ra}{dec}"
FITS_OUTNAME = "{outdir}/{objID}_{filter}_{obsid}_{filetype_ext}.{ext}"
LOG_OUTNAME = "{outdir}/{objID}.{ext}"
BASE_OUTNAME = "{objID}"
BASEDIR_OUTNAME = "{outdir}/{objID}"
FILETYPE_EXT = {'PSTH': 'psth', 'FLTD': 'fltd', 'CFLTD': 'cfltd', 'None': ''}
FITS_LC_OUTNAME = "{outdir}/lightcurve_{filter}_{filetype_ext}.{ext}"


def configure_logger(logger, MP=False, logfile=None, level=logging.NOTSET,
                     log_format=None, log_format_date=None):
    """
    Configure an existing logger with specified settings. Sets the format,
    logging level, and handlers for the given logger. If a logfile is provided,
    logs are written to both the console and the file with rotation. If no log
    format or date format is provided, default values are used.

    Parameters:
    - logger (logging.Logger): The logger to configure.
    - logfile (str, optional): Path to the log file. If `None`, logs to the console.
    - level (int): Logging level (e.g., `logging.INFO`, `logging.DEBUG`).
    - log_format (str, optional): Log message format (default is detailed format with function name).
    - log_format_date (str, optional): Date format for logs (default is `'%Y-%m-%d %H:%M:%S'`).
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

    # Update logger to see process inforation
    if MP is True:
        FORMAT = FORMAT.replace("[%(levelname)s]", "[%(levelname)s-%(processName)s]")

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
    if logger_has_stdout_handler(logger):
        print("Skipping adding stdout handler")
        return
    sh = logging.StreamHandler(sys.stdout)
    sh.setFormatter(formatter)
    sh.setLevel(level)
    handlers.append(sh)
    logger.addHandler(sh)
    return


def create_logger(logger=None, MP=False, logfile=None,
                  level=logging.NOTSET, log_format=None, log_format_date=None):
    """
    Configures and returns a logger with specified settings.
    Sets up logging based on provided level, format, and output file. Can be
    used for both `setup_logging` and other components.

    Parameters:
    - logger (logging.Logger, optional): The logger to configure. If `None`, a new logger
      is created.
    - logfile (str, optional): Path to the log file. If `None`, logs to the console.
    - level (int): Logging level (e.g., `logging.INFO`, `logging.DEBUG`).
    - log_format (str, optional): Format for log messages (e.g., `'%(asctime)s - %(message)s'`).
    - log_format_date (str, optional): Date format for logs (e.g., `'%Y-%m-%d %H:%M:%S'`).

    Returns:
    logging.Logger: The configured logger instance.

    Raises:
    - ValueError: If the log level or format is invalid.
    """

    if logger is None:
        logger = logging.getLogger(__name__)
    configure_logger(logger, MP=MP, logfile=logfile, level=level,
                     log_format=log_format, log_format_date=log_format_date)
    logging.basicConfig(handlers=logger.handlers, level=level)
    logger.propagate = False
    logger.info(f"Logging Created at level:{level}")
    return logger


def logger_has_stdout_handler(logger):
    """
    Check whether the given logger has a StreamHandler that
    writes to sys.stdout.
    Parameters:
     logger (logging.Logger): The logger instance to inspect.
    Returns:
     bool: True if a StreamHandler writing to sys.stdout is found,
           False otherwise.
    """
    for handler in logger.handlers:
        if isinstance(handler, logging.StreamHandler) and handler.stream == sys.stdout:
            return True
    return False


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


def update_wcs_matrix(header, x0, y0, proj='ZEA'):
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

    if proj == 'TAN':
        # Recompute CRVAL1/2 on the new center x0,y0
        CRVAL1, CRVAL2 = wcs.wcs_pix2world(x0, y0, 1)
        # Recast numpy objects as floats
        CRVAL1 = float(CRVAL1)
        CRVAL2 = float(CRVAL2)
        # Asign CRPIX1/2 on the new image
        CRPIX1 = 1
        CRPIX2 = 1
        # Update the values
        h['CRVAL1'] = CRVAL1
        h['CRVAL2'] = CRVAL2
        h['CRPIX1'] = CRPIX1
        h['CRPIX2'] = CRPIX2
        h['CTYPE1'] = 'RA---TAN'
        h['CTYPE2'] = 'DEC--TAN'
        # Delete some key that are not needed
        dkeys = ['PROJ', 'LONPOLE', 'LATPOLE', 'POLAR', 'ALPHA0', 'DELTA0', 'X0', 'Y0']
        for k in dkeys:
            h.delete(k)

    elif proj == 'ZEA':
        CRPIX1 = float(h['CRPIX1']) - x0
        CRPIX2 = float(h['CRPIX2']) - y0
        # Delete some key that are not needed
        dkeys = ['PROJ', 'POLAR', 'ALPHA0', 'DELTA0', 'X0', 'Y0']
        for k in dkeys:
            h.delete(k)
        h['CRPIX1'] = CRPIX1
        h['CRPIX2'] = CRPIX2
        LOGGER.debug(f"Update to CRPIX1:{CRPIX1}, CRPIX2:{CRPIX2}")

    else:
        raise NameError(f"Projection: {proj} not implemented")

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


def get_lightcurveFitsName(filter, filetype_ext,
                           refix=PREFIX, ext='fits', outdir=os.getcwd()):
    """ Common function to set the Fits lightcurve name """
    # Locals need to be captured at the end
    kw = locals()
    outname = FITS_LC_OUTNAME.format(**kw)
    return outname


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
    fits.close()
    return header, hdu


def get_NP(MP):

    """ Get the number of processors in the machine
    if MP == 0, use all available processor
    """
    # For it to be a integer
    MP = int(MP)
    if MP == 0:
        NP = int(multiprocessing.cpu_count())
    elif isinstance(MP, int):
        NP = MP
    else:
        raise ValueError('MP is wrong type: %s, integer type' % MP)
    return NP


def fitscutter(filename, ra, dec, cutout_names, rejected_names, lightcurve,
               objID=None, xsize=1.0, ysize=1.0, units='arcmin', get_lightcurve=False,
               prefix=PREFIX, outdir=None, clobber=True, logger=None, counter='',
               get_uniform_coverage=False, nofits=False,
               stage=False, stage_prefix='spt-dummy', obsid_names=False):

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
    rejected_names: dictionary
        Dict of rejected ids per filename with centers outside the FITS files
        or WGT=0
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
    get_uniform_coverage: bool
        Select only objects in the SPT uniform coverage

    Returns:
        cutout_names, rejected_names, lightcurve
    """

    # global timer for function
    t1 = time.time()
    if not logger:
        logger = LOGGER

    if not outdir:
        outdir = os.getcwd()

    # Check and fix inputs
    ra, dec, xsize, ysize, objID = check_inputs(ra, dec, xsize, ysize, objID)

    logger.info(f"Will cut: {len(ra)} stamps from FITS file: {filename} -- {counter}")

    # Check for the units
    if units == 'arcsec':
        scale = 1
    elif units == 'arcmin':
        scale = 60
    elif units == 'degree':
        scale = 3600
    else:
        raise Exception("ERROR: must define units as arcses/arcmin/degree only")

    # Stage if needed
    if stage:
        filename = stage_fitsfile(filename, stage_prefix=stage_prefix)

    # Get header/extensions/hdu
    t0 = time.time()
    header, hdunum = get_headers_hdus(filename)
    logger.debug(f"Done Getting header, hdus: {elapsed_time(t0)}")
    extnames = header.keys()  # Gets SCI and WGT
    logger.debug(f"Found EXTNAMES:{extnames}")

    # Fix DATE-END if set to None to avoid astropy warnings
    fix_date_END(header)

    # Get the dimensions of the parent image
    if 'EXTNAME' in header['SCI'] and header['SCI']['EXTNAME'].strip() == 'COMPRESSED_IMAGE':
        NAXIS1 = header['SCI']['ZNAXIS1']
        NAXIS2 = header['SCI']['ZNAXIS2']
    elif 'ZIMAGE' in header['SCI'] and header['SCI']['ZIMAGE'] is True:
        NAXIS1 = header['SCI']['ZNAXIS1']
        NAXIS2 = header['SCI']['ZNAXIS2']
    else:
        NAXIS1 = header['SCI']['NAXIS1']
        NAXIS2 = header['SCI']['NAXIS2']
    logger.debug(f"Found NAXIS1, NAXIS2: ({NAXIS1},{NAXIS2}) for: {filename}")

    # Get the pixel-scale of the input image
    pixelscale = astrometry.get_pixelscale(header['SCI'], units='arcsec')
    # Read in wcs with astropy
    wcs = WCS(header['SCI'])

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
    elif 'OBSID' in header['SCI']:
        obsid = str(header['SCI']['OBSID']).strip()
    else:
        raise Exception("ERROR: Cannot provide suitable OBS-ID from SCI header")

    # Extract FILETYPE from the header
    if 'FILETYPE' in header['SCI']:
        filetype = str(header['SCI']['FILETYPE']).strip()
    else:
        filetype = 'None'

    if 'DATE-BEG' in header['SCI']:
        date_beg = str(header['SCI']['DATE-BEG']).strip()
    else:
        raise Exception("ERROR: Cannot provide suitable DATE-BEG from SCI header")

    if 'DATE-END' in header['SCI']:
        date_end = str(header['SCI']['DATE-END']).strip()
    else:
        raise Exception("ERROR: Cannot provide suitable DATE-END from SCI header")

    # Get OBJECT, we will use as fieldname
    if 'FIELD' in header['SCI']:
        field = str(header['SCI']['FIELD']).strip()
    else:
        raise Exception("ERROR: Cannot provide suitable FIELD from SCI header")

    # # Check for object=None on yearly maps
    # if field == 'None' and obsid.find('yearly') != -1:
    # #    LOGGER.warning(f"Updating field to: {field}")

    # The extension to use for FILETYPE
    filetype_ext = FILETYPE_EXT[filetype]

    # Intitialize the FITS object
    t0 = time.time()
    ifits = fitsio.FITS(filename, 'r')
    logger.debug(f"Done loading fitsio.FITS({filename}): {elapsed_time(t0)}")

    if cutout_names is None:
        cutout_names = {}
    if rejected_names is None:
        rejected_names = {}
    if lightcurve is None:
        lightcurve = {}

    # Local lists/dicts
    outnames = []
    lc_local = {}
    rejected_ids = []

    # Define the ID for the lightcurve information from this filename
    if get_lightcurve:
        lcID = filename
        lc_local['DATE-BEG'] = date_beg
        lc_local['DATE-END'] = date_end
        lc_local['BAND'] = band
        lc_local['FILETYPE'] = filetype
        lc_local['OBSID'] = obsid

    ######################################
    # Loop over ra/dec and xsize,ysize
    for k in range(len(ra)):

        # The basename for the (ra,dec)
        if objID[k] is None:
            objID[k] = get_thumbBaseName(ra[k], dec[k], prefix=prefix)
        logger.debug(f"Defined objID: {objID[k]}")
        # image and header sections
        im_section = OrderedDict()
        h_section = OrderedDict()

        # Check if in field extent
        if get_uniform_coverage and not in_uniform_coverage(ra[k], dec[k], field):
            LOGGER.warning(f"position:{k} (RA,DEC):{ra[k]},{dec[k]} outside field: {field}")
            rejected_ids.append(objID[k])
            continue

        # Define the geometry of the thumbnail
        x0, y0 = wcs.wcs_world2pix(ra[k], dec[k], 0)
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
            LOGGER.warning(f"position:{k} (RA,DEC):{ra[k]},{dec[k]} outside {filename}")
            LOGGER.warning(f"position:{k} (x0,y0):{x0},{y0} > {NAXIS1},{NAXIS2}")
            rejected_ids.append(objID[k])
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

        LOGGER.debug(f"Working on position:{k} -- {objID[k]}")
        LOGGER.debug(f"Defined stamp naxis1,naxis2: {naxis1},{naxis2}")
        LOGGER.debug(f"Defined stamp x1,x2: {x1},{x2}")
        LOGGER.debug(f"Defined stamp y1,y2: {y1},{y2}")

        # Append data from (x0, y0) pixel for both extensions
        if get_lightcurve:
            HDU_SCI = hdunum['SCI']
            HDU_WGT = hdunum['WGT']
            try:
                data_WGT = float(ifits[HDU_WGT][int(y0), int(x0)][0][0])
                data_SCI = float(ifits[HDU_SCI][int(y0), int(x0)][0][0])
                if data_WGT != 0.0:
                    data_SCI = float(ifits[HDU_SCI][int(y0), int(x0)][0][0])
                else:
                    LOGGER.warning(f"position:{k} (RA,DEC):{ra[k]},{dec[k]} has zero flux weight: {data_WGT}")
                    rejected_ids.append(objID[k])
                    continue
            except Exception as e:
                logger.error(e)
                data_SCI = float("NaN")
                data_WGT = float("NaN")

            lc_local.setdefault('flux_WGT', []).append(data_WGT)
            lc_local.setdefault('flux_SCI', []).append(data_SCI)

            del data_SCI
            del data_WGT

        # Skip the fits part if notfits is true
        if nofits:
            LOGGER.debug(f"Skipping FITS file creation for objID:{objID[k]} (RA,DEC):{ra[k]},{dec[k]}")
            continue

        # Now we cut the fits stamp
        for EXTNAME in extnames:
            # The hdunum for that extname
            HDUNUM = hdunum[EXTNAME]
            # Create a canvas
            im_section[EXTNAME] = numpy.zeros((naxis1, naxis2))
            # Read in the image section we want for SCI/WGT
            im_section[EXTNAME] = ifits[HDUNUM][int(y1):int(y2), int(x1):int(x2)]
            # Correct NAXIS1 and NAXIS2
            naxis1 = numpy.shape(im_section[EXTNAME])[1]
            naxis2 = numpy.shape(im_section[EXTNAME])[0]
            # Update the WCS in the headers and make a copy
            h_section[EXTNAME] = update_wcs_matrix(header[EXTNAME], x1, y1)
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

        # Keep the objID information
        if obsid_names:
            if objID[k] not in cutout_names:
                cutout_names[objID[k]] = {}
            if band not in cutout_names[objID[k]].keys():
                cutout_names[objID[k]][band] = []
            cutout_names[objID[k]][band].append(outname)

        # Write out the file
        t0 = time.time()
        ofits = fitsio.FITS(outname, 'rw', clobber=clobber)
        for EXTNAME in extnames:
            ofits.write(im_section[EXTNAME], extname=EXTNAME, header=h_section[EXTNAME])
        ofits.close()
        logger.debug(f"Done writing position:{k} to {outname}: {elapsed_time(t0)}")

    ifits.close()
    logger.info(f"Done filename: {filename} in {elapsed_time(t1)} -- {counter}")

    # Assigning internal lists/dict to managed dictionaries
    if obsid_names is False:
        cutout_names[filename] = outnames
    rejected_names[filename] = rejected_ids

    if get_lightcurve:
        # Remove the rejected ids from objID list,
        # Otherwise index search will fail
        for id in rejected_ids:
            logger.debug(f"Removing rejected id:{id} from lightcurve[objID]")
            objID.remove(id)
        # We add the objID array after we pruned it from rejected ids
        lc_local['objID'] = objID
        lc_local['rejected_ids'] = rejected_ids
        lightcurve[lcID] = lc_local
    if len(rejected_ids) > 0:
        logger.info(f"Rejected {len(rejected_ids)} positions for {filename}")

    if stage:
        remove_staged_file(filename)

    # Clean up variables
    del ifits
    del outnames
    del lc_local
    del rejected_ids
    del im_section
    del h_section

    return cutout_names, rejected_names, lightcurve


def get_id_names(ra, dec, prefix):
    "Get the ID names associated with every position"
    names = []
    for k in range(len(ra)):
        names.append(get_thumbBaseName(ra[k], dec[k], prefix=prefix))
    return names


def get_size_on_disk(outdir, timeout=15):
    "Get the size of the outdir outputs"
    t0 = time.time()
    LOGGER.info(f"Getting size_on_disk with timeout={timeout}s.")
    try:
        size = subprocess.check_output(['du', '-sh', outdir], timeout=timeout).split()[0].decode('ascii')
    except subprocess.TimeoutExpired:
        LOGGER.warning(f"Cannot get_size_on_disk, timeout after {timeout}s.")
        size = f"Timed out: {timeout} sec, too large to compute"
    LOGGER.info(f"Done size_on_disk in: {elapsed_time(t0)}")
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

    LOGGER.info("Getting job metadata for manifest file")

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
    # args.size_on_disk = get_size_on_disk(args.outdir)
    args.size_on_disk = None
    args.files_on_disk = len(args.cutout_files)

    # Get the job information from k8s
    (args.JOB_ID, args.JOB_OUTPUT_DIR) = get_job_info(args)
    return args


def get_mean_date(date1, date2):
    """ Gets the mean date betwenn to timestamps"""
    # Need to try/except to catch dates for yearly maps
    try:
        D1 = pandas.to_datetime(date1)
        D2 = pandas.to_datetime(date2)
        date_mean = pandas.Timestamp((D1.value + D2.value)/2.).isoformat()
    except (TypeError, dateutil.parser._parser.ParserError, pandas._libs.tslibs.parsing.DateParseError):
        date_mean = date1
        # add a warning for getting yearly map if required.
        # This should not be in the light curve
        LOGGER.warning(f"Ran into yearly map: {date1}")
    return date_mean


def get_obs_dictionary(lightcurve):
    "Create a dictionary of obervations keyed to BAND and FILETYPE"

    LOGGER.info("Creating dictionary with observations")
    obs_dict = {}
    for obs in lightcurve:
        FILETYPE = lightcurve[obs]['FILETYPE']
        BAND = lightcurve[obs]['BAND']
        if BAND not in obs_dict:
            obs_dict[BAND] = {}
        if FILETYPE not in obs_dict[BAND]:
            obs_dict[BAND][FILETYPE] = []
        obs_dict[BAND][FILETYPE].append(obs)
    return obs_dict


def repack_lightcurve_band_filetype(lightcurve, BAND, FILETYPE, args):
    "Repack the lightcurve dictionary keyed by objID"

    t0 = time.time()
    LOGGER.info(f"Repacking lightcurve information for band: {BAND}, filetype: {FILETYPE}")
    LOGGER.debug(f"Memory: {psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3} Gb")
    process = psutil.Process(os.getpid())
    LOGGER.debug(f"Memory percent: {process.memory_percent()} %")

    # Select only the observation for the BAND/FILETYPE combination
    observations = args.obs_dict[BAND][FILETYPE]

    LC = {}
    for objID in args.id_names:

        dates_ave = []
        dates_beg = []
        dates_end = []
        flux_SCI = []
        flux_WGT = []
        obsids = []
        # Loop over the observations (OBS-ID + filetype)
        for obs in observations:

            if objID in lightcurve[obs]['rejected_ids']:
                LOGGER.debug(f"Ignoring {objID} for {obs} -- rejected")
                continue

            OBSID = lightcurve[obs]['OBSID']
            DATE_BEG = lightcurve[obs]['DATE-BEG']
            DATE_END = lightcurve[obs]['DATE-END']
            LOGGER.debug(f"{obs} DATE_BEG: {DATE_BEG}")
            LOGGER.debug(f"{obs} DATE_END: {DATE_END}")
            DATE_AVE = get_mean_date(DATE_BEG, DATE_END)
            if DATE_AVE.find('yearly') != -1:
                LOGGER.debug(f"Ignoring {objID} for {obs} -- yearly map")
                continue

            # Get the index for objID
            idx = lightcurve[obs]['objID'].index(objID)
            try:
                flux_wgt = lightcurve[obs]['flux_WGT'][idx]
                # Only store if flux is > 0
                if flux_wgt > 0:
                    flux_WGT.append(flux_wgt)
                    # storing dates
                    dates_ave.append(DATE_AVE)
                    dates_beg.append(DATE_BEG)
                    dates_end.append(DATE_END)
                    # storing obsid
                    obsids.append(OBSID)
                    # storing flux
                    flux_sci = lightcurve[obs]['flux_SCI'][idx]
                    flux_SCI.append(flux_sci)
            except KeyError:
                flux_wgt = None
                LOGGER.warning(f"NO flux_WGT - obs:{objID} date:{DATE_BEG} BAND:{BAND} FILETYPE: {FILETYPE}")

        # Put everything into a main dictionary, only if we get any hits
        # since now zero weights have been removed, need to be smarter to
        # include the objID into rejected for lightcurve.
        # Talk to Felipe about this.
        if len(flux_WGT) > 0:
            LC[objID] = {}
            LC[objID]['id'] = objID
            LC[objID]['dates_ave'] = Time(dates_ave).mjd     # converting the date array to mjd
            # LC[objID]['dates_beg'] = dates_beg
            # LC[objID]['dates_end'] = dates_end
            LC[objID]['obsids'] = obsids
            LC[objID]['flux_SCI'] = flux_SCI
            LC[objID]['flux_WGT'] = flux_WGT

    LOGGER.info(f"Done Re-packed lightcurve for {BAND}/{FILETYPE} in: {elapsed_time(t0)}")
    if len(LC) == 0:
        LOGGER.warning(f"Lightcurve for {BAND}/{FILETYPE} is empty -- will not write lightcurve table")
    else:
        write_lightcurve_band_filetype(LC, BAND, FILETYPE, args)
    del lightcurve
    return


def get_rejected_ids(args):
    "Parse the rejected_positions and get list of rejected ids"

    rejected_ids = []
    for key, list in args.rejected_positions.items():
        for item in list:
            id = item.split(', ')[2]
            if id not in rejected_ids:
                rejected_ids.append(id)
    LOGGER.info(f"Found {len(rejected_ids)} objID to reject")
    return rejected_ids


def write_lightcurve_band_filetype(lc, BAND, FILETYPE, args):

    t0 = time.time()
    max_epochs = 15000  # this has maximum number of epochs as 15k for fits table format
    fits_file = get_lightcurveFitsName(BAND, FILETYPE_EXT[FILETYPE], outdir=args.outdir)

    # LOGGER.info(f"Writing lightcurve to: {fits_file}")
    # Nested dictionaries cannot be sliced, so going through pandas route :(
    # as well as re-orienting
    df = pandas.DataFrame.from_dict(lc, orient='index')
    dict = df.to_dict()
    LOGGER.debug(f"Converted dictionary to pandas and back in: {elapsed_time(t0)}")
    col1 = fits.Column(name='id', format='30A', array=np.array(list(dict['id'].values()), dtype=object))
    col2 = fits.Column(name='dates_ave', format=f'PD({max_epochs})',
                       array=np.array(list(dict['dates_ave'].values()), dtype=object), unit='days, MJD')
    col3 = fits.Column(name='flux_SCI', format=f'PD({max_epochs})',
                       array=np.array(list(dict['flux_SCI'].values()), dtype=object), unit='mJy')
    col4 = fits.Column(name='flux_WGT', format=f'PD({max_epochs})',
                       array=np.array(list(dict['flux_WGT'].values()), dtype=object))
    col5 = fits.Column(name='obsids', format=f'PD({max_epochs})',
                       array=np.array(list(dict['obsids'].values()), dtype=object))
    hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4, col5])
    hdu.header.set('TELESCOP', 'South Pole Telescope')
    hdu.header.set('INSTRUME', 'SPT-3G')
    hdu.header.set('BAND', BAND)

    hdu.writeto(fits_file, overwrite=True)
    LOGGER.info(f"Wrote lightcurve file to: {fits_file} in: {elapsed_time(t0)}")


def write_lightcurve(args):

    t0 = time.time()
    d = datetime.datetime.today()
    date = d.isoformat('T', 'seconds')
    comment = f"# Lightcurve file created by: spt3g_cutter-{spt3g_cutter.__version__} on {date}\n"

    yaml_file = os.path.join(args.outdir, 'lightcurve.yaml')
    with open(yaml_file, 'w') as lightcurve_file:
        lightcurve_file.write(comment)
        yaml.dump(args.lc, lightcurve_file, sort_keys=False, default_flow_style=False)
    LOGGER.info(f"Wrote lightcurve file to: {yaml_file} in: {elapsed_time(t0)}")


def write_manifest(args):

    """Write file with files created and input options"""

    ordered = ['bands', 'date_start', 'date_end', 'tablename', 'dbname', 'np', 'outdir',
               'inputList', 'yearly', 'files', 'id_names', 'size_on_disk',
               'JOB_ID', 'JOB_OUTPUT_DIR', 'files_on_disk', 'cutout_files',
               'rejected_names', 'input_positions']
    manifest = {}

    t0 = time.time()

    dt = datetime.datetime.today()
    date = dt.isoformat('T', 'seconds')
    comment = f"Manifest file created by: spt3g_cutter-{spt3g_cutter.__version__} on {date}"
    d = args.__dict__
    manifest["metadata"] = comment
    for key in ordered:
        if isinstance(d[key], DictProxy):
            manifest[key] = d[key]._getvalue()
        else:
            manifest[key] = d[key]
    json_file = os.path.join(args.outdir, 'manifest.json')
    LOGGER.info(f"writing manifest to: {json_file}")
    with open(json_file, 'w') as manifest_file:
        json.dump(manifest, manifest_file, sort_keys=False, indent=6)
    LOGGER.info(f"Wrote manifest file to: {json_file} in: {elapsed_time(t0)}")


def in_uniform_coverage(ra, dec, field):
    """Returns True/False if a (ra,dec) pair is in an spt3f field"""

    # Get the extent, using default padding
    ra_range, dec_range = fields.get_field_extent(field, ra_pad=0, dec_pad=0)

    ra_range = list(ra_range)
    dec_range = list(dec_range)

    LOGGER.debug(f"RA:{ra}, DEC:{dec}")
    LOGGER.debug(f"range: {ra_range},{dec_range}")

    # Check to see if it crosses RA=0
    if ra_range[0] > ra_range[1]:
        crossRA0 = True
    else:
        crossRA0 = False

    if crossRA0 and ra_range[0] > 180:
        ra_range[0] = ra_range[0] - 360
    if crossRA0 and ra > 180:
        ra = ra - 360

    if ra_range[0] < ra < ra_range[1] and dec_range[0] < dec < dec_range[1]:
        in_field = True
    else:
        in_field = False

    LOGGER.debug(f"crossRA0: {crossRA0}")
    LOGGER.debug(f"field:{field}, (ra/dec) range: {ra_range},{dec_range}")
    LOGGER.debug(f"RA:{ra}, DEC:{dec} in_field:{field} -- {in_field}")

    return in_field


def create_dir(dirname):
    "Safely attempt to create a folder"
    if not os.path.isdir(dirname):
        LOGGER.info(f"Creating directory {dirname}")
        try:
            os.makedirs(dirname, mode=0o755, exist_ok=True)
        except OSError as e:
            if e.errno != errno.EEXIST:
                LOGGER.warning(f"Problem creating {dirname} -- proceeding with trepidation")


def stage_fitsfile(fitsfile, stage_prefix="spt", use_cp=False):
    """
    Stage input fitsfile to the stage directory
    """
    tmp_dir = mkdtemp(prefix=stage_prefix)
    fitsfile_copy = os.path.join(tmp_dir, os.path.basename(fitsfile))
    LOGGER.info(f"Will stage: {fitsfile} --> {fitsfile_copy}")
    # Make sure that the folder exists:
    create_dir(os.path.dirname(fitsfile_copy))
    if use_cp:
        LOGGER.warning("Will use system copy call")
        cmd = f"cp -pv {fitsfile} {fitsfile_copy}"
        os.system(cmd)
    else:
        shutil.copy2(fitsfile, fitsfile_copy)
    return fitsfile_copy


def remove_staged_file(fitsfile):
    LOGGER.info(f"Removing: {fitsfile}")
    os.remove(fitsfile)
    tmp_dir = os.path.dirname(fitsfile)
    LOGGER.info(f"Removing tmp dir: {tmp_dir}")
    shutil.rmtree(tmp_dir)


def fix_date_END(header):
    extnames = header.keys()  # Gets SCI and WGT
    for EXT in extnames:
        if 'DATE-END' in header[EXT] and header[EXT]['DATE-END'] is None:
            header[EXT]["DATE-END"] = header[EXT]["DATE-BEG"]
            LOGGER.debug(f"Updating DATE-END to DATE-BEG for {EXT}")


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
