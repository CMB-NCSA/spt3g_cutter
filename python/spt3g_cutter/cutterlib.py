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
import numpy as np
import pandas
import dateutil
from tempfile import mkdtemp
import errno
import shutil
import psutil
from astropy.io import fits
from astropy.time import Time

core_G3Units_deg = 0.017453292519943295
core_G3Units_rad = 1

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


def fitscutter(filename, ra, dec, cutout_names, rejected_positions, lightcurve,
               objID=None, xsize=1.0, ysize=1.0, units='arcmin', get_lightcurve=False,
               prefix=PREFIX, outdir=None, clobber=True, logger=None, counter='',
               get_uniform_coverage=False, nofits=False,
               stage=False, stage_prefix='spt-dummy'):

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
    get_uniform_coverage: bool
        Select only objects in the SPT uniform coverage

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
        raise Exception("ERROR: Cannot provide suitable DATE-BEG from SCI header")

    if 'DATE-END' in header['SCI']:
        date_end = str(header['SCI']['DATE-END']).strip()
    else:
        raise Exception("ERROR: Cannot provide suitable DATE-END from SCI header")

    # Get OBJECT, we will use as fieldname
    if 'OBJECT' in header['SCI']:
        object = str(header['SCI']['OBJECT']).strip()
    else:
        raise Exception("ERROR: Cannot provide suitable OBJECT from SCI header")
    # Check for object=None on yearly maps
    if object == 'None' and obsid.find('yearly') != -1:
        object = 'yearly'
        LOGGER.warning(f"Updating field to: {object}")

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
    rejected_ids = []
    
    # Define the ID for the lightcurve information from this filename
    if get_lightcurve:
        lcID = filename
        lc_local['DATE-BEG'] = date_beg
        lc_local['DATE-END'] = date_end
        lc_local['BAND'] = band
        lc_local['FILETYPE'] = filetype

    ######################################
    # Loop over ra/dec and xsize,ysize
    for k in range(len(ra)):

        # The basename for the (ra,dec)
        if objID[k] is None:
            objID[k] = get_thumbBaseName(ra[k], dec[k], prefix=prefix)

        # image and header sections
        im_section = OrderedDict()
        h_section = OrderedDict()

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

        # Check if in field extent
        if get_uniform_coverage and not in_uniform_coverage(ra[k], dec[k], object):
            LOGGER.debug(f"(RA,DEC):{ra[k]},{dec[k]} outside field extent")
            rejected.append(f"{ra[k]}, {dec[k]}, {objID[k]}")
            rejected_ids.append(objID[k])
            continue

        # Make sure the (x0,y0) is contained within the image
        if x0 < 0 or y0 < 0 or x0 > NAXIS1 or y0 > NAXIS2:
            LOGGER.debug(f"(RA,DEC):{ra[k]},{dec[k]} outside {filename}")
            LOGGER.debug(f"(x0,y0):{x0},{y0} > {NAXIS1},{NAXIS2}")
            rejected.append(f"{ra[k]}, {dec[k]}, {objID[k]}")
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

        LOGGER.debug(f"Working on object:{k} -- {objID[k]}")
        LOGGER.debug(f"Found naxis1,naxis2: {naxis1},{naxis2}")
        LOGGER.debug(f"Found x1,x2: {x1},{x2}")
        LOGGER.debug(f"Found y1,y2: {y1},{y2}")
        
        # Append data from (x0, y0) pixel from EXTNAME
        if get_lightcurve:
            HDUNUMW = hdunum['WGT']
            HDUNUMS = hdunum['SCI']
            try:
                data_extname = float(ifits[HDUNUMW][int(y0), int(x0)][0][0])
                if data_extname != 0.0:
                    lc_local.setdefault(f'flux_WGT', []).append(data_extname)
                    lc_local.setdefault(f'flux_SCI', []).append(float(ifits[HDUNUMS][int(y0), int(x0)][0][0]))
                else:
                    LOGGER.debug(f"(RA,DEC):{ra[k]},{dec[k]} zero flux weight")
                    rejected.append(f"{ra[k]}, {dec[k]}, {objID[k]}")
                    rejected_ids.append(objID[k])
                    continue
            except Exception as e:
                logger.error(e)
                data_extname = float("NaN")
            
            del data_extname

        for EXTNAME in extnames:
            # The hdunum for that extname
            HDUNUM = hdunum[EXTNAME]
            # Append data from (x0, y0) pixel from EXTNAME
            # Getting this out of the loop so we can do flux_wgt cut before it is added to the dictionnary for speed-up
            #if get_lightcurve:
            #    try:
            #        data_extname = float(ifits[HDUNUM][int(y0), int(x0)][0][0])
            #    except Exception as e:
            #        logger.error(e)
            #        data_extname = float("NaN")
            #    lc_local.setdefault(f'flux_{EXTNAME}', []).append(data_extname)
            #    del data_extname

            # Skip the fits part if notfits is true
            if nofits:
                continue
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

        # Skip the fits part if notfits is true
        if nofits:
            LOGGER.debug(f"Skipping FITS file creation for objID:{objID[k]} (RA,DEC):{ra[k]},{dec[k]}")
            continue

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
    logger.info(f"Done filename: {filename} in {elapsed_time(t1)} -- {counter}")

    # Assigning internal lists/dict to managed dictionaries
    cutout_names[filename] = outnames
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

    if len(rejected) > 0:
        rejected_positions[filename] = rejected
        logger.info(f"Rejected {len(rejected)} positions for {filename}")

    if stage:
        remove_staged_file(filename)

    # Clean up variables
    del ifits
    del outnames
    del rejected
    del lc_local
    del rejected_ids
    del im_section
    del h_section

    return cutout_names, rejected_positions, lightcurve


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
    #args.size_on_disk = get_size_on_disk(args.outdir)
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
    except (TypeError, dateutil.parser._parser.ParserError):
        date_mean = date1
        # add a warning for getting yearly map if required.
        # This should not be in the light curve
        LOGGER.debug(f"Ran into yearly  map: {date1}")
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
    LOGGER.info(f"Repacking lightcurve information for {BAND}, {FILETYPE}")
    LOGGER.info(f"Memory: {psutil.Process(os.getpid()).memory_info().rss / 1024 ** 3} Gb")
    process = psutil.Process(os.getpid())
    LOGGER.info(f"Memory percent: {process.memory_percent()} %")

    # Select only the observation for the BAND/FILETYPE combination
    observations = args.obs_dict[BAND][FILETYPE]

    LC = {}
    for objID in args.id_names:

        dates_ave = []
        dates_beg = []
        dates_end = []
        flux_SCI = []
        flux_WGT = []
        # Loop over the observations (OBS-ID + filetype)
        for obs in observations:

            if objID in lightcurve[obs]['rejected_ids']:
                LOGGER.debug(f"Ignoring {objID} for {obs} -- rejected")
                continue

            DATE_BEG = lightcurve[obs]['DATE-BEG']
            DATE_END = lightcurve[obs]['DATE-END']
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
            LC[objID]['flux_SCI'] = flux_SCI
            LC[objID]['flux_WGT'] = flux_WGT

    LOGGER.info(f"Done Re-packed lightcurve for {BAND}/{FILETYPE} in: {elapsed_time(t0)}")
    write_lightcurve_band_filetype(LC, BAND, FILETYPE, args)
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
    fits_file = os.path.join(args.outdir, f"lightcurve_{BAND}_{FILETYPE}.fits")
    LOGGER.info(f"Writing lightcurve to: {fits_file}")
    # Nested dictionaries cannot be sliced, so going through pandas route :(
    # as well as re-orienting
    df = pandas.DataFrame.from_dict(lc, orient='index')
    dict = df.to_dict()
    LOGGER.debug(f"Converted dictionary to pandas and back in: {elapsed_time(t0)}")
    col1 = fits.Column(name='id', format='30A', array=np.array(list(dict['id'].values()), dtype=object))
    col2 = fits.Column(name='dates_ave', format=f'PD({max_epochs})',
                       array=np.array(list(dict['dates_ave'].values()), dtype=object), unit='d, MJD')
    col3 = fits.Column(name='flux_SCI', format=f'PD({max_epochs})',
                       array=np.array(list(dict['flux_SCI'].values()), dtype=object), unit='mJy')
    col4 = fits.Column(name='flux_WGT', format=f'PD({max_epochs})',
                       array=np.array(list(dict['flux_WGT'].values()), dtype=object))
    hdu = fits.BinTableHDU.from_columns([col1, col2, col3, col4])
    hdu.header.set('TELESCOP', 'South Pole Telescope')
    hdu.header.set('INSTRUME', 'SPT-3G')
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

    """Write YAML file with files created and input options"""

    ordered = ['bands', 'date_start', 'date_end', 'tablename', 'dbname', 'np', 'outdir',
               'inputList', 'yearly', 'files', 'id_names', 'size_on_disk',
               'JOB_ID', 'JOB_OUTPUT_DIR', 'files_on_disk', 'cutout_files',
               'rejected_positions', 'input_positions']
    manifest = {}

    t0 = time.time()
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
    LOGGER.info(f"Wrote manifest file to: {yaml_file} in: {elapsed_time(t0)}")


def in_uniform_coverage(ra, dec, field):
    """Returns True/False if a (ra,dec) pair is in an spt3f field"""

    # Get the extent, using default padding
    ra_range, dec_range = get_field_extent(field, ra_pad=0, dec_pad=0)

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

    LOGGER.debug(f"RA:{ra}, DEC:{dec}")
    LOGGER.debug(f"range: {ra_range},{dec_range}")
    LOGGER.debug(f"in_field:{in_field}")

    return in_field


# Function lifted from:
# https://github.com/SouthPoleTelescope/spt3g_software/blob/master/sources/source_utils.py
def get_field_extent(
    field,
    ra_pad=3 * core_G3Units_deg,
    dec_pad=2 * core_G3Units_deg,
    sky_pad=True,
):
    """
    Get the extent of the given field.
    RA angles are always given between 0 and 360 degrees, with the left edge of
    the field given first. Dec angles are given between -90 and 90 degrees,with
    the lower edge of the field given first.
    For example, the SPT-3G winter field "ra0hdec-44.75" has ra range (310, 50)
    degrees, and dec range (-47.5, -42) degrees.
    The default padding parameters ensure that the extent of the SPT-3G fields
    covers the area that any bolometer on the focal plane may touch.  Setting
    the padding parameters to 0 instead returns the part of the field that has
    uniformly weighted coverage.  Padding in ra is defined in terms of degrees
    on the sky if ``sky_pad`` is True.
    The same padding considerations are not taken into account for the SPT-SZ
    and SPTpol fields.
    Arguments
    ---------
    field : str
        Name of the field
    ra_pad : float
        Padding to apply to extent in ra, in G3Units.  Values > 0 extend the
        field outward, and values < 0 shrink the field inward.  Set to 0 to
        select the uniformly weighted coverage region of the field.  If
        ``sky_pad`` is True, the padding is determined in terms of true degrees
        on the sky.
    dec_pad : float
        Padding to apply to extent in dec, in G3Units. Values > 0 extend the
        field outward, and values < 0 shrink the field inward.  Set to 0 to
        select the uniformly weighted coverage region of the field.
    sky_pad : bool
        If True, ``ra_pad`` is in terms of true degrees on the sky.  Otherwise,
        it is applied directly to the ra extent without correction.
    Returns
    -------
    ra_range : 2-tuple of floats
        Tuple of (min, max), the field extent in ra, in G3Units.
    dec_range : 2-tuple of floats, in G3Units
        Tuple of (min, max), the field extent in dec, in G3Units.
    """
    field = get_field_name(field)

    extents = {
        # sptsz fields
        "ra21hdec-42.5": ((300, 330), (-45, -40)),
        "ra21hdec-50": ((300, 330), (-55, -45)),
        "ra21hdec-60": ((300, 330), (-65, -55)),
        "ra22h30dec-55": ((330, 345), (-60, -50)),
        "ra23hdec-45": ((330, 360), (-50, -40)),
        "ra23hdec-62.5": ((330, 360), (-65, -60)),
        "ra23h30dec-55": ((345, 360), (-60, -50)),
        "ra0h50dec-50": ((0, 25), (-55, -45)),
        "ra1hdec-42.5": ((0, 30), (-45, -40)),
        "ra1hdec-60": ((0, 30), (-65, -55)),
        "ra2h30dec-50": ((25, 50), (-55, -45)),
        "ra3h30dec-42.5": ((30, 75), (-45, -40)),
        "ra3h30dec-60": ((30, 75), (-65, -55)),
        "ra4h10dec-50": ((50, 75), (-55, -45)),
        "ra5h30dec-45": ((75, 90), (-50, -40)),
        "ra5h30dec-55": ((75, 90), (-60, -50)),
        "ra6hdec-62.5": ((75, 105), (-65, -60)),
        "ra6h30dec-45": ((90, 105), (-50, -40)),
        "ra6h30dec-55": ((90, 105), (-60, -50)),
        "sptsz": ((300, 105), (-65, -40)),

        # sptpol fields
        "sptpol-100d": ((345, 360), (-60, -50)),
        "sptpol-500d": ((330, 30), (-65, -50)),
        "sptpol-kids": ((330, 52.5), (-36, -26)),
        "ra23hdec-25": ((330, 360), (-30, -20)),
        "ra23hdec-35": ((330, 360), (-40, -30)),
        "ra1hdec-25": ((0, 30), (-30, -20)),
        "ra1hdec-35": ((0, 30), (-40, -30)),
        "ra3hdec-25": ((30, 60), (-30, -20)),
        "ra3hdec-35": ((30, 60), (-40, -30)),
        "ra5hdec-25": ((60, 90), (-30, -20)),
        "ra5hdec-35": ((60, 90), (-40, -30)),
        "sptpol-ecs": ((330, 90), (-40, -20)),
        "ra11hdec-25": ((150, 180), (-30, -20)),
        "ra13hdec-25": ((180, 210), (-30, -20)),
        "sptpol-ecs-back": ((150, 210), (-30, -20)),

        # spt3g fields
        "ra0hdec-44.75": ((310, 50), (-47.5, -42)),
        "ra0hdec-52.25": ((310, 50), (-55, -49.5)),
        "ra0hdec-59.75": ((310, 50), (-62.5, -57)),
        "ra0hdec-67.25": ((310, 50), (-70, -64.5)),
        "spt3g-winter": ((310, 50), (-70, -42)),
        "ra5hdec-24.5": ((50, 100), (-27, -22)),
        "ra5hdec-31.5": ((50, 100), (-34, -29)),
        "ra5hdec-38.5": ((50, 100), (-41, -36)),
        "ra5hdec-45.5": ((50, 100), (-48, -43)),
        "ra5hdec-52.5": ((50, 100), (-55, -50)),
        "ra5hdec-59.5": ((50, 100), (-62, -57)),
        "ra5hdec-29.75": ((50, 100), (-30.5, -29)),
        "ra5hdec-33.25": ((50, 100), (-34, -32.5)),
        "ra5hdec-36.75": ((50, 100), (-37.5, -36)),
        "ra5hdec-40.25": ((50, 100), (-41, -39.5)),
        "spt3g-summera": ((50, 100), (-62, -29)),
        "ra1h40dec-29.75": ((0, 50), (-30.5, -29)),
        "ra1h40dec-33.25": ((0, 50), (-34, -32.5)),
        "ra1h40dec-36.75": ((0, 50), (-37.5, -36)),
        "ra1h40dec-40.25": ((0, 50), (-41, -39.5)),
        "spt3g-summerb": ((0, 50), (-41, -29)),
        "ra12h30dec-29.75": ((150, 225), (-30.5, -29)),
        "ra12h30dec-33.25": ((150, 225), (-34, -32.5)),
        "ra12h30dec-36.75": ((150, 225), (-37.5, -36)),
        "ra12h30dec-40.25": ((150, 225), (-41, -39.5)),
        "spt3g-summerc": ((150, 225), (-41, -29)),

        # Extra custom field for yearly
        "yearly": ((306.25, 53.25), (-72, -40.25)),
    }

    ra, dec = extents[field]

    # padding in true degrees ra
    deg = core_G3Units_deg
    if sky_pad:
        ra_pad /= np.cos(np.mean(dec) * deg / core_G3Units_rad)

    # apply padding
    ra = (
        (ra[0] * deg - ra_pad) % (360 * deg),
        (ra[1] * deg + ra_pad) % (360 * deg),
    )
    dec = (dec[0] * deg - dec_pad, dec[1] * deg + dec_pad)

    # Convert back to degrees
    ra = (ra[0]/deg, ra[1]/deg)
    dec = (dec[0]/deg, dec[1]/deg)

    return ra, dec


# Function lifted from:
# https://github.com/SouthPoleTelescope/spt3g_software/blob/master/sources/source_utils.py
def get_field_name(field):
    """
    Return the standard name for the given observing field.
    Arguments
    ---------
    field : str
        Field name, may be an alias
    Returns
    -------
    field : str
        Standard field name
    """
    field = field.lower()
    if field in ["spt3g", "winter"]:
        field = "spt3g-winter"
    elif field in ["summera", "summer-a", "summer", "spt3g-summer"]:
        field = "spt3g-summera"
    elif field in ["summerb", "summer-b"]:
        field = "spt3g-summerb"
    elif field in ["summerc", "summer-c"]:
        field = "spt3g-summerc"
    elif field in ["100d", "ra23h30dec-55"]:
        field = "sptpol-100d"
    elif field in ["500d", "ra0hdec-57.5"]:
        field = "sptpol-500d"
    elif field in ["kids", "ra0p75hdec-31"]:
        field = "sptpol-kids"
    elif field in ["ecs", "sptpol-summer"]:
        field = "sptpol-ecs"
    elif field in ["ecs-back", "sptpol-summer-back"]:
        field = "sptpol-ecs-back"
    return field


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
