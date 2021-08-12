#!/usr/bin/env python

"""
A set of simple proto-function to make postage stamps using fitsio
Based on desthumbs code circa 2015
"""

import fitsio
import os
import sys
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

# To avoid header warning from astropy
warnings.filterwarnings('ignore', category=AstropyWarning, append=True)

# Logger
LOGGER = logging.getLogger(__name__)

# Naming template
PREFIX = 'SPT3G'
FITS_OUTNAME = "{outdir}/{prefix}J{ra}{dec}_{filter}_{obsid}.{ext}"
LOG_OUTNAME = "{outdir}/{prefix}J{ra}{dec}.{ext}"
BASE_OUTNAME = "{prefix}J{ra}{dec}"
BASEDIR_OUTNAME = "{outdir}/{prefix}J{ra}{dec}"


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


def check_inputs(ra, dec, xsize, ysize):

    """ Check and fix inputs for cutout"""
    # Make sure that RA,DEC are the same type
    if type(ra) != type(dec):
        raise TypeError('RA and DEC need to be the same type()')
    # Make sure that XSIZE, YSIZE are same type
    if type(xsize) != type(ysize):
        raise TypeError('XSIZE and YSIZE need to be the same type()')
    # Make them iterable and proper length
    if hasattr(ra, '__iter__') is False and hasattr(dec, '__iter__') is False:
        ra = [ra]
        dec = [dec]
    if hasattr(xsize, '__iter__') is False and hasattr(ysize, '__iter__') is False:
        xsize = [xsize]*len(ra)
        ysize = [ysize]*len(ra)
    # Make sure they are all of the same length
    if len(ra) != len(dec):
        raise TypeError('RA and DEC need to be the same length')
    if len(xsize) != len(ysize):
        raise TypeError('XSIZE and YSIZE need to be the same length')
    if (len(ra) != len(xsize)) or (len(ra) != len(ysize)):
        raise TypeError('RA, DEC and XSIZE and YSIZE need to be the same length')
    return ra, dec, xsize, ysize


def get_thumbFitsName(ra, dec, filter, obsid, prefix=PREFIX, ext='fits', outdir=os.getcwd()):
    """ Common function to set the Fits thumbnail name """
    ra = astrometry.dec2deg(ra/15., sep="", plussign=False)
    dec = astrometry.dec2deg(dec, sep="", plussign=True)
    kw = locals()
    outname = FITS_OUTNAME.format(**kw)
    return outname


def get_thumbBaseDirName(ra, dec, prefix=PREFIX, outdir=os.getcwd()):
    """ Common function to set the Fits thumbnail name """
    ra = astrometry.dec2deg(ra/15., sep="", plussign=False)
    dec = astrometry.dec2deg(dec, sep="", plussign=True)
    kw = locals()
    basedir = BASEDIR_OUTNAME.format(**kw)
    return basedir


def get_thumbLogName(ra, dec, prefix=PREFIX, ext='log', outdir=os.getcwd()):
    """ Common function to set the Fits thumbnail name """
    ra = astrometry.dec2deg(ra/15., sep="", plussign=False)
    dec = astrometry.dec2deg(dec, sep="", plussign=True)
    kw = locals()
    outname = LOG_OUTNAME.format(**kw)
    return outname


def get_thumbBaseName(ra, dec, prefix=PREFIX):
    """ Common function to set the Fits thumbnail name """
    ra = astrometry.dec2deg(ra/15., sep="", plussign=False)
    dec = astrometry.dec2deg(dec, sep="", plussign=True)
    kw = locals()
    outname = BASE_OUTNAME.format(**kw)
    return outname


def get_headers_hdus(filename):

    header = OrderedDict()
    hdu = OrderedDict()

    # Case 1 -- for well-defined fitsfiles with EXTNAME
    with fitsio.FITS(filename) as fits:
        for k in range(len(fits)):
            h = fits[k].read_header()
            # Make sure that we can get the EXTNAME
            if not h.get('EXTNAME'):
                continue
            extname = h['EXTNAME'].strip()
            if extname == 'COMPRESSED_IMAGE':
                continue
            header[extname] = h
            hdu[extname] = k

    # Case 2 -- older DESDM files without EXTNAME
    if len(header) < 1:
        sci_hdu, wgt_hdu = get_fits_hdu_extensions_byfilename(filename)
        fits = fitsio.FITS(filename)
        header['SCI'] = fits[sci_hdu].read_header()
        header['WGT'] = fits[wgt_hdu].read_header()
        hdu['SCI'] = sci_hdu
        hdu['WGT'] = wgt_hdu

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


def fitscutter(filename, ra, dec, xsize=1.0, ysize=1.0, units='arcmin', prefix=PREFIX,
               outdir=None, clobber=True, logger=None, counter=''):

    """
    Makes cutouts around ra, dec for a give xsize and ysize
    ra,dec can be scalars or lists/arrays
    """

    t0 = time.time()
    if not logger:
        logger = LOGGER

    if not outdir:
        outdir = os.getcwd()

    # Check and fix inputs
    ra, dec, xsize, ysize = check_inputs(ra, dec, xsize, ysize)

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
    header, hdunum = get_headers_hdus(filename)
    extnames = header.keys()  # Gets SCI and WGT
    logger.debug(f"Found EXTNAMES:{extnames}")

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
    else:
        raise Exception("ERROR: Cannot provide suitable OBS-ID from SCI header")

    # Intitialize the FITS object
    ifits = fitsio.FITS(filename, 'r')

    ######################################
    # Loop over ra/dec and xsize,ysize
    for k in range(len(ra)):

        # Define the geometry of the thumbnail
        x0, y0 = wcs.wcs_world2pix(ra[k], dec[k], 1)
        yL = 10000
        xL = 10000
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

        if y1 < 0:
            y1 = 0
            dy = y0
            y2 = y0 + dy
        if y2 > yL:
            y2 = yL
            dy = yL - y0
            y1 = y0-dy

        if x1 < 0:
            x1 = 0
            dx = x0
            x2 = x0 + dx
        if x2 > xL:
            x2 = xL
            dx = xL - x0
            x1 = x0 - dx

        im_section = OrderedDict()
        h_section = OrderedDict()
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
            h_section[EXTNAME] = update_wcs_matrix(header[EXTNAME], x0, y0, naxis1, naxis2, ra[k], dec[k])

        # Get the basedir
        basedir = get_thumbBaseDirName(ra[k], dec[k], prefix=prefix, outdir=outdir)
        if not os.path.exists(basedir):
            os.makedirs(basedir)

        # Construct the name of the Thumbmail using BAND/FILTER/prefix/etc
        outname = get_thumbFitsName(ra[k], dec[k], band, obsid, prefix=prefix, outdir=basedir)

        # Write out the file
        ofits = fitsio.FITS(outname, 'rw', clobber=clobber)
        for EXTNAME in extnames:
            ofits.write(im_section[EXTNAME], extname=EXTNAME, header=h_section[EXTNAME])

        ofits.close()
        logger.debug(f"Wrote: {outname}")

    logger.info(f"Done {filename} in {elapsed_time(t0)} -- {counter}")
    return


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
