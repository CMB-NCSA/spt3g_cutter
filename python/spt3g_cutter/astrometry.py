"""
A collection of useful functions in astrometry. The functions ported
here correspond to a subset of inhereted from Felipe Menanteau's
astrometry.py old library. Using wcs/header transformations from astropy.wcs

The functions will:
     - format decimal <---> DDMMSS/HHMMMSS
     - greater circle distance(ra, dec)
     - area in polygon

 Requires:
     numpy

 Felipe Menanteau, Apr/Oct 2014
"""

import math
import numpy


def circle_distance(ra1, dec1, ra2, dec2, units='deg'):
    """
    Calculates great-circle distances between the two points that is,
    the shortest distance over the earth's surface using the Haversine
    formula, see http://www.movable-type.co.uk/scripts/latlong.html

    Parameters:
        ra1: float
            RA of the first point
        dec1: float
            DEC of the first point
        ra2: float
            RA of the second point
        dec2: float
            DEC of the secons point
        units: str
            The units of the coordinates. If 'deg' (default) then input and
            output are in degrees. Anything else is assumed to be radians.

    Returns:
        float of the great circle distance in the same units as the inputs.
    """

    cos = numpy.cos
    sin = numpy.sin
    acos = numpy.arccos

    if units == 'deg':
        ra1 = ra1 * math.pi / 180.
        ra2 = ra2 * math.pi / 180.
        dec1 = dec1 * math.pi / 180.
        dec2 = dec2 * math.pi / 180.

    x = sin(dec1) * sin(dec2) + cos(dec1) * cos(dec2) * cos(ra2 - ra1)
    x = numpy.where(x > 1.0, 1, x)  # Avoid x>1.0 values

    d = acos(x)
    if units == 'deg':
        d = d * 180.0 / math.pi
    return d


def deg2dec(deg, sep=":"):
    """
    Degrees to decimal, one element or list/array object.
    Parameters:
        deg: str
            The sexagesimal formatted string, must contain entries for degrees,
            minutes, and seconds. Can also be an array/list of strings
        sep: str
            The separator of the different elements. Default is ':'

    Returns:
        float representation of the input(s) in degrees. The return will be in
        the same format as the input (a single input will return a single
        output, an array of inputs will return an array of the same length)

    """
    if isinstance(deg, (list, numpy.ndarray, numpy.recarray)):
        return [deg2dec_one(d, sep=sep) for d in deg]

    return deg2dec_one(deg, sep=sep)


def deg2dec_one(deg, sep=":"):
    """
    Degrees to decimal, one element only.

    Parameters:
        deg: str
            The sexagesimal formatted string, must contain entries for degrees,
            minutes, and seconds.
        sep: str
            The separator of the different elements. Default is ':'

    Returns:
        float representation of the input in degrees.
    """
    vals = deg.split(sep)
    dd = float(vals[0])
    mm = float(vals[1]) / 60.
    ss = float(vals[2]) / 3600.
    if dd < 0 or vals[0] == '-00' or vals[0] == '-0':
        mm = -mm
        ss = -ss
    return dd + mm + ss


def dec2deg(dec, sep=":", plussign=False, short=False, sectol=1e-3):
    """
    From decimal to degrees, array or scalar

    Parameters:
        dec: float
            Decimal representation of the coordinate. Can be an array/list of
            floats.
        sep: str
            The separator to use when formatting the string. Default is ':'
        plussign: bool
            Whether to prepend a plus sign to the output for positive values.
            Default is False (do not prepend)
        short: misc
            What format to return the value(s) in. 'ra' will format it for
            right ascension, any other value that evaluates to True will only
            return degrees and minutes. Default is False.
        sectol: float
            The tolerance to round the seconds to the next minute.
            Default if 1e-3

    Returns:
        Sexagesimal formatted representation of the input. The type of output
        will match that of the input (e.g. a single input will have a single
        output, a list of inputs will have a list of outputs with the same
        length)
    """

    # Make it a numpy object if iterable
    if isinstance(dec, (list, numpy.ndarray, numpy.recarray)):
        dec = numpy.asarray(dec)
        # Keep the sign for later
        sig = numpy.where(dec < 0, -1, +1)
        dec = abs(dec)
        dd = abs(dec.astype("Int32"))
        mm = (abs(dec - dd) * 60).astype("Int32")
        ss = (abs(dec - dd) * 60 - mm) * 60
        # Truncating ss < 0.001
        ids = numpy.where(abs(ss - 60.) <= sectol)
        ss[ids] = 0.0
        mm[ids] = mm[ids] + 1

        # Make an numpy array structures -- probably unnecessary
        x = numpy.concatenate((sig, dd, mm, ss))
        # re-shape
        x = numpy.resize(x, (4, len(dec)))
        x = numpy.swapaxes(x, 0, 1)
        return [format_deg(element, short=short, sep=sep, plussign=plussign) for element in x]

    sig = math.copysign(1, dec)
    dd = int(dec)
    mm = int(abs(dec - dd) * 60.)
    ss = (abs(dec - dd) * 60 - mm) * 60

    if float(abs(ss - 60.)) < sectol:
        ss = 0.0
        mm = mm + 1
    return format_deg((sig, dd, mm, ss), short=short, sep=sep, plussign=plussign)


def format_deg(x, short=False, sep=":", plussign=False):
    """
    From decimal to degrees

    Parameters:
        x: float
            Decimal representation of the coordinate. Can be an array/list of
            floats.
        short: misc
            What format to return the value(s) in. 'ra' will format it for
            right ascension, any other value that evaluates to True will only
            return degrees and minutes. Default is False.
        sep: str
            The separator to use when formatting the string. Default is ':'
        plussign: bool
            Whether to prepend a plus sign to the output for positive values.
            Default is False (do not prepend)

    Returns:
        Sexagesimal formatted representation of the input.
    """

    sign, dd, mm, ss = x
    dd = int(dd)   # make sure these are ints
    mm = int(mm)
    if sign < 0:
        sig = "-"
    else:
        if plussign:
            sig = "+"
        else:
            sig = ""

    # Zero padded formatting for fields
    f1 = "{:02d}"
    f2 = sep + "{:02d}"
    f3 = sep + "{:04.1f}"

    if short == 'ra':
        fmt = sig + f1 + f2 + sep + "{:04.1f}"
        return fmt.format(abs(dd), mm, float(ss / 6))

    if short:
        fmt = sig + f1 + f2
        return fmt.format(abs(dd), mm)

    fmt = sig + f1 + f2 + f3
    return fmt.format(abs(dd), mm, ss)


def sky_area(ra, dec, units='degrees'):
    """
    Based on: 'Computing the Area of a Spherical Polygon" by Robert D. Miller
    in "Graphics Gems IV', Academic Press, 1994

    http://users.soe.ucsc.edu/~pang/160/f98/Gems/GemsIV/sph_poly.c

    Translated from J.D. Smith's IDL routine spherical_poly_area
    Returns the area (solid angle) projected in the sky by a polygon with
    vertices (ra, dec) in degrees

    Doesn't work well on wide range of RA's

    Parameters:
        ra: numpy.ndarray
            Array of the RA points of the polygon in degrees, should be the
            same length as dec.
        dec: numpy.ndarray
            Array of the DEC points of the polygon in degrees, should be the
            same length as ra.
        units: str
            The units of the output (not the input). Default is 'degrees',
            meaning square degrees. Any other value is interpreted as
            steradians.

    Returns:
        float of the requested area.
    """

    sterad2degsq = (180 / math.pi) ** 2  # sterad --> deg2
    RADEG = 180.0 / math.pi              # rad    --> degrees
    HalfPi = math.pi / 2.
    lam1 = ra / RADEG
    lam2 = numpy.roll(lam1, 1)
    beta1 = dec / RADEG
    beta2 = numpy.roll(beta1, 1)
    cbeta1 = numpy.cos(beta1)
    cbeta2 = numpy.roll(cbeta1, 1)

    HavA = numpy.sin((beta2 - beta1) / 2.0) ** 2 + cbeta1 * cbeta2 * numpy.sin((lam2 - lam1) / 2.0) ** 2

    A = 2.0 * numpy.arcsin(numpy.sqrt(HavA))
    B = HalfPi - beta2
    C = HalfPi - beta1
    S = 0.5 * (A + B + C)
    T = numpy.tan(S / 2.0) * numpy.tan((S - A) / 2.0) * numpy.tan((S - B) / 2.0) * numpy.tan((S - C) / 2.0)

    lam = (lam2 - lam1) + 2.0 * math.pi * (numpy.where(lam1 >= lam2, 1, 0))

    Excess = numpy.abs(4.0 * numpy.arctan(numpy.sqrt(numpy.abs(T)))) * \
                      (1.0 - 2.0 * (numpy.where(lam > math.pi, 1.0, 0.0)))

    area = abs((Excess * numpy.where(lam2 != lam1, 1, 0)).sum())
    if units == 'degrees':
        area = area * sterad2degsq

    return area


def get_pixelscale(header, units='arcsec'):
    """
    Returns the pixel-scale from the CDX_X matrix in an WCS-compiant header

    Parameters:
    header: fits style header
        The header to work with
    units: str
        The units to return the value in. Possibilities are: 'degree',
        'arcmin', 'arcsec' (the default).

    Returns:
        float of the requested value
    """
    if units == 'arcsec':
        scale = 3600
    elif units == 'arcmin':
        scale = 60
    elif units == 'degree':
        scale = 1
    else:
        raise ValueError("must define units as arcses/arcmin/degree only")

    try:
        CD1_1 = header['CD1_1']
        CD1_2 = header['CD1_2']
        CD2_1 = header['CD2_1']
        CD2_2 = header['CD2_2']
        pixscale = scale * math.sqrt(abs(CD1_1 * CD2_2 - CD1_2 * CD2_1))
    except KeyError:
        CDELT1 = header['CDELT1']
        CDELT2 = header['CDELT2']
        pixscale = scale * (abs(CDELT1) + abs(CDELT2))/2.

    return pixscale
