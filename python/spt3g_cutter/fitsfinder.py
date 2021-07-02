import numpy
import logging
import sqlite3

XSIZE_default = 10.0
YSIZE_default = 10.0

logger = logging.getLogger(__name__)


def check_xysize(df, xsize=None, ysize=None):

    """
    Check if  xsize,ysize are set from command-line or read from csv file
    """
    nobj = len(df.RA.values)
    if xsize:
        xsize = numpy.array([xsize]*nobj)
    else:
        try:
            xsize = df.XSIZE.values
        except Exception:
            xsize = numpy.array([XSIZE_default]*nobj)

    if ysize:
        ysize = numpy.array([ysize]*nobj)
    else:
        try:
            ysize = df.YSIZE.values
        except Exception:
            ysize = numpy.array([YSIZE_default]*nobj)

    return xsize, ysize


def get_query(tablename, bands=None, date_start=None, date_end=None, yearly=None):
    """Format query template"""

    query_files_template = """
    select ID, FILEPATH || '/' || FILENAME as FILE, BAND, DATE_BEG from {tablename}
      {where}
       {and_bands}
       {and_dates}
    """

    # Formatting bands
    if bands:
        in_bands = ','.join("\'{}\'".format(s) for s in bands)
        and_bands = f"BAND in ({in_bands})"
    else:
        and_bands = ''

    # Formatting dates
    if isinstance(date_start, str) and isinstance(date_end, str):
        and_dates = f"DATE_BEG between '{date_start}' and '{date_end}'"
        if isinstance(yearly, str):
            and_dates = f"{and_dates} or OBS_ID == '{yearly}'"
        if len(and_bands) > 1:
            and_dates = f"and ({and_dates})"
    else:
        and_dates = ''

    # Adding a where if needed
    if and_dates or and_bands:
        where = 'where'
    else:
        where = ''

    kw = {where: 'where', 'tablename': tablename, 'and_bands': and_bands, 'and_dates': and_dates}
    return query_files_template.format(**kw)


def connect_db(dbname):
    """Establish connection to DB"""
    logger.info(f"Establishing DB connection to: {dbname}")
    # SQLlite DB lives in a file
    con = sqlite3.connect(dbname)
    return con


def query2rec(query, dbhandle):
    """
    Queries DB and returns results as a numpy recarray.
    """
    # Get the cursor from the DB handle
    cur = dbhandle.cursor()
    # Execute
    cur.execute(query)
    tuples = cur.fetchall()

    # Return rec array
    if tuples:
        names = [d[0] for d in cur.description]
        return numpy.rec.array(tuples, names=names)
    else:
        logger.warning("# WARNING DB Query in query2rec() returned no results")
    return False
