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


def get_query(tablename, bands=None, filetypes=None, date_start=None, date_end=None, yearly=None):
    """Format query template"""

    query_files_template = """
    select ID, FILEPATH || '/' || FILENAME as FILE, BAND, DATE_BEG from {tablename}
      {where}
       {and_bands}
       {and_dates}
       {and_filetypes}
    """

    # Formatting bands
    if bands:
        in_bands = ','.join("\'{}\'".format(s) for s in bands)
        and_bands = f"BAND in ({in_bands})"
    else:
        and_bands = ''

    # Formatting filetypes
    if filetypes:
        in_filetypes = ','.join("\'{}\'".format(s) for s in filetypes)
        and_filetypes = f"FILETYPE in ({in_filetypes})"
        if bands is not None:
            and_filetypes = f"and ({and_filetypes})"
    else:
        and_filetypes = ''

    # Formatting dates
    if isinstance(date_start, str) and isinstance(date_end, str):
        and_dates = f"DATE_BEG between '{date_start}' and '{date_end}'"
        and_dates_or = ' or '
    else:
        and_dates = ''
        and_dates_or = ''
    if yearly is not None:
        in_yearly = ','.join("\'{}\'".format(s) for s in yearly)
        and_dates = f"{and_dates}{and_dates_or}OBS_ID in ({in_yearly})"
    if bands is not None or filetypes is not None:
        and_dates = f"and ({and_dates})"

    # Adding a where if needed
    if and_dates or and_bands or and_filetypes:
        where = 'where'
    else:
        where = ''

    kw = {where: 'where', 'tablename': tablename,
          'and_bands': and_bands,
          'and_filetypes': and_filetypes,
          'and_dates': and_dates}
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
        logger.error("# DB Query in query2rec() returned no results")
        msg = f"# Error with query:{query}"
        raise RuntimeError(msg)
    return False
