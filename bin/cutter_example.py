#!/usr/bin/env python

import spt3g_cutter.cutterlib as cutterlib
import logging
import time
import sys

if __name__ == "__main__":

    # Example of inputs:
    # ra,dec can be list or scalars
    filename = sys.argv[1]
    outdir = sys.argv[2]
    ra = [358.3406816, 355.8253677]
    dec = [-58.9660379, -57.1623017]

    xsize = [10]*len(ra)
    ysize = [10]*len(ra)

    # Create logger
    cutterlib.create_logger()
    logger = logging.getLogger(__name__)

    t0 = time.time()
    cutterlib.fitscutter(filename, ra, dec, xsize=xsize, ysize=ysize,
                         units='arcmin', prefix='spt3g', outdir=outdir)
    logger.info(f"Done: {cutterlib.elapsed_time(t0)}")
