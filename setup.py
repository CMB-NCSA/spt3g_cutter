#!/usr/bin/env python3
from distutils.core import setup
import glob

bin_files = glob.glob("bin/*")
data_files = [("", ["setpath.sh"])]

setup(
    name='spt3g_ingest',
    version='0.4.4',
    license="GPL",
    description="Thumbnail cutter for SPT3G",
    author="Felipe Menanteau",
    author_email="felipe@illinois.edu",
    packages=['spt3g_cutter'],
    package_dir={'': 'python'},
    scripts=bin_files,
    package_data={'': ['LICENSE']},
    data_files=data_files,
)
