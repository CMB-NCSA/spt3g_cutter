from setuptools import setup
import glob

bin_files = glob.glob("bin/*")
data_files = [("", ["setpath.sh"])]

setup(
    name='spt3g_ingest',
    version='0.4.0',
    license="GPL",
    description="Ingesting for SPT3G",
    author="Felipe Menanteau",
    author_email="felipe@illinois.edu",
    packages=['spt3g_ingest'],
    package_dir={'': 'python'},
    scripts=bin_files,
    package_data={'': ['LICENSE']},
    data_files=data_files,
)
