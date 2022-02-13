from distutils.core import setup
import glob

# Get the scripts/bin files
bin_files = glob.glob("bin/*")

# Build the structure for etc folder
etc_dirs = ['etc']
data_files = [("", ["setpath.sh"])]
for edir in etc_dirs:
    data_files.append((edir, glob.glob("{}/*".format(edir))))


# The main call
setup(name='spt3g_cutter',
      version='0.3.4',
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
