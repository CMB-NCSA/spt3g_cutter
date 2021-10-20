#!/usr/bin/env bash

if ! [ -n "${SPT3G_CUTTER_VERSION+1}" ]; then
  export SPT3G_CUTTER_VERSION=0.2.1
  echo "SPT3G_CUTTER_VERSION is undefined will set to: $SPT3G_CUTTER_VERSION"
else
  echo "SPT3G_CUTTER_VERSION is set to: $SPT3G_CUTTER_VERSION"
fi

conda-build . -c conda-forge -c menanteau

# Find the tar.bz2 file
tarbz=`ls $CONDA_PREFIX/conda-bld/*/spt3g_cutter-$SPT3G_CUTTER_VERSION*.tar.bz2`

echo " "
echo "--------------------"
echo "Now Please upload:"
echo " anaconda upload --force -u menanteau $tarbz"
echo " "
echo "--------------------"
