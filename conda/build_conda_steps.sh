export SPT3G_CUTTER_VERSION=0.1.0
conda-build . -c conda-forge -c menanteau

# Find the tar.bz2 file
tarbz=`ls $CONDA_PREFIX/conda-bld/*/spt3g_cutter-*.tar.bz2`

echo " "
echo "--------------------"
echo "Now Please upload:"
echo " anaconda upload --force -u menanteau $tarbz"
echo " "
echo "--------------------"
