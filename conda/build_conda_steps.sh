export SPT3G_CUTTER_VERSION=0.1.0
conda-build . -c conda-forge -c menanteau
echo " "
echo "--------------------"
echo "Now Please upload:"
echo " anaconda upload --force -u lsstts $CONDA_PREFIX/conda-bld/linux-64/spt3g_cutter-$SPT3G_CUTTER_VERSION-py38_0.tar.bz2"
echo " "
echo "--------------------"
