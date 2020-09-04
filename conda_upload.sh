# Only need to change these two variables
PKG_NAME=fastmlst
USER=enzoandree

mkdir ~/conda-bld
conda config --set anaconda_upload no
export CONDA_BLD_PATH=~/conda-bld
conda-build .
ls -liar $CONDA_BLD_PATH
anaconda -t $CONDA_UPLOAD_TOKEN upload -u $USER $CONDA_BLD_PATH/noarch/$PKG_NAME*.tar.bz2 --force
