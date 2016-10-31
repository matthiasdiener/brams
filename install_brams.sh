#!/bin/bash

sudo apt-get update

while ! sudo DEBIAN_FRONTEND=noninteractive apt-get -y --force-yes install make openmpi-bin libopenmpi-dev gfortran g++ unzip libc-dev zlib1g-dev
 do
	echo "Failed. Trying again..." >&2
	sleep 1
done

LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:$LIBRARY_PATH
export LIBRARY_PATH

############## Install mpich2 #################

#cd ~
#wget -c https://www.dropbox.com/s/wrk3x9ckkvvrbwf/mpich2-1.5.tar.gz
#tar -xzf mpich2-1.5.tar.gz
#cd mpich2-1.5/
#./configure --disable-fast
#make
#sudo make install

############### Install hdf5 ###################

cd ~
wget -c https://www.dropbox.com/s/uilwdfwecfyusni/hdf5-1.8.13.tar.gz
tar -xzf hdf5-1.8.13.tar.gz
cd ~/hdf5-1.8.13
./configure --enable-fortran --with-zlib=/usr/include/ --enable-parallel --enable-filters=all --enable-build-all --prefix=/opt/hdf5-gfortran
make
sudo make install
sudo ln -s /opt/hdf5-gfortran /opt/hdf5

############## Install BRAMS 5 #################

export LD_LIBRARY_PATH=/usr/local/lib/:$LIBRARY_PATH
export PATH=$PATH:/opt/hdf5/bin/

cd ~
wget -c https://www.dropbox.com/s/jafjvgh1odkev00/ediazc-brams-5.0-4d8df4e59dfc.zip
unzip ediazc-brams-5.0-4d8df4e59dfc.zip
mv ediazc-brams-5.0-4d8df4e59dfc BRAMS-5.0

cd ~/BRAMS-5.0/src/jules/LIB/
make -f Makefile_gfortran_x86

cd ~/BRAMS-5.0/build/bin
make -f Make_model OPT=opt.gfortran-x86 CHEM=CB07_TUV

mkdir ~/BRAMS5
cd ~/BRAMS5
cp ~/BRAMS-5.0/build/jules3.0-ccatt-brams-5.0-opt.gfortran-x86-CB07_TUV .

sudo mkdir /brams
sudo chown `whoami`:`whoami` /brams
ln -s /brams/RAMSIN RAMSIN
ln -s /brams/tables tables

cd ~
wget -c https://www.dropbox.com/s/t4ymsdlqaxk9946/light1gr.tar.gz
tar -xzf light1gr.tar.gz
cd light1gr/
ln -s /home/`whoami`/BRAMS5/jules3.0-ccatt-brams-5.0-opt.gfortran-x86-CB07_TUV jules3.0-ccatt-brams-5.0-opt.gfortran-x86-CB07_TUV

################################################

cd ~

sudo rm -f mpich2-1.5.tar.gz
sudo rm -f -f -R mpich2-1.5

sudo rm -f hdf5-1.8.13.tar.gz
sudo rm -f -R hdf5-1.8.13

sudo rm -f ediazc-brams-5.0-4d8df4e59dfc.zip
sudo rm -f -R BRAMS-5.0

sudo rm -f light1gr.tar.gz
