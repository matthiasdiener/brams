#! /bin/bash +x

if [ -z $1 ]; then
   echo
   echo "Falta mecanismo quimico".; echo
   echo "Use: comp.bash <Mecanismo> [clean]"
   echo "Onde, <Mecanismo>:"
   echo "    CO2_JULES"
   echo "    CB07_TUV"
   echo "    RACM_TUV"
   echo "    RELACS_TUV"
   echo "  ou variacoes "
   echo ""
   echo "O segundo argumento 'clean' eh opcional";echo
   echo "Exemplo: comp.bash RELACS_TUV "
   exit
fi
if [ "$1" = "clean" ]
then
   echo
   echo "Mecanismo quimico inadequado: " $1; echo
   echo "Use: comp.bash <Mecanismo> [clean]"
   echo "Onde, <Mecanismo>:"
   echo "    CO2_JULES"
   echo "    CB07_TUV"
   echo "    RACM_TUV"
   echo "    RELACS_TUV"
   echo "  ou variacoes "
   echo ""
   echo "O segundo argumento 'clean' eh opcional";echo
   echo "Exemplo: comp.bash RELACS_TUV clean"
   exit
fi

DIR=$(pwd)
rm -f ../jules3.0-brams5.0.x

cd ../../src/jules/LIB/
rm -f jules_fast.exe
make $2
if [ ! -s jules_fast.exe ] && [ "$2" != "clean" ]; then
   echo;echo "Erro na compilacao do JULES"
   exit
fi

cd $DIR

. /opt/modules/default/etc/modules.sh
#module load netcdf/3.6.2
module load PrgEnv-pgi
module load xt-mpich2
module load hdf5-parallel
module list

if [ "$2" = "clean" ]; then
   rm -f *.o *.oo *.a *.mod ../*.a ../jules3.0-ccatt-brams-5.0-opt.pgi-cray-$1
else   
   make -f Make_model OPT=opt.pgi-cray CHEM=$1
fi

