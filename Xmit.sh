#!/bin/bash
# 
# script to run BRAMS on pure MPI at tupa
#
# usage: Xmit.sh TotMpi 
#    where:
#    TotMpi:  total number of MPI ranks
#
#
#
# script arguments
#
export TotMpi=$1
export OmpMpi=1
export RamsIn=RAMSIN_DINAMICA_PURA
#
# pbs arguments
#
if ((${TotMpi} < 24))
then {
export Width=24
export MpiNod=${TotMpi}
export Nppn=24
}
else {
export Width=${TotMpi}
export MpiNod=24
export Nppn=${MpiNod}
}
fi
#
# job name 
#
RunName=R_${TotMpi}_${MpiNod}
OutName=Out_${TotMpi}M_${MpiNod}MpNode
DumpOut=Dump${OutName}
#
# directories
# executable full path
#
export DirBase=`pwd`
cd ${DirBase}
export executable=${DirBase}/brams.exe

#
# starts producing queue script file qsub.sh
#
cat <<EOF0> qsub.sh
#!/bin/bash
#PBS -A CPTEC
#PBS -l walltime=01:00:00
#PBS -l mppwidth=${Width}
#PBS -l mppdepth=${OmpMpi}
#PBS -l mppnppn=${Nppn}
#PBS -N ${RunName}
#PBS -j oe
#PBS -o ${DirBase}/${OutName}
#PBS -q pesq

cd ${DirBase}

export OMP_NUM_THREADS=${OmpMpi}

ulimit -m unlimited
ulimit -s unlimited
ulimit -v unlimited
#ulimit -c unlimited

time aprun -b -n ${TotMpi} -d ${OmpMpi} -N ${MpiNod} ${executable} -f ${RamsIn}
EOF0
#
# finishes producing file qsub.sh and moves to executable directory
#
chmod +x qsub.sh
#
# qsub with variable # PEs per node
#
qsub qsub.sh 

exit
