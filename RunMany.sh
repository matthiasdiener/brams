#!/bin/bash
export FIRST=$1
export LAST=$2
for (( proc=$FIRST ; proc <= $LAST ; proc++ ))
do
rm -rf ANL/* HIS/* POS/* Out* Dump*
echo "proc=" $proc
if  ./XBlock.sh $proc 
then
echo "diff -r --brief ANL ANL_1Proc_Curto"
if diff -r --brief ANL ANL_1Proc_Curto
#echo "diff -r --brief ANL ANL_1Proc_Curto_NoOPT"
#if diff -r --brief ANL ANL_1Proc_Curto_NoOPT
then
grep "Time integration" Out_${proc}M_* 
else
exit -1
fi
else
exit -1
fi
done
exit 0

