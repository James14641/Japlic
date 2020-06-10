#!/bin/bash -e

rm -rf constant/polyMesh [0-9]*

blockMesh
mkdir 0
cp constant/Tsave constant/T_analytic_init 
setAnalyticTracerField
setVelocityField

cp 0/T_analytic 0/T
cp constant/Tsave 0/KK
cp constant/Tsave 0/KK2

#gmtFoam -time 0 UT
#gv 0/UT.pdf &

## JimpExpEulerFoam -FEBE #>& log & sleep 0.01; tail -f log
JimpExpEulerFoam -SSP2BE #>& log & sleep 0.01; tail -f log

time=600
#ln -s ../0/CoImp $time
#ln -s ../0/U $time
#gmtFoam -time $time T
#gv $time/T.pdf &




