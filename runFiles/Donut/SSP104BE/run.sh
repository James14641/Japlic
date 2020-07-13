#!/bin/bash -e

rm -rf constant/polyMesh [0-9]* constant/T_analytic_init constant/Co* \
       constant/phi constant/U* log

blockMesh
mkdir 0
cp constant/Tsave constant/T_analytic_init 
setAnalyticTracerField
setVelocityField
mv 0/phi 0/U 0/Uf constant
cp 0/T_analytic 0/T

JimpExpEulerFoam -SSP104BE >& log & sleep 0.01; tail -f log

