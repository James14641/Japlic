#!/bin/bash -e
cd ./3200
rm -rf constant/polyMesh [0-9]*
./run.sh
cd ../100
rm -rf constant/polyMesh [0-9]*
./run.sh
cd ../200
rm -rf constant/polyMesh [0-9]*
./run.sh
cd ../400
rm -rf constant/polyMesh [0-9]*
./run.sh
cd ../800
rm -rf constant/polyMesh [0-9]*
./run.sh
cd ../1600
rm -rf constant/polyMesh [0-9]*
./run.sh
cd ..
python3 readDatFile.py
