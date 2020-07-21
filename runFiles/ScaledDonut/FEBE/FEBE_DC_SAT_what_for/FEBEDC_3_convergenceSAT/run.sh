#!/bin/bash -e

cd ./32
rm -rf constant/polyMesh [0-9]*
./run.sh
cd ../64
./run.sh
cd ../128
./run.sh
cd ../256
./run.sh
cd ../512
./run.sh
cd ..
python3 readDatFile.py


