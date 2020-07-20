#!/bin/bash -e
cd ./3200
./run.sh
cd ../100
./run.sh
cd ../200
./run.sh
cd ../400
./run.sh
cd ../800
./run.sh
cd ../1600
./run.sh
cd ..
python3 readDatFile.py
