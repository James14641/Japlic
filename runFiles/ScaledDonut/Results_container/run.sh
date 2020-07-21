#!/bin/bash -e

cd /home/james/OpenFOAM/james-7/Japplications/runFiles/ScaledDonut/FEBE/Tests_of_worth/FEBE_Upwind_Relative_T_Conv
./run.sh

cd /home/james/OpenFOAM/james-7/Japplications/runFiles/ScaledDonut/SSP104BE/tests_of_worth/SSP104BE_Relative_Time_Convergence
./run.sh
 ## these ones takes a while to run so it is read in rather than run .sh again.
cd /home/james/OpenFOAM/james-7/Japplications/runFiles/ScaledDonut/SSP104BE/tests_of_worth/sandt/SSP104_convergenceSAT
python3 readDatFile.py

cd /home/james/OpenFOAM/james-7/Japplications/runFiles/ScaledDonut/SSP104BE/tests_of_worth/sandt/SSP104BE_convergenceSAT
python3 readDatFile.py

cd /home/james/OpenFOAM/james-7/Japplications/runFiles/ScaledDonut/SSP104BE/tests_of_worth/sandt/SSP104BEDC_3it_convergenceSAT
python3 readDatFile.py


