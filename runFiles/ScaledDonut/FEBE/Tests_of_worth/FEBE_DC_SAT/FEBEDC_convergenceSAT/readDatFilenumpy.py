import numpy as np
from pathlib import Path
import pandas as pd

pathArray = [Path(r'./FEBEDC_32/globalSumT.dat'),
             Path(r'./FEBEDC_64/globalSumT.dat'),
             Path(r'./FEBEDC_128/globalSumT.dat'),
             Path(r'./FEBEDC_256/globalSumT.dat'),
             Path(r'./FEBEDC_512/globalSumT.dat')]


d1 = Path(r'./FEBEDC_32/globalSumT.dat')
r1 = open(d1)
print(r1.read())
array1 = np.loadtxt(d1)
print(array1)

arraytot = np.zeros(len(pathArray), len(array1))
for i in len(pathArray):
    arraytot[i,:] =  np.loadtxt(pathArray[i])





#fullDataFrame.to_csv(r'TestResults.csv', index = False)