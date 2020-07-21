import numpy as np
from pathlib import Path
import pandas as pd
import matplotlib.pyplot as plt

pathArrayT = [Path(r'./32/globalSumT.dat'),
             Path(r'./64/globalSumT.dat'),
             Path(r'./128/globalSumT.dat'),
             Path(r'./256/globalSumT.dat'),
             Path(r'./512/globalSumT.dat')]

pathArrayTdiff = [Path(r'./32/globalSumTdiff.dat'),
             Path(r'./64/globalSumTdiff.dat'),
             Path(r'./128/globalSumTdiff.dat'),
             Path(r'./256/globalSumTdiff.dat'),
             Path(r'./512/globalSumTdiff.dat')]

d1 = Path(r'./32/globalSumT.dat')
r1 = open(d1)
print(r1.read())
'''
array1 = np.loadtxt(d1)
print(array1)
'''

def createDataFrameFromPaths(pathArray):
    listOfDataFrames = []
    for path in pathArray:
        df = pd.read_csv(path, sep="\s+", index_col=None, header=0)
        listOfDataFrames.append(df)
    return pd.concat(listOfDataFrames, axis=0, ignore_index=True)

fullDataFrameT = createDataFrameFromPaths(pathArrayT)
fullDataFrameTdiff = createDataFrameFromPaths(pathArrayTdiff)
print(fullDataFrameT)
print(fullDataFrameTdiff)

#fullDataFrame.to_csv(r'TestResults.csv', index = False)#saves table

Tmin = fullDataFrameT['min'].to_numpy()
Tmax = fullDataFrameT['max'].to_numpy()
Tsum = fullDataFrameT['sum'].to_numpy()
Tvar = fullDataFrameT['variance'].to_numpy()
Tl1 = fullDataFrameTdiff['mag'].to_numpy()
Tl2 = fullDataFrameTdiff['RMS'].to_numpy()
Tlinf = fullDataFrameTdiff['inf'].to_numpy()
nx = 2**fullDataFrameT.index.values*32
nt = 2**fullDataFrameT.index.values*100


Table = fullDataFrameTdiff
Table.drop(columns = ['#time','variance','sum','variance','min','max'],inplace = True)
Table.rename(columns = {"mag":"l1", "RMS":"l2","inf":"linf"},inplace = True)
Table.insert(0,"nx",nx),Table.insert(0,"nt",nt);Table.insert(4,"Tmin",Tmin);
Table.insert(5,"Tmax",Tmax);
Table.insert(6,"Tsum",Tsum);
print(Table)
print(Table.to_latex(index=False))

#### Plotting the errors and their convergence 

#### I believe that the issue here is the test being carried out
'''
plt.subplot(121); 
plt.plot(nx,Tl1,'-o',label = "l1");
plt.plot(nx,Tl2,'-o',label = "l2");
plt.plot(nx,Tlinf,'-o',label = "linfty");
plt.title(r'$l^1, l^2,l^{\infty}$ error plot');plt.xlabel(r'$nx$');
plt.ylabel('error magnitude');plt.legend();
plt.subplot(122);
plt.loglog(nx,Tl1,'-o',label = "l1");
plt.loglog(nx,Tl2,'-o',label = "l2");
plt.loglog(nx,Tlinf,'-o',label = "linfty");
plt.loglog(nx,nx**(-1/2),'k',label = r'$O(nx^{-0.5})$');
plt.loglog(nx,(nx*1.0)**(-1),'k',label = r'$O(nx^{-1})$');
plt.title('loglogplot');plt.xlabel(r'$\log(nx)$');
plt.ylabel(r'$\log(err)$');plt.legend();
plt.suptitle("convergence")

plt.savefig('FEBEDC_convergence_SAT.png')

'''

plt.plot();
plt.loglog(nx,Tl1,'-o',label = "l1");
plt.loglog(nx,Tl2,'-o',label = "l2");
plt.loglog(nx,Tlinf,'-o',label = "linfty");
plt.loglog(nx,nx**(-1/2),'k',label = r'$O(nx^{-0.5})$');
plt.loglog(nx,(nx*1.0)**(-1),'k',label = r'$O(nx^{-1})$');
plt.title(r'$l^1, l^2,l^{\infty}$ log error plot');plt.xlabel(r'$\log(nx)$');
plt.ylabel(r'$\log(err)$');plt.legend();
plt.suptitle("FEBEDC_3_SAT")


plt.savefig('FEBEDC_3_convergence_SAT.png')



