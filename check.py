from netCDF4 import Dataset
import numpy as np

mycdf = Dataset('new_ctm.bpch.nc')
idlcdf = Dataset('bpch2nc_output.20040101.000000.nc')

for idlk, myk in [('BIOGSRCE__ISOP', 'BIOGSRCE_ISOP')]:
    myv = mycdf.variables[myk]
    idlv = idlcdf.variables[idlk]
    print 'my', myv.shape, np.percentile(myv[:], [0, 5, 10, 25, 50, 75, 90, 95, 100])
    print 'idl', idlv.shape, np.percentile(idlv[:], [0, 5, 10, 25, 50, 75, 90, 95, 100])
    all = myv[:] - idlv[:]
    