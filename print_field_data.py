import numpy as np
import glob2
import os
from netCDF4 import Dataset
import pandas as pd
from pandas import ExcelWriter
import time
start = time.time()
# sort input files in alphabetical order
ifiles= glob2.glob('era5_ua_*_timeseries_1800.nc')
ifiles.sort()
print ifiles
#start at the 0th entry in the list, ie the 1st one  because Python is daft
count = 0
models_field = []
mnames_list = []
for model in ifiles: 
    # extract model names
    mname = model.split('_')[2]
    print mname
    # extract data
    f=Dataset(model, mode='r')
    field=f.variables['u'][:]
    # convert precip units to mm/day
    #field=field*86400
    #convert temp units to Celsius
    #field=field-273.15
    # remove annoying square brackets
    field = ( " ".join(str(i) for i in field))
    field= str(field).replace('[',' ').replace(']]',',')
    # add next model's data onto end of list 
    models_field=np.append(models_field,field)
    mnames_list=np.append(mnames_list,mname)
    
models_field=models_field.tolist()
df=pd.DataFrame([mnames_list[i], models_field[i]] for i in np.arange(0,len(models_field)))
writer = ExcelWriter('era5_1800_uwind.xlsx')
df.to_excel(writer,'Sheet1', index=True)
writer.save()
print 'Spreadsheets!'
print 'Execution time =',time.time()-start, 'seconds' 
