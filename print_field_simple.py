import pandas as pd
from pandas import ExcelWriter
from pandas import ExcelFile
from netCDF4 import Dataset

ifile='merra2_omega_kenya_400_timeseries.nc'
f=Dataset(ifile, mode='r')
field=f.variables['wap'][:,:]
#field=field*86400
mname='omega'
print field
#print str(field).replace('[','').replace(']',' ')
#field=str(field).replace('[',' ').replace(']]',',')
field = ( " ".join(str(i) for i in field))
field= str(field).replace('[',' ').replace(']]]',',')

df=pd.DataFrame({mname:[field]})
writer = ExcelWriter('kenya_omega_monthly_means.xlsx')
df.to_excel(writer,'Sheet1', index=True)
writer.save()
