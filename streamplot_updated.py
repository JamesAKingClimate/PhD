# -*- coding: utf-8 -*-
# Script to make longtitude by height streamplots  from netCDF data. James A. King, University of Oxford, 2018.
import numpy as np
import matplotlib.pyplot as plt
import time
import glob2
import os
import scipy
from netCDF4 import Dataset

start=time.time()
fig=plt.figure(figsize=(16,9))
# sort input files by time they were made 
ifiles= glob2.glob('merra2_*_u_omega_1980-2008.nc')
ifiles.sort(key=os.path.getmtime)
#sort input files in alphabetical order
#ifiles.sort()
# start at the 0th entry in the list, ie the 1st one  because Python is daft
count = 0
print  len(ifiles)

for model in ifiles: 
    # extract model names for subplot titles
    mname = model.split('_')[1]
    print mname
    
    f=Dataset(model, mode='r')
    lons=f.variables['lon'][:]
    lats=f.variables['lat'][:]
    levels=f.variables['plev'][:]
    levels=levels/100
    omega=f.variables['wap'][0,:]
    uwind=f.variables['ua'][0,:]
    f.close()
    # define subplot no. of rows and columns, and run through files by adding 1 to count each time
    plots=plt.subplot(2,3,count+1)
    print count
    # remove extra time dimension from data for xsections
    omega = np.squeeze(omega)
    omega_stream = omega*1000
    uwind = np.squeeze(uwind)
    # create lon/level or lon/lat grid

    [lonsall, levelsall] = np.meshgrid(lons, levels)
    
    print omega.shape
    print uwind.shape 
    print lons.shape
    print lats.shape
    # define contour levels and plot filled contours
    #statistics for dataset - useful for analysis and for defining contour dimensions
    #print np.min(omega),np.max(omega),np.mean(omega),np.std(omega)     
    #print np.min(uwind), np.max(uwind), np.mean(uwind), np.std(uwind)
    #print np.min(omega_stream),np.max(omega_stream),np.mean(omega_stream),np.std(omega_stream)
    field_range=np.arange(-0.08,0.08,0.002) 
    # landcontour=plt.contour(lonsall,levelsall, omega, levels=[1e+20], colors='k',linestyles=('-'),linewidths=(2)) #zero contour
    mymap=plt.contourf(lonsall, levelsall, omega,field_range,cmap=plt.cm.seismic,extend='both')
    plt.streamplot(lonsall,levelsall,uwind, omega_stream,color='k', linewidth=0.75, density=2, arrowsize=0.75, arrowstyle='fancy') 
    #planline=plt.contour(lonsall,levelsall, levelsall, levels=[400], colors='k',linestyles=('-'),linewidths=(3))
    count +=1
    #axis labels and title
    plt.xlim(15,130)
    plt.ylim(200,1000)
    plt.title(mname, loc= 'left',fontsize=10)
    #if mname in ['ACCESS1-3', 'CMCC-CM', 'HadGEM2', 'MIROC5']:
    plt.ylabel('Pressure level (hPa)', rotation='vertical', labelpad=0.5, fontsize=10)
    plt.yticks(np.arange(200,1000,100), fontsize=7)
    #else:
    	#plt.yticks(np.arange(200,1000,100), [])

    #if mname in ['MIROC5', 'MPI-ESM-MR', 'MPI-ESM-LR','NorESM1-M','IPSL-CM5B-LR']:
    plt.xlabel('Longitude ($^\circ$)',labelpad=0.3, fontsize=10)
    plt.xticks(np.arange(15,130,10), fontsize=7)
    #else:
	#plt.xticks(np.arange(15,130,10), []) 
    plt.tight_layout()
    fig.subplots_adjust(left=0.05, right=0.98, bottom=0.38, top=0.94, wspace=0.12, hspace=0.16) 
    #invert y axis for xsections
    plt.gca().invert_yaxis() 
#plt.tight_layout()     
#colour bar
cbar_axis = fig.add_axes([0.21, 0.3, 0.6, 0.03])
cbar=plt.colorbar(mymap, cbar_axis, orientation='horizontal')
cbar.set_label('Omega (Pa/s)', rotation=0, fontsize=10)
cbar.ax.tick_params(labelsize=10,length=0)       

#fig.suptitle('CMIP5 dry models 1975-2005', y=0.98, fontsize=12)  
#save plot to png file
plt.savefig('King_figure_4.png',orientation='portrait',format='png', dpi=300)
plt.close()
plt.clf()
print 'Execution time =',time.time()-start, 'seconds'    
