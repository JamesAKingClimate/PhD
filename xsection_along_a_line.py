import matplotlib.pyplot as plt
import time
import numpy as np
import glob2
import os
from netCDF4 import Dataset
from scipy import interpolate
from sys import exit

start=time.time()
fig=plt.figure(figsize=(12,6))
#define a line - lons, lats, and number of steps (arbitrary)
lonline = np.linspace(31,40, 100)
latline = np.linspace(-4,8, 100)

# sort input files by time they were made 
ifiles= glob2.glob('jra55_*_windspeed.nc')
ifiles.sort(key=os.path.getmtime)

# start at the 0th entry in the list, ie the 1st one  because Python is daft
count = 0
print ifiles

for model in ifiles:
    # extract model names for subplot titles
    mname = model.split('_')[1]
    print mname

    f=Dataset(model, mode='r')
    lons=f.variables['g4_lon_3'][:]
    lats=f.variables['g4_lat_2'][:]
    levels=f.variables['lv_HYBL1'][:]
    #levels=levels/100
    field=f.variables['UGRD_GDS4_HYBL_S123'][0,:]
    field = np.squeeze(field)
    f.close()
    plots=plt.subplot(3,4,count+1)
    print count
    # generate an empty array of the same dimensions as your line
    interpolated_stuff=np.zeros((len(lonline),len(levels)))    
    [lonsall, levelsall] = np.meshgrid(lons, levels) 
    print field.shape  

    # define contour levels and plot filled contours
    #statistics for dataset - useful for analysis and for defining contour dimensions
    print np.min(field),np.max(field),np.mean(field),np.std(field)   

    #for each level, interpolate the field data onto the line

    for i,l in enumerate(levels):
        inter=interpolate.interp2d(lons, lats,field[i])
        interpolated_stuff[:,i]=[inter(lonline[j],latline[j]) for j in range(100)]
    count +=1
        
    
    
    field_range=np.arange(0,16,1)
    plt.ylim(1,19)
    plt.ylabel('Hybrid level', fontsize=7)
    plt.xlabel('$^\circ$E', fontsize=7, labelpad=0.15)
    plt.yticks(np.arange(1,19,1),fontsize=7)
    #plt.gca().invert_yaxis() 
    plt.title(mname, loc= 'left',fontsize=7.5)
    mymap=plt.contourf(lonline, levels, interpolated_stuff.T,field_range,cmap=plt.cm.RdYlGn_r,extend='max')
    fig.subplots_adjust(wspace=0.325, hspace=0.295)
    #plt.tight_layout()
cbar_axis = fig.add_axes([0.915, 0.2, 0.02, 0.6])
cbar=plt.colorbar(mymap,cbar_axis, orientation='vertical')
cbar.set_label('Windspeed (m/s)', rotation=90, fontsize=10)
cbar.ax.tick_params(labelsize=10,length=0)

plt.savefig('turkana_xsection_JRA55_windspeed.png', dpi=300)
print 'Execution time =',time.time()-start, 'seconds' 


