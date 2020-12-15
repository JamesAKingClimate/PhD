# -*- coding: utf-8 -*-
# Script to make maps from netCDF data. James A. King, University of Oxford, 2017.
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import time
import glob2
import os
from netCDF4 import Dataset
from matplotlib.ticker import MaxNLocator
from matplotlib.patches import Polygon

start=time.time()  

# Define boxes
def draw_screen_poly1(lons1, lats1, m):
     x1,y1 = m(lons1, lats1)
     x1y1 = zip(x1,y1)
     poly1 = Polygon(x1y1,edgecolor='k', fill='no', linewidth=2, alpha=1, facecolor='none')
     plt.gca().add_patch(poly1)

def draw_screen_poly2(lons2, lats2, m):
    x2,y2 = m(lons2, lats2)
    x2y2 = zip(x2,y2)
    poly2 = Polygon(x2y2,edgecolor='green', fill='no', linewidth=2, alpha=1, facecolor='none')
    plt.gca().add_patch(poly2)    
    
lons1 = [34, 42, 42, 34] 
lats1 = [-5, -5, 5, 5]

lons2 = [95, 105, 105, 95]
lats2=[-5, -5, 5, 5]

fig=plt.figure(figsize=(15,6))

#for multiple files:
# sort input files by time they were made 
ifiles= glob2.glob('merra2_*_400_omega_EA_IO_1980-2008.nc')
ifiles.sort(key=os.path.getmtime) #sort by time files were made
#ifiles.sort() # sort files in alphabetical order
# start at the 0th entry in the list, ie the 1st one  because Python is daft
count = 0
print ifiles

for model in ifiles: 
    # extract model names for subplot titles
    mname = model.split('_')[1]
    print mname
    f=Dataset(model, mode='r')
    lons=f.variables['lon'][:]
    lats=f.variables['lat'][:]
    field=f.variables['wap'][0,:]
    field=np.squeeze(field)
    f.close()
    # define subplot no. of rows and columns, and run through files by adding 1 to count each time
#for a single file with multiple timesteps:    
# for i in np.arange(12):  
#     ifiles='merra2_relative_humidity_turkana_850.nc'
#     f=Dataset(ifiles, mode='r')
#     field=f.variables['hur'][i,:,:]
#     lons=f.variables['lon'][:]
#     lats=f.variables['lat'][:]
#     f.close()   
    plots=plt.subplot(2,3,count+1)
    #plt.tight_layout(5,5)
    print count
    # create lat/lon grids
    [lonall,latall]=np.meshgrid(lons,lats)
    
    # define map - arange creates range of numbers from -180 to 180 in steps of 30 - used to draw lat/long grid on a map (vary according to scale)
    m=bm.Basemap(projection='cyl', resolution ='l',llcrnrlat=-25,llcrnrlon=15,urcrnrlat=20,urcrnrlon=130)
    m.drawcoastlines(linewidth=1.75)
    m.drawcountries(linewidth=1)
    m.drawmapboundary()
    m.drawmeridians(np.arange(-180,180,20),labels=[0,0,0,1],fontsize=7,linewidth=0.5)
    m.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0],fontsize=7,linewidth=0.5)
    draw_screen_poly1(lons1, lats1, m)
    draw_screen_poly2(lons2, lats2, m)
    # define contour levels and plot filled contours
    #statistics for dataset - useful for analysis and for defining contour dimensions
    print np.min(field),np.max(field),np.mean(field),np.std(field)  
    #print np.shape(lons)   
    print np.shape(lonall)
    #print np.shape(lats)
    print np.shape(latall)
    field=np.squeeze(field)
    print np.shape(field)
    field_range=np.arange(-0.08,0.08,0.002)
    #contourlines=MaxNLocator(nbins=25).tick_values(40,100)
    mymap=plt.contourf(lonall, latall, field,field_range,cmap=plt.cm.seismic,extend='both')
    count +=1
    plt.title(mname, loc= 'left', fontsize=10)    
     
plt.tight_layout() 
#subplot location - 5x5   
#fig.subplots_adjust(left=0.05, right=0.98, bottom=0.38, top=0.94, wspace=0.2, hspace=0.3)
#subplotlocation = 3x2
fig.subplots_adjust(left=0.05, right=0.98, bottom=0.175, top=0.92, wspace=0.15, hspace=0.2)
#colour bar - vertical
#cbar_axis = fig.add_axes([0.99999, 0.36, 0.03, 0.6])
#cbar=plt.colorbar(mymap,cbar_axis, orientation='vertical')
#cbar.set_label('(m/s)', rotation=90, fontsize=1
#colour bar - horizontal
cbar_axis = fig.add_axes([0.21, 0.075, 0.6, 0.03])
cbar=plt.colorbar(mymap, cbar_axis, orientation='horizontal')
cbar.set_label('Omega (Pa/s)', rotation=0, fontsize=10)
cbar.ax.tick_params(labelsize=10,length=0) 

# axis labels and title
#fig.suptitle('',fontsize=12,y=1)

# save plot to png file


plt.savefig('paper_1_fig_3.png',orientation='portrait',format='png')
plt.close()
plt.clf()
print 'Execution time =',time.time()-start, 'seconds' 
