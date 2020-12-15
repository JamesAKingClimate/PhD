# -*- coding: utf-8 -*-
# Script to make trend maps from netCDF data
import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.basemap as bm
import time
import glob2
import os
from netCDF4 import Dataset
from matplotlib.ticker import MaxNLocator
from numpy import *
from scipy import stats

def l_trend(var,lon,lat,time,sig=False):
    nlon=len(lon)
    nlat=len(lat)
    nt=len(time)
    vart=zeros(nlat*nlon)
    varp=zeros(nlat*nlon)
    
    if len(var.shape)== 3:        
        var=reshape(var,(nt,nlat*nlon)) 
        print('l_trend: assuming variable as 3D [time,lat,lon]')
        for i in range(nlat*nlon):
            v=var[:,i]  
            vart[i], intercept, r_value, varp[i], std_err=stats.linregress(time,v)
            
        vart=reshape(vart,(nlat,nlon))
        varp=reshape(varp,(nlat,nlon))
        return (vart,varp)
        
    else:
        raise ValueError('Variable shape is not 2D or 3D. Function requires variable in format var[time,lat,lon] or var[time,lon*lat]')
    
    if sig==False:
        return (vart, varp) 
    else:
        for i in range(nlat):
            for j in range (nlon):
                if varp[i,j]>sig:
                  vart[i,j]=nan
        return (vart, varp)

fig=plt.figure(figsize=(15,6))

#for multiple files:
# sort input files by time they were made 

ifiles = glob2.glob('jra55_rainy_seasons/jra55_*_*_mslp_1981-2016.nc')
ifiles.sort()
count = 0
print ifiles

for model in ifiles:
	#extract model name for subplot titles
	mname = model.split('_')[4]
	print mname
	f = Dataset(model, mode = 'r')
	lons = f.variables['g4_lon_2'][:]
	lats = f.variables['g4_lat_1'][:]
	field = f.variables['PRMSL_GDS4_MSL_S113'][:,:,:]
	field=field/100
	time = f.variables['initial_time0_hours'][:]
	f.close()
	plots=plt.subplot(2,3,count+1)
	nlon = len(lons)
	nlat = len(lats)
	ntime = len(time)
        year = linspace(1981, 2016, 36)
        #if mname in  ['EC-EARTH', 'HadGEM2-AO']: 
		#year = linspace(2020, 2099, 80)
	#else:
		#year = linspace(2020, 2100, 81)
        print np.shape(field)

	field_trend, field_p = l_trend(field, lons, lats, year)

	[lonall,latall]=np.meshgrid(lons,lats)
    
	# define map - arange creates range of numbers from -180 to 180 in steps of 30 - used to draw lat/long grid on a map (vary according to scale)
	m=bm.Basemap(projection='cyl', resolution ='l',llcrnrlat=-25,llcrnrlon=15,urcrnrlat=20,urcrnrlon=130)
	m.drawcoastlines(linewidth=1.75)
	m.drawcountries(linewidth=1)
	m.drawmapboundary()
	m.drawmeridians(np.arange(-180,180,20),labels=[0,0,0,1],fontsize=7,linewidth=0.5)
	m.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0],fontsize=7,linewidth=0.5)
 
	# define contour levels and plot filled contours
	#statistics for dataset - useful for analysis and for defining contour dimensions
	print np.min(field_trend),np.max(field_trend),np.mean(field_trend),np.std(field_trend)  
	print np.shape(field_trend)
	field_range=np.arange(-0.1,0.1,0.005)
	mymap=plt.contourf(lonall, latall, field_trend,field_range,cmap=plt.cm.RdBu_r,extend='both')
	#add stippling for significance test based on calculated trend p value
	intervals = [0,0.1,1]
	sig_plot = plt.contourf(lonall,latall,field_p,intervals,extend='neither', colors='none', edgecolor='black', hatches=['.....', None])
	count +=1
	plt.title(mname, loc= 'left', fontsize=10) 	

plt.tight_layout()    
fig.subplots_adjust(left=0.05, right=0.98, bottom=0.175, top=0.92, wspace=0.15, hspace=0.2)
#colour bar
cbar_axis = fig.add_axes([0.21, 0.075, 0.6, 0.03])
cbar=plt.colorbar(mymap, cbar_axis, orientation = 'horizontal')
cbar.set_label('SLP trend (hPa month$^-$$^1$ year$^-$$^1$)', fontsize=10)
cbar.ax.tick_params(labelsize=10,length=0)
# axis labels and title
#fig.suptitle('Dec',fontsize=12,y=1)

# save plot to png file
plt.savefig('JRA55_EA_SLP_trends_scipy_90_confidence.png',orientation='portrait',format='png', bbox_inches='tight', dpi=300)
plt.close()
plt.clf()


