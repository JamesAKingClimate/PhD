#Script to compute correlations in space against a timeseries
import numpy as np
from netCDF4 import Dataset
import scipy.stats as sps
import mpl_toolkits.basemap as bm
import matplotlib.pyplot as plt
import time
import glob2
import os


start = time.time()
#get files
precipfiles = glob2.glob('chirps_kenya_timeseries_monthly_years_*.nc')
omegafiles = glob2.glob('merra2_omega_monthly_400_EA_IO_*.nc')
precipfiles.sort(key=os.path.getmtime)
omegafiles.sort(key=os.path.getmtime)
count = 0 
print precipfiles
print omegafiles
fig = plt.figure(figsize=(15,6))

for precip,omega in zip(precipfiles,omegafiles):
	mname = precip.split('_')[5]
    	print mname
	##extract data from files - make sure number of timesteps are the same
	f = Dataset(precip, mode='r') 
	precip_timeseries=f.variables['precip'][:,0,0]
	print precip_timeseries.shape
	timesteps = precip_timeseries.shape[0]
	g = Dataset(omega, mode='r')
	omega_timeseries = g.variables['wap'][:,0,:,:]
	print omega_timeseries.shape
	lons = g.variables['lon'][:]
	lonsize = lons.shape[0]
	print lons.shape
	lats = g.variables['lat'][:]
	latsize=lats.shape[0]
	print lats.shape
	ngridcells=lonsize*latsize
	print ngridcells
	#convert spatial data to an array with dimensions (number of timesteps, number of grid cells)
	two_d_array = np.reshape(omega_timeseries,  (timesteps, ngridcells))
	plots=plt.subplot(2,3,count+1)
	print count
	#print precip_timeseries.shape
	#print omega_timeseries.shape
	#print two_d_array.shape
	#make some empty arrays to put things in
	list_1=[]
	t_test_results = []
	#loop through each grid point in gridded data and calculate Pearson coefficient
	for i in np.arange(0, ngridcells, 1):
        	pearson_coeff = sps.pearsonr((two_d_array[:,i]), precip_timeseries)
        	##print pearson_coeff[0]
        	list_1 = np.append(list_1, pearson_coeff[0])
        	t_test_results = np.append(t_test_results, pearson_coeff[1])

	list_array = np.asarray(list_1)
	#reshape output so it will fit on a map of the domain
	pearson_corr_array = np.reshape(list_array, (latsize, lonsize))
	#print pearson_corr_array.shape
	#print np.min(pearson_corr_array)
	#print np.max(pearson_corr_array)
	#make a map and put data on it


	m = bm.Basemap(projection='cyl', resolution = 'i', llcrnrlat=-25, llcrnrlon=15, urcrnrlat=20, urcrnrlon=130)
	m.drawcoastlines(linewidth=1.75)
	m.drawcountries(linewidth=1)
	m.drawmapboundary()
	m.drawmeridians(np.arange(-180,180,20),labels=[0,0,0,1],fontsize=7,linewidth=0.5)
	m.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0],fontsize=7,linewidth=0.5)

	[lonall,latall] = np.meshgrid(lons,lats)
	print lonall.shape
	print latall.shape
	#x,y=m(*grid)
	#print np.shape(x)
	#print np.shape(y)
	m_data = pearson_corr_array
	print np.shape(m_data)
	contour_levels = np.arange(-1,1,0.1)
	contour_plot = plt.contourf(lonall,latall,m_data,contour_levels, cmap=plt.cm.seismic, extend='both')

	#add t test results 
	t_test_data = np.asarray(t_test_results)
	t_test_data = np.reshape(t_test_data, (latsize, lonsize))
	intervals = [0,0.1,1]
	t_test_plot=plt.contourf(lonall,latall,t_test_data, intervals, extend = 'neither', colors ='none', edgecolor ='black', hatches = ['//', None])
	count +=1
	plt.title(mname, loc= 'left', fontsize=10) 
fig.subplots_adjust(left=0.05, right=0.98, bottom=0.175, top=0.92, wspace=0.15, hspace=0.2)
cbar_axis = fig.add_axes([0.21, 0.075, 0.6, 0.03])
cbar=plt.colorbar(contour_plot, cbar_axis, orientation = 'horizontal')
cbar.set_label('R', fontsize=10)
cbar.ax.tick_params(labelsize=10,length=0)
print 'Execution time =',time.time()-start, 'seconds'  
#plt.title('All Months')
plt.savefig('King_figure_5.png',orientation='portrait',format='png', dpi=300)
plt.close()
plt.clf()

#print 'Execution time =',time.time()-start, 'seconds'   
