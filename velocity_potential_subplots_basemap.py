 import  mpl_toolkits.basemap as bm
from mpl_toolkits.basemap import addcyclic
import matplotlib as mpl
import matplotlib.pyplot as plt
from netCDF4 import Dataset, date2index
import numpy as np
from windspharm.standard import VectorWind
from windspharm.tools import prep_data, recover_data, order_latdim
from windspharm.examples import example_data_path
import time
import glob2
import os
import sys as sys
#sys.path.append('/lustre/soge1/users/sedm4616/dphil_work/pyclimate')
#from pyclimate import diffoperators as diff
start = time.time()

mpl.rcParams['mathtext.default'] = 'regular'

fig=plt.figure(figsize=(15,6))

#for multiple files:
# sort input files by time they were made 
ifiles= glob2.glob('merra2_u_v_gaussian_200_*_.nc')
ifiles.sort(key=os.path.getmtime) #sort by time files were made
#ifiles.sort() # sort files in alphabetical order
# start at the 0th entry in the list, ie the 1st one  because Python is daft
count = 0
print ifiles

for model in ifiles: 
   
    # extract model names for subplot titles
    mname = model.split('_')[5]
    print mname
    # extract data
    f=Dataset(model, mode='r')
    uwnd = f.variables['ua'][0,0,:,:]
    vwnd = f.variables['va'][0,0,:,:]
    uwnd2 = np.array(uwnd)
    vwnd2 = np.array(vwnd)
    lons = f.variables['lon'][:]
    lats = f.variables['lat'][:]
    f.close()
    plots=plt.subplot(2,3, count+1)
    print count
    # velocity potential calculation - requires some data formatting using Windspharm's tools
    uwnd, uwnd_info = prep_data(uwnd, 'yx')
    vwnd, vwnd_info = prep_data(vwnd, 'yx')
    lats, uwnd, vwnd = order_latdim(lats, uwnd, vwnd)
    w = VectorWind(uwnd, vwnd, gridtype='gaussian')
    vp = w.velocitypotential()
    vp = recover_data(vp, uwnd_info)
    vp = vp*1e-06
    # divergence calculation
    divg = w.divergence()
    divg = recover_data(divg, uwnd_info)
    divg = divg*1e6
    #divergent wind compotnents calculation
    uchi,vchi = w.irrotationalcomponent()
    uchi = np.squeeze(uchi)
    vchi = np.squeeze(vchi)
    # print some data dimensions - useful for colourbar range etc
    print np.shape(lons)
    print np.shape(lats)
    print np.shape(uchi)
    print np.shape(vchi)
    print np.mean(uchi)
    print np.mean(vchi)
    # create lat/lon grids
    #[lonall,latall]=np.meshgrid(lons,lats)
    # plot velocity potential 
    m=bm.Basemap(projection='cyl', resolution ='l',llcrnrlat=-25,llcrnrlon=10,urcrnrlat=20,urcrnrlon=130)
    m.drawcoastlines(linewidth=1.75)
    m.drawcountries(linewidth=1)
    m.drawmapboundary()
    m.drawmeridians(np.arange(-180,180,20),labels=[0,0,0,1],fontsize=7,linewidth=0.5)
    m.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0],fontsize=7,linewidth=0.5)
    field_range=np.arange(-12,12,0.5)
    mymap=plt.contourf(lons, lats, vp,field_range,cmap=plt.cm.seismic,extend='both')
    # add divergence as open contour lines
    #contour = plt.contour(lons, lats, vp, levels=10, colors='k', extend = 'both')
    # add divergent wind as vector arrows
    xx, yy = m(*np.meshgrid(lons,lats))
    u_out, v_out = m.rotate_vector(uchi, vchi, lons, lats)
    Q = m.quiver(xx[::7,::7], yy[::7,::7], u_out[::7,::7], v_out[::7,::7], width=0.00325,scale=None,color ='lime')
    #plt.quiverkey(Q, X = 0.9, Y = -0.075, U = 8, label = '8 m $s^-1$', labelpos='W',fontproperties={'size':'7'})
   
    plt.title(mname, loc= 'left', fontsize=10)     
    count +=1
    
plt.tight_layout()
fig.subplots_adjust(left=0.05, right=0.98, bottom=0.175, top=0.92, wspace=0.15, hspace=0.2)
#colour bar - horizontal
cbar_axis = fig.add_axes([0.21, 0.075, 0.6, 0.03])
cbar=plt.colorbar(mymap, cbar_axis, orientation='horizontal')
cbar.set_label('Velocity Potential ($10^6$m$^2$s$^{-1}$)', rotation=0, fontsize=10)
#cbar.set_label('Div ($s^{-1}$)', rotation = 0, fontsize = 10)
cbar.ax.tick_params(labelsize=10,length=0) 

plt.savefig('MERRA2_200hPa_vp_divergentwindarrows_months.png',orientation='portrait',format='png',dpi=300)
plt.close()
plt.clf()
print 'Execution time =',time.time()-start, 'seconds' 
