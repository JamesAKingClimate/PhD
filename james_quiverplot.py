# Script to plot moisture flux vectors overlayed on a field. Adapted from Amy Creese by James A. King, University of Oxford, 2017. 
import sys as sys
from mpl_toolkits.basemap import Basemap, addcyclic
from netCDF4 import Dataset, date2index
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import time
sys.path.append('Z:/Project_stuff/project_code/pyclimate')
from pyclimate import diffoperators as diff
start =time.time()
# Define paths and files
obswinddir='Z:/Project_stuff/project_code/reanalysis/'
obsraindir='Z:/Project_stuff/project_code/Precip/monthly_composites/'
savedir='Z:/Project_stuff/project_code/reanalysis/'

# Define domain - East Africa/Indian Ocean Basin
lon0=15
lon1=130
lat0=-25
lat1=20


# Read in data 
ncu = Dataset(obswinddir+'erai_mm_u_850_Aug_global.nc','r')
ncv = Dataset(obswinddir+'erai_mm_v_850_Aug_global.nc','r')
ncq=Dataset(obswinddir+'erai_mm_q_850_Aug_global.nc','r')
#ncpre = Dataset(obsraindir+'chirps-v2.0_Aug_EA_mean.nc','r')
u = np.array(ncu.variables['u'][0,0,:-1,:])
v = np.array(ncv.variables['v'][0,0,:-1,:])
q = np.array(ncq.variables['q'][0,0,:-1,:])
#pre = np.array(ncpre.variables['precip'][0,:])

# Generate grid for precip
#lonx1 =np.array(ncpre.variables['longitude'][:])
#latx1 = np.array(ncpre.variables['latitude'][:])
#ncpre.close()

#Generate grid for vectors
lonx2 =np.array(ncu.variables['longitude'][:])
latx2 = np.array(ncu.variables['latitude'][:-1])
ncu.close()

#[lonall1, latall1] = np.meshgrid(lonx1, latx1)

[lonall2, latall2] = np.meshgrid(lonx2, latx2)

#[lonall3, latall3] = np.meshgrid(lonx3, latx3)

#add MFC!
#Calculate qu and qv
qu=q[:]*u[:]
qv=q[:]*v[:]


#Calculate convergence/divergence
Grd = diff.HGRADIENT(latx2, lonx2)
[dqu_dx, dqu_dy] = Grd.hgradient(u)
[dqv_dy, dqv_dx] = Grd.hgradient(v)
divg = dqu_dx + dqv_dy
divg = divg*1e6
print divg.shape


# Basemap features
plt.figure(figsize=(14,8))
m=Basemap(projection = 'cyl', resolution = 'l', llcrnrlat = lat0, urcrnrlat = lat1, llcrnrlon = lon0, urcrnrlon = lon1)
m.drawmapboundary(fill_color='1')
m.drawcoastlines(linewidth=1.5)
m.drawcountries(linewidth=1)
m.drawmeridians(np.arange(-180,180,10),labels=[0,0,0,1],fontsize=7,linewidth=0.5)
m.drawparallels(np.arange(-90,90,10),labels=[1,0,0,0],fontsize=7,linewidth=0.5)
levels=np.arange(-10,10,0.1)
cmap1=plt.cm.coolwarm

# Plot anything you want underneath as a colour map, e.g. q or precip
mymapf = plt.contourf(lonall2, latall2, divg, levels, cmap=cmap1, extend='both')

# This is important for plotting vectors but I can't remember why...
xx, yy = m(*np.meshgrid(lonx2, latx2))
u_out, v_out= m.rotate_vector(qu,qv, lonx2, latx2)

# Plotting vectors
# xx[::5,::5], yy[::5,::5], u_out[::5,::5], v_out[::5,::5] - this part is the spacing, this means 1 in every 5 vectors will be plotted (otherwise it can be too many to see). Note that more vectors need to be plotted for lower-res data e.g. NCEP
Q = m.quiver(xx[::4,::4], yy[::4,::4], u_out[::4,::4], v_out[::4,::4], width=0.0025,scale=0.85,color ='k', pivot ='mid') # width and scale are the size and length of vectors - play around with them to get it right
plt.quiverkey(Q, 0.9, -0.075, 0.05,  '0.05 kg m-1 s-1', labelpos='W',fontproperties={'size':'7'}) #refers to length of each arrow representing units of moisture flux


# Define position of colour bar on plot 
cbar=plt.colorbar(mymapf,orientation='horizontal', shrink=0.7, pad=0.05, extend='both')
cbar.solids.set_edgecolor("face")
cbar.set_label('qDiv', fontsize=8, labelpad=-0.75)
cbar.ax.tick_params(labelsize=7,length=0)

# Plot title
plt.title('ERAI mean Aug 850hPa moisture flux (convergence)',fontsize=12)
plt.tight_layout(pad=2.5,)
plt.savefig(savedir+'erai_mm_qdivg_Aug_850_vectors.png', orientation='landscape')
plt.close()

print 'Execution time =',time.time()-start,'seconds'