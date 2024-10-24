import netCDF4 as nc
import matplotlib.pyplot as plt
import numpy as np
import glob
import os
from scipy.ndimage import gaussian_filter

from skimage.feature import canny

# setup run info
run = 'PBL_control_DT3'
#run = 'PBL_control'
run = 'PBL_20210402_raw_noWind'
#run = 'PBL_20210402_aiqc_noWind'
#run = 'PBL_20210402_raw'
run = 'PBL_20210402_raw_chem_sfc'
run = 'PBL_20210402_GE_hetero'
run = 'PBL_S1_GE_hetero'
path = f'/data/huanghuai/VVM/DATA/{run}'
#path = '/data/wudingrong/VVM/DATA/pbl_control'
topo_file = f'{path}/TOPO.nc'

try: os.makedirs(f'{run}')
except: pass

therm = np.sort(glob.glob(f'{path}/archive/*.Thermodynamic*.nc'))
dyn = np.sort(glob.glob(f'{path}/archive/*.Dynamic*.nc'))
sfc = np.sort(glob.glob(f'{path}/archive/*.Surface*.nc'))
rad = np.sort(glob.glob(f'{path}/archive/*.Radiation*.nc'))
land = np.sort(glob.glob(f'{path}/archive/*.Landsurface*.nc'))
tracer = np.sort(glob.glob(f'{path}/archive/*.Tracer*.nc'))
chem = np.sort(glob.glob(f'{path}/archive/*.Chem*.nc'))


# read coordinate info
bar = np.fromfile(f'{path}/bar.dat', dtype=np.float32).reshape((-1,50))
P = bar[0]
RHO = bar[2]
datT = nc.Dataset(therm[1])
x = datT['xc'][:]/1000
y = datT['yc'][:]/1000
z = datT['zc'][:]/1000
dt = datT['time'][:] * 60

nt = len(therm)
nz = len(z)
ny = len(y)
nx = len(x)


wth_G = np.zeros((nt, nz)) + np.nan
wth_E = np.zeros((nt, nz))+ np.nan
wth_GE = np.zeros((nt, nz))+ np.nan
times = np.zeros(nt)


for t in range(nt):
    print(t)
    datT = nc.Dataset(therm[t])
    datD = nc.Dataset(dyn[t])
    datS = nc.Dataset(sfc[t])
    
    time = datT['time'][0] + 5*60 # init time
    day = int(time//1440)
    hr = int((time//60)%24)
    mn = int(time%60)
    
    W_ = datD['w'][0]
    W = np.zeros(W_.shape)
    W[1:] = 0.5*(W_[1:]+W_[:-1])
    W[0] = 0
    W = W-np.mean(W, axis=(1,2)).reshape((-1,1,1))
    
    TH = datT['th'][0]
    THbar = np.mean(TH, axis=(1,2))
    THp = TH-THbar.reshape((-1,1,1))
    
    wth = W*THp
    wth[0] = datS['wth']/RHO[0]
    
    times[t] = time
    wth_G[t] = np.mean(wth[:,:,:int(nx/2)], axis=(1,2))
    wth_E[t] = np.mean(wth[:,:,int(nx/2):], axis=(1,2))
    wth_GE[t] = np.mean(wth[:,:,:], axis=(1,2))
    
    '''
    f, ax = plt.subplots()
    c = ax.pcolormesh(x, z, np.mean(gaussian_filter(wth, (0,2,2)), axis=1),
                  vmin=-0.1, vmax=0.1, cmap='seismic')
    plt.colorbar(c, ticks=np.arange(-0.1,0.101,0.05), extend='both')
    ax.contour(x, z, np.mean(TH-TH[0], axis=1), levels=[0.5], colors='k', linewidths=3)
    ax.set_xticks(np.arange(0,25.7,6.4)) 
    ax.set_yticks(np.arange(0,2.1,0.5)) 
    ax.set_xlabel('x (km)') 
    ax.set_ylabel('z (km)') 
    ax.set_xlim(0,25.6)
    ax.set_title(f'{run}\nt = day{day} {"%02d"%hr}:{"%02d"%mn}, y-mean', fontsize=10, fontproperties='monospace', loc='right')
    ax.set_title("$w'\\theta '$ (K m/s)", fontsize=16, loc='left')
    #plt.savefig(f'{run}/wth_mean/{"%06d"%t}', dpi=250, bbox_inches='tight')

    plt.show()
    '''

#np.save(f'{run}_wth', [wth_G, wth_E, wth_GE])

#%%

f, axs = plt.subplots(1,3, figsize=(8,3), sharey='all')
c = axs[0].pcolormesh(times/60, z, wth_G.T, vmin=-0.03, vmax=0.03, cmap='RdBu_r')
c = axs[1].pcolormesh(times/60, z, wth_GE.T, vmin=-0.03, vmax=0.03, cmap='RdBu_r')
c = axs[2].pcolormesh(times/60, z, wth_E.T, vmin=-0.03, vmax=0.03, cmap='RdBu_r')

plt.colorbar(c, ax=axs, orientation='horizontal', aspect=50, extend='both', pad=0.21, label="$w'\\theta '$ (K m/s)")
[ax.set_xlim(5,22) for ax in axs]
[ax.grid() for ax in axs]
axs[1].set_xlabel('Time (hr)')
axs[0].set_ylabel('Height (km)')
axs[0].set_title('G')
axs[1].set_title('GE')
axs[2].set_title('E')
f.suptitle(run, y=1.05, fontsize=14)
#plt.savefig(f'{run}_wth_hov', dpi=300, bbox_inches='tight')




#%%

f, axs = plt.subplots(1,3, figsize=(8,3), sharey='all')

axs[0].plot(wth_G[times/60==10][0], z, c='y')
axs[0].plot(wth_E[times/60==10][0], z, c='forestgreen')
axs[0].plot(wth_GE[times/60==10][0], z, c='k', ls='--')

axs[0].set_title('09:00')

axs[1].plot(wth_G[times/60==12][0], z, c='y')
axs[1].plot(wth_E[times/60==12][0], z, c='forestgreen')
axs[1].plot(wth_GE[times/60==12][0], z, c='k', ls='--')
axs[1].set_title('12:00')

axs[2].plot(wth_G[times/60==15][0], z, c='y')
axs[2].plot(wth_E[times/60==15][0], z, c='forestgreen')
axs[2].plot(wth_GE[times/60==15][0], z, c='k', ls='--')
axs[2].set_title('15:00')


[ax.axis([-0.05,0.12,0,2]) for ax in axs]
[ax.grid() for ax in axs]
[ax.set_xlabel("$w'\\theta '$ (K m/s)") for ax in axs]
[ax.set_yticks(np.arange(0,2.1,0.5))  for ax in axs]
axs[0].set_ylabel('Height (km)')
axs[1].legend(['G', 'E', 'GE'], ncol=3, loc=(-0.1,-0.35))
f.suptitle(run, y=1.05, fontsize=14)

#plt.savefig(f'{run}_wth_profile', dpi=250, bbox_inches='tight')



