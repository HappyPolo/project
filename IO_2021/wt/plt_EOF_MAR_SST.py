'''
Author       : ChenKt
Date         : 2021-08-29 16:51:01
LastEditors  : ChenKt
LastEditTime : 2021-08-29 16:58:47
FilePath     : /project/IO_2021/plt_EOF_MAR_SST.py
Aim          : 
Mission      : 
'''
# %%
import cmaps
from mpl_toolkits.axes_grid1 import AxesGrid
from cartopy.mpl.ticker import LongitudeFormatter, LatitudeFormatter
from cartopy.mpl.geoaxes import GeoAxes
import xarray as xr
import numpy as np
import pickle
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
import matplotlib.gridspec as gridspec


def mapart(ax):
  '''
  添加地图元素
  '''
  proj = ccrs.PlateCarree()
  ax.coastlines(color='k', lw=1.2)
  ax.add_feature(cfeature.LAND, facecolor='white')
  # 设置地图范围
  ax.set_extent(extents, crs=proj)
  # 设置经纬度标签
  xticks = np.arange(extents[0], extents[1]+1, 20)
  yticks = np.arange(extents[2], extents[3]+1, 15)
  ax.set_xticks(xticks, crs=proj)
  ax.set_yticks(yticks, crs=proj)
  lon_formatter = LongitudeFormatter(zero_direction_label=True)
  lat_formatter = LatitudeFormatter()
  ax.xaxis.set_major_formatter(lon_formatter)
  ax.yaxis.set_major_formatter(lat_formatter)
  ax.tick_params(axis='both', which='major', labelsize=13, direction='out',
                 length=5, width=0.9, pad=0.2, top=True, right=True)
  ax.minorticks_on()
  ax.tick_params(axis='both', which='minor', direction='out',
                 width=0.6, top=True, right=True)


def barart(ax):
  ax.set_ylim(-3., 3.)
  ax.tick_params(axis='both', which='major', labelsize=11, direction='out',
                 length=2.8, width=0.9, pad=0.2, top=False, right=False)
  ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
  # ax.xaxis.set_tick_params(which='minor', bottom=False)
  ax.tick_params(axis='both', which='major', labelsize=13, direction='out',
                 length=5, width=0.9, pad=0.2, top=True, right=True)
  ax.tick_params(axis='both', which='minor', direction='out',
                 width=0.6, top=True, right=True)
  ax.axhline(y=0, c="black", ls="-", lw=0.5)
  ax.legend(loc='upper left', frameon=False, fontsize=5, ncol=4)
  ax.tick_params(labelsize=14)


def plt_eofs(lp, le, frac, title):
  plt.close
  fig = plt.figure(figsize=(12, 12))
  proj = ccrs.PlateCarree()
  lon = lp.lon
  lat = lp.lat
  x = range(1982, 2022)
  ax1 = fig.add_subplot(3, 2, 1, projection=proj)
  mapart(ax1)
  p1 = ax1.contourf(
      lon, lat, lp[0], transfrom=proj, extend='both', cmap=cmaps.ViBlGrWhYeOrRe,
      alpha=0.8, levels=np.arange(-1., 1.01, 0.1), vmin=-1., vmax=1., zorder=0)
  ax1.set_title(
      title+' MAR EOF1 (%s' % (round(frac[0], 2))+'%)', fontdict=font, loc='left')
  ax2 = fig.add_subplot(3, 2, 2)
  # ax2.plot(x, le[0], )
  b1 = ax2.bar(x, le[0], color='r')
  for bar, height in zip(b1, le[0]):
    if height < 0:
      bar.set(color='blue')
  barart(ax2)
  ax2.set_title('PC1' % (round(frac[0], 2)), fontdict=font, loc='left'
                )

  ax3 = fig.add_subplot(3, 2, 3, projection=proj)
  mapart(ax3)
  p2 = ax3.contourf(
      lon, lat, lp[1], transfrom=proj, extend='both', cmap=cmaps.ViBlGrWhYeOrRe,
      alpha=0.8, levels=np.arange(-1., 1.01, 0.1), vmin=-1., vmax=1., zorder=0)
  ax3.set_title(
      title+' MAR EOF2 (%s' % (round(frac[1], 2))+"%)", fontdict=font, loc='left')

  ax4 = fig.add_subplot(3, 2, 4)
  b2 = ax4.bar(x, le[1], color='r')
  for bar, height in zip(b2, le[1]):
    if height < 0:
      bar.set(color='blue')
  barart(ax4)
  ax4.set_title('PC2' % (round(frac[1], 2)), fontdict=font, loc='left')

  ax5 = fig.add_subplot(3, 2, 5, projection=proj)
  mapart(ax5)
  p3 = ax5.contourf(
      lon, lat, lp[2], transfrom=proj, extend='both', cmap=cmaps.ViBlGrWhYeOrRe,
      alpha=0.8, levels=np.arange(-1., 1.01, 0.1), vmin=-1., vmax=1., zorder=0)
  ax5.set_title(
      title+' MAR EOF3 (%s' % (round(frac[2], 2))+'%)', fontdict=font, loc='left')

  ax6 = fig.add_subplot(3, 2, 6)
  b3 = ax6.bar(x, le[2], color='r')
  for bar, height in zip(b3, le[2]):
    if height < 0:
      bar.set(color='blue')
  barart(ax6)
  ax6.set_title('PC3' % (round(frac[2], 2)), fontdict=font, loc='left')

  fig.tight_layout()

  fig.subplots_adjust(bottom=0.1)
  l = 0.25
  b = 0.04
  w = 0.6
  h = 0.015

  rect = [l, b, w, h]
  cbar_ax = fig.add_axes(rect)

  c = plt.colorbar(p1, cax=cbar_ax, orientation='horizontal',
                   aspect=20, pad=0.1)
  c.ax.tick_params(labelsize=14)

  # plt.subplots_adjust(wspace=-0.2, hspace=0.3)
  plt.savefig('../../plot/IO_2021/EOF/MAR/EOF_MAR_SST_'+title+'.jpg', dpi=300, format='jpg',
              bbox_inches='tight', transparent=True, pad_inches=0)


# %%
path = '../../data/data/'
LlatS = [-60.,   0., -15., -60., -30., -60.]  # -30.,
LlatN = [30.,  30.,  15.,   0.,   0., -30.]  # 30.
# LlonL = [  50.,  50.,  50., 50.,   50.]
# LlonR = [ 110., 110., 110.,110.,  110.]
LlonL, LlonR = 50., 110.
titles = ['BIO', 'NIO', 'TIO', 'SIO1', 'SIO2', 'SIO3']
font = {'family': 'Times New Roman',
        'weight': 'normal',
        'size': 14,
        }


# %%
for i in range(6):
  extents = [LlonL, LlonR, LlatS[i], LlatN[i]]
  dir = path+'EOF/EOF_MAR_'+titles[i]+'_pickle.dat'
  with open(dir, 'rb') as f:
    eofs = pickle.load(f)
  lp, le, frac = eofs[0], eofs[1], eofs[2]
  plt_eofs(lp, le, frac.data*100, titles[i])
