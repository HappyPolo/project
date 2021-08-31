'''
Author       : ChenKt
Date         : 2021-08-28 19:29:25
LastEditors  : ChenKt
LastEditTime : 2021-08-30 18:05:56
FilePath     : /project/IO_2021/plt_MCA_MAY_SST_PRC_UV.py
Aim          : 
Mission      : 
'''
# %%from matplotlib.patches import Rectangle
from cartopy.util import add_cyclic_point
from cartopy.mpl.ticker import LatitudeFormatter, LongitudeFormatter
import matplotlib.gridspec as gridspec
import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import matplotlib.colors as colors
import cmaps
import pickle
import xarray as xr
import numpy as np
from matplotlib.patches import Rectangle


def mapart1(ax, extents):
  proj = ccrs.PlateCarree()
  ax.coastlines(color='k', lw=.7)
  ax.add_feature(cfeature.LAND, facecolor='white')
  # 设置地图范围
  ax.set_extent(extents, crs=proj)
  xticks = np.arange(extents[0], extents[1]+1, 20)
  yticks = np.arange(extents[2], extents[3]+1, 15)
  ax.set_xticks(xticks, crs=proj)
  ax.set_yticks(yticks, crs=proj)
  lon_formatter = LongitudeFormatter(zero_direction_label=True)
  lat_formatter = LatitudeFormatter()
  ax.xaxis.set_major_formatter(lon_formatter)
  ax.yaxis.set_major_formatter(lat_formatter)
  ax.tick_params(axis='both', which='major', labelsize=7, direction='out',
                 length=5, width=0.9, pad=0.2, top=True, right=True)
  ax.minorticks_on()
  ax.tick_params(axis='both', which='minor', direction='out',
                 width=0.6, top=True, right=True)


def mapart2(ax, extents):
  proj = ccrs.PlateCarree()
  ax.coastlines(color='k', lw=.7)
  ax.add_feature(cfeature.LAND, facecolor='white')

  # 设置地图范围

  ax.set_extent(extents, crs=proj)
  xticks = np.arange(extents[0], extents[1]+1, 40)
  yticks = np.arange(extents[2], extents[3]+1, 30)
  ax.set_xticks(xticks, crs=proj)
  ax.set_yticks(yticks, crs=proj)
  lon_formatter = LongitudeFormatter(zero_direction_label=True)
  lat_formatter = LatitudeFormatter()
  ax.xaxis.set_major_formatter(lon_formatter)
  ax.yaxis.set_major_formatter(lat_formatter)
  ax.tick_params(axis='both', which='major', labelsize=7, direction='out',
                 length=5, width=0.9, pad=0.2, top=True, right=True)
  ax.minorticks_on()
  ax.tick_params(axis='both', which='minor', direction='out',
                 width=0.6, top=True, right=True)


def lineart(ax):
  ax.set_ylim(-3., 3.)
  ax.tick_params(axis='both', which='major', labelsize=11, direction='out',
                 length=2.8, width=0.9, pad=0.2, top=False, right=False)
  ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
  ax.xaxis.set_tick_params(which='minor', bottom=False)
  ax.xaxis.set_tick_params(which='major', top=False)
  ax.axhline(y=0, c="black", ls="-", lw=0.5)
  ax.legend(loc='upper left', frameon=False, fontsize=5, ncol=4)
  ax.tick_params(labelsize=6)
  labels = ax.get_xticklabels() + ax.get_yticklabels()
  [label.set_fontname('Times New Roman') for label in labels]


def cbart(cb):
  cb.set_label('', fontdict=font)
  cb.ax.tick_params(labelsize=6)
  cb.ax.tick_params(which='major', direction='out',
                    labelsize=6, length=2)  # 主刻度设置


def plt_sig(da, ax, n, i, area):
  da_cyc, lon_cyc = add_cyclic_point(da[i, ::n, ::n], coord=da.lon[::n])
  nx, ny = np.meshgrid(lon_cyc, da.lat[::n])
  sig = ax.scatter(
      nx[area], ny[area], marker='.', s=1.5, c='red', alpha=0.6, transform=proj
  )


def plt_mcas_BIO(s, u, v, extents):
  for i in range(N):
    # plt.close()
    fig = plt.figure()
    gs = fig.add_gridspec(14, 13, wspace=2, hspace=.1)
    ax1 = fig.add_subplot(gs[4:14, 0:6], projection=proj)
    ax2 = fig.add_subplot(gs[0:14, 7:12], projection=proj)
    ax3 = fig.add_subplot(gs[11:14, 7:13])

    mapart1(ax1, extents)

    im = ax1.contourf(
        s[0].lon, s[0].lat, s[0][i], transform=proj, extend='neither', cmap=cmaps.ViBlGrWhYeOrRe, alpha=0.8, levels=np.arange(-1., 1.09, 0.1), zorder=0
    )
    ax1.set_title('(a) '+MON+' mode'+str(i+1), fontdict=font, loc='left')
    ax1.set_title('{:.2f}'.format(
        s[4][i].values)+'%, R='+'{:.2f}'.format(s[7][i].values), fontdict=font, loc='right')

    cbposition = fig.add_axes([0.45, 0.12, 0.015, 0.55])

    cb1 = fig.colorbar(im, cax=cbposition, orientation='vertical',
                       spacing='proportional', format='%.1f', extend='both')
    cbart(cb1)
    area = np.where(s[5][i, ::3, ::3].data < 0.05)
    n = 3
    plt_sig(s[5], ax1, n, i, area)

    extents2 = [30, 180, -30, 60]
    mapart2(ax2, extents2)
    color_map1 = ['#af6c29', '#c68741', '#dba767', '#f0ce87', '#f2e0a1', '#fdf5bd',
                  '#ffffff', '#ffffff', '#ffffff', '#ffffff',
                  '#d6fecf', '#aaefb4', '#6abda0', '#7aabc8', '#2b71b0', '#103557']
    color_map = colors.ListedColormap(color_map1)
    rm = ax2.contourf(
        s[1].lon, s[1].lat, s[1][i], transform=proj, extend='neither', cmap=color_map, alpha=0.8, levels=np.arange(-1., 1.09, 0.1), zorder=0
    )
    ax2.set_title('(b) '+MON+' BIO PRC', fontdict=font, loc='left')
    cbposition = fig.add_axes([0.85, 0.38, 0.015, 0.25])
    cb2 = fig.colorbar(rm, cax=cbposition, orientation='vertical',
                       spacing='proportional', format='%.1f', extend='both')

    cbart(cb2)

    area = np.where(s[6][i, ::6, ::6].data < 0.05)
    n = 6
    plt_sig(s[6], ax2, n, i, area)
    Q = ax2.quiver(
        u[1].lon[::6].data, u[1].lat[::6].data,
        u[1][i, ::6, ::6].data, v[1][i, ::6, ::6].data, transform=proj, zorder=2,
        units='xy', angles='uv', scale=0.1, scale_units='xy',
        width=0.8
    )
    w, h = 0.12, 0.12
    rect = Rectangle(
        (1 - w, 0), w, h, transform=ax2.transAxes,
        fc='white', ec='k', lw=0.5, zorder=1.1
    )
    ax2.add_patch(rect)
    # 添加quiverkey.
    # U指定风箭头对应的速度.
    qk = ax2.quiverkey(
        Q, X=1-w/2, Y=0.7*h, U=1,
        label=f'{1}', labelpos='S', labelsep=0.02,
        fontproperties={'size': 4}
    )

    x = s[2].time

    ax3.plot(x, s[2][i], ls='-', c='b', label='sst')
    ax3.plot(x, s[3][i], ls='-', c='r', label='prc')
    lineart(ax3)
    ax3.set_title('(c) ', fontdict=font, loc='left')
    import os
    savedir = '../../plot/IO_2021/MCA/MAY/MAY_BIO_MCA'+str(i+1)+'.png'
    if os.path.isfile(savedir):
        os.remove(savedir)
    fig.savefig(savedir, dpi=600,
                bbox_inches='tight', pad_inches=0.0, facecolor='w', format='png')  # %%


def plt_mcas_NIO(s, u, v, extents):
  for i in range(N):
    # plt.close()
    fig = plt.figure()
    gs = fig.add_gridspec(14, 13, wspace=2, hspace=.1)
    ax1 = fig.add_subplot(gs[4:14, 0:6], projection=proj)
    ax2 = fig.add_subplot(gs[0:14, 7:12], projection=proj)
    ax3 = fig.add_subplot(gs[11:14, 7:13])

    mapart1(ax1, extents)

    im = ax1.contourf(
        s[0].lon, s[0].lat, s[0][i], transform=proj, extend='neither', cmap=cmaps.ViBlGrWhYeOrRe, alpha=0.8, levels=np.arange(-1., 1.09, 0.1), zorder=0)
    ax1.set_title('(a) '+MON+' mode'+str(i+1), fontdict=font, loc='left')
    ax1.set_title('{:.2f}'.format(
        s[4][i].values)+'%, R='+'{:.2f}'.format(s[7][i].values), fontdict=font, loc='right')

    cbposition = fig.add_axes([0.12, 0.14, 0.35, 0.03])

    cb1 = fig.colorbar(im, cax=cbposition, orientation='horizontal',
                       spacing='proportional', format='%.1f', extend='both')
    cbart(cb1)
    area = np.where(s[5][i, ::3, ::3].data < 0.05)
    n = 3
    plt_sig(s[5], ax1, n, i, area)

    extents2 = [30, 180, -30, 60]
    mapart2(ax2, extents2)
    color_map1 = ['#af6c29', '#c68741', '#dba767', '#f0ce87', '#f2e0a1', '#fdf5bd',
                  '#ffffff', '#ffffff', '#ffffff', '#ffffff',
                  '#d6fecf', '#aaefb4', '#6abda0', '#7aabc8', '#2b71b0', '#103557']
    color_map = colors.ListedColormap(color_map1)
    rm = ax2.contourf(
        s[1].lon, s[1].lat, s[1][i], transform=proj, extend='neither', cmap=color_map, alpha=0.8, levels=np.arange(-1., 1.09, 0.1), zorder=0
    )
    ax2.set_title('(b) '+MON+' NIO PRC', fontdict=font, loc='left')
    cbposition = fig.add_axes([0.85, 0.38, 0.015, 0.25])
    cb2 = fig.colorbar(rm, cax=cbposition, orientation='vertical',
                       spacing='proportional', format='%.1f', extend='both')
    cbart(cb2)

    area = np.where(s[6][i, ::6, ::6].data < 0.05)
    n = 6
    plt_sig(s[6], ax2, n, i, area)
    Q = ax2.quiver(
        u[1].lon[::6].data, u[1].lat[::6].data,
        u[1][i, ::6, ::6].data, v[1][i, ::6, ::6].data, transform=proj, zorder=2,
        units='xy', angles='uv', scale=0.1, scale_units='xy',
        width=0.8
    )
    w, h = 0.12, 0.12
    rect = Rectangle(
        (1 - w, 0), w, h, transform=ax2.transAxes,
        fc='white', ec='k', lw=0.5, zorder=1.1
    )
    ax2.add_patch(rect)
    # 添加quiverkey.
    # U指定风箭头对应的速度.
    qk = ax2.quiverkey(
        Q, X=1-w/2, Y=0.7*h, U=1,
        label=f'{1}', labelpos='S', labelsep=0.02,
        fontproperties={'size': 4}
    )

    x = s[2].time

    ax3.plot(x, s[2][i], ls='-', c='b', label='sst')
    ax3.plot(x, s[3][i], ls='-', c='r', label='prc')
    lineart(ax3)
    ax3.set_title('(c) ', fontdict=font, loc='left')
    import os
    savedir = '../../plot/IO_2021/MCA/MAY/MAY_NIO_MCA'+str(i+1)+'.png'
    if os.path.isfile(savedir):
        os.remove(savedir)
    fig.savefig(savedir, dpi=600,
                bbox_inches='tight', pad_inches=0.0, facecolor='w', format='png')  # %%


def plt_mcas_TIO(s, u, v, extents):
  for i in range(N):
    # plt.close()
    fig = plt.figure()
    gs = fig.add_gridspec(14, 13, wspace=2, hspace=.1)
    ax1 = fig.add_subplot(gs[4:14, 0:6], projection=proj)
    ax2 = fig.add_subplot(gs[0:14, 7:12], projection=proj)
    ax3 = fig.add_subplot(gs[11:14, 7:13])

    mapart1(ax1, extents)

    im = ax1.contourf(
        s[0].lon, s[0].lat, s[0][i], transform=proj, extend='neither', cmap=cmaps.ViBlGrWhYeOrRe, alpha=0.8, levels=np.arange(-1., 1.09, 0.1), zorder=0
    )
    ax1.set_title('(a) '+MON+' mode'+str(i+1), fontdict=font, loc='left')
    ax1.set_title('{:.2f}'.format(
        s[4][i].values)+'%, R='+'{:.2f}'.format(s[7][i].values), fontdict=font, loc='right')

    cbposition = fig.add_axes([0.12, 0.14, 0.35, 0.03])

    cb1 = fig.colorbar(im, cax=cbposition, orientation='horizontal',
                       spacing='proportional', format='%.1f', extend='both')
    cbart(cb1)
    area = np.where(s[5][i, ::3, ::3].data < 0.05)
    n = 3
    plt_sig(s[5], ax1, n, i, area)

    extents2 = [30, 180, -30, 60]
    mapart2(ax2, extents2)
    color_map1 = ['#af6c29', '#c68741', '#dba767', '#f0ce87', '#f2e0a1', '#fdf5bd',
                  '#ffffff', '#ffffff', '#ffffff', '#ffffff',
                  '#d6fecf', '#aaefb4', '#6abda0', '#7aabc8', '#2b71b0', '#103557']
    color_map = colors.ListedColormap(color_map1)
    rm = ax2.contourf(
        s[1].lon, s[1].lat, s[1][i], transform=proj, extend='neither', cmap=color_map, alpha=0.8, levels=np.arange(-1., 1.09, 0.1), zorder=0
    )
    ax2.set_title('(b) '+MON+' TIO PRC', fontdict=font, loc='left')
    cbposition = fig.add_axes([0.85, 0.38, 0.015, 0.25])
    cb2 = fig.colorbar(rm, cax=cbposition, orientation='vertical',
                       spacing='proportional', format='%.1f', extend='both')
    cbart(cb2)

    area = np.where(s[6][i, ::6, ::6].data < 0.05)
    n = 6
    plt_sig(s[6], ax2, n, i, area)
    Q = ax2.quiver(
        u[1].lon[::6].data, u[1].lat[::6].data,
        u[1][i, ::6, ::6].data, v[1][i, ::6, ::6].data, transform=proj, zorder=2,
        units='xy', angles='uv', scale=0.1, scale_units='xy',
        width=0.8
    )
    w, h = 0.12, 0.12
    rect = Rectangle(
        (1 - w, 0), w, h, transform=ax2.transAxes,
        fc='white', ec='k', lw=0.5, zorder=1.1
    )
    ax2.add_patch(rect)
    # 添加quiverkey.
    # U指定风箭头对应的速度.
    qk = ax2.quiverkey(
        Q, X=1-w/2, Y=0.7*h, U=1,
        label=f'{1}', labelpos='S', labelsep=0.02,
        fontproperties={'size': 4}
    )

    x = s[2].time

    ax3.plot(x, s[2][i], ls='-', c='b', label='sst')
    ax3.plot(x, s[3][i], ls='-', c='r', label='prc')
    lineart(ax3)
    ax3.set_title('(c) ', fontdict=font, loc='left')
    import os
    savedir = '../../plot/IO_2021/MCA/MAY/MAY_TIO_MCA'+str(i+1)+'.png'
    if os.path.isfile(savedir):
      os.remove(savedir)
    fig.savefig(savedir, dpi=600,
                bbox_inches='tight', pad_inches=0.0, facecolor='w', format='png')  # %%


def plt_mcas_SIO2(s, u, v, extents):
  for i in range(N):
    # plt.close()
    fig = plt.figure()
    gs = fig.add_gridspec(14, 13, wspace=2, hspace=.1)
    ax1 = fig.add_subplot(gs[4:14, 0:6], projection=proj)
    ax2 = fig.add_subplot(gs[0:14, 7:12], projection=proj)
    ax3 = fig.add_subplot(gs[11:14, 7:13])

    mapart1(ax1, extents)

    im = ax1.contourf(
        s[0].lon, s[0].lat, s[0][i], transform=proj, extend='neither', cmap=cmaps.ViBlGrWhYeOrRe, alpha=0.8, levels=np.arange(-1., 1.09, 0.1), zorder=0
    )
    ax1.set_title('(a) '+MON+' mode'+str(i+1), fontdict=font, loc='left')
    ax1.set_title('{:.2f}'.format(
        s[4][i].values)+'%, R='+'{:.2f}'.format(s[7][i].values), fontdict=font, loc='right')

    cbposition = fig.add_axes([0.12, 0.14, 0.35, 0.03])

    cb1 = fig.colorbar(im, cax=cbposition, orientation='horizontal',
                       spacing='proportional', format='%.1f', extend='both')
    cbart(cb1)
    area = np.where(s[5][i, ::3, ::3].data < 0.05)
    n = 3
    plt_sig(s[5], ax1, n, i, area)

    extents2 = [30, 180, -30, 60]
    mapart2(ax2, extents2)
    color_map1 = ['#af6c29', '#c68741', '#dba767', '#f0ce87', '#f2e0a1', '#fdf5bd',
                  '#ffffff', '#ffffff', '#ffffff', '#ffffff',
                  '#d6fecf', '#aaefb4', '#6abda0', '#7aabc8', '#2b71b0', '#103557']
    color_map = colors.ListedColormap(color_map1)
    rm = ax2.contourf(
        s[1].lon, s[1].lat, s[1][i], transform=proj, extend='neither', cmap=color_map, alpha=0.8, levels=np.arange(-1., 1.09, 0.1), zorder=0
    )
    ax2.set_title('(b) '+MON+' SIO1 PRC', fontdict=font, loc='left')
    cbposition = fig.add_axes([0.85, 0.38, 0.015, 0.25])
    cb2 = fig.colorbar(rm, cax=cbposition, orientation='vertical',
                       spacing='proportional', format='%.1f', extend='both')
    cbart(cb2)

    area = np.where(s[6][i, ::6, ::6].data < 0.05)
    n = 6
    plt_sig(s[6], ax2, n, i, area)
    Q = ax2.quiver(
        u[1].lon[::6].data, u[1].lat[::6].data,
        u[1][i, ::6, ::6].data, v[1][i, ::6, ::6].data, transform=proj, zorder=2,
        units='xy', angles='uv', scale=0.1, scale_units='xy',
        width=0.8
    )
    w, h = 0.12, 0.12
    rect = Rectangle(
        (1 - w, 0), w, h, transform=ax2.transAxes,
        fc='white', ec='k', lw=0.5, zorder=1.1
    )
    ax2.add_patch(rect)
    # 添加quiverkey.
    # U指定风箭头对应的速度.
    qk = ax2.quiverkey(
        Q, X=1-w/2, Y=0.7*h, U=1,
        label=f'{1}', labelpos='S', labelsep=0.02,
        fontproperties={'size': 4}
    )

    x = s[2].time

    ax3.plot(x, s[2][i], ls='-', c='b', label='sst')
    ax3.plot(x, s[3][i], ls='-', c='r', label='prc')
    lineart(ax3)
    ax3.set_title('(c) ', fontdict=font, loc='left')
    savedir = '../../plot/IO_2021/MCA/MAY/MAY_SIO2_MCA'+str(i+1)+'.png'
    import os
    if os.path.isfile(savedir):
      os.remove(savedir)
    fig.savefig(savedir, dpi=600,
                bbox_inches='tight', pad_inches=0.0, facecolor='w', format='png')  # %%


def plt_mcas_SIO3(s, u, v, extents):
  for i in range(N):
    # plt.close()
    fig = plt.figure()
    gs = fig.add_gridspec(14, 13, wspace=2, hspace=.1)
    ax1 = fig.add_subplot(gs[4:14, 0:6], projection=proj)
    ax2 = fig.add_subplot(gs[0:14, 7:12], projection=proj)
    ax3 = fig.add_subplot(gs[11:14, 7:13])

    mapart1(ax1, extents)

    im = ax1.contourf(
        s[0].lon, s[0].lat, s[0][i], transform=proj, extend='neither', cmap=cmaps.ViBlGrWhYeOrRe, alpha=0.8, levels=np.arange(-1., 1.09, 0.1), zorder=0
    )
    ax1.set_title('(a) '+MON+' mode'+str(i+1), fontdict=font, loc='left')
    ax1.set_title('{:.2f}'.format(
        s[4][i].values)+'%, R='+'{:.2f}'.format(s[7][i].values), fontdict=font, loc='right')

    cbposition = fig.add_axes([0.12, 0.14, 0.35, 0.03])

    cb1 = fig.colorbar(im, cax=cbposition, orientation='horizontal',
                       spacing='proportional', format='%.1f', extend='both')
    cbart(cb1)
    area = np.where(s[5][i, ::3, ::3].data < 0.05)
    n = 3
    plt_sig(s[5], ax1, n, i, area)

    extents2 = [30, 180, -30, 60]
    mapart2(ax2, extents2)
    color_map1 = ['#af6c29', '#c68741', '#dba767', '#f0ce87', '#f2e0a1', '#fdf5bd',
                  '#ffffff', '#ffffff', '#ffffff', '#ffffff',
                  '#d6fecf', '#aaefb4', '#6abda0', '#7aabc8', '#2b71b0', '#103557']
    color_map = colors.ListedColormap(color_map1)
    rm = ax2.contourf(
        s[1].lon, s[1].lat, s[1][i], transform=proj, extend='neither', cmap=color_map, alpha=0.8, levels=np.arange(-1., 1.09, 0.1), zorder=0
    )
    ax2.set_title('(b) '+MON+' SIO2 PRC', fontdict=font, loc='left')
    cbposition = fig.add_axes([0.85, 0.38, 0.015, 0.25])
    cb2 = fig.colorbar(rm, cax=cbposition, orientation='vertical',
                       spacing='proportional', format='%.1f', extend='both')
    cbart(cb2)

    area = np.where(s[6][i, ::6, ::6].data < 0.05)
    n = 6
    plt_sig(s[6], ax2, n, i, area)
    Q = ax2.quiver(
        u[1].lon[::6].data, u[1].lat[::6].data,
        u[1][i, ::6, ::6].data, v[1][i, ::6, ::6].data, transform=proj, zorder=2,
        units='xy', angles='uv', scale=0.1, scale_units='xy',
        width=0.8
    )
    w, h = 0.12, 0.12
    rect = Rectangle(
        (1 - w, 0), w, h, transform=ax2.transAxes,
        fc='white', ec='k', lw=0.5, zorder=1.1
    )
    ax2.add_patch(rect)
    # 添加quiverkey.
    # U指定风箭头对应的速度.
    qk = ax2.quiverkey(
        Q, X=1-w/2, Y=0.7*h, U=1,
        label=f'{1}', labelpos='S', labelsep=0.02,
        fontproperties={'size': 4}
    )

    x = s[2].time

    ax3.plot(x, s[2][i], ls='-', c='b', label='sst')
    ax3.plot(x, s[3][i], ls='-', c='r', label='prc')
    lineart(ax3)
    ax3.set_title('(c) ', fontdict=font, loc='left')
    import os
    savedir = '../../plot/IO_2021/MCA/MAY/MAY_SIO3_MCA'+str(i+1)+'.png'
    if os.path.isfile(savedir):
      os.remove(savedir)
    fig.savefig(savedir, dpi=600,
                bbox_inches='tight', pad_inches=0.0, facecolor='w', format='png')  # %%


def plt_mcas_SIO1(s, u, v, extents):
  for i in range(N):
    # plt.close()
    fig = plt.figure()
    gs = fig.add_gridspec(14, 13, wspace=2, hspace=.1)
    ax1 = fig.add_subplot(gs[2:14, 0:6], projection=proj)
    ax2 = fig.add_subplot(gs[0:14, 7:12], projection=proj)
    ax3 = fig.add_subplot(gs[11:14, 7:13])

    mapart1(ax1, extents)

    im = ax1.contourf(
        s[0].lon, s[0].lat, s[0][i], transform=proj, extend='neither', cmap=cmaps.ViBlGrWhYeOrRe, alpha=0.8, levels=np.arange(-1., 1.09, 0.1), zorder=0
    )
    ax1.set_title('(a) '+MON+' mode'+str(i+1), fontdict=font, loc='left')
    ax1.set_title('{:.2f}'.format(
        s[4][i].values)+'%, R='+'{:.2f}'.format(s[7][i].values), fontdict=font, loc='right')

    cbposition = fig.add_axes([0.12, 0.1, 0.35, 0.03])

    cb1 = fig.colorbar(im, cax=cbposition, orientation='horizontal',
                       spacing='proportional', format='%.1f', extend='both')
    cbart(cb1)
    area = np.where(s[5][i, ::3, ::3].data < 0.05)
    n = 3
    plt_sig(s[5], ax1, n, i, area)

    extents2 = [30, 180, -30, 60]
    mapart2(ax2, extents2)
    color_map1 = ['#af6c29', '#c68741', '#dba767', '#f0ce87', '#f2e0a1', '#fdf5bd',
                  '#ffffff', '#ffffff', '#ffffff', '#ffffff',
                  '#d6fecf', '#aaefb4', '#6abda0', '#7aabc8', '#2b71b0', '#103557']
    color_map = colors.ListedColormap(color_map1)
    rm = ax2.contourf(
        s[1].lon, s[1].lat, s[1][i], transform=proj, extend='neither', cmap=color_map, alpha=0.8, levels=np.arange(-1., 1.09, 0.1), zorder=0
    )
    ax2.set_title('(b) '+MON+' SIO3 PRC', fontdict=font, loc='left')
    cbposition = fig.add_axes([0.85, 0.38, 0.015, 0.25])
    cb2 = fig.colorbar(rm, cax=cbposition, orientation='vertical',
                       spacing='proportional', format='%.1f', extend='both')
    cbart(cb2)

    area = np.where(s[6][i, ::6, ::6].data < 0.05)
    n = 6
    plt_sig(s[6], ax2, n, i, area)
    Q = ax2.quiver(
        u[1].lon[::6].data, u[1].lat[::6].data,
        u[1][i, ::6, ::6].data, v[1][i, ::6, ::6].data, transform=proj, zorder=2,
        units='xy', angles='uv', scale=0.1, scale_units='xy',
        width=0.8
    )
    w, h = 0.12, 0.12
    rect = Rectangle(
        (1 - w, 0), w, h, transform=ax2.transAxes,
        fc='white', ec='k', lw=0.5, zorder=1.1
    )
    ax2.add_patch(rect)
    # 添加quiverkey.
    # U指定风箭头对应的速度.
    qk = ax2.quiverkey(
        Q, X=1-w/2, Y=0.7*h, U=1,
        label=f'{1}', labelpos='S', labelsep=0.02,
        fontproperties={'size': 4}
    )

    x = s[2].time

    ax3.plot(x, s[2][i], ls='-', c='b', label='sst')
    ax3.plot(x, s[3][i], ls='-', c='r', label='prc')
    lineart(ax3)
    ax3.set_title('(c) ', fontdict=font, loc='left')
    import os
    savedir = '../../plot/IO_2021/MCA/MAY/MAY_SIO1_MCA'+str(i+1)+'.png'
    if os.path.isfile(savedir):
      os.remove(savedir)
    fig.savefig(savedir, dpi=600,
                bbox_inches='tight', pad_inches=0.0, facecolor='w', format='png')  # %%


# %%
N = 3
proj = ccrs.PlateCarree()
MON = 'MAY'
font = {'family': 'Times New Roman',
        'weight': 'normal',
        'size': 6,
        }


# %%
# note:for BIO
dir = '../../data/data/MCA/MCA_MAY_BIO_pickle.dat'
with open(dir, 'rb') as f:
  mcas = pickle.load(f)
extents = [50, 110, -60, 30]
s, u, v = mcas[0], mcas[1], mcas[2]
extents = [50, 110, -60, 30]
plt_mcas_BIO(s, u, v, extents)

# %%
# note:for NIO
dir = '../../data/data/MCA/MCA_MAY_NIO_pickle.dat'
with open(dir, 'rb') as f:
  mcas = pickle.load(f)
s, u, v = mcas[0], mcas[1], mcas[2]
extents = [50, 110, 0, 30]
plt_mcas_NIO(s, u, v, extents)
# %%
# note:for TIO
dir = '../../data/data/MCA/MCA_MAY_TIO_pickle.dat'
with open(dir, 'rb') as f:
  mcas = pickle.load(f)
s, u, v = mcas[0], mcas[1], mcas[2]
extents = [50, 110, -15, 15]
plt_mcas_TIO(s, u, v, extents)

# %%
# note:for SIO1

dir = '../../data/data/MCA/MCA_MAY_SIO1_pickle.dat'
with open(dir, 'rb') as f:
  mcas = pickle.load(f)
s, u, v = mcas[0], mcas[1], mcas[2]
extents = [50, 110, -60, 0]
plt_mcas_SIO1(s, u, v, extents)

# %%
# note:for SIO2

dir = '../../data/data/MCA/MCA_MAY_SIO2_pickle.dat'
with open(dir, 'rb') as f:
  mcas = pickle.load(f)
s, u, v = mcas[0], mcas[1], mcas[2]
extents = [50, 110, -30, 0]
plt_mcas_SIO2(s, u, v, extents)
# %%
# note:for SIO3

dir = '../../data/data/MCA/MCA_MAY_SIO3_pickle.dat'
with open(dir, 'rb') as f:
  mcas = pickle.load(f)
s, u, v = mcas[0], mcas[1], mcas[2]
extents = [50, 110, -60, -30]
plt_mcas_SIO3(s, u, v, extents)
