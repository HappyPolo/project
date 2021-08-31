'''
Author       : ChenKt
Date         : 2021-08-29 17:05:56
LastEditors  : ChenKt
LastEditTime : 2021-08-30 18:34:13
FilePath     : /project/IO_2021/plt_MI_PCs_APR.py
Aim          : 
Mission      : 
'''
# %%
import os
import xarray as xr
from xMCA import xMCA
import pickle
import matplotlib.pyplot as plt
import numpy as np
import matplotlib.ticker as ticker
# %%


def lineart(ax):
  ax.set_ylim(-3., 3.)
  ax.set_xlim(1980, 2025)
  ax.tick_params(axis='both', which='major', labelsize=11, direction='out',
                 length=2.8, width=0.9, pad=0.2, top=False, right=False)
  ax.yaxis.set_minor_locator(ticker.AutoMinorLocator())
  ax.xaxis.set_tick_params(which='minor', bottom=False)
  ax.xaxis.set_tick_params(which='major', top=False)
  ax.axhline(y=0, c="black", ls="-", lw=0.5)
  ax.legend(loc='lower right', frameon=False, fontsize=9, ncol=4)
  ax.tick_params(labelsize=8)
  labels = ax.get_xticklabels() + ax.get_yticklabels()
  [label.set_fontname('Times New Roman') for label in labels]


def plt_xy(MI, PC, tles, savedir):
  fig, ax = plt.subplots(4, 1, sharex='all', sharey='all', figsize=(8, 5))
  fig.subplots_adjust(hspace=0)
  b = PC
  x = np.arange(1982, 2022, 1)
  for i in range(4):
    a = MI[i]
    a['time'] = b['time']
    R = xr.corr(a, b, dim='time')
    if abs(R.values) >= 0.312:
      lbl = a.name + ' R = %s' % (np.around(R, 3).values)+'${*}$'
    else:
      lbl = a.name + ' R = %s' % (np.around(R, 3).values)+''
    ax[i].plot(x, MIs[i], color='black', ls='-', label=lbl, lw=1.5)
    br = ax[i].bar(x, b, color='tomato')
    for bar, height in zip(br, b):
      if height < 0:
        bar.set(color='slateblue')
    lineart(ax[i])
    ax[0].set_title(tles, loc='left')

  if os.path.isfile(savedir):
    os.remove(savedir)
  fig.savefig(savedir, dpi=600,
              bbox_inches='tight', pad_inches=0.0, facecolor='w', format='png')  # %%


# %%
path = '../../data/data/'
titles = ['BIO', 'NIO', 'TIO', 'SIO1', 'SIO2', 'SIO3']
# %%
dir1 = path+'MIs_JJA.dat'
with open(dir1, 'rb') as f:
  MIs = pickle.load(f)
for j in range(6):
  for k in range(3):
    dir = path+'EOF/EOF_APR_'+titles[j]+'_pickle.dat'
    with open(dir, 'rb') as f:
      eofs = pickle.load(f)
    lt = eofs[1]
    title = 'APR '+titles[j]+' PC'+str(k+1)
    savedir = '../../plot/IO_2021/MIs/APR/APR_'+titles[j]+'_PC'+str(k+1)+'.png'
    plt_xy(MIs, lt[k], title, savedir)
