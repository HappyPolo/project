'''
 ┌─────────────────────────────────────────────────────────────┐
 │┌───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┬───┐│
 ││Esc│!1 │@2 │#3 │$4 │%5 │^6 │&7 │*8 │(9 │)0 │_- │+= │|\ │`~ ││
 │├───┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴───┤│
 ││ Tab │ Q │ W │ E │ R │ T │ Y │ U │ I │ O │ P │{[ │}] │ BS  ││
 │├─────┴┬──┴┬──┴┬──┴┬──┴┬──┴┬──┴┬──┴┬──┴┬──┴┬──┴┬──┴┬──┴─────┤│
 ││ Ctrl │ A │ S │ D │ F │ G │ H │ J │ K │ L │: ;│" '│ Enter  ││
 │├──────┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴─┬─┴────┬───┤│
 ││ Shift  │ Z │ X │ C │ V │ B │ N │ M │< ,│> .│? /│Shift │Fn ││
 │└─────┬──┴┬──┴──┬┴───┴───┴───┴───┴───┴──┬┴───┴┬──┴┬─────┴───┘│
 │      │Fn │ Alt │         Space         │ Alt │Win│   HHKB   │
 │      └───┴─────┴───────────────────────┴─────┴───┘          │
 └─────────────────────────────────────────────────────────────┘
'''

'''
Author       : ChenKt
Date         : 2021-08-28 20:27:16
LastEditors  : ChenKt
LastEditTime : 2021-08-28 20:35:11
FilePath     : /project/IO_2021/MIs.py
Aim          : with long term trend
Mission      :
'''
# %%
# %%




import matplotlib.ticker as ticker
import matplotlib.pyplot as plt
import pickle
import numpy as np
import xarray as xr
import scipy.stats as stats
def standardize(x):
  return (x - x.mean()) / x.std()


def SAM(v):
  # V850 - V200
  lon = v.lon
  lat = v.lat
  lon_range = lon[(lon >= 70.) & (lon <= 110.)]
  lat_range = lat[(lat >= 10.) & (lat <= 30.)]
  v850 = v.sel(level=850, lon=lon_range,
               lat=lat_range).mean(dim=['lat', 'lon']).drop('level')
  v200 = v.sel(level=200, lon=lon_range,
               lat=lat_range).mean(dim=['lat', 'lon']).drop('level')
  sam = standardize(v850-v200)
  # sam = v850-v200
  return sam


def WY(u):
  # U850 - U200
  lon = u.lon
  lat = u.lat
  lon_range = lon[(lon >= 40.) & (lon <= 110.)]
  lat_range = lat[(lat >= 5.) & (lat <= 20.)]
  u850 = u.sel(level=850, lon=lon_range,
               lat=lat_range).mean(dim=['lat', 'lon']).drop('level')
  u200 = u.sel(level=200, lon=lon_range,
               lat=lat_range).mean(dim=['lat', 'lon']).drop('level')
  # wy = u850-u200
  wy = standardize(u850-u200)
  return wy


def SEAM(u):
  lon = u.lon
  lat = u.lat
  lon_range1 = lon[(lon >= 90.) & (lon <= 130.)]
  lat_range1 = lat[(lat >= 5.) & (lat <= 15.)]
  lon_range2 = lon[(lon >= 110.) & (lon <= 140.)]
  lat_range2 = lat[(lat >= 22.) & (lat <= 33.)]
  u1 = u.sel(level=850, lon=lon_range1,
             lat=lat_range1).mean(dim=['lat', 'lon']).drop('level')
  u2 = u.sel(level=850, lon=lon_range2,
             lat=lat_range2).mean(dim=['lat', 'lon']).drop('level')
  # seam = u1-u2
  seam = standardize(u1-u2)
  return seam


def EAM(u):
  lon = u.lon
  lat = u.lat
  lon_range1 = lon[(lon >= 110.) & (lon <= 150.)]
  lat_range1 = lat[(lat >= 25.) & (lat <= 35.)]
  lon_range2 = lon[(lon >= 110.) & (lon <= 150.)]
  lat_range2 = lat[(lat >= 40.) & (lat <= 50.)]
  u1 = u.sel(level=200, lon=lon_range1,
             lat=lat_range1).mean(dim=['lat', 'lon']).drop('level')
  u2 = u.sel(level=200, lon=lon_range2,
             lat=lat_range2).mean(dim=['lat', 'lon']).drop('level')
  # eam = u1-u2
  eam = standardize(u1-u2)
  return eam


def month_to_season(xMon, season):
  """ This function takes an xarray dataset containing monthly data spanning years and
      returns a dataset with one sample per year, for a specified three-month season.

      Time stamps are centered on the season, e.g. seasons='DJF' returns January timestamps.

      If a calculated season's timestamp falls outside the original range of monthly values, then the calculated mean
      is dropped.  For example, if the monthly data's time range is [Jan-2000, Dec-2003] and the season is "DJF", the
      seasonal mean computed from the single month of Dec-2003 is dropped.
  """
  startDate = xMon.time[0]
  endDate = xMon.time[-1]
  seasons_pd = {
      'DJF': ('QS-DEC', 1),
      'JFM': ('QS-JAN', 2),
      'FMA': ('QS-FEB', 3),
      'MAM': ('QS-MAR', 4),
      'AMJ': ('QS-APR', 5),
      'MJJ': ('QS-MAY', 6),
      'JJA': ('QS-JUN', 7),
      'JAS': ('QS-JUL', 8),
      'ASO': ('QS-AUG', 9),
      'SON': ('QS-SEP', 10),
      'OND': ('QS-OCT', 11),
      'NDJ': ('QS-NOV', 12)
  }
  try:
    (season_pd, season_sel) = seasons_pd[season]
  except KeyError:
    raise ValueError("contributed: month_to_season: bad season: SEASON = " +
                     season)

  # Compute the three-month means, moving time labels ahead to the middle month.
  month_offset = 'MS'
  xSeasons = xMon.resample(time=season_pd, loffset=month_offset).mean()

  # Filter just the desired season, and trim to the desired time range.
  xSea = xSeasons.sel(time=xSeasons.time.dt.month == season_sel)
  xSea = xSea.sel(time=slice(startDate, endDate))
  return xSea


def filplonlat(ds):
  # To facilitate data subsetting
  # print(da.attrs)
  '''
  print(
      f'\n\nBefore flip, lon range is [{ds["lon"].min().data}, {ds["lon"].max().data}].'
  )
  ds["lon"] = ((ds["lon"] + 180) % 360) - 180
  # Sort lons, so that subset operations end up being simpler.
  ds = ds.sortby("lon")
  '''

  ds = ds.sortby("lat", ascending=True)
  # print(ds.attrs)
  print('\n\nAfter sorting lat values, ds["lat"] is:')
  print(ds["lat"])
  return ds


def detrend_dim(da, dim, deg=1):
  # detrend along a single dimension
  p = da.polyfit(dim=dim, deg=deg)
  fit = xr.polyval(da[dim], p.polyfit_coefficients)
  return da - fit


def selYear(da, startYear, endYear):
  startDate = da.sel(time=da.time.dt.year == startYear).time[0]
  endDate = da.sel(time=da.time.dt.year == endYear).time[-1]
  da = da.sel(time=slice(startDate, endDate))
  return da


# %%
path = '../../../data/data/'
fp = 'precip.mon.mean.nc'
fu = 'uwnd.mon.mean.nc'
fv = 'vwnd.mon.mean.nc'
with xr.open_dataset(path+fp) as f:
  PRC = filplonlat(f['precip'])
with xr.open_dataset(path+fu) as f:
  UWD = filplonlat(f['uwnd'])
with xr.open_dataset(path+fv) as f:
  VWD = filplonlat(f['vwnd'])

PRC = selYear(PRC, 1982, 2021)
U = selYear(UWD, 1982, 2021)
V = selYear(VWD, 1982, 2021)

# prc_JJA = standardize(detrend_dim(month_to_season(PRC, 'JJA'), 'time'))
# u_JJA = standardize(detrend_dim(month_to_season(U, 'JJA'), 'time'))
# v_JJA = standardize(detrend_dim(month_to_season(V, 'JJA'), 'time'))
prc_JJA = standardize(month_to_season(PRC, 'JJA'))
u_JJA = standardize(month_to_season(U, 'JJA'))
v_JJA = standardize(month_to_season(V, 'JJA'))
# %%
SAMIs_JJA = SAM(v_JJA)
WYIs_JJA = WY(u_JJA)
SEAMIs_JJA = SEAM(u_JJA)
EAMIs_JJA = EAM(u_JJA)
SAMIs_JJA.name = 'SAM'
WYIs_JJA.name = 'WY'
SEAMIs_JJA.name = 'SEAM'
EAMIs_JJA.name = 'EAM'
# %%
MIs = [SAMIs_JJA, WYIs_JJA, SEAMIs_JJA, EAMIs_JJA]
with open(path+'wt/MIs_JJA.dat', 'wb') as f:
  pkl = pickle.dump(MIs, f, protocol=-1)

# %%
SAMIs_JJA.plot()
# %%

fig = plt.figure(figsize=(12, 4))
ax = plt.subplot(111)
EAMIs_JJA.plot(ax=ax, label='EAM', lw=3)
SAMIs_JJA.plot(ax=ax, label='SAM', lw=3)
SEAMIs_JJA.plot(ax=ax, label='SEAM', lw=3)
WYIs_JJA.plot(ax=ax, label='WY', lw=3)
# %%


EAMIs_JJA.name


# %%
SAMIs_JJA.plot()
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


def plt_xy(MI):
  fig, ax = plt.subplots(4, 1, sharex='all', sharey='all', figsize=(8, 5))
  fig.subplots_adjust(hspace=0)
  # b = PC
  x = np.arange(1982, 2022, 1)
  for i in range(4):
    a = MI[i]
    # a['time'] = b['time']
    # R = xr.corr(a, b, dim='time')
    # if abs(R.values) >= 0.312:
    lbl = a.name
    # else:
    print(a.name)
    # lbl = a.name + ' R = %s' % (np.around(R, 3).values)+''
    ax[i].plot(x, MIs[i], color='black', ls='-', label=a.name, lw=1.5)
    # br = ax[i].bar(x, b, color='tomato')
    # for bar, height in zip(br, b):
    #   if height < 0:
    #     bar.set(color='slateblue')
    lineart(ax[i])
    ax[0].set_title('JJA', loc='left')

  savedir = '../../../plot/IO_2021/MIs/wt_MI_JJA.png'
  if os.path.isfile(savedir):
    os.remove(savedir)
  fig.savefig(savedir, dpi=600,
              bbox_inches='tight', pad_inches=0.0, facecolor='w', format='png')  # %%


# %%
plt_xy(MIs)

# %%
standardize(EAMIs_JJA).plot()
# %%
