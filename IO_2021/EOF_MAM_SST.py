'''
Author       : ChenKt
Date         : 2021-08-29 19:50:40
LastEditors  : ChenKt
LastEditTime : 2021-08-30 21:51:37
FilePath     : /project/IO_2021/EOF_MAM_SST.py
Aim          : 
Mission      : 
'''
# %%
from eofs.xarray import Eof
import scipy.stats as stats
import xarray as xr
from xMCA import xMCA
import pickle
import matplotlib.pyplot as plt
import numpy as np


def standardize(x):
  return (x - x.mean()) / x.std()


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
  # print(
  #     f'\n\nAfter flip, lon range is [{ds["lon"].min().data}, {ds["lon"].max().data}].'
  # )
  # To facilitate data subsetting

  ds = ds.sortby("lat", ascending=True)
  # print(ds.attrs)
  print('\n\nAfter sorting lat values, ds["lat"] is:')
  print(ds["lat"])
  return ds


def weightslat(ds):
  deg2rad = np.pi/180.
  clat = ds['lat'].astype(np.float64)
  clat = np.sqrt(np.cos(deg2rad * clat))
  print('\n\nclat:\n')
  print(clat)
  # Xarray will apply lat-based weights to all lons and timesteps automatically.
  # This is called "broadcasting".
  wds = ds
  wds.attrs = ds.attrs
  ds = clat * ds
  wds.attrs['long_name'] = 'Wgt: '+wds.attrs['long_name']
  return wds


def EOF(da, N):
  coslat = np.cos(np.deg2rad(da.coords['lat'].values)).clip(0., 1.)
  wgts = np.sqrt(coslat)[..., np.newaxis]
  solver = Eof(da, weights=wgts)
  eof = solver.eofsAsCovariance(neofs=N, pcscaling=1)
  pc = solver.pcs(npcs=N, pcscaling=1)
  pc = pc.transpose('mode', 'time')
  vF = solver.varianceFraction(neigs=3)
  return eof, pc, vF


def selYear(da, startYear, endYear):
  startDate = da.sel(time=da.time.dt.year == startYear).time[0]
  endDate = da.sel(time=da.time.dt.year == endYear).time[-1]
  da = da.sel(time=slice(startDate, endDate))
  return da


def selMon(da, Mon):
  return da.sel(time=da.time.dt.month == Mon)


def lsmask(ds, lsdir, label='ocean'):
  with xr.open_dataset(lsdir) as f:
    dd = f['mask'][0]
  landsea = filplonlat(dd)
  ds.coords['mask'] = (('lat', 'lon'), landsea.values)
  if label == 'land':
    ds = ds.where(ds.mask < 1)
  elif label == 'ocean':
    ds = ds.where(ds.mask > 0)
  del ds['mask']
  return ds


def detrend_dim(da, dim, deg=1):
  # detrend along a single dimension
  p = da.polyfit(dim=dim, deg=deg)
  fit = xr.polyval(da[dim], p.polyfit_coefficients)
  return da - fit


def test_ng(eof, pc, N):
  fig, ax1 = plt.subplots(N, 2, figsize=(9, 8))
  for i in range(N):
    eof[i].plot(ax=ax1[i][0], cmap='bwr')
    pc[i].plot(ax=ax1[i][1])


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


# %%
path = '../../data/data/'
fs = 'sst.mnmean.nc'
fls = 'lsmask.nc'
with xr.open_dataset(path+fs) as f:
  SST = filplonlat(f["sst"])
N = 3
LlatS = [-60.,   0., -15., -60., -30., -60.]  # -30.,
LlatN = [30.,  30.,  15.,   0.,   0., -30.]  # 30.

LlonL, LlonR = 50., 110.

times = ['MAM', 'MAR', 'APR', 'MAM']
# %%
# note: MAM
SST = selYear(SST, 1982, 2021)
# sst_MAM = selMon(SST, 5)
sst_MAM = month_to_season(SST, "MAM")

# %%
# note: BIO
I = 0
da = standardize(detrend_dim(sst_MAM, 'time'))
da = lsmask(da, path+fls).sel(
    lat=slice(LlatS[I], LlatN[I]), lon=slice(LlonL, LlonR))
da.name = SST.name
# %%
lp, le, frac = EOF(da, N)

# test_ng(lp,le,N)

# %%
# lp[0], le[0] = -lp[0], -le[0]
lp[1], le[1] = -lp[1], -le[1]
lp[2], le[2] = -lp[2], -le[2]
test_ng(lp, le, N)
# %%
with open(path+'EOF/EOF_MAM_BIO_pickle.dat', 'wb') as f:
  pkl = pickle.dump([lp, le, frac], f, protocol=-1)
# %%
# note: NIO
I = 1
da = standardize(detrend_dim(sst_MAM, 'time'))
da = lsmask(da, path+fls).sel(
    lat=slice(LlatS[I], LlatN[I]), lon=slice(LlonL, LlonR))
da.name = SST.name
# %%
lp, le, frac = EOF(da, N)

# test_ng(lp, le, N)
# %%
lp[0], le[0] = -lp[0], -le[0]
lp[1], le[1] = -lp[1], -le[1]
# lp[2], le[2] = -lp[2], -le[2]
test_ng(lp, le, N)
# %%
with open(path+'EOF/EOF_MAM_NIO_pickle.dat', 'wb') as f:
  pkl = pickle.dump([lp, le, frac], f, protocol=-1)
# %%
# note: TIO
I = 2
del da
# da = standardize(detrend_dim(sst_MAM, 'time'))
da = standardize(sst_MAM)
da = standardize(detrend_dim(sst_MAM, 'time'))
da = lsmask(da, path+fls).sel(
    lat=slice(LlatS[I], LlatN[I]), lon=slice(LlonL, LlonR))

da.name = SST.name
# %%
lp, le, frac = EOF(da, N)
print(frac)
# test_ng(lp, le, N)
# %%
# lp[0], le[0] = -lp[0], -le[0]
lp[1], le[1] = -lp[1], -le[1]
# lp[2], le[2] = -lp[2], -le[2]
test_ng(lp, le, N)
# %%
with open(path+'EOF/EOF_MAM_TIO_pickle.dat', 'wb') as f:
  pkl = pickle.dump([lp, le, frac], f, protocol=-1)
# %%
# note: SIO1
I = 3
da = standardize(detrend_dim(sst_MAM, 'time'))
da = lsmask(da, path+fls).sel(
    lat=slice(LlatS[I], LlatN[I]), lon=slice(LlonL, LlonR))

da.name = SST.name
# %%
lp, le, frac = EOF(da, N)

# test_ng(lp, le, N)
# %%
# lp[0], le[0] = -lp[0], -le[0]
lp[1], le[1] = -lp[1], -le[1]
# lp[2], le[2] = -lp[2], -le[2]
test_ng(lp, le, N)
# %%
with open(path+'EOF/EOF_MAM_SIO1_pickle.dat', 'wb') as f:
  pkl = pickle.dump([lp, le, frac], f, protocol=-1)
# %%
# note: SIO2
I = 4
da = standardize(detrend_dim(sst_MAM, 'time'))
da = lsmask(da, path+fls).sel(
    lat=slice(LlatS[I], LlatN[I]), lon=slice(LlonL, LlonR))

da.name = SST.name
# %%
lp, le, frac = EOF(da, N)

# test_ng(lp, le, N)
# %%
# lp[0], le[0] = -lp[0], -le[0]
# lp[1], le[1] = -lp[1], -le[1]
# lp[2], le[2] = -lp[2], -le[2]
test_ng(lp, le, N)
# %%
with open(path+'EOF/EOF_MAM_SIO2_pickle.dat', 'wb') as f:
  pkl = pickle.dump([lp, le, frac], f, protocol=-1)
# %%
# note: SIO3
I = 5
da = standardize(detrend_dim(sst_MAM, 'time'))
da = lsmask(da, path+fls).sel(
    lat=slice(LlatS[I], LlatN[I]), lon=slice(LlonL, LlonR))

da.name = SST.name
# %%
lp, le, frac = EOF(da, N)

# test_ng(lp, le, N)
# %%
lp[0], le[0] = -lp[0], -le[0]
lp[1], le[1] = -lp[1], -le[1]
lp[2], le[2] = -lp[2], -le[2]
test_ng(lp, le, N)
# %%
with open(path+'EOF/EOF_MAM_SIO3_pickle.dat', 'wb') as f:
  pkl = pickle.dump([lp, le, frac], f, protocol=-1)
# %%
frac
# %%
