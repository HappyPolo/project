'''
                               |~~~~~~~|
                               |       |
                               |       |
                               |       |
                               |       |
                               |       |
    |~.\\\_\~~~~~~~~~~~~~~xx~~~         ~~~~~~~~~~~~~~~~~~~~~/_//;~|
    |  \  o \_         ,XXXXX),                         _..-~ o /  |
    |    ~~\  ~-.     XXXXX`)))),                 _.--~~   .-~~~   |
     ~~~~~~~`\   ~\~~~XXX' _/ ';))     |~~~~~~..-~     _.-~ ~~~~~~~
              `\   ~~--`_\~\, ;;;\)__.---.~~~      _.-~
                ~-.       `:;;/;; \          _..-~~
                   ~-._      `''        /-~-~
                       `\              /  /
                         |         ,   | |
                          |  '        /  |
                           \/;          |
                            ;;          |
                            `;   .       |
                            |~~~-----.....|
                           | \             \
                          | /\~~--...__    |
                          (|  `\       __-\|
                          ||    \_   /~    |
                          |)     \~-'      |
                           |      | \      '
                           |      |  \    :
                            \     |  |    |
                             |    )  (    )
                              \  /;  /\  |
                              |    |/   |
                              |    |   |
                               \  .'  ||
                               |  |  | |
                               (  | |  |
                               |   \ \ |
                               || o `.)|
                               |`\\) |
                               |       |
                               |       |
'''

'''
Author       : ChenKt
Date         : 2021-08-25 19:27:22
LastEditors  : ChenKt
LastEditTime : 2021-08-25 20:40:35
FilePath     : /project/IO_2021/EOF_MAY_SST.py
Aim          : 对MAY IO SST进行EOF并保存数据
Mission      :
'''

# %%




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
    sst_ts = xMCA(da, da.rename(da.name+'_copy'))
    sst_ts.solver()
    lp, _ = sst_ts.patterns(n=N)
    le, _ = sst_ts.expansionCoefs(n=N)
    frac = sst_ts.covFracs(n=N)
    return lp, le, frac


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


# %%
path = '../../data/data/'
fs = 'sst.mnmean.nc'
fls = 'lsmask.nc'
with xr.open_dataset(path+fs) as f:
    SST = weightslat(filplonlat(f["sst"]))
N = 3
LlatS = [-60.,   0., -15., -60., -30., -60.]  # -30.,
LlatN = [30.,  30.,  15.,   0.,   0., -30.]  # 30.

LlonL, LlonR = 50., 110.

times = ['MAM', 'MAR', 'APR', 'MAY']
# %%
# note: MAY
SST = selYear(SST, 1982, 2021)
sst_MAY = selMon(SST, 5)
# %%
# note: BIO
I = 0
da = standardize(detrend_dim(sst_MAY, 'time'))
da = lsmask(da, path+fls).sel(
    lat=slice(LlatS[I], LlatN[I]), lon=slice(LlonL, LlonR))
da.name = SST.name
# %%
lp, le, frac = EOF(da, N)

# test_ng(lp,le,N)

# %%
# lp[0], le[0] = -lp[0],-le[0]
# lp[1], le[1] = -lp[1], -le[1]
# lp[2], le[2] = -lp[2], -le[2]
test_ng(lp, le, N)
# %%
with open(path+'EOF/EOF_MAY_BIO_pickle.dat', 'wb') as f:
    pkl = pickle.dump([lp, le, frac], f, protocol=-1)
# %%
# note: NIO
I = 1
da = standardize(detrend_dim(sst_MAY, 'time'))
da = lsmask(da, path+fls).sel(
    lat=slice(LlatS[I], LlatN[I]), lon=slice(LlonL, LlonR))
da.name = SST.name
# %%
lp, le, frac = EOF(da, N)

# test_ng(lp, le, N)
# %%
lp[0], le[0] = -lp[0], -le[0]
# lp[1], le[1] = -lp[1], -le[1]
# lp[2], le[2] = -lp[2], -le[2]
test_ng(lp, le, N)
# %%
with open(path+'EOF/EOF_MAY_NIO_pickle.dat', 'wb') as f:
    pkl = pickle.dump([lp, le, frac], f, protocol=-1)
# %%
# note: TIO
I = 2
da = standardize(detrend_dim(sst_MAY, 'time'))
da = lsmask(da, path+fls).sel(
    lat=slice(LlatS[I], LlatN[I]), lon=slice(LlonL, LlonR))

da.name = SST.name
# %%
lp, le, frac = EOF(da, N)

# test_ng(lp, le, N)
# %%
# lp[0], le[0] = -lp[0],-le[0]
# lp[1], le[1] = -lp[1], -le[1]
# lp[2], le[2] = -lp[2], -le[2]
test_ng(lp, le, N)
# %%
with open(path+'EOF/EOF_MAY_TIO_pickle.dat', 'wb') as f:
    pkl = pickle.dump([lp, le, frac], f, protocol=-1)
# %%
# note: SIO1
I = 3
da = standardize(detrend_dim(sst_MAY, 'time'))
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
with open(path+'EOF/EOF_MAY_SIO1_pickle.dat', 'wb') as f:
    pkl = pickle.dump([lp, le, frac], f, protocol=-1)
# %%
# note: SIO2
I = 4
da = standardize(detrend_dim(sst_MAY, 'time'))
da = lsmask(da, path+fls).sel(
    lat=slice(LlatS[I], LlatN[I]), lon=slice(LlonL, LlonR))

da.name = SST.name
# %%
lp, le, frac = EOF(da, N)

# test_ng(lp, le, N)
# %%
lp[0], le[0] = -lp[0], -le[0]
# lp[1], le[1] = -lp[1], -le[1]
# lp[2], le[2] = -lp[2], -le[2]
test_ng(lp, le, N)
# %%
with open(path+'EOF/EOF_MAY_SIO2_pickle.dat', 'wb') as f:
    pkl = pickle.dump([lp, le, frac], f, protocol=-1)
# %%
# note: SIO3
I = 5
da = standardize(detrend_dim(sst_MAY, 'time'))
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
with open(path+'EOF/EOF_MAY_SIO3_pickle.dat', 'wb') as f:
    pkl = pickle.dump([lp, le, frac], f, protocol=-1)
