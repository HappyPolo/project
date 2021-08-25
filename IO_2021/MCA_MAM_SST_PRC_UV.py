'''
                  江城子 . 程序员之歌

              十年生死两茫茫，写程序，到天亮。
                  千行代码，Bug何处藏。
              纵使上线又怎样，朝令改，夕断肠。

              领导每天新想法，天天改，日日忙。
                  相顾无言，惟有泪千行。
              每晚灯火阑珊处，夜难寐，加班狂。

'''

'''
Author       : ChenKt
Date         : 2021-08-24 21:14:02
LastEditors  : ChenKt
LastEditTime : 2021-08-25 10:14:02
FilePath     : /kt/project/IO_2021/MCA_SST_PRC_UV.py
Aim          : 对春季IO(BIO/NIO/SIO) SST 与 降水/风场进行MCA，并保存数据
Mission      : 
  1.读取数据并进行预处理（包括提取所需时间段数据，去趋势，标准化，对维度进行加权，海温还要进行掩膜）
  2.进行MCA计算，并使得对应模态同号，保存数据
  3.区域依次进行
'''
# %%
# %%




import xarray as xr
import numpy as np
import pickle
from xMCA import xMCA
import matplotlib.pyplot as plt
def detrend_dim(da, dim, deg=1):
    # detrend along a single dimension
    p = da.polyfit(dim=dim, deg=deg)
    fit = xr.polyval(da[dim], p.polyfit_coefficients)
    return da - fit


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


def SVD(L, R, N):
    sst_ts = xMCA(L, R)
    sst_ts.solver()
    lp, rp = sst_ts.patterns(n=N)
    lt, rt = sst_ts.expansionCoefs(n=N)
    le, re, lphet, rphet = sst_ts.heterogeneousPatterns(
        n=N, statistical_test=True)
    frac = sst_ts.covFracs(n=N)
    Frac = frac * 100.
    Corr = xr.corr(lt, rt, dim='time')
    return le, re, lt, rt, Frac, lphet, rphet, Corr


def test_pn(le, ule, vle):
    for i in range(3):
        fig, ax1 = plt.subplots(3, figsize=(3, 5))
        le[i].plot(ax=ax1[0], cmap='bwr')
        ule[i].plot(ax=ax1[1], cmap='bwr')
        vle[i].plot(ax=ax1[2], cmap='bwr')


def selYear(da, startYear, endYear):
    startDate = da.sel(time=da.time.dt.year == startYear).time[0]
    endDate = da.sel(time=da.time.dt.year == endYear).time[-1]
    da = da.sel(time=slice(startDate, endDate))
    return da


def selMon(da, Mon):
    return da.sel(time=da.time.dt.month == Mon)


def lsmask(ds, lsdir, label='ocean'):
    landsea = filplonlat(xr.open_dataset(lsdir)['mask'][0])
    ds.coords['mask'] = (('lat', 'lon'), landsea.values)
    if label == 'land':
        ds = ds.where(ds.mask < 1)
    elif label == 'ocean':
        ds = ds.where(ds.mask > 0)
    del ds['mask']
    return ds


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
fp = 'precip.mon.mean.1x1.nc'
fu = 'uwnd.mon.mean.1x1.nc'
fv = 'vwnd.mon.mean.1x1.nc'
fls = 'lsmask.nc'

with xr.open_dataset(path+fs) as f:
    SST = weightslat(filplonlat(f["sst"]))
with xr.open_dataset(path+fp) as f:
    PRC = weightslat(filplonlat(f['precip']))
with xr.open_dataset(path+fu) as f:
    UWD = weightslat(filplonlat(f['uwnd']))
with xr.open_dataset(path+fv) as f:
    VWD = weightslat(filplonlat(f['vwnd']))

N = 3

LlatS = [-60.,   0., -15., -60., -30., -60.]  # -30.,
LlatN = [30.,  30.,  15.,   0.,   0., -30.]  # 30.
# LlonL = [  50.,  50.,  50., 50.,   50.]
# LlonR = [ 110., 110., 110.,110.,  110.]
LlonL, LlonR = 50., 110.

RlatS, RlatN, RlonL, RlonR = -30., 60., 30., 180.

times = ['MAM', 'MAR', 'APR', 'MAY']

# %%
SST = selYear(SST, 1982, 2021)
PRC = selYear(PRC, 1982, 2021)
U = selYear(UWD, 1982, 2021)
V = selYear(VWD, 1982, 2021)
# %%
# note : MAM
sst_MAM = standardize(detrend_dim(
    lsmask(month_to_season(SST, 'MAM'), path+fls), 'time'))
prc_MAM = standardize(detrend_dim(month_to_season(PRC, 'MAM'), 'time'))
u_MAM = standardize(detrend_dim(month_to_season(U, 'MAM'), 'time'))
v_MAM = standardize(detrend_dim(month_to_season(V, 'MAM'), 'time'))
# %%
# note : BIO
I = 0
L = sst_MAM.sel(lat=slice(LlatS[I], LlatN[I]), lon=slice(LlonL, LlonR))

Rp = prc_MAM.sel(lat=slice(RlatS, RlatN), lon=slice(RlonL, RlonR))
Ru = u_MAM.sel(lat=slice(RlatS, RlatN), lon=slice(RlonL, RlonR), level=850)
Rv = v_MAM.sel(lat=slice(RlatS, RlatN), lon=slice(RlonL, RlonR), level=850)
x = Rp.time
# %%
le, re, lt, rt, Frac, lphet, rphet, Corr = SVD(L, Rp, N)
ule, ure, ult, urt, ulphet, Fracu, rphetu, Corru = SVD(L, Ru, N)
vle, vre, vlt, vrt, vlphet, Fracv, rphetv, Corrv = SVD(L, Rv, N)
ure = ure.where(rphetu <= 0.1)
vre = vre.where(rphetv <= 0.1)
# %%
# test_pn(le, ule, vle)

# del L, Rp, Ru, Rv, x
# %%
le[0], re[0], lt[0], rt[0] = -le[0], -re[0], -lt[0], -rt[0]
# ule[0], ure[0], ult[0], urt[0] = -ule[0], -ure[0], -ult[0], -urt[0]
vle[0], vre[0], vlt[0], vrt[0] = -vle[0], -vre[0], -vlt[0], -vrt[0]
# le[1], re[1], lt[1], rt[1] = -le[1], -re[1], -lt[1], -rt[1]
# ule[1], ure[1], ult[1], urt[1] = -ule[1], -ure[1], -ult[1], -urt[1]
vle[1], vre[1], vlt[1], vrt[1] = -vle[1], -vre[1], -vlt[1], -vrt[1]
# le[2], re[2], lt[2], rt[2] = -le[2], -re[2], -lt[2], -rt[2]
# ule[2], ure[2], ult[2], urt[2] = -ule[2], -ure[2], -ult[2], -urt[2]
# vle[2], vre[2], vlt[2], vrt[2] = -vle[2], -vre[2], -vlt[2], -vrt[2]
test_pn(le, ule, vle)
# %%
s_out = [le, re, lt, rt, Frac, lphet, rphet, Corr]
u_out = [ule, ure, ult, urt, ulphet, Fracu, rphetu, Corru]
v_out = [vle, vre, vlt, vrt, vlphet, Fracv, rphetv, Corrv]
svd_out = [s_out, u_out, v_out]
with open(path+'MCA/MCA_BIO_pickle.dat', 'wb') as f:
    pkl = pickle.dump([s_out, u_out, v_out], f, protocol=-1)
#
# %%
# note : for NIO ----------------------------------------------------------------------
I = 1

L = sst_MAM.sel(lat=slice(LlatS[I], LlatN[I]), lon=slice(LlonL, LlonR))

Rp = prc_MAM.sel(lat=slice(RlatS, RlatN), lon=slice(RlonL, RlonR))
Ru = u_MAM.sel(lat=slice(RlatS, RlatN), lon=slice(RlonL, RlonR), level=850)
Rv = v_MAM.sel(lat=slice(RlatS, RlatN), lon=slice(RlonL, RlonR), level=850)
x = Rp.time

# %%
le, re, lt, rt, Frac, lphet, rphet, Corr = SVD(L, Rp, N)
ule, ure, ult, urt, ulphet, Fracu, rphetu, Corru = SVD(L, Ru, N)
vle, vre, vlt, vrt, vlphet, Fracv, rphetv, Corrv = SVD(L, Rv, N)
ure = ure.where(rphetu <= 0.1)
vre = vre.where(rphetv <= 0.1)
# %%
# test_pn(le, ule, vle)

# %%
le[0], re[0], lt[0], rt[0] = -le[0], -re[0], -lt[0], -rt[0]
ule[0], ure[0], ult[0], urt[0] = -ule[0], -ure[0], -ult[0], -urt[0]
# vle[0], vre[0], vlt[0], vrt[0] = -vle[0], -vre[0], -vlt[0], -vrt[0]
le[1], re[1], lt[1], rt[1] = -le[1], -re[1], -lt[1], -rt[1]
# ule[1], ure[1], ult[1], urt[1] = -ule[1], -ure[1], -ult[1], -urt[1]
# vle[1], vre[1], vlt[1], vrt[1] = -vle[1], -vre[1], -vlt[1], -vrt[1]
# le[2], re[2], lt[2], rt[2] = -le[2], -re[2], -lt[2], -rt[2]
# ule[2], ure[2], ult[2], urt[2] = -ule[2], -ure[2], -ult[2], -urt[2]
# vle[2], vre[2], vlt[2], vrt[2] = -vle[2], -vre[2], -vlt[2], -vrt[2]
test_pn(le, ule, vle)

# %%
s_out = [le, re, lt, rt, Frac, lphet, rphet, Corr]
u_out = [ule, ure, ult, urt, ulphet, Fracu, rphetu, Corru]
v_out = [vle, vre, vlt, vrt, vlphet, Fracv, rphetv, Corrv]
svd_out = [s_out, u_out, v_out]
with open(path+'MCA/MCA_NIO_pickle.dat', 'wb') as f:
    pkl = pickle.dump([s_out, u_out, v_out], f, protocol=-1)

# %%
# note : TIO
I = 2
L = sst_MAM.sel(lat=slice(LlatS[I], LlatN[I]), lon=slice(LlonL, LlonR))

Rp = prc_MAM.sel(lat=slice(RlatS, RlatN), lon=slice(RlonL, RlonR))
Ru = u_MAM.sel(lat=slice(RlatS, RlatN), lon=slice(RlonL, RlonR), level=850)
Rv = v_MAM.sel(lat=slice(RlatS, RlatN), lon=slice(RlonL, RlonR), level=850)
x = Rp.time
# %%
le, re, lt, rt, Frac, lphet, rphet, Corr = SVD(L, Rp, N)
ule, ure, ult, urt, ulphet, Fracu, rphetu, Corru = SVD(L, Ru, N)
vle, vre, vlt, vrt, vlphet, Fracv, rphetv, Corrv = SVD(L, Rv, N)
ure = ure.where(rphetu <= 0.1)
vre = vre.where(rphetv <= 0.1)

# %%
# test_pn(le, ule, vle)

# %%
# le[0], re[0], lt[0], rt[0] = -le[0], -re[0], -lt[0], -rt[0]
ule[0], ure[0], ult[0], urt[0] = -ule[0], -ure[0], -ult[0], -urt[0]
# vle[0], vre[0], vlt[0], vrt[0] = -vle[0], -vre[0], -vlt[0], -vrt[0]
# le[1], re[1], lt[1], rt[1] = -le[1], -re[1], -lt[1], -rt[1]
ule[1], ure[1], ult[1], urt[1] = -ule[1], -ure[1], -ult[1], -urt[1]
# vle[1], vre[1], vlt[1], vrt[1] = -vle[1], -vre[1], -vlt[1], -vrt[1]
le[2], re[2], lt[2], rt[2] = -le[2], -re[2], -lt[2], -rt[2]
# ule[2], ure[2], ult[2], urt[2] = -ule[2], -ure[2], -ult[2], -urt[2]
# vle[2], vre[2], vlt[2], vrt[2] = -vle[2], -vre[2], -vlt[2], -vrt[2]
test_pn(le, ule, vle)

# %%
s_out = [le, re, lt, rt, Frac, lphet, rphet, Corr]
u_out = [ule, ure, ult, urt, ulphet, Fracu, rphetu, Corru]
v_out = [vle, vre, vlt, vrt, vlphet, Fracv, rphetv, Corrv]
svd_out = [s_out, u_out, v_out]
with open(path+'MCA/MCA_TIO_pickle.dat', 'wb') as f:
    pkl = pickle.dump([s_out, u_out, v_out], f, protocol=-1)

# %%
# note : SIO1
I = 3
L = sst_MAM.sel(lat=slice(LlatS[I], LlatN[I]), lon=slice(LlonL, LlonR))

Rp = prc_MAM.sel(lat=slice(RlatS, RlatN), lon=slice(RlonL, RlonR))
Ru = u_MAM.sel(lat=slice(RlatS, RlatN), lon=slice(RlonL, RlonR), level=850)
Rv = v_MAM.sel(lat=slice(RlatS, RlatN), lon=slice(RlonL, RlonR), level=850)
x = Rp.time
# %%
le, re, lt, rt, Frac, lphet, rphet, Corr = SVD(L, Rp, N)
ule, ure, ult, urt, ulphet, Fracu, rphetu, Corru = SVD(L, Ru, N)
vle, vre, vlt, vrt, vlphet, Fracv, rphetv, Corrv = SVD(L, Rv, N)
ure = ure.where(rphetu <= 0.1)
vre = vre.where(rphetv <= 0.1)
# %%
# test_pn(le, ule, vle)

# %%
le[0], re[0], lt[0], rt[0] = -le[0], -re[0], -lt[0], -rt[0]
# ule[0], ure[0], ult[0], urt[0] = -ule[0], -ure[0], -ult[0], -urt[0]
vle[0], vre[0], vlt[0], vrt[0] = -vle[0], -vre[0], -vlt[0], -vrt[0]
# le[1], re[1], lt[1], rt[1] = -le[1], -re[1], -lt[1], -rt[1]
# ule[1], ure[1], ult[1], urt[1] = -ule[1], -ure[1], -ult[1], -urt[1]
vle[1], vre[1], vlt[1], vrt[1] = -vle[1], -vre[1], -vlt[1], -vrt[1]
le[2], re[2], lt[2], rt[2] = -le[2], -re[2], -lt[2], -rt[2]
ule[2], ure[2], ult[2], urt[2] = -ule[2], -ure[2], -ult[2], -urt[2]
# vle[2], vre[2], vlt[2], vrt[2] = -vle[2], -vre[2], -vlt[2], -vrt[2]
test_pn(le, ule, vle)

# %%
s_out = [le, re, lt, rt, Frac, lphet, rphet, Corr]
u_out = [ule, ure, ult, urt, ulphet, Fracu, rphetu, Corru]
v_out = [vle, vre, vlt, vrt, vlphet, Fracv, rphetv, Corrv]
svd_out = [s_out, u_out, v_out]
with open(path+'MCA/MCA_SIO1_pickle.dat', 'wb') as f:
    pkl = pickle.dump([s_out, u_out, v_out], f, protocol=-1)

# %%
# note : SIO2
I = 4
L = sst_MAM.sel(lat=slice(LlatS[I], LlatN[I]), lon=slice(LlonL, LlonR))

Rp = prc_MAM.sel(lat=slice(RlatS, RlatN), lon=slice(RlonL, RlonR))
Ru = u_MAM.sel(lat=slice(RlatS, RlatN), lon=slice(RlonL, RlonR), level=850)
Rv = v_MAM.sel(lat=slice(RlatS, RlatN), lon=slice(RlonL, RlonR), level=850)
x = Rp.time
# %%
le, re, lt, rt, Frac, lphet, rphet, Corr = SVD(L, Rp, N)
ule, ure, ult, urt, ulphet, Fracu, rphetu, Corru = SVD(L, Ru, N)
vle, vre, vlt, vrt, vlphet, Fracv, rphetv, Corrv = SVD(L, Rv, N)
ure = ure.where(rphetu <= 0.1)
vre = vre.where(rphetv <= 0.1)
# %%
# test_pn(le, ule, vle)

# %%
# le[0], re[0], lt[0], rt[0] = -le[0], -re[0], -lt[0], -rt[0]
ule[0], ure[0], ult[0], urt[0] = -ule[0], -ure[0], -ult[0], -urt[0]
# vle[0], vre[0], vlt[0], vrt[0] = -vle[0], -vre[0], -vlt[0], -vrt[0]
# le[1], re[1], lt[1], rt[1] = -le[1], -re[1], -lt[1], -rt[1]
# ule[1], ure[1], ult[1], urt[1] = -ule[1], -ure[1], -ult[1], -urt[1]
vle[1], vre[1], vlt[1], vrt[1] = -vle[1], -vre[1], -vlt[1], -vrt[1]
# le[2], re[2], lt[2], rt[2] = -le[2], -re[2], -lt[2], -rt[2]
# ule[2], ure[2], ult[2], urt[2] = -ule[2], -ure[2], -ult[2], -urt[2]
vle[2], vre[2], vlt[2], vrt[2] = -vle[2], -vre[2], -vlt[2], -vrt[2]
test_pn(le, ule, vle)

# %%
s_out = [le, re, lt, rt, Frac, lphet, rphet, Corr]
u_out = [ule, ure, ult, urt, ulphet, Fracu, rphetu, Corru]
v_out = [vle, vre, vlt, vrt, vlphet, Fracv, rphetv, Corrv]
svd_out = [s_out, u_out, v_out]
with open(path+'MCA/MCA_SIO2_pickle.dat', 'wb') as f:
    pkl = pickle.dump([s_out, u_out, v_out], f, protocol=-1)

# %%
# note : SIO3
I = 5
L = sst_MAM.sel(lat=slice(LlatS[I], LlatN[I]), lon=slice(LlonL, LlonR))

Rp = prc_MAM.sel(lat=slice(RlatS, RlatN), lon=slice(RlonL, RlonR))
Ru = u_MAM.sel(lat=slice(RlatS, RlatN), lon=slice(RlonL, RlonR), level=850)
Rv = v_MAM.sel(lat=slice(RlatS, RlatN), lon=slice(RlonL, RlonR), level=850)
x = Rp.time
# %%
le, re, lt, rt, Frac, lphet, rphet, Corr = SVD(L, Rp, N)
ule, ure, ult, urt, ulphet, Fracu, rphetu, Corru = SVD(L, Ru, N)
vle, vre, vlt, vrt, vlphet, Fracv, rphetv, Corrv = SVD(L, Rv, N)
ure = ure.where(rphetu <= 0.1)
vre = vre.where(rphetv <= 0.1)
# %%
# test_pn(le, ule, vle)

# %%
# %%
le[0], re[0], lt[0], rt[0] = -le[0], -re[0], -lt[0], -rt[0]
# ule[0], ure[0], ult[0], urt[0] = -ule[0], -ure[0], -ult[0], -urt[0]
vle[0], vre[0], vlt[0], vrt[0] = -vle[0], -vre[0], -vlt[0], -vrt[0]
le[1], re[1], lt[1], rt[1] = -le[1], -re[1], -lt[1], -rt[1]
# ule[1], ure[1], ult[1], urt[1] = -ule[1], -ure[1], -ult[1], -urt[1]
# vle[1], vre[1], vlt[1], vrt[1] = -vle[1], -vre[1], -vlt[1], -vrt[1]
# le[2], re[2], lt[2], rt[2] = -le[2], -re[2], -lt[2], -rt[2]
ule[2], ure[2], ult[2], urt[2] = -ule[2], -ure[2], -ult[2], -urt[2]
# vle[2], vre[2], vlt[2], vrt[2] = -vle[2], -vre[2], -vlt[2], -vrt[2]
test_pn(le, ule, vle)

# %%
s_out = [le, re, lt, rt, Frac, lphet, rphet, Corr]
u_out = [ule, ure, ult, urt, ulphet, Fracu, rphetu, Corru]
v_out = [vle, vre, vlt, vrt, vlphet, Fracv, rphetv, Corrv]
svd_out = [s_out, u_out, v_out]
with open(path+'MCA/MCA_SIO3_pickle.dat', 'wb') as f:
    pkl = pickle.dump([s_out, u_out, v_out], f, protocol=-1)

# %%
