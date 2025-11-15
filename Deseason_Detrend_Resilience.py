#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import sys
import numpy as np
import pandas as pd
import geopandas as gpd

import scipy.stats
from scipy.signal import argrelmin
from scipy.optimize import curve_fit
import scipy.stats as st

from multiprocessing import Pool
import itertools, time, os
import traceback, shutil

import numpy.ma as ma
import glob

import warnings
warnings.filterwarnings('ignore')

#%%
def prep_file(fid, eco=False):
    data = pd.read_csv(fid)
    gb = data.groupby('fid')
    sers = {}
    for g in gb:
        if not eco:
            idnr = g[0]
        else:
            idn, eid = g[0]
            idnr = str(eid) + '_' + str(idn) 
        date = g[1].date
        val = g[1]['mean'].values
        date = [pd.Timestamp(x) for x in date]
        ser = pd.Series(val, index=date)
        sers[idnr] = ser
    return sers

def calc_ar1(x):
    return ma.corrcoef(ma.masked_invalid(x[:-1]), ma.masked_invalid(x[1:]))[0,1]

def compute_lam(x, dt=1):
    dx = (x[1:] - x[:-1]) / dt
    x0 = x[:-1]
    mask = ~np.isnan(x0) & ~np.isnan(dx)
    return st.linregress(x0[mask], dx[mask])[0]

def compute_sigma(x, dt=1):
    dx = (x[1:] - x[:-1]) / dt
    lamb = compute_lam(x, dt)
    diff = dx - lamb * x[:-1]
    return np.nanstd(diff) * np.sqrt(dt)

def harmonic_fit(ser, order=3):
    import statsmodels.api as sm
    harm_freq = list(range(1, order+1))
    
    x, y = ser.index, ser.values
    
    x = np.array([(x - pd.Timestamp('1970-01-01')).days for x in x]) / 365.25
    x_rad = x * 2 * np.pi #Convert days to radians for harmonic fitting
    
    #Create empty array to hold the independents
    nr_indep = order*2 + 2
    indep = np.empty((y.shape[0], nr_indep))
    
    #Add constant for intercept and then time
    indep[:,0] = 1
    indep[:,1] = x_rad
    
    #Now create the harmonic variables
    i = 2
    for freq in harm_freq:
        cos = np.cos(x_rad * freq)
        sin = np.sin(x_rad * freq)
        indep[:,i] = cos
        i = i + 1
        indep[:,i] = sin
        i = i + 1
        
    model = sm.OLS(y, indep, missing='drop').fit()
    coefs = model.params
    fitted = []
    for t in range(x_rad.shape[0]):
        data = indep[t,:]
        harm_term = np.nansum(coefs*data)
        fitted.append(harm_term)
    fitted = np.array(fitted)
    return pd.Series(ser.values - fitted, index=ser.index)

def runmean(x, w):
    n = x.shape[0]
    xs = np.zeros_like(x)
    for i in range(w // 2):
        xs[i] = np.nanmean(x[: i + w // 2 + 1])
    for i in range(n - w // 2, n):
        xs[i] = np.nanmean(x[i - w // 2 + 1:])
    for i in range(w // 2, n - w // 2):
        xs[i] = np.nanmean(x[i - w // 2 : i + w // 2 + 1])
    return x - xs

def deseason_detrend(ser, yrs=5, yl=23):
    rm_offline = pd.Series(runmean(ser.values, yrs*yl), index=ser.index)
    deseason_rolling = harmonic_fit(rm_offline, order=3)
    return deseason_rolling

def proc_ts(rint, ts, res):
    outdf = pd.DataFrame()
    outdf['rint'] = [rint]
    outdf['fullmn'] = [np.nanmean(ts.values)]
    
    car_raw = calc_ar1(ts.values)
    car_detrend = calc_ar1(res.values)
    if np.isnan(car_raw):
        car_raw = np.nan
    if np.isnan(car_detrend):
        car_detrend = np.nan
    outdf['ar1_raw'] = [car_raw]
    outdf['ar1_resid'] = [car_detrend]
    lg = np.log(car_detrend)
    if np.isnan(lg):
        lg = np.nan
    outdf['lambda_ar1'] = [lg]
    
    outdf['variance_raw'] = [np.nanvar(ts.values)]
    var = np.nanvar(res.values)
    outdf['variance_resid'] = [var]
    
    sigma = compute_sigma(res.values)
    outdf['sigma_lamb'] = [sigma]
    rr_var = 0.5 * np.log(1-sigma**2 / var)
    if np.isnan(rr_var):
        rr_var = np.nan
    outdf['lambda_variance'] = [rr_var]
    
    return outdf

def run_data(fid):
    ser_dict = prep_file(fid)
    
    output = []
    for idnr in ser_dict.keys():
        ser = ser_dict[idnr]
        res = deseason_detrend(ser, yrs=5, yl=23) #Hard coded for MODIS NDVI/EVI/kNDVI, change for LAI/GPP
        try:
            outdf = proc_ts(idnr, ser, res) #Get the resilience metrics, as well as a mean
            output.append(outdf)
        except:
            pass
    try:
        full = pd.concat(output).reset_index(drop=True)
        return full
    except:
        pass
         

#%%
ty = 'EVI'
res = '5k'

#Get a list of all the individual blocks of time series created via Earth Engine
fidlist = glob.glob(ty + '_MODLC_' + res + '/*StratSamp*')
compiled_fid = ty + '_' + res + '_MODLCStratSamp_Trends.csv'

pool = Pool(processes=n_proc)
results = list(pool.imap_unordered(run_data, fidlist))
pool.close()
pool.join()

df = pd.concat(results).reset_index(drop=True)
df.to_csv(compiled_fid)
del df