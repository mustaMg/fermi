#%%
from decimal import Decimal
from astropy.io import fits
from numpy import sqrt
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
# reading fits file
hdul = fits.open('/Volumes/GoogleDrive/My Drive/Python/Fermi/swift/heasoft/data/200303onesec_4ch.lc')

# creating colom names
e_min = (pd.DataFrame(hdul[2].data))['E_MIN'].to_list()
e_min = ([int(i) for i in e_min])
e_max = (pd.DataFrame(hdul[2].data))['E_MAX'].to_list()
e_max = ([int(i) for i in e_max])
column = [f'{e_min[i]}-{e_max[i]}'for i in range(len(e_min))]

NGOODPIX = hdul[1].header['NGOODPIX']   # value that will be used later

rate = pd.DataFrame(columns=column)     # count rate 
for i, ii in enumerate(column):         # creating dataframe
    rate[str(ii)] = [x[1][i]*NGOODPIX for x in hdul[1].data]
rate['15-50'] = rate['15-25'] + rate['25-50']

time = pd.DataFrame([hdul[1].data[i][0] for i in range(len(hdul[1].data))])
TRIGTIME = hdul[1].header['TRIGTIME']
time = time - TRIGTIME
rate.index = time[0]

rate = rate[(-100 < rate.index) & (rate.index < 350)]

error = pd.DataFrame([hdul[1].data[i][2] for i in range(len(hdul[1].data))],\
    index=time[0], columns=column)
print(error)
error = error * NGOODPIX

error = error[(-100 < error.index) & (error.index < 350)]
error['15-50'] = sqrt((error['15-25'])**2 + (error['25-50'])**2)

tbin = time[0][1]- time[0][0]

#######################################################################
# UNWEIGHTED
path = '/Users/mustafagumustas/Downloads/Swift_BAT/bat_data/grb200303A/LC/unweighted.lc'

hhdul = fits.open(path)
bg = pd.DataFrame(columns=column)       # unweighted count rate
for i, ii in enumerate(column):
    bg[str(ii)] = [x[1][i] for x in hhdul[1].data]
bg['15-50'] = bg['15-25'] + bg['25-50'] # creating new data 

time = pd.DataFrame([hhdul[1].data[i][0] for i in range(len(hhdul[1].data))])
TRIGTIME = hhdul[1].header['TRIGTIME']
time = time - TRIGTIME
bg.index = time[0]

bg = bg[(-100 < bg.index) & (bg.index < 350)]

dbg = sqrt(bg * tbin) / tbin    # error of unweighted data
dbg['15-50'] = (dbg['15-25'] + dbg['25-50']) / 2


#######################################################################
# REBIN DATA
# selecting data between 5-350 secs

def rebin_not(df, n):
    tablo = (df.groupby(df.index // n).keys)
    unique, counts = np.unique(tablo, return_counts=True)
    sonuc = (dict(zip(unique, counts)))
    last = {
    k: v
    for k, v in sonuc.items()
    if v != n}
    df = (df.groupby(df.index // n).sum())
    df = df.drop(list(last.keys()))
    return df

def rebin_v4_not(df, n):
    remain = (len(df) % n)
    df = df[:df.index[-(remain+1)]]
    b_df = [sum(i) for i in np.array_split( df, (len(df)/n) ) if sum(i) != 0]
    return pd.DataFrame(b_df)

def rebin_v5(df, n):
    a_df = np.array(df)
    remain = (len(df) % n)
    a_df = a_df[:-remain]
    return [sum(i) for i in a_df.reshape(len(df)//n,n)]

def nD_bin(df, n):
    ndf = pd.DataFrame()
    for i in df:
        ndf[i] = rebin_v5(df[i],n)
    return ndf
   
rate2 = rate[(5 < rate.index) & (rate.index < 350)]
binned_rate2  = nD_bin(rate2, 4)

error2 = error[(5 < error.index) & (error.index < 350)]
bin_error2 = nD_bin(error2, 4)

bbg = bg[(5 < bg.index) & (bg.index < 350)]
bin_bbg = nD_bin(bbg, 4)

dbbg = dbg[(5 < dbg.index) & (dbg.index < 350)]
bin_dbbg = nD_bin(dbbg, 4)

####################
# SIGNAL TO NOISE
def div_nD_df(df1, df2):
    # it has to be the same size and same column names
    result = pd.DataFrame()
    for i in df1:
        result[i] = [x/y for x, y in zip(df1[i], df2[i])]
    return result

def sum_nD_df(df1, df2):
    # it has to be the same size and same column names
    result = pd.DataFrame()
    for i in df1:
        result[i] = [x + y for x, y in zip(np.array(df1[i]), np.array(df2[i]))]

def mult_nD_df(df1, df2):
    # it has to be the same size and same column names
    result = pd.DataFrame()
    for i in df1:
        result[i] = [float(round(Decimal(x*y),4)) for x, y in zip(np.array(df1[i]), np.array(df2[i]))]

    return result

def extr_nD_df(df1, df2):
        # it has to be the same size and same column names
    result = pd.DataFrame()
    for i in df1:
        result[i] = [x-abs(y) for x, y in zip(df1[i], df2[i])]
    return result

def SNR(live, phot, dphot, cr, dcr):
    counts = cr * live
    dcounts = dcr * live

    sigma = sqrt( phot * live ) 
    dsigma = ( 1/2 * div_nD_df(dphot,phot) ) 
    dsigma.index = sigma.index
    dsigma = mult_nD_df(dsigma,sigma)
    # SQRT((dcounts / counts)^2 + (dsig02 / sig02)^2) * snr02
    snr = div_nD_df(counts,sigma)
    dsnr = sqrt( (div_nD_df(dcounts,counts)**2) + (div_nD_df(dsigma,sigma)**2)) * snr
    return snr, dsnr


# SIGNAL TO NOISE   To-100 : To+5
snr, dsnr= SNR(4, bin_bbg, bin_dbbg, binned_rate2, bin_error2)
# SIGNAL TO NOISE   To-100 : To+5s
rate3 = rate[rate.index < 5]
rate3= rate3[rate3 != 0]

error3 = error[error.index < 5]
error3= error3[error3 != 0]
print(error3)
print(rate3)
bbg2 = bg[bg.index < 5]
bbg2= bbg2[bbg2 != 0]

dbbg2 = dbg[dbg.index < 5]
dbbg2= dbbg2[dbbg2 != 0]

# print(rate2)
# print(rate3)
snr2, dsnr2 = SNR(1, bbg2, dbbg2, rate3, error3)
 
# ------------------ #
def snr_diff(col1, col2):
    df = pd.DataFrame()
    df['snr'] = col1
    df['dsnr'] = col2
    df['EE'] = df['snr']- abs(df['dsnr'])
    df['index'] = df.index
    gt = []
    values = []
    for i in range(len(df)-1):
        if (df['EE'].iloc[i]) > 1.5:
            # if (df['index'].iloc[i+1]) - (df['index'].iloc[i]) == 1:
            values.append(df['index'].iloc[i])
        elif (df['EE'].iloc[i]) < 1.5 and len(values) >3 and values not in gt:
            gt.append(values)
        else:
            values = []
    # for i, ii in zip(df.index, df['EE']):
        # print(i, ii)
    return gt


def EE_finder(snr, dsnr):
    gt = {}
    for i in snr:
        gt[i] = (snr_diff(snr[i], dsnr[i]))
    return gt


# print(EE_finder(snr2,dsnr2))
# print(snr2)
# print(rate3)
# print(bbg2)
plt.step(snr2.index, (snr2['100-150']-dsnr2['100-150']), zorder=1)
# plt.step(snr2.index, rate3['100-150'])
for keys, values in EE_finder(snr2,dsnr2).items():
    if keys == '100-150':
        (keys, values)
        xx = max(values[0])
        y = min(values[0]) -1
        print(y, xx)
        plt.axvline(x=xx, color='red')
        plt.axvline(x=y, color='red')
        plt.axvspan(xx, y, alpha=0.1, color='red', hatch='/')

yy = [i-(ii) for i, ii in zip(snr2['100-150'], dsnr2['100-150'])]
xx = [i-0.5 for i in snr2.index]
plt.errorbar(x=xx, y=yy,ls='none', ecolor='black', yerr=dsnr2['100-150'], zorder=2)#     
plt.axhline(y=1.5, color='red')
# plt.show()

#  15-25       25-50      50-100    100-150       15-50

# rate          count rate
# error         error of count
# bg            unwe
# ighted data count rate
# dbg           error of unweighted
# rate2         new selection of rate, limiting time between 5350
# error2        limiting error same above
# bbg           limiting unweighted with new time interval
# dbbg          same above for unweighted error
# binned_rate2  binning rate data with 4 sec
# binned_error2 binning error with 4 sec bins
# binned        binning unweighted count
# dbinned       binning error of unweighted
#%%