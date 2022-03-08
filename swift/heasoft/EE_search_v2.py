import math
from astropy.io import fits
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt

# reading weighted data 
hdul = fits.open('/Volumes/GoogleDrive/My Drive/Python/Fermi/swift/heasoft/data/200410onesec_4ch.lc')

#creating tiem index values
time = pd.DataFrame([i[0] for i in hdul[1].data])
TRIGTIME = hdul[1].header['TRIGTIME']
time = time - TRIGTIME

# creating count and error values in DF
count = pd.DataFrame([i[1] for i in hdul[1].data], index=time[0] ) 

error = pd.DataFrame([i[2] for i in hdul[1].data], index=time[0])

# multiplying NGOODPIX with every count value
NGOODPIX = hdul[1].header['NGOODPIX']
count = count * NGOODPIX
error = error * NGOODPIX

# filtering data between -100 sec and 250
count = count[(count.index > -100) & (count.index < 350)]
error = error[(error.index > -100) & (error.index < 350)]

# creating 15-50 energy channel count and error values
count[4] = count[0] + count[1] 
error[4] = np.sqrt( (error[0]**2) + (error[1]**2) )

################################################################################
#                                                                              #
#                                UNWEIGHTED                                    #
#                                                                              #
################################################################################

hdul2 = fits.open('/Users/mustafagumustas/Downloads/Swift_BAT/bat_data/grb200410A/LC/unweighted.lc')

bg = pd.DataFrame([i[1] for i in hdul2[1].data], index=time[0])
bg = bg[(bg.index > -100) & (bg.index < 350)]
dbg = np.sqrt(bg)

bg[4] = bg[0] + bg[1]
dbg[4] = (dbg[0] + dbg[1]) / 2

################################################################################
#                                                                              #
#                                REBIN DATA                                    #
#                                                                              #
################################################################################

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

# bin the datas after 5sec
rcount = count[(count.index) > 5]
rerror = error[(error.index) > 5]
rbg    = bg[(bg.index) > 5]
rdbg   = dbg[(dbg.index) > 5]

rcount = nD_bin(rcount, 4)
rerror = nD_bin(rerror, 4)
rbg    = nD_bin(rbg, 4)
rdbg   = nD_bin(rdbg, 4)

# fixing their index 
rcount.index = [i for i in range(5,349, 4)]
rerror.index = [i for i in range(5,349, 4)]
rbg.index    = [i for i in range(5,349, 4)]
rdbg.index   = [i for i in range(5,349, 4)]


################################################################################
#                                                                              #
#                              SIGNAL TO NOISE                                 #
#                                                                              #
################################################################################

def SNR(live, phot, dphot, cr, dcr):
    counts = cr * live
    dcounts = dcr * live 

    phot = phot * live
    sigma = np.sqrt(phot)
    dsigma = ( (1/2) * (dphot/phot) ) * sigma

    snr = counts / sigma
    dsnr = ( np.sqrt( ((dcounts/phot)**2) + ((dsigma/sigma)**2) ) ) * snr
    return snr , dsnr

# snr of 5 sec and after
snr ,dsnr = SNR(4, rbg, rdbg, rcount, rerror)

# snr of 5 sec and before
bg2    = bg[(bg.index > -100) & (bg.index < 5)]
dbg2   = dbg[(dbg.index > -100) & (dbg.index < 5)]
count2 = count[(count.index > -100) & (count.index < 5)]
error2  = error[(error.index > -100) & (error.index < 5)]

snr2, dsnr2 = SNR(4, bg2, dbg2, count2, error2)

snr2.index = [i for i in range(-100, 5)]
dsnr2.index = [i for i in range(-100, 5)]


################################################################################
#                                                                              #
#                           FIND EXTENDED EMISSION                             #
#                                                                              #
################################################################################



def snr_diff(col1, col2):
    df = pd.DataFrame()
    df['snr'] = col1; df['dsnr'] = col2; df['EE'] = df['snr']- (abs(df['dsnr']))
    df['index'] = df.index
    gt = [];    values = []
    for i in range(len(df)-1):
        if (df['EE'].iloc[i]) > 1.5:
            values.append(df['index'].iloc[i])
        elif (df['EE'].iloc[i]) < 1.5 and len(values) >=3 \
            and values not in gt:
            gt.append(values)
            values = []
        else:
            values = []
    return gt
def EE_finder(snr, dsnr):
    gt = {}
    for i in snr:
        gt[i] = (snr_diff(snr[i], dsnr[i]))
    return gt

# plt.step(snr.index, snr[0], zorder=1, where='mid')
# y = [i-(abs(ii)) for i,ii in zip(snr[0], dsnr[0]) ]
# plt.axhline(y=1.5, color='red')
# 

fig, axs = plt.subplots(math.ceil(len(snr.columns)))
for keys, values in EE_finder(snr,dsnr).items():
    merged = pd.concat([snr2[keys], snr[keys]])
    axs[keys].step(merged.index, merged, where='mid')
    # 1.5 line
    axs[keys].axhline(y=1.5, color='red')

    # snr values that greater than 1.5 and their min max
    for i in range(len(values)):
        xx = max(values[i]) + 2
        y = min(values[i]) - 2
        axs[keys].axvline(x=xx, color='red')
        axs[keys].axvline(x=y, color='red')
        # EE filled 
        axs[keys].axvspan(xx, y, alpha=0.1, color='red', hatch='/')

    xx = [i-2 for i in snr.index]
    # axs[keys].scatter(xx, (snr[keys]-(2*abs(dsnr[keys]))), color='red')
    # axs[keys].scatter(snr.index, snr[keys])

    # error for -100 5
    yy = [i-(abs(ii)) for i, ii in zip(snr2[keys], dsnr2[keys])]
    xx = [i-0.5 for i in snr2.index]
    axs[keys].errorbar(snr.index, snr[keys], yerr=dsnr[keys], ls='none', zorder=2, ecolor='black')

    # error for 5 350
    yy = [i-(abs(ii)) for i, ii in zip(snr[keys], dsnr[keys])][1:]
    xx = [i-2 for i in snr.index][1:]
    axs[keys].errorbar(snr.index, snr[keys], yerr=dsnr[keys], ls='none', zorder=2, ecolor='black')

    
    # axs[keys].legend()
print(EE_finder(snr,dsnr))
# plt.show()