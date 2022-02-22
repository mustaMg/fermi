#%%
from astropy.io import fits
from numpy import sqrt
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
hdul = fits.open('/Volumes/GoogleDrive/My Drive/Python/Fermi/swift/heasoft/data/onesec_4ch.lc')

e_min = (pd.DataFrame(hdul[2].data))['E_MIN'].to_list()
e_min = ([int(i) for i in e_min])
e_max = (pd.DataFrame(hdul[2].data))['E_MAX'].to_list()
e_max = ([int(i) for i in e_max])
column = [f'{e_min[i]}-{e_max[i]}'for i in range(len(e_min))]

NGOODPIX = hdul[1].header['NGOODPIX']

rate = pd.DataFrame(columns=column)     # count rate 
for i, ii in enumerate(column):
    rate[str(ii)] = [x[1][i]*NGOODPIX for x in hdul[1].data]
rate['15-50'] = rate['15-25'] + rate['25-50']

time = pd.DataFrame([hdul[1].data[i][0] for i in range(len(hdul[1].data))])
TRIGTIME = hdul[1].header['TRIGTIME']
time = time - TRIGTIME
rate.index = time[0]

rate = rate[(-100 < rate.index) & (rate.index < 350)]

error = pd.DataFrame([hdul[1].data[i][2] for i in range(len(hdul[1].data))],\
    index=time[0], columns=column)
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
    bg[str(ii)] = [x[1][i]*NGOODPIX for x in hhdul[1].data]
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


# in order to bin the data we need a little modification to index
# because python doesn't have a function to that we use groupby and sum them
# it takes the index and look how many n in that index and groups the similars
# then sums in one row. 
   
rate2 = rate[(5 < rate.index) & (rate.index < 350)]
rate2.index = rate2.index - 1   # we have to do this and use 4th s later filter
binned_rate2  = (rate2.groupby(rate2.index // 4).sum())

error2 = error[(5 < error.index) & (error.index < 350)]
error2.index = error2.index - 1
bin_error2 = (error2.groupby(error2.index // 4).sum())

bbg = bg[(5 < bg.index) & (bg.index < 350)]
bbg.index = bbg.index - 1
bin_bbg = (bbg.groupby(bbg.index // 4).sum())

dbbg = dbg[(5 < dbg.index) & (dbg.index < 350)]
dbbg.index = dbbg.index - 1
bin_dbbg = (dbbg.groupby(dbbg.index // 4).sum())

# y0    count
# dy0   error 
# bg0    unweighted
# dbg0  unweighted error
# yy0   count, time between 5-350
# dyy0  error, time between 5-350
# bbg0  unweighted, time between 5-350
# dbbg0 unweighted error, time between 5-350
# ryy0  rebin of count 
# rbg0  rebin of unweighted

###
# there is a problem here
# this works great but the first elemnt is not binned with 4
 
####################
# SIGNAL TO NOISE

def SNR(live, phot, dphot, cr, dcr):
    # live = 4
    # phot = bin_bbg
    # dphot = bin_dbbg
    # cr = binned_rate2
    # dcr = bin_error2
    counts = cr * live
    dcounts = dcr * live

    sigma = sqrt( bin_bbg * live ) 
    dsigma = ( 1/2 * (dphot/phot) ) * sigma

    snr = counts / sigma
    dsnr = sqrt( (dcounts / counts)**2 + (dsigma / sigma)**2 ) * snr
    return snr, dsnr

# SIGNAL TO NOISE   To-100 : To+5
snr, dsnr= SNR(4, bin_bbg, bin_dbbg, binned_rate2, bin_error2)  

# SIGNAL TO NOISE   To-100 : To+5s
rate3 = rate[rate.index < 5]
error3 = error[error.index < 5]

bbg2 = bg[bg.index < 5]
dbbg2 = dbg[dbg.index < 5]

snr2, dsnr2 = SNR(1, bbg2, dbbg2, rate3, error3)

#   ERROR   #
# every values is nan because of division, counts is ok but sigma is not


live = 1
phot = bbg2
dphot = dbbg2
cr = rate3
dcr = error3

counts = cr * live
dcounts = dcr * live
sigma = sqrt( bin_bbg * live ) 
dsigma = ( 1/2 * (dphot/phot) ) * sigma

snr = counts / sigma
dsnr = sqrt( (dcounts / counts)**2 + (dsigma / sigma)**2 ) * snr
a = (((1/2) *(dphot/phot)))
counts = counts[counts.index > -96]
print(sqrt( bin_bbg * live ) )
print(( 1/2 * (dphot/phot) ))
print(sqrt( bin_bbg * live ))
print(12976.095561/-304.859408)
print(0.000098 * 12976.095561)

# ------------------ #
# snr22, dsnr22 = 







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