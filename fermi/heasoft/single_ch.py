#%%
import pandas as pd
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt

hdul = fits.open('msec64.lc')

# need to create energy ranges to give column names
e_min = (pd.DataFrame(hdul[2].data))['E_MIN'].to_list()
e_min = ([int(i) for i in e_min])
e_max = (pd.DataFrame(hdul[2].data))['E_MAX'].to_list()
e_max = ([int(i) for i in e_max])
column = [f'{e_min[i]}-{e_max[i]}'for i in range(len(e_min))]

# time data
time = pd.DataFrame([hdul[1].data[i][0] for i in range(len(hdul[1].data))])
TRIGTIME = hdul[1].header['TRIGTIME']
# we have to substracy trig time from all time values to get 0 at real trigger time
time = time - TRIGTIME

# count data
count = pd.DataFrame([hdul[1].data[i][1] for i in range(len(hdul[1].data))], columns=column, index=time[0])
NGOODPIX = hdul[1].header['NGOODPIX']
# we have to multiply count values with NGOODPIX to get real count numbers
count = count * NGOODPIX

# error = pd.DataFrame([hdul[1].data[i][2] for i in range(len(hdul[1].data))])
max_c = max(count['15-150'])
sec5_count = count[(-0.1 < count.index) & (count.index < 5)]
thrty_p = max(sec5_count['15-150'])* 0.3
ind = count.index
for i, photon in (zip(sec5_count.index, sec5_count['15-150'])):
    if photon > thrty_p:
        print(i, photon , max_c, thrty_p)

#%%

###################################################################################################
# plotting begins

f, (ax1, ax2) = plt.subplots(2,1 ,sharex=False, sharey=False)
ax1.plot(count['15-150'])
ax1.set_xlim([-100,350])
ax1.axvline(x = 5, color='b', linestyle='--', label='5.th sec')
ax1.set_ylabel('Count rate/sec')

# SECOND MORPHOLOGICAL CRITERIA
# if count rate is below 11k look under %40 else %30
thrty_p = max_c * 0.4 if max_c < 11000 else max_c * 0.3


# sec5_count = count[(-0.5 < count.index) & (count.index < 6)]


# ax2 is zoomed one
ax2.plot(count['15-150'])
ax2.set_xlim([-10,30])
ax2.axhline(y = thrty_p, color='r', label='%30 count')
ax2.axhline(y = 0, color='r', linestyle='--', label='background level')
ax2.axvline(x = 5, color='b', linestyle='--', label='5.th sec')
params = count['15-150'].iloc[(count['15-150']-thrty_p).abs().argsort()[:1]]
params = params[(params.index>0)]
ax2.axvline(x = params.index[0], color='g', linestyle='--', label='%30 index')
ax2.legend(loc="upper right")
ax2.set_ylabel('Count rate/sec')
ax2.set_xlabel('Time since the trigger (sec)')


plt.plot()
plt.show()

#%%
