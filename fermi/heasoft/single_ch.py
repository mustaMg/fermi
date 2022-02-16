#%%
import pandas as pd
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 200

hdul = fits.open('msec64.lc')

# need to create energy ranges to give column names
e_min = (pd.DataFrame(hdul[2].data))['E_MIN'].to_list()
e_min = ([int(i) for i in e_min])
e_max = (pd.DataFrame(hdul[2].data))['E_MAX'].to_list()
e_max = ([int(i) for i in e_max])
column = [f'{e_min[i]}-{e_max[i]}'for i in range(len(e_min))]

# TIME DATA
time = pd.DataFrame([hdul[1].data[i][0] for i in range(len(hdul[1].data))])
TRIGTIME = hdul[1].header['TRIGTIME']
# we have to substracy trig time from all time values to get 0 at
# real trigger time
time = time - TRIGTIME

# COUNT DATA
count = pd.DataFrame([hdul[1].data[i][1] for i in range(len(hdul[1].data))],\
    columns=column, index=time[0])
NGOODPIX = hdul[1].header['NGOODPIX']
# we have to multiply count values with NGOODPIX to get real count numbers
count = count * NGOODPIX

# ERROR DATA
# error = pd.DataFrame([hdul[1].data[i][2] for i in range(len(hdul[1].data))])

# between -1 and 5th sec
sec5_count = count[(-1 < count.index) & (count.index < 5)] 
max_c = max(sec5_count['15-150'])   # max photon count value

max_i = sec5_count.loc[sec5_count['15-150'] == max_c].index[0]

# FIRST MORPHOLOGICAL CRITERIA
# if the max count rate after the fifth second thats a long GRB
while True:
    if max_i < 5:
        break 
    else:
        exit()


# SECOND MORPHOLOGICAL CRITERIA
# if count rate is below 11k look under %40 else %30
thrty_p = max_c * 0.4 if max_c < 11000 else max_c * 0.3

#  data set tht ones after %30 count values
prcnt30_count = sec5_count[max_i < sec5_count.index]

params = prcnt30_count['15-150']-thrty_p
params = (params.abs()).sort_values()

################################################################################
# plotting begins

f, (ax1, ax2) = plt.subplots(2,1 ,sharex=False, sharey=False)
ax1.plot(count['15-150'])
ax1.set_xlim([-100,350])
ax1.axvline(x = 5, color='b', linestyle='--', label='5.th sec')
ax1.set_ylabel('Count rate/sec')

# ax2 is zoomed one
ax2.plot(count['15-150'])
ax2.set_xlim([-5,10])
ax2.axhline(y = thrty_p, color='r', label='%30 count')
ax2.axhline(y = 0, color='r', linestyle='--', label='background level')
ax2.axvline(x = 5, color='b', linestyle='--', label='5.th sec')
ax2.axvline(x = params.index[0], color='g', linestyle='--', label='%30 index')
ax2.set_ylabel('Count rate/sec')
ax2.set_xlabel('Time since the trigger (sec)')


plt.show()
#%%
