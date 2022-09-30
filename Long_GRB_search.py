import pandas as pd
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 180

data = pd.read_csv('/Volumes/GoogleDrive/My Drive/Python/Fermi/swift/heasoft/grb_table.txt', sep='\t')
os.chdir('/Users/mustafagumustas/Downloads/Swift_BAT/bat_data')
results = pd.read_csv('/Users/mustafagumustas/Downloads/Swift_BAT/sample_list.csv', sep=';')
last_results = results
results['T_90'] = results['T_90'].astype(float)
grbs = [i for i in os.listdir() if i != '.DS_Store' and 'pdf' not in i]

Longs = []

for i in grbs:
    file_name = f'/Users/mustafagumustas/Downloads/Swift_BAT/bat_data/{i}/LC/msec64.lc'
    hdul = fits.open(file_name)

    # TIME DATA
    time = pd.DataFrame([i[0] for i in hdul[1].data])
    TRIGTIME = hdul[1].header['TRIGTIME']
    time = time - TRIGTIME

    # COUNT DATA
    count = pd.DataFrame([i[1] for i in hdul[1].data], index=time[0])
    NGOODPIX = hdul[1].header['NGOODPIX']
    count = count * NGOODPIX

    # ERROR DATA
    error = pd.DataFrame([i[2] for i in hdul[1].data], index=time[0])
    error = error * NGOODPIX


    # between -1 and 5th sec
    sec5_count = count[(-1 < count.index) & (count.index < 5)]
    max_c = max(count[0])   # max photon count value
    max_i = count.loc[count[0] == max_c].index[0]
    
    results['LongCriteria'] = False

    # FILTER OUT SHORT GRBS
    results = results[results['T_90']>5]

    # FIRST MORPHOLOGICAL CRITERIA
    if -1 <max_i and max_i < 5:

        # SECOND MORPHOLOGICAL CRITERIA
        # if count rate is below 11k look under %40 else %30
        thrty_p = max_c * 0.4 if max_c < 11000 else max_c * 0.3

        # we need to look after the peak point, data elimination
        after_max = sec5_count[max_i < sec5_count.index]

        # data elimination, only taking count values under %30/40
        below_thrty = after_max[after_max.values < thrty_p]
        bel_c, col = (below_thrty.shape)

        # some GRBs not drops under %30
        if bel_c > 0:
            # vertical bar that shows the first value that under %30/40 of max count rate
            thry_bar = min(below_thrty.index)

            # data after vertical bar
            # we need total number after that bar and number that below horizontal bar
            # if %50 of them under h bar its good to go
            after_thrty = after_max[after_max.index > thry_bar]

            if len(below_thrty) > (len(after_thrty)/2):

                Longs.append(i)

                f, (ax1, ax2) = plt.subplots(2,1 ,sharex=False, sharey=False)
                count = count[(-50 < count.index) & (count.index < 350)] 
                ax1.set_title(f'{i[:3].upper()} {i[3:]}', loc='center')
                ax1.step(count[0].index, count[0].values)
                ax1.set_ylim([0,None])
                ax1.axvline(x = 5, color='b', linestyle='--', label='5.th sec')
                ax1.set_ylabel('Count rate/sec')

                # setting T90 value on plt
                T90 = results['T_90'].loc[results['GRB Name'] == i[3:]].values
                ax1.set_title(f'T_90: {T90}', loc='right')

                # ax2 is zoomed one
                count = count[(-6 < count.index) & (count.index < 10)] 
                error = error[(-6 < error.index) & (error.index < 10)] 
                error = error/2
                ax2.step(count[0].index, count[0].values, lw=0.5)
                ax2.set_ylim([0,None])
                ax2.axhline(y = thrty_p, color='r', label='%30 count')
                ax2.axhline(y = 0, color='r', linestyle='--', label='background level')
                ax2.axvline(x = 5, color='b', linestyle='--', label='5.th sec')
                ax2.axvline(x = thry_bar, color='g', linestyle='--', label='%30 index')
                ax2.set_ylabel('Count rate/sec')
                ax2.set_xlabel('Time since the trigger (sec)')
                xdata = count.index
                xdata = xdata - 0.03
                ax2.errorbar(xdata, count.values, ls='none', yerr=error[0], 
                ecolor='black', elinewidth=0.5)

                # plt.savefig(f'/nextq/mustafa/bat_data/{grb}/LC/{file_name}.pdf')
                # plt.savefig(f'/Users/mustafagumustas/Downloads/{grb[:3].upper()}{grb[3:]}_{file_name}.pdf')
    else:
        continue
for i in results['GRB Name']:
    if f'grb{i}' in Longs:
        last_results.loc[last_results['GRB Name'] == i, 'LongCriteria'] = True

last_results.to_csv('/Users/mustafagumustas/Downloads/Swift_BAT/sample_list.csv', sep=';', index=False)
print(last_results)
plt.show()
