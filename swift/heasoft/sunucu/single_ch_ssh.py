#%%
import pandas as pd
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 180

data = pd.read_csv('/Volumes/GoogleDrive/My Drive/Python/Fermi/swift/heasoft/grb_table3.txt', sep='\t')
os.chdir('/Users/mustafagumustas/Downloads/Swift_BAT/bat_data')
results = pd.read_csv('/Users/mustafagumustas/Downloads/Swift_BAT/sample_list.csv', sep=';')


# to get T90 values we nee to read a catalogue
# data = pd.read_csv('/nextq/mustafa/grb_table.txt', sep='\t')
# os.chdir('nextq/mustafa/bat_data')
# results = pd.read_excel("/nextq/mustafa/sample_list.xlsx")

grbs = [i for i in os.listdir() if i != '.DS_Store']

def single_ch(grb, file_name):
    hdul = fits.open(f'{file_name}.lc')

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

    # count data
    count = pd.DataFrame([hdul[1].data[i][1] for i in range(len(hdul[1].data))],\
        columns=column, index=time[0])
    NGOODPIX = hdul[1].header['NGOODPIX']
    # we have to multiply count values with NGOODPIX to get real count numbers
    count = count * NGOODPIX

    # ERROR DATA
    error = pd.DataFrame([hdul[1].data[i][2] for i in range(len(hdul[1].data))]\
    ,index=time[0])
    error = error * NGOODPIX
    # between -1 and 5th sec
    sec5_count = count[(-1 < count.index) & (count.index < 5)]
    max_c = max(count['15-150'])   # max photon count value
    max_i = count.loc[count['15-150'] == max_c].index[0]

    # FIRST MORPHOLOGICAL CRITERIA
    # if the max count rate after the fifth second thats a long GRB
    T90 = data['T90'].loc[data.GRB == grb[3:]].values[0]
    if float(T90) > 5:
        if -1 <max_i and max_i < 5:
            # SECOND MORPHOLOGICAL CRITERIA
            # if count rate is below 11k look under %40 else %30
            thrty_p = max_c * 0.4 if max_c < 11000 else max_c * 0.3

            # we need to look after the peak point, data elimination
            after_max = sec5_count[max_i < sec5_count.index]

            # data elimination, only taking count values under %30/40
            below_thrty = after_max[after_max.values < thrty_p]

            # vertical bar that shows the first value that under %30/40 of max count rate
            thry_bar = min(below_thrty.index)

            # data after vertical bar
            # we need total number after that bar and number that below horizontal bar
            # if %50 of them under h bar its good to go
            after_thrty = after_max[after_max.index > thry_bar]
            

            if len(below_thrty) > (len(after_thrty)/2):
                # results.loc[results["GRB Name"] == grb, "LongCriteria"] = 'True'
                f, (ax1, ax2) = plt.subplots(2,1 ,sharex=False, sharey=False)
                count = count[(-50 < count.index) & (count.index < 350)] 
                ax1.set_title(f'{grb[:3].upper()} {grb[3:]}', loc='center')
                ax1.step(count['15-150'].index, count['15-150'].values)
                ax1.set_ylim([0,None])
                ax1.axvline(x = 5, color='b', linestyle='--', label='5.th sec')
                ax1.set_ylabel('Count rate/sec')

                # setting T90 value on plt
                ax1.set_title(f'T_90: {T90}', loc='right')

                # ax2 is zoomed one
                count = count[(-6 < count.index) & (count.index < 10)] 
                error = error[(-6 < error.index) & (error.index < 10)] 
                ax2.step(count['15-150'].index, count['15-150'].values, lw=0.5)
                ax2.set_ylim([0,None])
                ax2.axhline(y = thrty_p, color='r', label='%30 count')
                ax2.axhline(y = 0, color='r', linestyle='--', label='background level')
                ax2.axvline(x = 5, color='b', linestyle='--', label='5.th sec')
                ax2.axvline(x = thry_bar, color='g', linestyle='--', label='%30 index')
                ax2.set_ylabel('Count rate/sec')
                ax2.set_xlabel('Time since the trigger (sec)')
                xdata = count.index
                xdata = xdata - 0.03
                ax2.errorbar(xdata, count.values, ls='none', yerr=error[0], ecolor='black',
                elinewidth=0.5)

                # plt.savefig(f'/nextq/mustafa/bat_data/{grb}/LC/{file_name}.pdf')
                # plt.savefig(f'/Users/mustafagumustas/Downloads/{grb[:3].upper()}{grb[3:]}_{file_name}.pdf')
    else:
        results.loc[results["GRB Name"] == grb, "LongCriteria"] = 'False'


def lc_picker(grb):
    # gets the path to lc files
    given_path = ('/Users/mustafagumustas/Downloads/Swift_BAT/bat_data')
    # given_path = ('nextq/mustafa/bat_data')
    os.chdir(f'{given_path}/{grb}')
    try:
        # os.chdir(f'{given_path}/{grb}')
        file = os.listdir()
        # There might be errors in future, diff os systems creates diff names!!!
        if 'LC' in file:
            os.chdir(f'{given_path}/{grb}/LC')
            for i in os.listdir():
                if i == 'msec64.lc':
                    file_name = i.split('.')[0]
                    single_ch(grb, file_name)
                else:
                    continue
            return 0
    except:
        # print(f'Error, {given_path}/{grb}/{file}/bat/event check if this dir exists')
        return 1

print([lc_picker(grb) for grb in grbs if 'pdf' not in grb])
print(results)
plt.show()
#%%