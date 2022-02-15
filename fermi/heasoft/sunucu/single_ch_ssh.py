#%%
import pandas as pd
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
import os
import matplotlib as mpl
mpl.rcParams['figure.dpi'] = 200

def get_path():
    try: 
        ask_path = input('Enter a path: ')
        #given_path = '/Users/mustafagumustas/Downloads/Swift_BAT/bat_data'
        return ask_path
    except:
        print('Path error!')
    
if __name__ == '__main__':
    # os.chdir("/Users/mustafagumustas/Downloads/Swift_BAT/bat_data")
    # gonan use below later as main
    given_path = get_path()
    os.chdir(given_path)

grbs = [i for i in os.listdir() if i != '.DS_Store']

def single_ch(grb, file_name):
    hdul = fits.open(f'{file_name}.lc')

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

    ###################################################################################################
    # plotting begins

    f, (ax1, ax2) = plt.subplots(2,1 ,sharex=False, sharey=False)
    ax1.set_title(f'{grb}')
    ax1.plot(count['15-150'])
    ax1.set_xlim([-100,350])
    ax1.axvline(x = 5, color='b', linestyle='--', label='5.th sec')
    ax1.set_ylabel('Count rate/sec')

    # data between -1 and 5th sec
    sec5_count = count[(-1 < count.index) & (count.index < 5)] 
    max_c = max(sec5_count['15-150'])   # max photon count value

    # SECOND MORPHOLOGICAL CRITERIA
    # if count rate is below 11k look under %40 else %30
    thrty_p = max_c * 0.4 if max_c < 11000 else max_c * 0.3

    max_i = [i for i, ct in zip(sec5_count.index, sec5_count['15-150']) if ct == max_c][0]

    prcnt30_count = sec5_count[max_i < sec5_count.index]
    params = prcnt30_count['15-150']-thrty_p
    params = (params.abs()).sort_values()

    # ax2 is zoomed one
    ax2.plot(count['15-150'])
    ax2.set_xlim([-5,10])
    ax2.axhline(y = thrty_p, color='r')
    ax2.axhline(y = 0, color='r', linestyle='--')
    ax2.axvline(x = 5, color='b', linestyle='--')
    ax2.axvline(x = params.index[0], color='g', linestyle='--')
    ax2.set_ylabel('Count rate/sec')
    ax2.set_xlabel('Time since the trigger (sec)')

    print(f'{given_path}/{grb}/LC/{file_name}.pdf')
    # plt.savefig(f'{given_path}/{grb}/LC/{file_name}.pdf')


def lc_picker(grb):
    # gets the path to lc files
    try:
        os.chdir(f'{given_path}/{grb}')
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
        
        else:
            print('lc file not found!')
            return 1
    except:
        # print(f'Error, {given_path}/{grb}/{file}/bat/event check if this dir exists')
        return 1

print([lc_picker(grb) for grb in grbs])
plt.show()
#%%