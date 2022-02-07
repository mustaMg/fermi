# in this file we are gonna use 10 different GRB events files
# and heasoft to create some LC etc.

import os
import pandas as pd

#ask_path = input('grb dosyalarinin bulundugu path i giriniz')
try:    
    given_path = '/Users/mustafagumustas/Downloads/Swift_BAT/bat_data'
except:
    print('Path error!')

if __name__ == '__main__':
    os.chdir("/Users/mustafagumustas/Downloads/Swift_BAT/bat_data")

    # gonan use below lateras main
    #os.chdir(ask_path)

    # os.system("pwd")

# creating list of names 
grb = (os.listdir())


for i, grb in enumerate(grb):
    # i put if here because there are systme related files
    if 'grb' in grb:
        try:
            os.chdir(given_path+f'/{grb}')
            file = os.listdir()
            # There might be errors in future, diff os systems creates diff names!!!
            file = [i for i in file if i != '.DS_Store'][0]
            os.chdir(given_path+f'/{grb}/'+file+'/bat/event')
            event_files = os.listdir()
            infile = [i for i in event_files if 'bevshsp' in i][0]
            detmask = [i for i in event_files if 'bevtr' in i][0] 
        except:
            print(f'Error, {given_path}+/{grb}/+{file}/bat/event check if this dir exists')


# /Swift_BAT/bat_data/grb200303A/00959431000/bat/event