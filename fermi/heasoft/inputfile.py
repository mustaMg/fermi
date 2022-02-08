# in this file we are gonna use 10 different GRB events files
# and heasoft to create some LC etc.
import os
import pandas as pd

def get_path():
    try: 
        ask_path = input('Enter a path: ')
        #given_path = '/Users/mustafagumustas/Downloads/Swift_BAT/bat_data'
    except:
        print('Path error!')
    return ask_path

if __name__ == '__main__':
    # os.chdir("/Users/mustafagumustas/Downloads/Swift_BAT/bat_data")
    # gonan use below later as main
    given_path = get_path()
    os.chdir(given_path)
    
# creating list of names 
grbs = [i for i in os.listdir() if i != '.DS_Store']

def get_file(grb):
    try:
        os.chdir(given_path+f'/{grb}')
        file = os.listdir()
        # There might be errors in future, diff os systems creates diff names!!!
        file = [i for i in file if i != '.DS_Store'][0]
        os.chdir(given_path+f'/{grb}/'+file+'/bat/event')
        event_files = os.listdir()
        infile = [i for i in event_files if 'bevshsp' in i][0]
        detmask = [i for i in event_files if 'bevtr' in i][0] 
        return infile, detmask
    except:
        print(f'Error, {given_path}+/{grb}/+{file}/bat/event check if this dir exists')
        return 1

for grb in grbs:
    print(get_file(grb))

# /Swift_BAT/bat_data/grb200303A/00959431000/bat/event
# /Users/mustafagumustas/Downloads/Swift_BAT/bat_data