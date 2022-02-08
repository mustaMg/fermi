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

def onesec_4ch(grb, infile, detmask):
    os.system(f'batbinevt infile={infile} outfile={given_path}+/{grb}/LC/onesec_4ch.lc outtype=LC timedel=1.0\
    timebinalg=u energybins=15-25,25-50,50-100,100-150 detmask={detmask} clobber=YES')
    return 0

def onesec(grb, infile, detmask):
    os.system(f'batbinevt infile={infile} outfile={given_path}+/{grb}/LC/onesec.lc outtype=LC timedel=1.0\
    timebinalg=u energybins=15-150 detmask={detmask} clobber=YES')
    return 0

def msec64_4ch(grb, infile, detmask):
    os.system(f'batbinevt infile={infile} outfile={given_path}+/{grb}/LC/msec64_4ch.lc outtype=LC timedel=0.064\
    timebinalg=u energybins=15-25,25-50,50-100,100-150 detmask={detmask} clobber=YES')
    return 0

def msec64(grb, infile, detmask):
    os.system(f'batbinevt infile={infile} outfile={given_path}+/{grb}/LC/msec64.lc outtype=LC timedel=0.064\
    timebinalg=u energybins=15-150 detmask={detmask} clobber=YES')
    return 0

def unweighted_onesec(grb, infile):
    os.system(f'batbinevt infile={infile} outfile={given_path}+/{grb}/LC/total.dpi \
        outtype=DPI timedel=1.0 timebinalg=u energybins=15-25,25-50,50-100,100-150 \
        weighted=NO outunits=COUNTS clobber=YES')


def lc_creator(grb):
    try:
        os.chdir(given_path+f'/{grb}')
        file = os.listdir()
        # There might be errors in future, diff os systems creates diff names!!!
        file = [i for i in file if i != '.DS_Store'][0]
        os.chdir(given_path+f'/{grb}/'+file+'/bat/event')
        event_files = os.listdir()
        infile = [i for i in event_files if 'bevshsp' in i][0]
        detmask = [i for i in event_files if 'bevtr' in i][0]
        #onesec_4ch(grb, infile, detmask)
        onesec(grb, infile, detmask)
        #msec64_4ch(grb, infile, detmask)
        #msec64(grb, infile, detmask)
        return infile, detmask

    except:
        print(f'Error, {given_path}+/{grb}/+{file}/bat/event check if this dir exists')
        return 1


names = pd.DataFrame([lc_creator(grb) for grb in grbs])

#names.to_csv('/nextq/mustafa//deneme.txt', sep='\t') # for ssh
#names.to_csv('/Volumes/GoogleDrive/My\ Drive/Python/Fermi/fermi/heasoft/deneme.txt', sep='\t') # local

# /Swift_BAT/bat_data/grb200303A/00959431000/bat/event
# /Users/mustafagumustas/Downloads/Swift_BAT/bat_data