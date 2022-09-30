import math
from astropy.io import fits
from matplotlib import patches
import matplotlib
import pandas as pd 
import numpy as np
import matplotlib.pyplot as plt
matplotlib.rcParams['figure.dpi'] = 180
results = pd.read_csv('/Users/mustafagumustas/Downloads/Swift_BAT/sample_list.csv', sep=';')
for i in results['GRB Name']:
    if results['LongCriteria'][results['GRB Name'] ==i].values:
        # reading weighted data 
        hdul = fits.open(f'/Users/mustafagumustas/Downloads/Swift_BAT/bat_data/grb{i}/LC/onesec_4ch.lc')
        count = pd.DataFrame([i[1] for i in hdul[1].data])
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
        error[4] = ( ( error[0].pow(2) ).add( error[1].pow(2) ) ).pow(1./2) 
        ################################################################################
        #                                                                              #
        #                                UNWEIGHTED                                    #
        #                                                                              #
        ################################################################################

        hdul2 = fits.open(f'/Users/mustafagumustas/Downloads/Swift_BAT/bat_data/grb{i}/LC/unweighted.lc')

        bg = pd.DataFrame([i[1] for i in hdul2[1].data], index=time[0])
        bg = bg[(bg.index > -100) & (bg.index < 350)]
        dbg = bg.pow(1./2)
        bg[4] = bg[0] + bg[1]

        dbg[4] = (dbg[0].add(dbg[1])).div(2)
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
        rcount.index = [i for i in range(5,349, 4)];rerror.index = [i for i in range(5,349, 4)]
        rbg.index    = [i for i in range(5,349, 4)];rdbg.index   = [i for i in range(5,349, 4)]


        ################################################################################
        #                                                                              #
        #                              SIGNAL TO NOISE                                 #
        #                                                                              #
        ################################################################################

        def SNR(live, phot, dphot, cr, dcr):
            counts = cr * live
            dcounts = dcr * live 

            sigma = (phot * live)**(1/2)

            dsigma = ( (dphot/phot) * (1/2) ) * sigma

            snr = counts/sigma
            dsnr = ( ( ((dcounts / counts)**2) + ((dsigma / sigma)**2))**(1/2) ) * snr

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
            df['snr'] = col1; df['dsnr'] = col2/2; df['EE'] = df['snr'] - abs(df['dsnr'])
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

        fig, axs = plt.subplots(math.ceil(len(snr.columns)), sharex=True)
        axs[0].set_title(f'GRB {i}')
        for keys, values in EE_finder(snr,dsnr).items():
            energy_ch = ['15-25', '25-50', '50-100', '100-150', '15-50']
            merged = pd.concat([snr2[keys], snr[keys]])
            axs[keys].step(merged.index, merged, where='mid', linewidth=0.5)
            lbl = patches.Patch(label=f'{energy_ch[keys]}')
            axs[keys].legend(handles=[lbl], prop={'size': 6})
            # 1.5 line
            axs[keys].axhline(y=1.5, color='red', linestyle='--', linewidth=0.5)

            # snr values that greater than 1.5 and their min max
            for i in range(len(values)):
                xx = max(values[i]) + 2
                y = min(values[i]) - 2
                axs[keys].axvline(x=xx, color='red', linestyle='--', linewidth=0.5)
                axs[keys].axvline(x=y, color='red', linestyle='--', linewidth=0.5)
                # EE filled 
                axs[keys].axvspan(xx, y, alpha=0.1, color='red', hatch='/')

            xx = [i-2 for i in snr.index]

            # error for -100 5
            yy = [i-(abs(ii)) for i, ii in zip(snr2[keys], dsnr2[keys])]
            xx = [i-0.5 for i in snr2.index]
            y_er = dsnr2[keys]/2
            axs[keys].errorbar(snr2.index, snr2[keys], yerr=y_er, ls='none', zorder=2, ecolor='black', elinewidth=0.5)

            # error for 5 350
            yy = [i-(abs(ii)) for i, ii in zip(snr[keys], dsnr[keys])][1:]
            y_er = dsnr[keys]/2
            axs[keys].errorbar(snr.index, snr[keys], yerr=y_er, ls='none', zorder=2, ecolor='black', elinewidth=0.5)


plt.legend()
plt.show()
