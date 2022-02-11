from isort import place_module
import pandas as pd
from astropy.io import fits
import numpy as np
import matplotlib.pyplot as plt
hdul = fits.open('onesec_4ch.lc')
#info = (hdul.info())

# loop both rate and count they r 4d, so you must create seperate dfs then merge them

e_min = (pd.DataFrame(hdul[2].data))['E_MIN'].to_list()
e_min = ([int(i) for i in e_min])
e_max = (pd.DataFrame(hdul[2].data))['E_MAX'].to_list()
e_max = ([int(i) for i in e_max])
column = [f'{e_min[i]}-{e_max[i]}'for i in range(len(e_min))]

count = pd.DataFrame([hdul[1].data[i][1] for i in range(len(hdul[1].data))], columns=column)

error = pd.DataFrame([hdul[1].data[i][2] for i in range(len(hdul[1].data))])

time = pd.DataFrame([hdul[1].data[i][0] for i in range(len(hdul[1].data))])

NGOODPIX = hdul[1].header['NGOODPIX']
TRIGTIME = hdul[1].header['TRIGTIME']

count = count * NGOODPIX
time = time - TRIGTIME
# deneme

max_c = (max(count['15-25'])*0.30)
f, (ax1, ax2) = plt.subplots(2, 1, sharex=False, sharey=True)

ax1.plot(time[0], count['15-25'])
# ax1.axhline(y = max_c, color='r', label='max count')
# ax1.axhline(y = 0, color='r', linestyle='--', label='background level')
# ax1.axvline(x = 5, color='b', linestyle='--', label='5.th sec')
# ax1.legend(loc="upper right")
ax1.set_xlim([-100,350])


ax2.plot(time[0], count['15-25'])
ax2.axhline(y = max_c, color='r', label='max count')
ax2.axhline(y = 0, color='r', linestyle='--', label='background level')
ax2.axvline(x = 5, color='b', linestyle='--', label='5.th sec')
ax2.legend(loc="upper right")
ax2.set_xlim([-20,60])
plt.show()

# No.    Name      Ver    Type      Cards   Dimensions   Format
#   0  PRIMARY       1 PrimaryHDU     108   ()      
#   1  RATE          1 BinTableHDU    306   1202R x 5C   [D, 4D, 4D, J, D]   
#   2  EBOUNDS       1 BinTableHDU    288   4R x 3C   [I, E, E]   
#   3  STDGTI        1 BinTableHDU    144   1R x 2C   [D, D] 



# Filename: onesec.lc
# No.    Name      Ver    Type      Cards   Dimensions   Format
#   0  PRIMARY       1 PrimaryHDU     108   ()      
#   1  RATE          1 BinTableHDU    299   1202R x 5C   [D, D, D, J, D]   
#   2  EBOUNDS       1 BinTableHDU    288   1R x 3C   [I, E, E]   
#   3  STDGTI        1 BinTableHDU    144   1R x 2C   [D, D]  