#%%
from astropy.io import fits         #reading fits files
import matplotlib.pyplot as plt     #plotting the graph
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

#hdul = fits.open('fermi/T90/fermi_all_error.fits')
hdul = fits.open('fermi_all_error.fits')
data = pd.DataFrame(hdul[1].data)

# bin degerlerini olusturma
bin = np.array(np.logspace(-2,3, num= 50), dtype=np.float64)
data = np.array(data['T90'].to_list(), dtype=np.float64)

plt.figure()

freq= plt.hist(data, bins=bin, histtype='step', zorder=1)
plt.xscale('log')
ydata = np.append(freq[0], 0)
xdata = freq[1]

def gauss(x, H, A, x0, sigma):
    return H + A * np.exp(-(x - x0) ** 2 / (2 * sigma ** 2))

def gauss_fit(x, y):
    mean = sum(x * y) / sum(y)
    sigma = np.sqrt(sum(y * (x - mean) ** 2) / sum(y))
    popt, pcov = curve_fit(gauss, x, y, maxfev=10**6,p0=[min(y), max(y), mean, sigma])
    return popt

# for error plot
bincenters = 0.5*(xdata[1:]+xdata[:-1])
menStd     = np.sqrt(ydata)
width      = 0.001
bincenters = np.append(bincenters, 0) # len of these data do not match, so I add zero to and
x =18
gauss_1 = gauss(xdata[:x], *gauss_fit(xdata[:x], ydata[:x]))
plt.plot(xdata[:x], gauss_1, 'r', label='fit', zorder=3)

gauss_2 = gauss(xdata[x-1:x+14], *gauss_fit(xdata[x-1:x+14], ydata[x-1:x+14]))
gauss_3 = gauss(xdata[x+13:x+18], *gauss_fit(xdata[x+13:x+18], ydata[x+13:x+18]))
gauss_4 = gauss(xdata[x+16:], *gauss_fit(xdata[x+16:], ydata[x+16:]))
gauss_1[-1] = gauss_2[0]
gauss_2[-1] = gauss_3[0]
gauss_3[-1] = gauss_4[0]

# in order to take symmetri of plots i have to do these:
#       add x to index and start to plot inverse y

i_ydata_1 = ydata[:x:]

i_gauss_1 = gauss(xdata[x-1:x+17], *gauss_fit(xdata[x-1:x+17], i_ydata_1))
plt.plot(xdata[x-1:x+17], i_gauss_1[::-1], 'black', linestyle='--', label='i_fit', zorder=3)

i_ydata_2 = ydata[x+16:]
i_xdata_2 = xdata[x+16::-1]
i_xdata_2 = i_xdata_2[:len(i_ydata_2)]

plt.plot(i_xdata_2, gauss_4[::], 'black', linestyle='--', label='fit', zorder=5)

# plt.plot(xdata[x-1:x+14], gauss_2, 'r', label='fit', zorder=4)
# plt.plot(xdata[x+13:x+18], gauss_3, 'r', label='fit', zorder=5)
plt.plot(xdata[x+16:], gauss_4, 'r', label='fit', zorder=5)

plt.bar(x=bincenters, height = ydata, width=width, color='w', yerr=menStd, zorder=2)
# color is white because i coulndt find a way with out plotting data
# plt.legend()
# plt.gca().set_xscale("log")
plt.xlabel('T90 sec')
plt.ylabel('Number of Burst')
plt.show()
#%%