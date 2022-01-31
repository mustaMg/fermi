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
    popt, pcov = curve_fit(gauss, x, y, p0=[min(y), max(y), mean, sigma])
    return popt

# for error plot
bincenters = 0.5*(xdata[1:]+xdata[:-1])
menStd     = np.sqrt(ydata)
width      = 0.001
bincenters = np.append(bincenters, 0) # len of these data do not match, so I add zero to and

gauss_1 = gauss(xdata[:25], *gauss_fit(xdata[:25], ydata[:25]))
gauss_2 = gauss(xdata, *gauss_fit(xdata, ydata))[24:]
gauss_1[-1] = gauss_2[0]

plt.plot(xdata[:25], gauss_1, 'r', label='fit', zorder=3)
plt.plot(xdata[24:], gauss_2, 'r', label='fit', zorder=4)

plt.bar(x=bincenters, height = ydata, width=width, color='w', yerr=menStd, zorder=2)
# color is white because i coulndt find a way with out plotting data
# plt.legend()
# plt.gca().set_xscale("log")
plt.xlabel('T90 sec')
plt.ylabel('Number of Burst')
plt.show()