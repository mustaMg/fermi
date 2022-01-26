#%%
#https://stackoverflow.com/questions/48874733/fit-triple-gauss-to-data-python

# https://stackoverflow.com/questions/15556930/turn-scatter-data-into-binned-data-with-errors-bars-equal-to-standard-deviation
# yukaridaki linkte binlenmis datanin sigmasini almayi anlatiyor
from astropy.io import fits         #reading fits files
import matplotlib.pyplot as plt     #plotting the graph
import numpy as np
import pandas as pd
from scipy.optimize import curve_fit

#hdul = fits.open('fermi_all.fits')    
hdul = fits.open('fermi_all_error.fits')
# our interested data is in table 1
data = pd.DataFrame(hdul[1].data)

# we are interested only with t90's
#data = (data['T90'])
bin = np.array(np.logspace(-2,3, num= 50), dtype=np.float64)
x = np.array(data['T90'].to_list(), dtype=np.float64)
error = np.array(data['T90_ERROR'].to_list(), dtype=np.float64)
frequancy = ([len(data[data['T90'].between(bin[i-1], bin[i])]) for i in range(0, len(bin-1))])
# frequancy plot ederken kullanilmiyor, amaci gauss fit cizerken sigma verisini elde etmek

#   frequancy = ([len(data['T90'][(data['T90'] > bin[i]) & (data['T90'] < bin[i+1])]) for i in range(0, len(bin)-1)])
# GAUSS FIT
import pylab as plb
from math import e
from numpy import exp

n = len(bin)                          #the number of data
mean = sum(bin*frequancy)/n                   #note this correction
sigma = sum(frequancy*(bin-mean)**2)/n        #note this correction

def gaus(x,a,x0,sigma):
    return a*exp(-(x-x0)**2/(2*sigma**2))

popt,pcov = curve_fit(gaus,frequancy,bin,p0=[1,mean,sigma])
# burada plt.hist degistirldi, hesaplamalari kolaylastirmak icin
plt.figure()
plt.hist(x, bins=bin, histtype='step')
#plt.plot(bin,gaus(frequancy,*popt),'ro:',label='fit')

# for error plot
y,binEdges = np.histogram(data['T90'],bins=bin)
bincenters = 0.5*(binEdges[1:]+binEdges[:-1])
menStd     = np.sqrt(y)
width      = 0.008
plt.bar(x=bincenters, height = y, width=width, color='w', yerr=menStd)
# color is white because i coulndt find a way with out plotting data
#plt.legend()
plt.gca().set_xscale("log")
plt.xlabel("T90 sec.")
plt.ylabel("Number of Burst")
plt.show()
#%%

# ok i understand it now 
# I did it xd