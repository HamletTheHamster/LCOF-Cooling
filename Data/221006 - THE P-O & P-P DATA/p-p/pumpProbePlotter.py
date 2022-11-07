#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 13:46:13 2022

@author: joel
"""

import pandas as pd
import matplotlib.pyplot as plt
#import numpy as np
from scipy.optimize import curve_fit

# import data
path = "/Users/joel/Research/LCOF Cooling/Data/221006 - THE P-O & P-P DATA/p-p/"

powers = ['0', '55', '110', '165']
file = ['bas', 'ras', 'bs', 'rs']
files = 5

df = {}
for pow in powers:
  for f in file:
    for n in range(files):
      df[pow+f+str(n)] = pd.read_csv(path+pow+"/"+f+str(n)+".csv", skiprows=1)
      df[pow+f+str(n)]['Lin'] = 10**6 * 10**(df[pow+f+str(n)]['Amp']/10)
      
      if n == files-1:
        
        # Average across Lin columns
        df[pow+f] = {}
        df[pow+f]['Sum'] = df[pow+f+'0']['Lin']
        
        for i in range(1,files):
          df[pow+f]['Sum'] = df[pow+f]['Sum'] + df[pow+f+str(i)]['Lin'] # Don't use +=
          
        df[pow+f]['Avg'] = df[pow+f]['Sum'] / float(files)
        
        df[pow+f]['Freq'] = df[pow+f+'0']['Freq']
        df[pow+f]['Sig'] = df[pow+f]['Avg']
        
        # Calculate sigma
        df[pow+f]['Σ(i-avg)^2'] = (df[pow+f+'0']['Lin'] - df[pow+f]['Avg'])**2
        for i in range(1,n):
          df[pow+f]['Σ(i-avg)^2'] += (df[pow+f+str(i)]['Lin'] - df[pow+f]['Avg'])**2
        
        df[pow+f]['σ'] = ((df[pow+f]['Σ(i-avg)^2']/float(files-1))**.5)/files**.5
        
        # Subtract
        if f == 'ras':
          df[pow+'aS'] = pd.DataFrame({
            'Freq': df[pow+f]['Freq']/1e9,
            'Sig': df[pow+'ras']['Sig'] - df[pow+'bas']['Sig'],
            'σ': (df[pow+'ras']['σ']**2 + df[pow+'bas']['σ']**2)**.5})
        if f == 'rs':
          df[pow+'s'] = pd.DataFrame({
            'Freq': df[pow+f]['Freq']/1e9,
            'sig': df[pow+'rs']['Sig'] - df[pow+'rs']['Sig'],
            'σ': (df[pow+'rs']['σ']**2 + df[pow+'bs']['σ']**2)**.5})
      
aS = dict(zip(powers, [df['0aS'], df['55aS'], df['110aS'], df['165aS']]))
s = dict(zip(powers, [df['0s'], df['55s'], df['110s'], df['165s']]))

# fit
a, w, c = 1.75, .1, 2.275

def f(x, amp, wid, cen):
  # fwhm wid lorentzian
  return amp * wid**2 / (wid**2 + ((x-cen))**2)

guess = a, w, c

# plot
plt.title("P-P anti-Stokes")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Spectral Density (uV)")

for pow in powers:
  popt, pcov = curve_fit(f, aS['0']['Freq'], aS[pow]['Sig'], guess, sigma=aS[pow]['σ'], absolute_sigma=True)
  amp, wid, cen = popt[0], 2*popt[1], popt[2]
  print(pow)
  print("Amp: " + str(amp) + " Wid: " + str(wid) + " Cen: " + str(cen))
  print(pcov)
  σAmp, σWid, σCen = pcov[0][0]**.5, pcov[1][1]**.5, pcov[2][2]**.5
  print("σAmp: " + str(σAmp) + " σWid: " + str(σWid) + " σCen: " + str(σCen))
  yfit = f(aS[pow]['Freq'], *popt)
  
  
  plt.scatter(aS[pow]['Freq'], aS[pow]['Sig'], marker=".")
  plt.errorbar(aS[pow]['Freq'], aS[pow]['Sig'], yerr=aS[pow]['σ'])
  
  plt.plot(aS[pow]['Freq'], yfit)