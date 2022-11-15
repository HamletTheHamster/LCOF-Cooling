#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct 28 13:46:13 2022

@author: joel
"""

import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import curve_fit as curveFit
import datetime
import os
import sys

timestamp = datetime.datetime.now()

if len(sys.argv) > 0:
    note = sys.argv[1]

# import data
path = "Data/221006 - THE P-O & P-P DATA/p-p/"
save = "Plots/" + timestamp.strftime("%Y-%b-%d") + "/P-P/" + timestamp.strftime("%X") + ": " + note + "/"
if not os.path.exists(save):
  os.makedirs(save)

powers = ['0', '55', '110', '165']
truePowers = []
for pow in powers:
    truePowers.append(float(int(pow) - 10)*.17**.5)

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
a, w, ce, c = 1.75, .1, 2.275, .06

def l(x, amp, wid, cen, c):
  # fwhm wid lorentzian
  return amp * wid**2 / (wid**2 + ((x-cen))**2) + c

def lin(x, m, b):
  return m*x+b

guess = a, w, ce, c

# plot fits w error
plt.figure(dpi=250)
plt.title("Pump-Probe anti-Stokes")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Spectral Density (μV)")
plt.xlim(2,2.5)
plt.ylim(0,2)
plt.minorticks_on()
plt.tick_params(which='both', direction='in', pad=5)
plt.tick_params(which='minor', axis='y', length=0)

paletteDict = {'0': (31/255, 211/255, 172/255), '55': (255/255, 122/255, 180/255), '110': (122/255, 156/255, 255/255), '165': (255/255, 182/255, 110/255)}

fwhm, fwhmσ = [], []
for (pow, truPow) in zip(powers, truePowers):
  popt, pcov = curveFit(l, aS['0']['Freq'], aS[pow]['Sig'], guess, sigma=aS[pow]['σ'], absolute_sigma=True)
  amp, wid, cen, c = popt[0], abs(2*popt[1]), popt[2], popt[3]
  fwhm.append(wid*1000)
  print(pow)
  print(f"Amp: {amp:.3f} μV \t Wid: {wid*1000:.3f} MHz \t Cen: {cen:.3f} GHz \t\t C: {c:.3f} μV")
  #print(pcov)
  σAmp, σWid, σCen, σC = pcov[0][0]**.5, pcov[1][1]**.5, pcov[2][2]**.5, pcov[3][3]**.5
  fwhmσ.append(σWid*1000)
  print(f"σAmp: {σAmp:.4f} μV \t σWid: {σWid*1000: .4f} MHz \t σCen: {σCen: .4f} GHz \t σC: {σC: .4f} μV")
  plt.errorbar(aS[pow]['Freq'], aS[pow]['Sig'], yerr=aS[pow]['σ'], fmt="None", elinewidth=.25, color=paletteDict[pow], alpha=.5, capsize=1, capthick=.25)
  plt.plot(aS[pow]['Freq'], l(aS[pow]['Freq'], *popt), color=paletteDict[pow], linewidth=1, label=f"{truPow: .2f}"+' mW')
  #plt.scatter(aS[pow]['Freq'], aS[pow]['Sig'], .5, marker=".", color=palette[pow], )#label=f"{truPow: .2f}"+' mW')

plt.legend()
plt.savefig(save + "P-P anti-Stokes Fits.pdf", format="pdf")
plt.savefig(save + "P-P anti-Stokes Fits.png", format="png")

# plot pow v wid
popt, pcov = curveFit(lin, truePowers, fwhm, [.1, 96.], sigma=fwhmσ, absolute_sigma=True)
m, b = popt[0], popt[1]
print(f"m: {m:.4f}, b: {b:.4f}")

paletteList = [(31/255, 211/255, 172/255), (255/255, 122/255, 180/255), (122/255, 156/255, 255/255), (255/255, 182/255, 110/255)]

plt.figure(dpi=250)
plt.title("Pump-Probe anti-Stokes Pow v Wid")
plt.xlabel("Pump Power (mW)")
plt.ylabel("fwhm (MHz)")
plt.xlim(-2,72)
#plt.ylim()
plt.minorticks_on()
plt.tick_params(which='both', direction='in', pad=5)

plt.errorbar(truePowers, fwhm, yerr=fwhmσ, fmt="None", elinewidth=.5, color='Gray', alpha=.5, capsize=1, capthick=.5)
plt.scatter(truePowers, fwhm, 5, color=paletteList)
plt.plot(np.array([-2,72]), lin(np.array([-2,72]), m, b), color="Black", linewidth=1)

plt.savefig(save + "P-P anti-Stokes Pow v Wid.pdf", format="pdf")
plt.savefig(save + "P-P anti-Stokes Pow v Wid.png", format="png")
