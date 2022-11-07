#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Nov  4 12:30:08 2022

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
path = "/Users/joel/Research/LCOF Cooling/Data/221006 - THE P-O & P-P DATA/p-o/"
base = "/Users/joel/Research/LCOF Cooling/Python/LCOF-Cooling/Plots/"
save = base + timestamp.strftime("%Y-%b-%d") + "/P-O/" + timestamp.strftime("%X") + ": " + note + "/"
if not os.path.exists(save):
  os.makedirs(save)

powers = ['10', '30', '50', '70', '90', '110', '130', '150', '170', '190',
          '210', '230', '250', '270', '290']
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
            'Sig': df[pow+'rs']['Sig'] - df[pow+'bs']['Sig'],
            'σ': (df[pow+'rs']['σ']**2 + df[pow+'bs']['σ']**2)**.5})

aS = dict(zip(powers, [df['10aS'], df['30aS'], df['50aS'], df['70aS'],
                       df['90aS'], df['110aS'], df['130aS'], df['150aS'],
                       df['170aS'], df['190aS'], df['210aS'], df['230aS'],
                       df['250aS'], df['270aS'], df['290aS']]))
s = dict(zip(powers, [df['10s'], df['30s'], df['50s'], df['70s'],
                      df['90s'], df['110s'], df['130s'], df['150s'],
                      df['170s'], df['190s'], df['210s'], df['230s'],
                      df['250s'], df['270s'], df['290s']]))

def l(x, amp, wid, cen, c):
  #fwhm wid lorentzian
  return amp * wid**2 / (wid**2 + ((x-cen))**2) + c

def lin(x, m, b):
  return m*x+b

# fit
a, w, ce, c = 2, .1, 2.275, 0
guess = a, w, ce, c

# plot anti-Stokes fits w error
plt.figure(dpi=250)
plt.title("Pump-Only anti-Stokes")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Spectral Density (μV)")
plt.xlim(2,2.5)
#plt.ylim()
plt.minorticks_on()
plt.tick_params(which='both', direction='in', pad=5)

print("anti-Stokes")
aSfwhm, aSfwhmσ = [], []
for pow in powers:
  popt, pcov = curveFit(l, aS[pow]['Freq'], aS[pow]['Sig'], guess, sigma=aS[pow]['σ'], absolute_sigma=True)
  amp, wid, cen, c = popt[0], abs(2*popt[1]), popt[2], popt[3]
  aSfwhm.append(wid*1000)
  print(pow)
  print(f"Amp: {amp:.3f} μV \t Wid: {wid*1000:.3f} MHz \t Cen: {cen:.3f} GHz \t\t C: {c:.3f} μV")
  #print(pcov)
  σAmp, σWid, σCen, σC = pcov[0][0]**.5, pcov[1][1]**.5, pcov[2][2]**.5, pcov[3][3]**.5
  aSfwhmσ.append(σWid*1000)
  print(f"σAmp: {σAmp:.4f} μV \t σWid: {σWid*1000: .4f} MHz \t σCen: {σCen: .4f} GHz \t σC: {σC: .4f} μV")
  plt.errorbar(aS[pow]['Freq'], aS[pow]['Sig'], yerr=aS[pow]['σ'], fmt="None", elinewidth=.25, alpha=.5, capsize=1, capthick=.25)
  plt.plot(aS[pow]['Freq'], l(aS[pow]['Freq'], *popt), linewidth=1, label=pow+' mW')

plt.legend(fontsize=7.5)
plt.savefig(save + "P-O anti-Stokes Fits.pdf", format="pdf")
plt.savefig(save + "P-O anti-Stokes Fits.png", format="png")

# plot pow v wid
npPowers = np.arange(10, 291, 20)
popt, pcov = curveFit(lin, npPowers, aSfwhm, [.1, 100], sigma=aSfwhmσ, absolute_sigma=True)
maS, baS = popt[0], popt[1]
print(f"m: {maS:.4f}, b: {baS:.4f}")

plt.figure(dpi=250)
plt.title("Pump-Only anti-Stokes Pow v Wid")
plt.xlabel("Pump Power (mW)")
plt.ylabel("fwhm (MHz)")
plt.xlim(0,300)
#plt.ylim()
plt.minorticks_on()
plt.tick_params(which='both', direction='in', pad=5)

plt.errorbar(npPowers, aSfwhm, yerr=aSfwhmσ, fmt="None", elinewidth=.5, color='Gray', alpha=.5, capsize=1, capthick=.5)
#plt.scatter(npPowers, aSfwhm, 1)
plt.plot(np.array([0,300]), lin(np.array([0,300]), maS, baS), color="Black", linewidth=1)

plt.savefig(save + "P-O anti-Stokes Pow v Wid.pdf", format="pdf")
plt.savefig(save + "P-O anti-Stokes Pow v Wid.png", format="png")


# plot Stokes fits w error
plt.figure(dpi=250)
plt.title("Pump-Only Stokes")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Spectral Density (μV)")
plt.xlim(2,2.5)
#plt.ylim()
plt.minorticks_on()
plt.tick_params(which='both', direction='in', pad=5)

print("Stokes")
sfwhm, sfwhmσ = [], []
for pow in powers:
  popt, pcov = curveFit(l, s[pow]['Freq'], s[pow]['Sig'], guess, sigma=s[pow]['σ'], absolute_sigma=True)
  amp, wid, cen, c = popt[0], abs(2*popt[1]), popt[2], popt[3]
  sfwhm.append(wid*1000)
  print(pow)
  print(f"Amp: {amp:.3f} μV \t Wid: {wid*1000:.3f} MHz \t Cen: {cen:.3f} GHz \t\t C: {c:.3f} μV")
  #print(pcov)
  σAmp, σWid, σCen, σC = pcov[0][0]**.5, pcov[1][1]**.5, pcov[2][2]**.5, pcov[3][3]**.5
  sfwhmσ.append(σWid*1000)
  print(f"σAmp: {σAmp:.4f} μV \t σWid: {σWid*1000: .4f} MHz \t σCen: {σCen: .4f} GHz \t σC: {σC: .4f} μV")
  plt.errorbar(s[pow]['Freq'], s[pow]['Sig'], yerr=s[pow]['σ'], fmt="None", elinewidth=.25, alpha=.5, capsize=1, capthick=.25)
  plt.plot(s[pow]['Freq'], l(s[pow]['Freq'], *popt), linewidth=1, label=pow+' mW')

plt.legend(fontsize=7.5)
plt.savefig(save + "P-O Stokes Fits.pdf", format="pdf")
plt.savefig(save + "P-O Stokes Fits.png", format="png")

# plot pow v wid
npPowers = np.arange(10, 291, 20)
popt, pcov = curveFit(lin, npPowers, sfwhm, [.1, 100], sigma=sfwhmσ, absolute_sigma=True)
ms, bs = popt[0], popt[1]
print(f"m: {ms:.4f}, b: {bs:.4f}")

plt.figure(dpi=250)
plt.title("Pump-Only Stokes Pow v Wid")
plt.xlabel("Pump Power (mW)")
plt.ylabel("fwhm (MHz)")
plt.xlim(0,300)
#plt.ylim()
plt.minorticks_on()
plt.tick_params(which='both', direction='in', pad=5)

plt.errorbar(npPowers, sfwhm, yerr=sfwhmσ, fmt="None", elinewidth=.5, color='Gray', alpha=.5, capsize=1, capthick=.5)
#plt.scatter(npPowers, sfwhm, 1)
plt.plot(np.array([0,300]), lin(np.array([0,300]), ms, bs), color="Black", linewidth=1)

plt.savefig(save + "P-O Stokes Pow v Wid.pdf", format="pdf")
plt.savefig(save + "P-O Stokes Pow v Wid.png", format="png")

# plot linewidths
plt.figure(dpi=250)
plt.title("Pump-Only Linewidths")
plt.xlabel("Pump Power (mW)")
plt.ylabel("fwhm (MHz)")
plt.xlim(0,300)
plt.minorticks_on()
plt.tick_params(which='both', direction='in', pad=5)

plt.errorbar(npPowers, aSfwhm, yerr=aSfwhmσ, fmt="None", elinewidth=.5, color='Blue', alpha=.5, capsize=1, capthick=.5)
plt.plot(np.array([0,300]), lin(np.array([0,300]), maS, baS), color="Blue", linewidth=1, label='anti-Stokes')
plt.errorbar(npPowers, sfwhm, yerr=sfwhmσ, fmt="None", elinewidth=.5, color='Red', alpha=.5, capsize=1, capthick=.5)
plt.plot(np.array([0,300]), lin(np.array([0,300]), ms, bs), color="Red", linewidth=1, label='Stokes')
# P-P m: 0.0912, b: 96.8517
plt.plot(np.array([0,300]), lin(np.array([0,300]), 0.0912, 96.8517), color="Black", linewidth=1, label='Pump-Probe anti-Stokes')
plt.legend()

plt.savefig(save + "P-O Linewidths.pdf", format="pdf")
plt.savefig(save + "P-O Linewidths.png", format="png")
