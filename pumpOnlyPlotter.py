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
path = "Data/221006 - THE P-O & P-P DATA/p-o/"
save = "Plots/" + timestamp.strftime("%Y-%b-%d") + "/P-O/" + timestamp.strftime("%X") + ": " + note + "/"
if not os.path.exists(save):
  os.makedirs(save)

powers = ['10', '30', '50', '70', '90', '110', '130', '150', '170', '190', '210', '230', '250', '270', '290']
truePowers = []
for pow in powers:
    truePowers.append(float(pow)*.17**.5)

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
plt.ylim(0,120)
plt.minorticks_on()
plt.tick_params(which='both', direction='in', pad=5)

blueGradient = []
for i in range(len(powers)):
    mod = (185/15)*i/255
    blueGradient.append((0+mod, 0+mod, 1))

#blueGradient.reverse()
blueGrad = dict(zip(powers, blueGradient))

print("anti-Stokes")
aSfwhm, aSfwhmσ = [], []
aSAmp, aSAmpσ = {}, {}
for (pow, truPow) in zip(powers, truePowers):
      popt, pcov = curveFit(l, aS[pow]['Freq'], aS[pow]['Sig'], guess, sigma=aS[pow]['σ'], absolute_sigma=True)
      amp, wid, cen, c = popt[0], abs(2*popt[1]), popt[2], popt[3]
      aSfwhm.append(wid*1000)
      aSAmp[pow] = amp
      print(f"pow: {pow} \t truPow: {truPow: .4f}")
      print(f"Amp: {amp:.3f} μV \t Wid: {wid*1000:.3f} MHz \t Cen: {cen:.3f} GHz \t\t C: {c:.3f} μV")
      #print(pcov)
      σAmp, σWid, σCen, σC = pcov[0][0]**.5, pcov[1][1]**.5, pcov[2][2]**.5, pcov[3][3]**.5
      aSfwhmσ.append(σWid*1000)
      aSAmpσ[pow] = σAmp
      print(f"σAmp: {σAmp:.4f} μV \t σWid: {σWid*1000: .4f} MHz \t σCen: {σCen: .4f} GHz \t σC: {σC: .4f} μV")
      plt.scatter(aS[pow]['Freq'], aS[pow]['Sig'], 1, color=blueGrad[pow], label=f"{truPow: .2f}"+" mW")
      #plt.errorbar(aS[pow]['Freq'], aS[pow]['Sig'], yerr=aS[pow]['σ'], fmt="None", elinewidth=.25, alpha=.5, capsize=1, capthick=.25)
      #plt.plot(aS[pow]['Freq'], l(aS[pow]['Freq'], *popt), color=blueGrad[pow], linewidth=1, label=f"{truPow: .2f}"+' mW')

plt.legend(fontsize=7.5)
plt.savefig(save + "P-O anti-Stokes Fits.pdf", format="pdf")
plt.savefig(save + "P-O anti-Stokes Fits.png", format="png")

# plot pow v wid
popt, pcov = curveFit(lin, truePowers, aSfwhm, [.1, 100], sigma=aSfwhmσ, absolute_sigma=True)
maS, baS = popt[0], popt[1]
print(f"m: {maS:.4f}, b: {baS:.4f}")

plt.figure(dpi=250)
plt.title("Pump-Only anti-Stokes Pow v Wid")
plt.xlabel("Pump Power (mW)")
plt.ylabel("fwhm (MHz)")
plt.xlim(0,125)
#plt.ylim()
plt.minorticks_on()
plt.tick_params(which='both', direction='in', pad=5)

plt.errorbar(truePowers, aSfwhm, yerr=aSfwhmσ, fmt="None", elinewidth=.5, color='Blue', alpha=.5, capsize=1, capthick=.5)
#plt.scatter(truePowers, aSfwhm, 1)
plt.plot(np.array([0,125]), lin(np.array([0,125]), maS, baS), color="Black", linewidth=1)

plt.savefig(save + "P-O anti-Stokes Pow v Wid.pdf", format="pdf")
plt.savefig(save + "P-O anti-Stokes Pow v Wid.png", format="png")


# plot Stokes fits w error
plt.figure(dpi=250)
plt.title("Pump-Only Stokes")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Spectral Density (μV)")
plt.xlim(2,2.5)
plt.ylim(0,120)
plt.minorticks_on()
plt.tick_params(which='both', direction='in', pad=5)

redGradient = []
for i in range(len(powers)):
    mod = (185/15)*i/255
    redGradient.append((1, 0+mod, 0+mod))

#redGradient.reverse()
redGrad = dict(zip(powers, redGradient))

print("Stokes")
sfwhm, sfwhmσ = [], []
sAmp, sAmpσ = {}, {}
for (pow, truPow) in zip(powers, truePowers):
  popt, pcov = curveFit(l, s[pow]['Freq'], s[pow]['Sig'], guess, sigma=s[pow]['σ'], absolute_sigma=True)
  amp, wid, cen, c = popt[0], abs(2*popt[1]), popt[2], popt[3]
  sfwhm.append(wid*1000)
  sAmp[pow] = amp
  print(f"pow: {pow} \t truPow: {truPow: .4f}")
  print(f"Amp: {amp:.3f} μV \t Wid: {wid*1000:.3f} MHz \t Cen: {cen:.3f} GHz \t\t C: {c:.3f} μV")
  #print(pcov)
  σAmp, σWid, σCen, σC = pcov[0][0]**.5, pcov[1][1]**.5, pcov[2][2]**.5, pcov[3][3]**.5
  sfwhmσ.append(σWid*1000)
  sAmpσ[pow] = σAmp
  print(f"σAmp: {σAmp:.4f} μV \t σWid: {σWid*1000: .4f} MHz \t σCen: {σCen: .4f} GHz \t σC: {σC: .4f} μV")
  plt.scatter(s[pow]['Freq'], s[pow]['Sig'], 1, color=redGrad[pow], label=f"{truPow: .2f}"+" mW")
  #plt.errorbar(s[pow]['Freq'], s[pow]['Sig'], yerr=s[pow]['σ'], fmt="None", elinewidth=.25, color='Red', alpha=.5, capsize=1, capthick=.25)
  #plt.plot(s[pow]['Freq'], l(s[pow]['Freq'], *popt), linewidth=1, label=f"{truPow: .2f}"+' mW')

plt.legend(fontsize=7.5)
plt.savefig(save + "P-O Stokes Fits.pdf", format="pdf")
plt.savefig(save + "P-O Stokes Fits.png", format="png")


# plot pow v wid
popt, pcov = curveFit(lin, truePowers, sfwhm, [.1, 100], sigma=sfwhmσ, absolute_sigma=True)
ms, bs = popt[0], popt[1]
print(f"m: {ms:.4f}, b: {bs:.4f}")

plt.figure(dpi=250)
plt.title("Pump-Only Stokes Pow v Wid")
plt.xlabel("Pump Power (mW)")
plt.ylabel("fwhm (MHz)")
plt.xlim(0,125)
#plt.ylim()
plt.minorticks_on()
plt.tick_params(which='both', direction='in', pad=5)

plt.errorbar(truePowers, sfwhm, yerr=sfwhmσ, fmt="None", elinewidth=.5, color='Red', alpha=.5, capsize=1, capthick=.5)
#plt.scatter(truePowers, sfwhm, 1)
plt.plot(np.array([0,125]), lin(np.array([0,125]), ms, bs), color="Black", linewidth=1)

plt.savefig(save + "P-O Stokes Pow v Wid.pdf", format="pdf")
plt.savefig(save + "P-O Stokes Pow v Wid.png", format="png")


# plot linewidths
plt.figure(dpi=250)
plt.title("Pump-Only Linewidths")
plt.xlabel("Pump Power (mW)")
plt.ylabel("fwhm (MHz)")
plt.xlim(0,125)
plt.ylim(90,110)
plt.minorticks_on()
plt.tick_params(which='both', direction='in', pad=5)

#deltaT = T_roomTemp*(gamma_eff - gamma)/gamma_eff
roomTemp = 293
gamma = 2*np.pi*baS
gammaEff = 2*np.pi*lin(300, maS, baS)
deltaTaS = 293*(gammaEff - gamma)/gammaEff

gamma = 2*np.pi*bs
gammaEff = 2*np.pi*lin(125, ms, bs)
deltaTS = 293*(gammaEff - gamma)/gammaEff
print(f"degrees cooled: {deltaTaS: .4f} \t\t degrees heated: {deltaTS: .4f}")

tempAxisMax = 293 + 293*(2*np.pi*110 - 2*np.pi*baS)/(2*np.pi*110)
tempAxisMin = 293 + 293*(2*np.pi*90 - 2*np.pi*baS)/(2*np.pi*90)

print(f"tempAxisMin: {tempAxisMin: .4f} \t\t tempAxisMax: {tempAxisMax: .4f}")

plt.scatter(truePowers, aSfwhm, 1, color="blue", label="anti-Stokes")
plt.scatter(truePowers, sfwhm, 1, color="red", label="Stokes")

#plt.errorbar(truePowers, aSfwhm, yerr=aSfwhmσ, fmt="None", elinewidth=.5, color='Blue', alpha=.5, capsize=1, capthick=.5)
plt.plot(np.array([0,125]), lin(np.array([0,125]), maS, baS), color="Blue", linewidth=1, label='anti-Stokes')
#plt.errorbar(truePowers, sfwhm, yerr=sfwhmσ, fmt="None", elinewidth=.5, color='Red', alpha=.5, capsize=1, capthick=.5)
plt.plot(np.array([0,125]), lin(np.array([0,125]), ms, bs), color="Red", linewidth=1, label='Stokes')
#ppm, ppb = 0.2213, 97.7640
#plt.plot(np.array([0,125]), lin(np.array([0,125]), ppm, ppb), color="Black", linewidth=1, label='Pump-Probe anti-Stokes')

# gamma_eff = 2pi*gamma_eff(1 + GPL/4) # G = 3.025
#gammaEff = 98.5985 + 98.5985*3.025*.119315*1.00/4
#plt.plot(np.array([0,125]), np.array([98.5985,gammaEff]), color="Black", linewidth=1, label='find gain')
plt.legend()

ax2 = plt.twinx()
ax2.set_ylabel('Temperature (K)')
ax2.set_ylim(tempAxisMin,tempAxisMax)

plt.savefig(save + "P-O Linewidths.pdf", format="pdf")
plt.savefig(save + "P-O Linewidths.png", format="png")


# plot height ratios
plt.figure(dpi=250)
plt.title("Pump-Only Height Ratios")
plt.xlabel("Pump Power (mW)")
plt.ylabel("Stokes/anti-Stokes Amplitude")
plt.xlim(0,125)
plt.minorticks_on()
plt.tick_params(which='both', direction='in', pad=5)

ratios, ampErr = [], []
for pow in powers:
    ratios.append(sAmp[pow]/aSAmp[pow])
    ampErr.append(((aSAmpσ[pow]/aSAmp[pow])**2 + (sAmpσ[pow]/sAmp[pow])**2)**.5)

popt, pcov = curveFit(lin, truePowers, ratios, [.1, 1], sigma=ampErr, absolute_sigma=True)
m, b = popt[0], popt[1]
print(f"m: {m: .6f} \t\t b: {b: .2f}")

plt.errorbar(truePowers, ratios, yerr=ampErr, fmt="None", elinewidth=.5, color='Green', alpha=.5, capsize=1, capthick=.5)
plt.plot(np.array([0,125]), lin(np.array([0,125]), m, b), color="Black", linewidth=1)

plt.savefig(save + "P-O Height Ratios.pdf", format="pdf")
plt.savefig(save + "P-O Height Ratios.png", format="png")


# plot peak height normalized by power vs power, s & aS same plot, pts only
plt.figure(dpi=250)
plt.title("Power-normalized Peak Height vs Power")
plt.xlabel("Pump Power (mW)")
plt.ylabel("Power-normalized Peak Height")
plt.xlim(0,125)
#plt.ylim()
plt.minorticks_on()
plt.tick_params(which='both', direction='in', pad=5)

aSAmpNorm, sAmpNorm = [], []
for (pow, truPow) in zip(powers, truePowers):
    aSAmpNorm.append(aSAmp[pow]/truPow)
    sAmpNorm.append(sAmp[pow]/truPow)

plt.scatter(truePowers, aSAmpNorm, 1, color="blue", label="anti-Stokes")
plt.scatter(truePowers, sAmpNorm, 1, color="red", label="Stokes")
plt.legend()
plt.savefig(f"{save} P-O Pow-norm Heights vs Pow.pdf", format="pdf")
plt.savefig(f"{save} P-O Pow-norm Heights vs Pow.png", format="png")


#Simulation data GB_797nm_LCOF (Nov 30, 2022 email from ryan)
simData = pd.read_csv("Simulations/simulated_GB_797nm_LCOF.csv", skiprows=5)
sim = pd.DataFrame({
    'Freq': simData['freq']/1e9,
    'Sig': simData['gb']
})

popt, pcov = curveFit(l, sim['Freq'], sim['Sig'], guess)
simAmp, wid, cen, c = popt[0], abs(2*popt[1]), popt[2], popt[3]

aS['10']['ScaledSig'] = aS['10']['Sig']*(simAmp/aSAmp['10'])

plt.figure(dpi=250)
plt.title("10 mW anti-Stokes Pump-Only Power")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Spectral Density (μV)")
plt.xlim(2,2.5)
#plt.ylim()
plt.minorticks_on()
plt.tick_params(which='both', direction='in', pad=5)

plt.scatter(aS['10']['Freq'], aS['10']['ScaledSig'], 1, label='Observed Data')
plt.scatter(sim['Freq'], sim['Sig'], 1, label='Simulation')

plt.legend(fontsize=7.5)
plt.savefig(f"{save} P-O Simulation.pdf", format="pdf")
plt.savefig(f"{save} P-O Simulation.png", format="png")
