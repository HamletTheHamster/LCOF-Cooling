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
from scipy import stats
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
    truePowers.append(float(int(pow))*.17**.5) # - 10 mW

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
        for i in range(1,files):
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
    return amp * wid**2 / (wid**2 + (x-cen)**2) + c

def lin(x, m, b):
    return m*x+b

def lnorm(x, gamma_0, gamma_eff, cen):
     return gamma_0**2 / (gamma_eff**2 + (x-cen)**2)

guess = a, w, ce, c

#------------------------------------------------------------------#

binGHz = 0.01
nBins = int((aS['0']['Freq'].iloc[-1] - aS['0']['Freq'].iloc[0]) / binGHz + 1)

bin = {}
for pow in powers:
    bound = aS[pow]['Freq'].iloc[0]

    binFreqs = []
    binSigs = []
    binErrs = []

    for n in range(nBins):
        bound += binGHz

        # Collect all Sig points in this bin
        sigsInBin = []
        for (freq, sig) in zip(aS[pow]['Freq'], aS[pow]['Sig']):
            if (freq < bound) and (freq > bound - binGHz):
                sigsInBin.append(sig)

        # Frequency for the bin’s center
        binFreqs.append(bound - binGHz / 2)

        # Mean value in this bin
        if len(sigsInBin) > 0:
            meanVal = np.mean(sigsInBin)
        else:
            meanVal = np.nan
        binSigs.append(meanVal)

        # Standard error of the mean for that bin
        if len(sigsInBin) > 1:
            semVal = np.std(sigsInBin, ddof=1) / np.sqrt(len(sigsInBin))
        else:
            semVal = np.nan  # Or 0, or however you want to handle empty bins
        binErrs.append(semVal)

    bin[pow] = pd.DataFrame({
        'Freq': binFreqs,
        'Sig': binSigs,
        'σ': binErrs,
    })

#--------------------------------------------------------------#

# plot fits w error
plt.figure(dpi=250)
plt.title("Experiment B: anti-Stokes Spectral Fits")
plt.xlabel("Frequency [(⍵ - ⍵$_{P}$)/2π] (GHz)")
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
  fwhmσ.append(2*σWid*1000)
  print(f"σAmp: {σAmp:.4f} μV \t σWid: {σWid*1000: .4f} MHz \t σCen: {σCen: .4f} GHz \t σC: {σC: .4f} μV")
  #plt.errorbar(bin[pow]['Freq'], bin[pow]['Sig'], yerr=bin[pow]['σ'], fmt="None", elinewidth=.25, color=paletteDict[pow], alpha=.25, capsize=1, capthick=.25)
  plt.plot(aS[pow]['Freq'], l(aS[pow]['Freq'], *popt), color=paletteDict[pow], linewidth=2, label=f"{truPow: .1f}")
  #plt.scatter(bin[pow]['Freq'], bin[pow]['Sig'], 40, edgecolors=paletteDict[pow], facecolors="none", marker="o", alpha=.5 )#label=f"{truPow: .2f}"+' mW')

plt.legend(title="Pump Power (mW)")
plt.savefig(f"{save}P-P anti-Stokes Fits.pdf", format="pdf")
plt.savefig(f"{save}P-P anti-Stokes Fits.png", format="png")

# Create lists to store results (if you still want to do a combined "pow v wid" plot later)
fwhm = []
fwhmσ = []
amps = []
ampsσ = []

# Loop over each power and produce *individual* fits & plots
for (pow, truPow) in zip(powers, truePowers):
    # 1. Fit this power's data
    popt, pcov = curveFit(l, aS[pow]['Freq'], aS[pow]['Sig'],
                          p0=guess,
                          sigma=aS[pow]['σ'],
                          absolute_sigma=True)
    amp, widHalf, cen, c = popt
    # Note: widHalf is the half-width in the fitted Lorentzian formula,
    #       so the full FWHM is 2*widHalf. We'll store that in "wid."
    wid = abs(2*widHalf)

    # Extract errors
    σAmp = np.sqrt(pcov[0,0])
    σWidHalf = np.sqrt(pcov[1,1])
    # Error in the full width = 2 * σWidHalf
    σWid = abs(2*σWidHalf)
    σCen = np.sqrt(pcov[2,2])
    σC   = np.sqrt(pcov[3,3])

    # Save or print out for your own records
    print(f"\n=== Power {pow} (true power {truPow:.2f} mW) ===")
    print(f"Amp: {amp:.3f} μV   FWHM: {wid*1000:.3f} MHz   Cen: {cen:.3f} GHz   C: {c:.3f} μV")
    print(f"σAmp: {σAmp:.4f}, σWid: {σWid*1000:.4f} MHz, σCen: {σCen:.4f} GHz, σC: {σC:.4f} μV")

    # 2. Store FWHM and its error for later "pow v wid" plot
    fwhm.append(wid*1000)       # store in MHz
    fwhmσ.append(σWid*1000)
    amps.append(amp)
    ampsσ.append(σAmp)

    # 3. Now create a *new figure* for each power
    plt.figure(dpi=250)
    plt.title(f"Experiment B: anti-Stokes Power Spectrum\nP$_P$ = {truPow:.1f} mW")
    plt.xlabel("Frequency [(⍵ - ⍵$_{P}$)/2π] (GHz)")
    plt.ylabel("Spectral Density (μV)")
    plt.xlim(2, 2.5)
    plt.ylim(0, 2)
    plt.minorticks_on()
    plt.tick_params(which='both', direction='in', pad=5)
    plt.tick_params(which='minor', axis='y', length=0)

    freq_array = bin[pow]['Freq']

    # Plot the *data points* for this single power
    plt.errorbar(bin[pow]['Freq'],
                bin[pow]['Sig'],
                yerr=bin[pow]['σ'],
                fmt="o",
                elinewidth=1,
                color=paletteDict[pow],
                ecolor=paletteDict[pow],
                alpha=.5,
                capsize=3,
                label="Observed Data"
                # capthick=.5,
    )

    # plt.scatter(freq_array,
    #             bin[pow]['Sig'],
    #             s=80,
    #             edgecolors=paletteDict[pow],
    #             facecolors="none",
    #             marker="o",
    #             #alpha=.5,
    #             label="Observed Data")

    # Plot the *fitted curve* for this single power
    plt.plot(aS[pow]['Freq'],
             l(aS[pow]['Freq'], *popt),
             color=paletteDict[pow],
             linewidth=3,
             label="Lorentzian Fit")

    plt.legend(title=f"$P_P$ = {truPow:.1f} mW")

    # 4. Save figure (your 'save' directory was already created above).
    plt.savefig(f"{save}P-P anti-Stokes Fit - {pow}mW.pdf", format="pdf")
    plt.savefig(f"{save}P-P anti-Stokes Fit - {pow}mW.png", format="png")

    # 5. (Optional) close the figure so we don't overwrite
    plt.close()

#-------------------------------------------------------------------#

# plot pow v wid
popt, pcov = curveFit(lin, truePowers, fwhm, [.1, 96.], sigma=fwhmσ, absolute_sigma=True)
m, b = popt[0], popt[1]
print(f"Linewidths | m: {m:.4f}, b: {b:.4f}")

paletteList = [(31/255, 211/255, 172/255), (255/255, 122/255, 180/255), (122/255, 156/255, 255/255), (255/255, 182/255, 110/255)]

plt.figure(dpi=250)
plt.title("Experiment B: anti-Stokes Linewidths vs Pump Power")
plt.xlabel("Pump Power (mW)")
plt.ylabel("Linewidth (MHz)")
plt.xlim(-2,72)
#plt.ylim(,)
plt.minorticks_on()
plt.tick_params(which='both', direction='in', pad=5)

plt.plot(np.array([-2,72]), lin(np.array([-2,72]), m, b), color="darkGray", linestyle='--', linewidth=3, label="Linear Fit")

for xval, yval, err, c in zip(truePowers, fwhm, fwhmσ, paletteList):
    plt.errorbar(
        xval, yval, yerr=err,
        fmt="o",            # marker style
        color=c,            # marker color
        ecolor=c,           # error bar color
        elinewidth=1,
        alpha=.5,
        capsize=3,
        label=f"{xval: .1f}"
    )

    # plt.scatter(
    #     xval, yval,
    #     100,
    #     edgecolors=c,
    #     facecolors="none",
    #     marker="o",
    #     label=f"{xval: .1f}"
    # )

plt.legend(title="Pump Power (mW)")

plt.savefig(f"{save}P-P anti-Stokes Wid v Pow.pdf", format="pdf")
plt.savefig(f"{save}P-P anti-Stokes Wid v Pow.png", format="png")

#-------------------------------------------------------------------#

# plot heights v wid

popt, pcov = curveFit(lin, truePowers, amps, [-.1, 2], sigma=ampsσ, absolute_sigma=True)
m, b = popt[0], popt[1]
print(f"Amplitudes | m: {m:.4f}, b: {b:.4f}")

plt.figure(dpi=250)
plt.title("Experiment B: anti-Stokes Peak Amplitudes vs Pump Power")
plt.xlabel("Pump Power (mW)")
plt.ylabel("Peak Spectral Density (μV)")
plt.xlim(-2, 72)
plt.ylim(1.25, 2)
plt.minorticks_on()
plt.tick_params(which='both', direction='in', pad=5)
plt.tick_params(which='minor', axis='y', length=0)

plt.plot(np.array([-2,72]), lin(np.array([-2,72]), m, b), color="darkGray", linestyle='--', linewidth=3, label="Linear Fit")

for xval, yval, err, c in zip(truePowers, amps, ampsσ, paletteList):
    plt.errorbar(
        xval, yval, yerr=err,
        fmt="o",            # marker style
        color=c,            # marker color
        ecolor=c,           # error bar color
        elinewidth=1,
        alpha=.5,
        capsize=3,
        label=f"{xval: .1f}"
    )

    # plt.scatter(
    #     xval, yval,
    #     100,
    #     edgecolors=c,
    #     facecolors="none",
    #     marker="o",
    #     label=f"{xval: .1f}"
    # )

plt.legend(title="Pump Power (mW)")

plt.savefig(f"{save}P-P anti-Stokes Height v Pow.pdf", format="pdf")
plt.savefig(f"{save}P-P anti-Stokes Height v Pow.png", format="png")

# #normalize
# aSnorm = {}
# for (pow, truPow) in zip(powers, truePowers):
#     aSnorm[pow] = pd.DataFrame({
#     'Freq': aS[pow]['Freq'],
#     'Sig': aS[pow]['Sig']/aS[pow]['Sig'].max(axis=0),
#     'σ': aS[pow]['σ']/aS[pow]['Sig'].max(axis=0)})
#
#
# #bin
# #binGHz = .007
# nBins = int((aSnorm['0']['Freq'][len(aSnorm['0']['Freq']) - 1] - aSnorm['0']['Freq'][0])/binGHz + 1)
#
# bin = {}
# for pow in powers:
#     bound = aSnorm[pow]['Freq'][0]
#
#     binFreqs = []
#     binSigs = []
#     for n in range(nBins):
#         bound += binGHz
#
#         sigsInBin = []
#         for (freq, sig) in zip(aSnorm[pow]['Freq'], aSnorm[pow]['Sig']):
#             if (freq < bound) and (freq > bound - binGHz):
#                 sigsInBin.append(sig)
#
#         binFreqs.append(bound - (binGHz/2))
#         binSigs.append(np.mean(sigsInBin))
#
#     bin[pow] = pd.DataFrame({
#         'Freq': binFreqs,
#         'Sig': binSigs,
#     })
#
# # normalized plots
# plt.figure(dpi=250)
# plt.title(f"Experiment B: anti-Stokes (Normalized by Power)")
# plt.xlabel("Frequency [(⍵ - ⍵$_{P}$)/2π] (GHz)")
# plt.ylabel("Spectral Density (μV)")
# plt.xlim(2,2.5)
# plt.ylim(0,1.1)
# plt.minorticks_on()
# plt.tick_params(which='both', direction='in', pad=5)
# plt.tick_params(which='minor', axis='y', length=0)
#
# blueGradient = []
# for i in range(len(powers)):
#     mod = (185/4)*i/255
#     blueGradient.append((0+mod, 0+mod, 1))
#
# blueGradient.reverse()
# blueGrad = dict(zip(powers, blueGradient))
#
# pumpOnlyGamma_eff = {'0': 97.5e-3/2, '55': 100.213e-3/2, '110': 101.828e-3/2, '165': 103.442e-3/2}
# for (pow, truPow) in zip(powers, truePowers):
#     plt.scatter(bin[pow]['Freq'], bin[pow]['Sig'], 1, color=paletteDict[pow])#, label=f"{truPow: .2f}"+" mW")
#     #plt.errorbar(aSnorm[pow]['Freq'], aSnorm[pow]['Sig'], yerr=aSnorm[pow]['σ'], fmt='none', elinewidth=.25, color=paletteDict[pow], alpha=.5, capsize=1, capthick=.25, label=f"{truPow: .2f} mW Pump-Probe Data")
#     plt.plot(aSnorm[pow]['Freq'], lnorm(aS[pow]['Freq'], 97.315e-3/2, pumpOnlyGamma_eff[pow], 2.269), color=paletteDict[pow], linewidth=1, label=f"{truPow: .2f} mW Pump-Only Synthesized")
#
# #plt.legend()
# plt.savefig(f"{save}P-P Normalized anti-Stokes Fits.pdf", format="pdf")
# plt.savefig(f"{save}P-P Normalized anti-Stokes Fits.png", format="png")
