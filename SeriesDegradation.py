import pandas as pd
from matplotlib import pyplot as plt
from scipy import constants as const
from scipy.io import loadmat
import os
import numpy as np
from datetime import timedelta
import json
from scipy.optimize import curve_fit

import sys
import plotHelperLatex
plotHelperLatex.setMatplotSettings()
#sys.path.insert(1, r"C:\Users\M\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm")
#from ShowDelayScan import fromFiles

from argparse import ArgumentParser
from fileParsingMethods import getTimes, parseTime, parseSummaryFileToArray, parseFilenames, removeBackground, parseSummaryFiletoRaw



import plotHelperLatex
plotHelperLatex.setMatplotSettings()


def autoCompensation(filenames, degConst, p_ratio, dirPath=r"", backgroundMean = 15, linearBg = False):
    '''assumes they are continuous measurements/saturation type file, will not work properly if they are not\\ 
    returns absorbance and transmittance arrays + correction factor array'''
    deltaTimes = parseTime(getTimes(filenames, dirPath))
    dArray, delay, backgroundParameters = parseSummaryFileToArray(filenames, dirPath, linearBackgroundSubtract=linearBg, backgroundLen=backgroundMean)
    dArray, aArray, _,  = removeBackground(dArray, backgroundMean)
    correctionFactors = degradationCompensation(degConst, deltaTimes, len(delay), p_ratio)
    return aArray, dArray, correctionFactors,delay, backgroundParameters


def degradationCompensation(degConstant, times, decaysteps: int, powerDensities=1):
    '''returns correction factor for each measurement\\
        not sure what to do with powerDensities, got to try a few things and see if they work but also read up on it
        '''

    meanTime = np.mean((times[1:]-times[:-1]))
    #print(decaysteps)
    Correction = np.zeros((len(times), decaysteps))

    if not hasattr(powerDensities, "len"):
        powerDensities = np.ones(np.shape(times)[0])
    #should also correct for within a measurement, though that part is more estimate (should be fine for long decay times)
    #Correction = np.zeros((len(times), decaysteps))
    for i in range(len(times)-1):
        Correction[i] = np.exp((times[i]+(times[i+1]-times[i])*np.linspace(0,1, decaysteps))*powerDensities[i]/degConstant)
    Correction[len(times)-1] = np.exp((times[len(times)-1]+meanTime*np.linspace(0,1, decaysteps))*powerDensities[i]/degConstant)
    return Correction



def plotTrend(filePaths, dirPath=r"", figsize_in=(8,4), symlog=False):
    #need to switch getTimes to the same system as parseSummary to allow any type of input
    times = getTimes(filePaths, dirPath)
    dArray, delay, _ = parseSummaryFileToArray(filePaths, dirPath)
    dArray, OD, _ = removeBackground(dArray, 10)
    #get delta t of times
    timesFromZero = np.array(parseTime(times), dtype=int)

    delays, measurementTimes = np.meshgrid(delay, timesFromZero)
    fig, ax = plt.subplots(1,1, figsize=figsize_in, dpi = 288)
    map = ax.pcolor(delays*1e-3, measurementTimes, OD, cmap="plasma")
    ax.set_xlabel('delay time / ps')
    ax.set_ylabel('total time elapsed at start / s')
    fig.colorbar(map, ax=ax, label=r"$\Delta A$ / mOD")
    ax.set_yticks(timesFromZero-(timesFromZero[1]-timesFromZero[0])/2, timesFromZero)
    if symlog==True:
        ax.set_xscale("symlog")
    plt.tight_layout()
    plt.show()



#degradationCompensation(1, [r"C:\Users\M\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm\Pump653Probe493_Degradation\saturation_2024-01-19_16-14.mat"], 1)
if __name__ == "__main__":
    """
    parser = ArgumentParser(description="Shows response trend over mulitple adjacent measurements")
    parser.add_argument("startFileNumber")
    parser.add_argument("endFileNumber")
    parser.add_argument("pathDir")
    
    args = parser.parse_args()

    start = int(args.startFileNumber)
    stop = int(args.endFileNumber)
    plotTrend(args.pathDir, start, stop)"""
    

def fitDegradation(inputArray, times, powerDensity=1, sliding_window_len = 5, t0=0):
    '''only works for single series without interruptions\n 
    returns popt [tau]'''
    #find the mean time to adjust for the final part of the measurement
    for i in range(len(times)-1):
        dtime = times[i+1]-times[i]
    dtime = dtime/(len(times)-1)
    times = np.array(times)

    working_array = np.zeros((len(times), np.shape(inputArray)[1]-sliding_window_len+1))
    for i in range(len(times)):
        working_array[i,:] = np.convolve(inputArray[i,:], np.ones(sliding_window_len)/sliding_window_len, mode = 'valid')


    def costFunction(times, tau):
        time_fit = lambda time: -time/tau
        residual = 0
        #cut away any data before temporal overlap, more sophisticated methods did not work
        for delayInd in range(int(t0),np.shape(working_array)[1]):
            ratio = abs((working_array[:, delayInd])/(working_array[0,delayInd]))
            residual += abs(np.sum(np.log(ratio[:])-time_fit(times)))**2

        return residual

    popt, pcov = curve_fit(costFunction, times, np.zeros(np.shape(working_array)[0],dtype=float), p0 = [1e4])

    return popt, pcov[0,0]





def ignore():
    #this is the constant pump power and probe power variation
    filenames=[]
    for x in range(6645,6654,1):
        filenames.append("TA_fourier_%4d" %(x))
    dirPath = r"C:\Users\M\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm\STD-680\2023.12.22_PumpPowerVar"


    dirPath = r"C:\Users\M\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm\Pump653Probe493_Degradation"
    filenames = [r"\saturation_2024-01-19_16-14"]
        
    dArray, delay = parseSummaryFileToArray(filenames, dirPath)

    for ind in range(len(filenames)):
        filenames[ind] = dirPath +r"\\"[0] + filenames[ind]
    dtime = getTimes(filenames, dirPath)
    dtime = parseTime(dtime)
    #dOD = powerVar['dOD / mOD']
    print(np.shape(dArray))
    dArray, aArray, _ = removeBackground(dArray, 20)
    popt, pcov = fitDegradation(aArray[:,:], dtime[:])

    print(popt)
    print(np.sqrt(pcov))



    file = open(r"C:\Users\M\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm\Pump653Probe493_Degradation\TAfitPump653Probe493_OUTPUT.JSON")
    entries = json.load(file)['entries']
    entries_files = []
    timeAt = np.array([0,1e3,5e3]) #fs
    dOD = np.zeros((len(entries), len(timeAt)))
    expDecay = lambda ampTau1, ampTau2, t: ampTau1[0]*np.exp(-t/ampTau1[1])+ampTau2[0]*np.exp(-t/ampTau2[1])

    fig, ax = plt.subplots(1,1, figsize=(plotHelperLatex.figSizer(2,3)), dpi = 144)

    indices = [20,25,40]
    print(np.shape(aArray))
    #print(aArray)
    for i in indices:
        ax.plot(dtime, aArray[:,i], label='std')
        ax.plot(dtime, aArray[:,i]*np.exp(dtime/popt))
    #print((1-degradationCompensation(popt, dtime, np.shape(aArray)[0])[:,0])/abs(aArray[:,0]))
    #print(np.exp(dtime/popt))
    #print(aArray)
    #ax.plot(dtime, aArray[:,50]*np.exp(dtime/popt))
    ax.legend()
    ax.set_xlabel('time / s')
    ax.set_ylabel('absorbance / mOD')

    plt.show()




