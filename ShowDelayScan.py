import numpy as np
from matplotlib import pyplot as plt
from scipy.io import loadmat
import os.path
from argparse import ArgumentParser
import sys
sys.path.insert(1, r"C:\Users\M\Documents\Books\masterprojectinformation")
from fileParsingMethods import parseSummaryFileToArray, removeBackground
import plotHelperLatex
plotHelperLatex.setMatplotSettings()



def fromFiles(filenames, dirPath =r"", fig_size=(8,3), showMean=True, showLegend=False):
    datArray, delay, _ = parseSummaryFileToArray(filenames, dirPath)
    nFiles = np.shape(datArray)[0]

    dArray, aArray, _ = removeBackground(datArray, 10)

    fig, ax = plt.subplots(1,1, figsize = fig_size, dpi=288)

    #colors according to colormap
    colors = plt.cm.cool(np.linspace(0,1,nFiles))

    for i in range(np.shape(aArray)[0]):
        ax.plot(delay*1e-3, aArray[i,:], color = colors[i], label = "%s" %filenames[i])
    #plot mean
    if showMean == True:
        ax.plot(delay*1e-3, np.mean(aArray, axis = 0), 'k')
    ax.set_ylabel(r'$\Delta A$ / mOD')
    ax.set_xlabel('time delay / ps')
    if showLegend == True:
        plt.legend()
    #plt.tight_layout()
    plt.show()

def fromRawData(data, delay, fig_size=(4,2), showMean = True):
    nFiles = np.shape(data)[0]

    #dArray, aArray, _ = removeBackground(data, 10)

    fig, ax = plt.subplots(1,1, figsize = fig_size, dpi=288)

    #colors according to colormap
    colors = plt.cm.cool(np.linspace(0,1,nFiles))
    if np.ndim(data) == 1:
        ax.plot(delay*1e-3, data)
    else:
        for i in range(np.shape(data)[0]):
            ax.plot(delay*1e-3, data[i,:], color = colors[i])
        #plot mean
        if showMean == True:
            ax.plot(delay*1e-3, np.mean(data, axis = 0), 'k')
    ax.set_ylabel(r'$\Delta A$ / mOD')
    ax.set_xlabel('time delay / ps')
    plt.tight_layout()
    plt.show()


def wav_scan(filenames, wavelengths, dirPath=r"", figsize=(8,3), sliding_window_len=3):
    return
    datArray, delay = parseSummaryFileToArray(filenames, dirPath)
    aArray = removeBackground(datArray, 15)

    working_array = np.zeros((len(times), np.shape(aArray)[1]-sliding_window_len+1))
    for i in range(len(times)):
        working_array[i,:] = np.convolve(aArray[i,:], np.ones(sliding_window_len)/sliding_window_len, mode = 'valid')

    fig, ax = plt.subplots(1,1, figsize=(8,4), dpi = 288)
    map = ax.pcolor(delay*1e-3, measurementTimes, CorrArray, cmap="plasma", snap=True)
    plt

if __name__ == "__main__":
    parser = ArgumentParser(description="Quickly shows individual measurements")
    parser.add_argument("filename")
    parser.add_argument("--showMean", type = str, help="If the mean should be shown or not (default True)")
    args = parser.parse_args()
    print(args.showMean)
    if args.showMean == "False":
        fromFiles(args.filename, showMean = False)
    else:
        fromFiles(args.filename)
    



