#this is used to set all plots the same
from matplotlib import pyplot as plt

def setMatplotSettings():
    plt.rcParams['text.usetex'] = True
    font = {'size': 11}
    plt.rc('font', **font)
    plt.rcParams["axes.formatter.limits"] = [-3,3]

def figSizer(fractionWidth, aspectRatio=1.618,  scale = 1):
    '''scale = 0.5 means that the plot is meant to be used at scale = 0.5 in latex\\
        fractionwidth means fraction of textwidth in latex'''
    latexwidth = 5.788 # inches
    targetWidth = latexwidth/(fractionWidth*scale)
    targetHeight = targetWidth/aspectRatio
    return (targetWidth, targetHeight)
