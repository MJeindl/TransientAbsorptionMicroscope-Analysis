from scipy.io import loadmat
import sys
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from mpl_toolkits.mplot3d import Axes3D 

plt.rcParams['text.usetex'] = True

r'''
#sumFilename = sys.argv[1]
#relPosBool = bool(sys.argv[2])
sumFilename = r"C:\Users\Maximilian\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm\SampleAreaImagingPump653Probe680\images_2024-01-05_04-09.mat"
#sumFilename = r"C:\Users\Maximilian\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm\SampleAreaImagingPump653Probe680\images_2024-01-05_08-38.mat"
relPosBool = False

xlabel_step = 10
ylabel_step = 10

data = loadmat(sumFilename)
filenames = data["filenames"][0]
delays = data['delays'][0]

x0 = data['x0_vec'][0][0]
y0 = data['y0_vec'][0][0]




#to normalise values
meas_number = len(filenames)

count = 0
for filename in filenames:
    if os.path.isfile(os.path.dirname(sumFilename) + "/" + filename[0] + ".mat") == False:
        #correct which files are there
        raise IOError
        delays.delete(count)
        count -= 1
        continue
    
    
    currentData = loadmat(os.path.dirname(sumFilename)+"/" + filename[0] + ".mat")
    d_map = currentData["I1_vec"]/currentData["I2_vec"]

    if count == 0:
        d_maps = np.zeros((meas_number, np.shape(d_map)[0], np.shape(d_map)[1]))
        x_vec = np.array(currentData["X"])
        y_vec = np.array(currentData["Y"])
        ratio_yx = (np.max(y_vec) - np.min(y_vec))/(np.max(x_vec) - np.min(x_vec))

    d_maps[count] = d_map

    count += 1


if relPosBool == 1:
    x_vec = x_vec - x0
    y_vec = y_vec - y0

for index in range(meas_number-1):
    fig, ax = plt.subplots(subplot_kw={"projection": "3d", "proj_type": "ortho"})
    ax.plot_surface(x_vec, y_vec, -1e3*np.log10(d_maps[index+1]- d_maps[0] + 1), cmap=cm.coolwarm)
    ax.view_init(90,0)
    ax.axis("equal")
    ax.set_xlim([np.min(x_vec), np.max(x_vec)])
    ax.set_ylim([np.min(y_vec), np.max(y_vec)])
    #plt.xticks(np.arange(len(x_vec))[::xlabel_step], np.round(x_vec[::xlabel_step], 4))
    #plt.yticks(np.arange(len(y_vec))[::ylabel_step], np.round(y_vec[::ylabel_step], 4))
    #plt.xlabel("x / mm")
    #plt.ylabel("y / mm")
    fig.suptitle("Delay: %.1d ps\n%s - %s" %(delays[index+1]*1e-3, filenames[index+1], filenames[0]))
    #plt.title("Delay: %.1d ps\n%s - %s" %(delays[index+1]*1e-3, filenames[index+1], filenames[0]))
plt.show()

pathing = r"C:\Users\M\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm\SampleAreaImagingPump653Probe680"
showImageFullPy(pathing, r"\image_0245", r"\image_0244", cbar_range=[-3,3])
'''

#decided to do a blatant copy of showImageFull


def showImageFullPy(pathing, filename, bgfilename, x0y0 = [0,0], cbar_range=[None,None]):
    fullPath = lambda fname: pathing + fname + r".mat"
    file = loadmat(fullPath(filename))
    bgfile = loadmat(fullPath(bgfilename))

    #if cbar_range == None:
    #    cbar_range = [None,None]

    T = file['I'] #transient transmissivity

    tempBack = np.zeros(np.shape(file["I"]))
    if np.shape(file["X"]) != np.shape(bgfile["X"]):
        #temporary fix, not fully implemented
        if np.mod(np.shape(bgfile["X"][0]), 2) == 1:
            #for odd
            timesX = (np.shape(bgfile['X'])[0]-1)/(np.shape(file['X'])[0]-1)
        else:
            #for even
            timesX = (np.shape(bgfile['X'])[0])/(np.shape(file['X'])[0])
        if np.mod(np.shape(bgfile['X'])[1], 2) == 1:
            #for odd
            timesY = (np.shape(bgfile['X'])[1]-1)/(np.shape(file['X'])[1]-1)
        else:
            #for even
            timesY = (np.shape(bgfile['X'])[1])/(np.shape(file['X'])[1])
    else:
        timesX = 1
        timesY = 1
    tempBack = bgfile["I"][0:np.shape(bgfile['I'])[0]:int(timesX), 0:np.shape(bgfile['I'])[1]:int(timesY)]
    T = T - tempBack + 1

    I = -np.log10(T)*1000

    X_arr = file['X'] - x0y0[0]
    Y_arr = file['Y'] - x0y0[1]

    fig, ax = plt.subplots(1,1, figsize=(5.78,5.78/3), dpi = 288)
    map = ax.pcolor(X_arr, Y_arr, I, cmap="plasma", vmin = cbar_range[0], vmax = cbar_range[1])
    plt.colorbar(map, ax = ax, label=r"$\Delta A$ / mOD")
    ax.set_ylabel('y / mm')
    ax.set_xlabel('x / mm')
    ax.set_aspect('equal')
    plt.show()

