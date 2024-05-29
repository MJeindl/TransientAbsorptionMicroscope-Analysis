#Version 1.3.2 from 2024.05.14

#XML should be working
#Introduced interpolation into getNumericCorrection, seems to be working as well as the fitted one for good overlap conditions
#note that with the camera pixel density it likely can still fluctuate strongly with less ideal conditions
#it is likely to be superior when beamshapes are distorted off diagonal though and should be used to check the analytic output

#idea
#use maximum gradient + orthogonal for fitting
#probably not needed with

###TODO

## Implement gradient and orthogonal for fitting -> seems not to be feasible as it makes analytical comparison impossible
## Implement self.scale DONE, Q: should it be automatically applied to overlap correction? I say no because power is not included so far



import numpy as np
import json
import os.path
from matplotlib import pyplot as plt
import warnings
import scipy.signal as signal
import scipy.optimize as optimize
from pandas import DataFrame
import pandas

from xml.dom import minidom
import base64


class ArtFit:
    
    #constructor
    def __init__(self, *inp):
        '''
        Inputs are (parenthesis only for visual reasons, do not use for input): \n
        -> Single filepath to quickly plot a square map \n
        -> [Path, anything] to plot with standard dimensions 1280x720\n
        -> [BackgroundJsonPath, probePath, width of average] average width default is 1\n
        -> [BackgroundJsonPath, probePath, pumpPath, width]
        '''
        self.scale = np.NaN
        self.repetitionRate = 4e4
        self.size = np.array([0,0], dtype=int)
        if len(inp) == 1:
            print("single input")
            self.construct_QuickPlot(inp[0], False)
        if len(inp) == 2:
            self.construct_QuickPlot(inp[0], True)
        elif len(inp) == 3:
            print("background and probe")
            backgroundJson = inp[0]
            probePath = inp[1]
            width = inp[2]
            self.construct_Single(backgroundJson, probePath, width)
        elif len(inp) == 4:
            print("bg, probe and pump")
            backgroundJson = inp[0]
            probePath = inp[1]
            pumpPath = inp[2]
            width = inp[3]
            self.construct_Dual(backgroundJson, probePath, pumpPath, width)



    #for two images
    def construct_Dual(self, backgroundJson, probePath, pumpPath, width = 1):
        #hmmm, the maximum value actually read is 4080 which is missing 4 bit on 4095...
        self.max_val = 4095
        
        if probePath.endswith('.txt'):
            #load byte array
            byteArr = self.loadByteArr(probePath)
            #convert to integer array
            self.probeMap = self.convertByteToInt(byteArr)
        elif probePath.endswith('.xml'):
            byteArr, hstart, hstop, vstart, vstop, self.date = self.loadXMLimage(probePath)
            self.x_pos = hstart-1
            self.y_pos = vstart -1
            #hstop is including hstop pixel
            self.size[0] = abs(hstop - hstart)+1
            self.size[1] = abs(vstop - vstart)+1
            self.probeMap = self.convertByteToInt(byteArr, self.size[0], self.size[1])
        elif (probePath.endswith('.tiff') or probePath.endswith('.TIFF') or probePath.endswith('.tif')):
            #additional tiff support now
            self.probeMap, self.size = self.loadTIFF(probePath)


        if pumpPath.endswith('.txt'):
            byteArr = self.loadByteArr(pumpPath)
            self.pumpMap = self.convertByteToInt(byteArr)
        elif pumpPath.endswith('.xml'):
            byteArr, hstart, hstop, vstart, vstop, date = self.loadXMLimage(pumpPath)
            if date != self.date:
                print("Pump date is not the same as probe date")
            if abs(hstop-hstart) != self.size[0] or abs(vstop - vstart) != self.size[1]:
                raise ImportError("Beware pump and probe map are not the same dimension")
            if (hstart != self.x_pos +1) or (vstart != self.y_pos +1):
                raise ImportError("Beware pump and probe ROI not identical")
            self.pumpMap = self.convertByteToInt(byteArr, self.size[0], self.size[1])
        elif (pumpPath.endswith('.tiff') or pumpPath.endswith('.TIFF')or pumpPath.endswith('.tif')):
            self.pumpMap, checkSize = self.loadTIFF(pumpPath)
            if np.any(checkSize != self.size):
                raise ImportError("Beware pump and probe map are not the same dimension")

        #would be cleaner to forward the actual size instead of the sqrt
        """ Not needed right now since it is actually checked for in convertByteToInt (suboptimal implementation)
        if self.probeSize != self.pumpSize:
            raise ImportError("The imported raws are not of the same size")
        """
        #doing this in case another camera comes around, to make it easier to find a place to start
        bgSize = np.array([1280, 720], dtype=int)
        #end with background
        bg_path = self.getJSON(backgroundJson)
        if bg_path.endswith('.txt'):
            bg_map = self.convertByteToInt(self.loadByteArr(bg_path), bgSize[0], bgSize[1])

        elif bg_path.endswith('.xml'):
            bg_map = self.convertByteToInt(self.loadXMLimage(bg_path)[0], bgSize[0], bgSize[1])
        elif bg_path.endswith('.TIFF') or bg_path.endswith('.tif') or bg_path.endswith('.tiff'):
            bg_map = self.convertByteToInt(self.loadTIFF(bg_path), bgSize[0], bgSize[1])
        else:
            raise ImportError("bg_map has unknown format")
        self.bgMap = self.cut_map(bg_map)
        #print(self.size)


        #note that the width is the average over how many lines one takes
        #supposed to help with worse beamshapes
        self.width = width
        self.status_identifier = 2
    
    def construct_Single(self, backgroundJson, probePath, width = 1):
        self.max_val = 4095
        if probePath.endswith('.txt'):
            #load byte array
            byteArr = self.loadByteArr(probePath)
            #convert to integer array
            self.probeMap = self.convertByteToInt(byteArr)
        elif probePath.endswith('.xml'):
            byteArr, hstart, hstop, vstart, vstop, self.date = self.loadXMLimage(probePath)
            self.x_pos = hstart -1
            self.y_pos = vstart -1
            #hstop is including hstop pixel
            self.size[0] = abs(hstop - hstart)+1
            self.size[1] = abs(vstop - vstart)+1
            self.probeMap = self.convertByteToInt(byteArr, self.size[0], self.size[1])

        #end with background
        bgSize = np.array([1280, 720], dtype=int)
        #end with background
        bg_path = self.getJSON(backgroundJson)
        if bg_path.endswith('.txt'):
            bg_map = self.convertByteToInt(self.loadByteArr(bg_path), bgSize[0], bgSize[1])
        elif bg_path.endswith('.xml'):
            bg_map = self.convertByteToInt(self.loadXMLimage(bg_path)[0], bgSize[0], bgSize[1])
        elif bg_path.endswith('.TIFF') or bg_path.endswith('.tif') or bg_path.endswith('.tiff'):
            bg_map = self.convertByteToInt(self.loadTIFF(bg_path), bgSize[0], bgSize[1])
        else:
            raise ImportError("bg_map has unknown format")
        self.bgMap = self.cut_map(bg_map)

        self.width = width
        self.status_identifier = 1

    #to just plot something quickly
    def construct_QuickPlot(self, probePath, sizeBool):
        '''
        sizeBool is for standard size = 1280*720 instead of square map
        '''
        if probePath.endswith('.txt'):
            #load byte array
            byteArr = self.loadByteArr(probePath)
            #convert to integer array
            if sizeBool == True:
                self.probeMap = self.convertByteToInt(byteArr, 1280, 720)
            else:
                self.probeMap = self.convertByteToInt(byteArr)
        elif probePath.endswith('.xml'):
            byteArr, hstart, hstop, vstart, vstop, self.date = self.loadXMLimage(probePath)
            self.x_pos = hstart-1
            self.y_pos = vstart -1
            #hstop is including hstop pixel
            self.size[0] = abs(hstop - hstart)+1
            self.size[1] = abs(vstop - vstart)+1
            self.probeMap = self.convertByteToInt(byteArr, self.size[0], self.size[1])
        
        
        self.plotMap(1, 2, False)
        
    def plainFit(self, id, width = None):
        '''
        id = 1 for probe / 2 for pump\n
        width is width of average used; default = 1
        '''
        if width == None:
            width = self.width
        #id corresponds to 1 (probe) or 2 (pump)
        #plainFit only does horizontal and vertical fitting
        #returns gaussian fit parameters
        if id == 1:
            current_map = self.probeMap - self.bgMap
        elif id == 2:
            if (hasattr(self, 'pumpMap')):
                current_map = self.pumpMap
            else:
                raise RuntimeError("No pumpMap loaded")
        else:
            raise RuntimeError("id choice for plainFit is 1 (probe) or 2 (pump)")

        #find maximum value position ignoring values under 10% threshold
        thresh_map = np.array(current_map)
        thresh_map[thresh_map < 0.08*4095] = 0
        #horizontal slice (sum over y)
        horX_slice = np.sum(thresh_map, axis = 0)
        #vertical slice (sum over x)
        verY_slice = np.sum(thresh_map, axis = 1)

        if width > 1:
            #find center positions

            x_Peak, x_heights = signal.find_peaks(horX_slice, height = (None, None))
            x_ind = np.argmax(x_heights["peak_heights"])
            x_Peak = x_Peak[x_ind]

            y_Peak, y_heights = signal.find_peaks(verY_slice, height = (None, None))
            y_ind = np.argmax(y_heights["peak_heights"])
            y_Peak = y_Peak[y_ind]


            
            if width%2 == 1:
                w_modifier = (width-1)/2
            else:
                w_modifier = width/2

            x_Peak_slice = np.arange(x_Peak - w_modifier, x_Peak + w_modifier, 1, dtype=int)
            y_Peak_slice = np.arange(y_Peak - w_modifier, y_Peak + w_modifier, 1, dtype=int)
            #get adjusted slices (hopefully correct slicing here)
            horX_slice = np.sum(current_map[y_Peak_slice, :], axis = 0)/width
            verY_slice = np.sum(current_map[:, x_Peak_slice], axis = 1)/width
        else:
            #single width slice
            x_Peak = np.argmax(horX_slice)
            y_Peak = np.argmax(verY_slice)
            #get adjusted slices (hopefully correct slicing here)
            horX_slice = current_map[y_Peak, :]
            verY_slice = current_map[:, x_Peak]




        #do the final fit
        #gaussBounds = [(0, self.size), (0, 3*self.size), (0, np.inf), (-np.inf, np.inf)]
        gaussBounds = ([0, 0, 0, -np.inf], [max(self.size), 3*max(self.size), np.inf, np.inf])
        #x
        #x_center, sig, Amplitude, offset
        xMean = np.mean(horX_slice)
        p0x = [x_Peak, 10, max(horX_slice)-xMean, xMean]
        #setting these points zero to avoid problems
        horX_slice[horX_slice<0] = 0
        [pOutX, covX] = optimize.curve_fit(self.gaussianSplit, np.arange(0, self.size[0], 1), horX_slice[:], p0x, bounds = gaussBounds)
        
        #y
        #y_center, sig, Amplitude, offset
        yMean = np.mean(verY_slice)
        p0y = [y_Peak, 10, max(verY_slice)-yMean, yMean]
        #setting these points zero to avoid problems
        verY_slice[verY_slice<0] = 0
        [pOutY, covY] = optimize.curve_fit(self.gaussianSplit, np.arange(0, self.size[1], 1), verY_slice[:], p0y, bounds = gaussBounds)

        if id == 1:
            self.XprobeCoeffs = pOutX
            self.YprobeCoeffs = pOutY
            self.XprobeCov = covX
            self.YprobeCov = covY
        elif id == 2:
            self.XpumpCoeffs = pOutX
            self.YpumpCoeffs = pOutY
            self.XpumpCov = covX
            self.YpumpCov = covY


        #also return directions, so the gauss plotting can be streamlined

    

    def plotGauss(self, printData = False, background = False):
        '''
        Plots Gaussian fits with contour levels to show how well the fit is\n
        Needs to already have the fit values via one of the fitting functions (plainFit)\n
        use printData = True to have fit data printed
        background = True subtracts background to show how well the subtractions is working
        '''
        #expect there at least to be the probe one (aka the first one)
        #got to do a translation for the axis direction (rotation)

        #number of contour lines
        n_levels = 4
        #contour points
        cont_points = 500
        
        #create figure
        fig, axs = plt.subplots(2, 2)
        #plot probe
        #x slice
        #setting scale 1 again, no point in not doing px maps I think
        #scale = self.getScale()
        scale = 1
        x = np.arange(0,scale*self.size[0],0.1)
        if hasattr(self, 'XprobeCoeffs'):
            axs[0,0].plot(x, self.gaussian(x, self.XprobeCoeffs))
        y = np.arange(0,scale*self.size[1],0.1)
        if hasattr(self, 'YprobeCoeffs'):
            axs[1,0].plot(y, self.gaussian(y, self.YprobeCoeffs))

        #check if coefficients for pump are there
        if hasattr(self, 'XpumpCoeffs'):
            axs[0,0].plot(x, self.gaussian(x, self.XpumpCoeffs), label = "2nd")
        if hasattr(self, 'YpumpCoeffs'):
            axs[1,0].plot(x, self.gaussian(y, self.YpumpCoeffs), label = "2nd")

        
        #plot equal height 
        #1st beam
        if background == False:
            axs[0,1].imshow(self.probeMap, vmin = 0, vmax = 4095)
        if background == True:
            axs[0,1].imshow(self.probeMap-self.bgMap, vmin = 0, vmax = 4095)
        

        #both need to be fitted to allow for this to plot
        if hasattr(self, 'XprobeCoeffs') and hasattr(self, 'YprobeCoeffs'):
            self.contourPlotter(axs[0,1], self.XprobeCoeffs, self.YprobeCoeffs, 4)
        

        #2nd beam
        if hasattr(self, 'pumpMap'):
            if background == False:
                axs[1,1].imshow(self.pumpMap, vmin = 0, vmax = 4095)
            elif background == True:
                axs[1,1].imshow(self.pumpMap-self.bgMap, vmin = 0, vmax = 4095)

        if hasattr(self, 'XpumpCoeffs') and hasattr(self, 'YpumpCoeffs'):
            self.contourPlotter(axs[1,1], self.XpumpCoeffs, self.YpumpCoeffs, 4)


        #titles
        axs[0,0].set_title('x-slice')
        axs[1,0].set_title('y-slice')
        axs[0,1].set_title('1st Beam contour')
        axs[1,1].set_title('2nd Beam contour')

        #print Fit Data
        self.readableFitData(printData)
        #not sure what info to put here yet
        plt.tight_layout()
        plt.show()
    
    def plotComp(self):
        '''shows overlap and gives overlap values\n
        should show maximum intensities as well'''
        #top left text only
        #top right is contour plots
        #bottom left is 1st (x) axis, right is 2nd axis (y)
        fig, ax = plt.subplots(2,2)
        plt.tight_layout()
        xaxis_vector = np.linspace(0,self.size[0], 1000)
        yaxis_vector = np.linspace(0,self.size[1], 1000)
        xOverlap = 0
        yOverlap = 0
        xOverlapCorr = 0
        yOverlapCorr = 0
        if hasattr(self, "XprobeCoeffs") and hasattr(self, "XpumpCoeffs"):
            xOverlapCorr = ArtFit.calculateOverlapCorrection(self.XpumpCoeffs, self.XprobeCoeffs)
            xOverlap = ArtFit.calculateOverlap(self.XpumpCoeffs, self.XprobeCoeffs)
            #plot xAxis as PDF
            ax[1,0].set_title("First Axis")
            XprobeCoeffs = self.XprobeCoeffs
            XprobeCoeffs[2:4] = [1,0]
            XpumpCoeffs = self.XpumpCoeffs
            XpumpCoeffs[2:4] = [1,0]
            ax[1,0].plot(xaxis_vector, self.gaussian(xaxis_vector, XprobeCoeffs))
            ax[1,0].plot(xaxis_vector, self.gaussian(xaxis_vector, XpumpCoeffs))
            
            
        if hasattr(self, "YprobeCoeffs") and hasattr(self, "YpumpCoeffs"):
            yOverlapCorr = ArtFit.calculateOverlapCorrection(self.YprobeCoeffs, self.YpumpCoeffs)
            yOverlap = ArtFit.calculateOverlap(self.YpumpCoeffs, self.YprobeCoeffs)
            #plot yAxis as PDF
            ax[1,1].set_title("Second Axis")
            YprobeCoeffs = self.YprobeCoeffs
            YprobeCoeffs[2:4] = [1,0]
            YpumpCoeffs = self.YpumpCoeffs
            YpumpCoeffs[2:4] = [1,0]
            ax[1,1].plot(yaxis_vector, self.gaussian(yaxis_vector, YprobeCoeffs))
            ax[1,1].plot(yaxis_vector, self.gaussian(yaxis_vector, YpumpCoeffs))


        #text only
        ax[0,0].set_axis_off()
        overlap_text = "X-Overlap is: %.4f\nY-Overlap is: %.4f\nTotal Overlap: %.4f\n" %(xOverlap, yOverlap, xOverlap*yOverlap)
        overlap_text += "X-Correction is: %.4f\nY-Correction is: %.4f\nTotal Correction is: %.4f" %(xOverlapCorr, yOverlapCorr, xOverlapCorr*yOverlapCorr)


        ax[0,0].text(0, 0.5, overlap_text)
        #contour plots
        ax[0,1].set_title("Beam shapes")
        #probe
        if hasattr(self, "XprobeCoeffs") and hasattr(self, "YprobeCoeffs"):
            self.contourPlotter(ax[0,1], self.XprobeCoeffs, self.YprobeCoeffs, 4, colormap = "Purples")
        #pump
        if hasattr(self, "XpumpCoeffs") and hasattr(self, "YpumpCoeffs"):
            self.contourPlotter(ax[0,1], self.XpumpCoeffs, self.YpumpCoeffs, 4, colormap = "Greens")


    def getNumericCorrection(self, threshold = 0.04, interpDensity = 32):
        #check if corresponding maps are there
        if hasattr(self, "pumpMap")  and hasattr(self, "probeMap"):
            #actually do this via the fit and a sigma related distance
            #check if the fits exist and if not try doing them (too lazy to implement a check for this
            #should be obvious if it doesn't work)
            if not hasattr(self, "XprobeCoeffs") or not hasattr(self, "YprobeCoeffs"):
                self.plainFit(1)
            if not hasattr(self, "XpumpCoeffs") or not hasattr(self, "YpumpCoeffs"):
                self.plainFit(2)

            cutPumpMap = self.pumpMap - self.bgMap
            cutProbeMap = self.probeMap - self.bgMap
            
            cutPumpMap[self.pumpMap < threshold*4095] = 0
            cutProbeMap[self.probeMap < threshold*4095] = 0
            

            #get pixels within "sigfactor" sigma distance of center
            #mainly done to remove noise
            cutProbeMap = self.cutMapSig(cutProbeMap, self.XpumpCoeffs, self.YpumpCoeffs, sigfactor = 3.5)
            cutPumpMap = self.cutMapSig(cutPumpMap, self.XpumpCoeffs, self.YpumpCoeffs, sigfactor = 3.5)
            
            #do interpolation 
            #need to make a map of the coordinates
            x_vec = np.arange(0, self.size[0],1)
            y_vec = np.arange(0, self.size[1],1)
            from scipy.interpolate import RegularGridInterpolator as Interpol
            cutProbeMapInterp = Interpol((x_vec, y_vec), cutProbeMap)
            cutPumpMapInterp = Interpol((x_vec, y_vec), cutPumpMap)
            
            #adding and subtracting one is a cheap trick to not have to think about this more
            x_interp = np.linspace(1, self.size[0]-1, int(self.size[0]/int(max([self.XpumpCoeffs[1], self.XprobeCoeffs[1]]))*interpDensity))
            y_interp = np.linspace(1, self.size[1]-1, int(self.size[1]/int(max([self.YpumpCoeffs[1], self.YprobeCoeffs[1]]))*interpDensity))
            x_interp, y_interp = np.meshgrid(x_interp, y_interp)
            x_interp = np.reshape(x_interp,(-1))
            y_interp = np.reshape(y_interp, (-1))
            full_interpol = list(zip(x_interp, y_interp))

            #basically calculate the same factor as for the normal correction factor, just with discrete pixels
            denominator = np.sum(np.multiply(cutProbeMapInterp(full_interpol), cutPumpMapInterp(full_interpol))[:])
            constantPumpIntensity = np.max(cutPumpMapInterp(full_interpol))
            enumerator = np.sum(constantPumpIntensity*cutProbeMapInterp(full_interpol)[:])
            #this is missing various constant factors, since this is not using the gaussian shape, but is discretized
            correction = enumerator/denominator
            
            #plot to check if cutMap makes sense
            fig, axs = plt.subplots(2,1)
            #probe 
            axs[0].imshow(cutProbeMap, vmin = 0, vmax = 4095)
            self.contourPlotter(axs[0], self.XprobeCoeffs, self.YprobeCoeffs, 6)
            #pump
            axs[1].imshow(cutPumpMap, vmin = 0, vmax = 4095)
            self.contourPlotter(axs[1], self.XpumpCoeffs, self.YpumpCoeffs, 6)
            #print(correction)
            plt.show()
            return correction
        

    @classmethod
    def cutMapSig(self, inputMap, XCoeffs, YCoeffs, sigfactor=3, thresh = 0.04):
        '''sets any value farther than 3sig from center point to 0 and returns cut map\\
            using the larger sigma value to define the radius\\
            setting any value under maximum value (4095)*thresh to 0'''
        inputMap[inputMap < thresh] = 0

        sigX = XCoeffs[1]
        muX = XCoeffs[0]
        sigY = YCoeffs[1]
        muY = YCoeffs[0]
        #max of the two is radius
        radius = max([sigX, sigY])
        for x_index in range(np.shape(inputMap)[0]):
            for y_index in range(np.shape(inputMap)[1]):
                if (x_index-muX)**2 + (y_index-muY)**2 > (sigfactor*radius)**2:
                    inputMap[y_index, x_index] = 0
        return inputMap



    def contourPlotter(self, axis_handle, coeffs1, coeffs2, lvl = 6, x_size = 0, y_size = 0, res = 0, colormap = "viridis"):
        '''plots the contour plot with lvl = number of equipotential lines'''
        if x_size == 0:
            x_size = self.size[0]
        if y_size == 0:
            y_size = self.size[1]
        if res == 0:
            res = 5

        #scale = self.getScale()
        #not showing this in px does not make much sense I think
        scale = 1

            
        X_contour = np.linspace(0,self.size[0]*scale, res*x_size)
        Y_contour = np.linspace(0, self.size[1]*scale, res*y_size)
        [X_contour, Y_contour] = np.meshgrid(X_contour, Y_contour)
        Z_contour = self.gaussian(X_contour, coeffs1)*self.gaussian(Y_contour, coeffs2)

        #adjust for correct height
        Z_contour = Z_contour/np.sqrt(coeffs1[2]*coeffs2[2])
        axis_handle.contour(X_contour, Y_contour, Z_contour, levels = lvl, cmap = colormap)
        axis_handle.set_aspect('equal')


    
    
    def titleFromFit(self, id):
        '''
        deprecated for now'''
        RuntimeWarning("titleFromFit is deprecated for now; use readableFitData instead")
        if id == 1:
            coeffs = self.XprobeCoeffs
            ycoeffs = self.YprobeCoeffs
            cov = self.XprobeCov
        FWHM_px = self.sigToFWHM(coeffs[1])
        sigFWHM_px = self.sigToFWHM(np.sqrt(cov[1,1]))
        titleString = "center: %.1d \u00B1 %.1e \nFWHM: (%.2d \u00B1 %.1e )px / %.1d \u00B1 %.1e µm\n" %(coeffs[0], np.sqrt(cov[0,0]), FWHM_px, sigFWHM_px, FWHM_px*self.scale, sigFWHM_px*self.scale)
        return titleString
    
    def readableFitData(self, printBool=False):
        '''
        aim is to have readable fitdata for either printing to console or file\n
        using pandas.DataFrame for that'''
        #iterate through possible data inputs and use bool array to find names

        scale, scale_bool = self.getScale(True)
        if scale_bool == False:
            scale = np.NaN
        
        listOfDicts = list()
        arrayOfNames = np.array(["1st-x", "1st-x Error", "2nd-x", "2nd-x Error", "1st-y", "1st-y Error", "2nd-y", "2nd-y Error"])
        bool_array = np.zeros(np.shape(arrayOfNames), dtype = bool)
        #probe x
        if hasattr(self, 'XprobeCoeffs'):
            listOfDicts.append(self.coeffsToDict(self.XprobeCoeffs, scale))
            listOfDicts.append(self.coeffsToDict(self.XprobeCov, scale))
            bool_array[0:2] = True

        #pump x
        if hasattr(self, 'XpumpCoeffs'):
            listOfDicts.append(self.coeffsToDict(self.XpumpCoeffs, scale))
            listOfDicts.append(self.coeffsToDict(self.XpumpCov, scale))
            bool_array[2:4] = True

        #probe y
        if hasattr(self, 'YprobeCoeffs'):
            listOfDicts.append(self.coeffsToDict(self.YprobeCoeffs, scale))
            listOfDicts.append(self.coeffsToDict(self.YprobeCov, scale))
            bool_array[4:6] = True

        #pump y
        if hasattr(self, 'YpumpCoeffs'):
            listOfDicts.append(self.coeffsToDict(self.YpumpCoeffs, scale))
            listOfDicts.append(self.coeffsToDict(self.YpumpCov, scale))
            bool_array[6:8] = True
        
        #set instead of bool
        x_correction = np.NaN
        y_correction = np.NaN
        
        #overlaps
        if hasattr(self, 'XprobeCoeffs') and hasattr(self, 'XpumpCoeffs'):
            x_correction = ArtFit.calculateOverlapCorrection(self.XpumpCoeffs, self.XprobeCoeffs)
        if hasattr(self, 'YprobeCoeffs') and hasattr(self, 'YpumpCoeffs'):
            y_correction = ArtFit.calculateOverlapCorrection(self.YpumpCoeffs, self.YprobeCoeffs)
        
        correction = {'x': x_correction, 'y': y_correction}
        #generate dataframe and then introduce a new column with the names
        fitData = DataFrame(listOfDicts)
        #cut down list of names and insert
        fitData.insert(0, "Names", arrayOfNames[bool_array])
        if printBool == True:
            print(fitData)
            print(correction)
        return fitData, correction
    

    @classmethod
    @staticmethod
    def calculateOverlapCorrection(pumpCoeff, probeCoeff, cutoff = 3):
        '''calculates correction factor for overlap of two sets of coefficients (pump and probe) and returns a single value with which to multiply the dOD value\n
        best possible overlap coefficient is 1: when probe is much smaller than pump and entire probe feels a constant pump intensity\n
        if the scaling factor is larger than a arbitrarily chosen factor NaN is returned
        analytical; not identical to Artray in C# where a sum is done'''
        #"coeffs = x_center, sig, Amplitude, offset"
        sig_squared = pumpCoeff[1]**2 + probeCoeff[1]**2
        c_Scaling = np.sqrt(sig_squared/pumpCoeff[1]**2)*np.exp((pumpCoeff[0]-probeCoeff[0])**2/(2*sig_squared))
        #this cutoff may need correcting to a lower value
        if c_Scaling > cutoff:
            c_Scaling = np.nan
        return c_Scaling
    
    @classmethod
    @staticmethod
    def parseFitFile(fitPath):
        '''Parses fitfile into coefficient arrays in a dict as used within ArtFit\\
        only use with strict probe (1st) pump (2nd) files for now \\
        coeffs = x_center, sig, Amplitude, offset
        '''
        #Names;µ / px;FWHM / px;Amplitude;Offset;µm/px
        #names given as "1st-x" and similar, grab first letter for number and last letter for x or y to parse
        fitnames = np.genfromtxt(fitPath, delimiter = ";", dtype=str, usecols=(0), skip_header = 1)
        fitfile = np.genfromtxt(fitPath, delimiter = ";", dtype=float, usecols=(1,2,3,4), skip_header = 1)
        

        #make a tuple array
        iterationTuple = list()
        beam_unique = []
        for row in range(0,int(np.shape(fitfile)[0]),2):
            if fitnames[row][0] == str(1): 
                pumpOrProbe = "probe"
            elif fitnames[row][0] == str(2):
                pumpOrProbe = "pump"

            xory = fitnames[row][-1]
            iterationTuple.append((row, pumpOrProbe, xory))
            if pumpOrProbe not in beam_unique:
                beam_unique.append(pumpOrProbe)
        
        #read names from first column and put into dictionary
        sorter = lambda array: array[:][1]
        iterationTuple.sort(key=sorter)
        #currently only if both probe and pump are given
        templateDict = lambda: {"x": np.zeros((4,2), dtype = float), "y": np.zeros((4,2), dtype = float)}
        fitdata = { "probe" : templateDict(), "pump": templateDict()}

        for tuple_loop in iterationTuple:
            fitdata[tuple_loop[1]][tuple_loop[2]][:,0] = fitfile[tuple_loop[0]]*np.array([1,1/(2*np.sqrt(2*np.log(2))),1,1])
            fitdata[tuple_loop[1]][tuple_loop[2]][:,1] = fitfile[tuple_loop[0]+1]*np.array([1,1/(2*np.sqrt(2*np.log(2))),1,1])
        return fitdata
        

    @classmethod
    @staticmethod
    def calculateOverlap(pumpCoeff, probeCoeff):
        '''calculates sum Overlap of two sets of coefficients (pump and probe)\n
        best possible overlap coefficient is 1 for when both coefficients are identical
        analytical; resembles Artray overlap(not identical to Artray in C# where a sum is done)'''
        #yet to be tested
        sig1 = pumpCoeff[1]
        sig2 = probeCoeff[1]
        sig_squared = sig1**2 + sig2**2
        return np.sqrt((2*sig1*sig2)/sig_squared)*np.exp(-(pumpCoeff[0]-probeCoeff[0])**2/(4*sig_squared))

    


    @staticmethod
    def coeffsToDict(coeffs, scale = np.NaN):
        '''coeffs may either be straight up coefficients or may be a nxn array of a cov matrix\n
        returns a dictionary with corresponding outputs name'''
        "coeffs = x_center, sig, Amplitude, offset"
        if np.shape(coeffs)[0] != np.size(coeffs):
            coeffs = np.sqrt(np.diag(coeffs))

        dict_out = {
            "µ / px": coeffs[0],
            "FWHM / px": ArtFit.sigToFWHM(coeffs[1]),
            "Amplitude": coeffs[2],
            "Offset": coeffs[3], 
            "µm/px": scale}
        
        return dict_out

        
    @staticmethod
    def loadByteArr(path):
        #not sure about the encoding
        print(path)
        file = open(path, 'rb')
        data = bytearray(file.read())
        return data
    
    @classmethod
    @staticmethod
    def loadXMLimage(path):
        #using Base64 encoding for raw data
        
        #print(event, node)
        document = minidom.parse(path)
        #I give up, it is time for simple brute force
        day = document.getElementsByTagName("day")[0].firstChild.nodeValue
        month = document.getElementsByTagName("month")[0].firstChild.nodeValue
        year = document.getElementsByTagName("year")[0].firstChild.nodeValue
        date = year+month+day
        #path(date)
        hstart = document.getElementsByTagName("hstart")[0].firstChild.nodeValue
        hstop = document.getElementsByTagName("hstop")[0].firstChild.nodeValue
        vstart = document.getElementsByTagName("vstart")[0].firstChild.nodeValue
        vstop = document.getElementsByTagName("vstop")[0].firstChild.nodeValue

        #now for the important thing:
        rawDatEncoded = document.getElementsByTagName("data")[0].firstChild.nodeValue
        print("rawDat from XML")
        #print(rawDatEncoded)
        binaryMap = base64.b64decode(rawDatEncoded, validate = True)
        return binaryMap, int(hstart), int(hstop), int(vstart), int(vstop), date



    
    @staticmethod
    def loadTIFF(pathImage):
        '''Tiff path to map and returns dimensions'''
        #importing here since never used usually
        from PIL import Image
        print(pathImage)
        tempImg = Image.open(pathImage)
        
        width, height = tempImg.size
        #print(width/height)
        mapImg = np.array(tempImg)
        #print(np.shape(mapImg))
        #plt.imshow(mapImg)
        #axes seem to be correct (comparing to loaded image)

        return mapImg, [width,height]


    def convertByteToInt(self, byte_arr, x_width = 0, y_height = 0):
        '''
        takes bytearray and optionally x_width and y_height dimensions\n
        returns intMap'''
        #checking for size consistency here
        #not sure why the warning gets printed for bg here, since it doesn't go into the if condition
        #print(x_width, y_height)
        #print((x_width == 0 and y_height == 0))
        if (x_width == 0 and y_height == 0):
            #print(x_width, y_height)
            #division by six (? do I need to do 8? aka 2*4(RGBA); apparently 6 is fine and A is not saved?)
            size = np.ones(2, dtype=int)*int(np.sqrt(len(byte_arr)/6))
            #np.any in the end so that it doesn't trigger on the first definition
            if (hasattr(self, 'size') and np.any(size != self.size) and not np.any(self.size)):
                warnings.warn("Dimensions of maps may not match up")
            if self.size[0] == 0 and self.size[1] == 0:
                self.size = size
            x_width = self.size[0]
            y_height = self.size[1]

        intMap = np.ndarray((y_height, x_width))
        #print(y_height, x_width)
        #print(np.size(intMap))
        #print(np.size(byte_arr))
        for y_indx in range(y_height):
            for x_indx in range(x_width):
                ind = x_indx + y_indx*x_width
    
                #fixed a flip in vertical axis, doesn't change anything but locations

                intMap[y_height-y_indx-1, x_indx] = ((byte_arr[(3*ind*2+1)] & 0b1111) << 8) | (byte_arr[3*ind*2])
        #plt.figure()
        #print(np.shape(intMap))
        #plt.imshow(intMap)
        return intMap
    
    def setScale(self, scale):
        '''setter for scale in µm/px (aka 10E-6 m / px unit)'''
        self.scale = scale



    def getScale(self, bool = False):
        if np.isnan(self.scale):
            scale = 1
            bool_out = False
        else:
            scale = self.scale
            bool_out = True
        if bool == True:
            return scale, bool_out
        else:
            return scale

    @classmethod
    @staticmethod
    def backgroundJSON_Gen(pathToSave, path, x_pos, y_pos, scale, date="?"):
        json_template = {
        "path": path,
        "date": date,
        "x_pos": x_pos,
        "y_pos": y_pos,
        "scale: µm/px": scale
        }
        file = open(pathToSave, 'x', encoding='utf-8')
        json.dump(json_template, file, ensure_ascii=False, indent = 4)


    def getJSON(self, path):
        '''takes background JSON path'''
        f_JSON = open(path)
        data_JSON = json.load(f_JSON)
        if not hasattr(self, "x_pos"):
            self.x_pos = data_JSON['x_pos']
        if not hasattr(self, "y_pos"):
            self.y_pos = data_JSON['y_pos']
        self.setScale(data_JSON["scale: Âµm/px"])
        return data_JSON['path']

    #cuts map (background) to size of input
    def cut_map(self, map):
        #not entirely sure if these are the right dimensions or flipped
        new_map = np.ndarray((self.size[1], self.size[0]))
        for y_indx in range(self.size[1]):
            for x_indx in range(self.size[0]):
                new_map[y_indx, x_indx] = map[self.y_pos + y_indx, self.x_pos + x_indx]
        return new_map

    #dim decides if it is a 2d colour plot or a 3d surface plot
    def plotMap(self, map_number, dim = 2, bg_sub = True, maxVal=True):
        if map_number == 0:
            map = self.bgMap
        elif map_number == 1:
            map = self.probeMap
        elif map_number == 2:
            map = self.pumpMap
        else:
            raise IndexError("map_number for plotMap has to be 0, 1 or 2")

        if bg_sub == True and map_number != 0:
            map = map - self.bgMap
        if maxVal==False:
            maximumValue = 4095
        else:
            maximumValue = np.max(map)
        
        if dim == 2:
            fig, axs = plt.subplots()
            ax_call = axs.imshow(map, vmin = 0, vmax = maximumValue)
        elif dim == 3:
            fig, axs = plt.subplots(subplot_kw={"projection": "3d"})
            x_map = np.arange(0, self.size[0])
            y_map = np.arange(0, self.size[1])
            [x_map, y_map] = np.meshgrid(x_map, y_map)
            fig.properties()
            ax_call = axs.plot_surface(x_map, y_map, map, vmin = 0, vmax = maximumValue)

        plt.colorbar(ax_call)
        plt.show()
        print("Input to end")
        input()
    
    @classmethod
    @staticmethod
    def gaussian(x, coeffs):
        "coeffs = x_center, sig, Amplitude, offset"
        return 1/(np.sqrt(2*np.pi)*coeffs[1])*coeffs[2]*np.exp(-(x-coeffs[0])**2/(2*coeffs[1]**2)) + coeffs[3]
    
    #seems I need this for curve_fit, wonde why I didn't do it like that in the first place (probably some reason for it)
    @classmethod
    @staticmethod
    def gaussianSplit(x, center, sig, Amp, offset):
        "coeffs = x_center, sig, Amplitude, offset"
        return 1/(np.sqrt(2*np.pi)*sig)*Amp*np.exp(-(x-center)**2/(2*sig**2)) + offset
    
    @classmethod
    @staticmethod
    def sigToFWHM(sig):
        #return 2*np.sqrt(2*np)
        return 2.354820045*sig

    def setPumpPower(self, pumpPower):
        '''Set pump power in W\n Setting probe power does not make sense as it has no influence on the absorbance (as long as it is low enough to not lead to non-linear power effects)'''
        self.pumpPower = pumpPower
    def setProbePower(self, probePower):
        self.probePower = probePower



    def getCorrectionFactor(self, limits = 3):
        '''Get correction factor as dictionary considering power density and overlap\n
        Normed due to powerdensity inclusion\n
        returns correctly written value for all combinations'''

        #check if a correctionfactor is to be calculated (meaning are there two beamshapes input?)
        if self.status_identifier != 2:
            Warning("No correction factor to be calculated. Only one beam input")
            return
        powerDensity, bool_scale = self.getPowerDensity(2)

        if hasattr(self, "XprobeCoeffs") and hasattr(self, "XpumpCoeffs"):
            x_correction = self.calculateOverlapCorrection(self.XpumpCoeffs, self.XprobeCoeffs, limits)
            if np.isnan(x_correction):
                Exception("getCorrectionFactor for X exceeds set limit %.1e" %limits)

        if hasattr(self, "YprobeCoeffs") and hasattr(self, "YpumpCoeffs"):
            y_correction = self.calculateOverlapCorrection(self.YpumpCoeffs, self.YprobeCoeffs, limits)
            if np.isnan(x_correction):
                Exception("getCorrectionFactor for Y exceeds set limit %.1e" %limits)
        
        if not(np.isnan(x_correction) and np.isnan(y_correction)):
            #divide or multiply with powerdensity?
            if bool_scale == False: 
                print("Warning")
                Warning("Power density per px, no scale given.\n powerDensity = %.2e W/px^2" %powerDensity)
            return x_correction*y_correction/powerDensity
        else:
            Warning("Incomplete values for CorrectionFactor. Not outputting anything for now.")
            return
        


        
    def getPowerDensity(self, pulse_id):
        '''pulse_id: 1 for probe; 2 for pump\n returns W/m^2 and scalebool\\(I am sticking with metric)        '''
        # according to liveFit
        # peakRadiance = power / (pGaussX(3) * pGaussY(3)) / (umPerPixel * 1e-4) ^ 2 / (2 * pi); (for radiance in W/cm^2)
        # 1e-6 m = 1 µm
        # 1e-2 m = 1 cm = 1e4 µm
        

        try:
            if (pulse_id == 1):
                sigX = self.XprobeCoeffs[1]
                sigY = self.YprobeCoeffs[1]
                power = self.probePower
            elif (pulse_id == 2):
                sigX = self.XpumpCoeffs[1]
                sigY = self.YpumpCoeffs[1]
                power = self.pumpPower
            else:
                UserWarning("Only pulse_id 1 and 2 allowed for getPowerDensity")
        except:
            UserWarning("Something went wrong in getPowerDensity\n Likely that there is no fit for the chosen pulse_id %d" %pulse_id)
        umPerPx, scaleBool = self.getScale(True)
        
        try: 
            if (pulse_id == 1):
                power = self.probePower
            else:
                power = self.pumpPower
        except:
            UserWarning("No power given for the chosen pulse_id %d.\nSetting power = 1 µW" %pulse_id)
            power = 1e-6
        
        peakRadiance = power / (sigX*sigY*umPerPx**2*1e-12*2*np.pi)

        outputString = ''
        outputString += "Peak radiance of "
        if pulse_id == 1:
            outputString += "Probe: "
        elif pulse_id == 2:
            outputString +="Pump: "
        outputString += "%.2e" %peakRadiance
        if scaleBool == True:
            outputString += "W/m^2\nat %.2e W" %power
        else:
            outputString += "W/px^2\nat %.2e W" %power

        #print(outputString)

        return peakRadiance, scaleBool
    

    @classmethod
    def rerunFromJSON(self, JSONpath, overwrite = False, showPlots = False, bgSubtract = False, maxCorr = 3):
        '''currently only for full pump and probe files, implement only single beam later
        currently always uses backgroundfile scale'''
        print(JSONpath)
        f_JSON = open(JSONpath, 'r')
        summary_JSON = json.load(f_JSON)
        f_JSON.close()
        Routine = ArtFit(summary_JSON["BackgroundJSON"], summary_JSON["pathProbe"], summary_JSON["pathPump"], summary_JSON["FitAvgWidth / px"])
        Routine.setPumpPower(summary_JSON["Pump power: W"])
        #Routine.setScale(summary_JSON["scale: Âµm/px"])

        Routine.plainFit(1)
        Routine.plainFit(2)

        if showPlots == True:
            Routine.plotComp()
            Routine.plotGauss(background=bgSubtract)
        
        print("old corr: %.3e \nnew corr: %.3e" %(summary_JSON["Correction Factor: per W/m^2 pump power"], Routine.getCorrectionFactor(5)))
        print("old power density: %.3e \nnew power density: %.3e" %(summary_JSON["Pump power density: W/m^2"], Routine.getPowerDensity(2)[0]))

        print("old fit values:\n")
        oldFit = pandas.read_csv(summary_JSON["FittingPath"], delimiter = ";")
        print(oldFit)
        print("new fit values:\n")
        Routine.readableFitData(True)
        if overwrite == False:
            #I think this option turned out more annoying than anything
            #owBool = input("Do you want to overwrite? yes/no? \n")
            owBool=False
            if owBool == "yes" or owBool == "y":
                overwrite = True
            else:
                print("Not saving")
        if overwrite == True:
            owBool = input("Do you want to overwrite? yes/no? \n")
            if owBool == "yes" or owBool == "y":
                print("Saving")
                Routine.readableFitData()[0].to_csv(summary_JSON["FittingPath"], sep = ";", index = False, mode = 'w', encoding="utf-8", decimal = '.')
                summary_JSON["scale: Âµm/px"] = Routine.getScale()
                summary_JSON["Correction Factor: per W/m^2 pump power"] = Routine.getCorrectionFactor(maxCorr)
                summary_JSON["Pump power density: W/m^2"] = Routine.getPowerDensity(2)[0]
                f_JSON = open(JSONpath, "w")
                json.dump(summary_JSON, f_JSON, ensure_ascii=False, indent = 4)
                f_JSON.truncate()
                f_JSON.close()
            else:
                print("Not saving")
        return Routine


        

        

        
#JSON_summary = r"C:\Users\Maximilian\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm\Pump653Probe493_ArtrayCompensationControl\SummaryPump653Probe493Caution.JSON"
#ArtFit.rerunFromJSON(JSON_summary, True, True, True)


r"""
Pump = r"C:\Users\Maximilian\Documents\Artray-notsynced-test_RAWS\31.10.2023\Pump653forProbe380.txt"
Probe = r"C:\Users\Maximilian\Documents\Artray-notsynced-test_RAWS\31.10.2023\Probe380.txt"
bgJSON = r"C:\Users\Maximilian\Documents\Artray-notsynced-test_RAWS\BackgroundJSON-24.10.2023"

#currentFit = ArtFit(bgJSON, Probe, 1)
currentFit = ArtFit(bgJSON, Probe, Pump, 1)
#currentFit.plotMap(1, dim=2, bg_sub = True)
currentFit.plainFit(1, 2)
currentFit.plainFit(2, 2)
currentFit.setScale(4)

#currentFit.plotComp()
#currentFit.plotGauss()
currentFit.setPumpPower(7.5e-6)
print("Correction factor %.3e per W/m^2" %currentFit.getCorrectionFactor())

bg_map = r"C:\Users\Maximilian\Documents\CameraData_07.12.2023\Artray\ArtrayBackground_10.10.2023_exp15ms.txt"
save_path = r"C:\Users\Maximilian\Documents\CameraData_07.12.2023\Artray\JSON_ArtrayBG_27.10.2023.txt"
ArtFit.backgroundJSON_Gen(save_path, bg_map, 780, 100, 3.5228, "27.10.2023")

#bgPlot = ArtFit(bg_map, True)

"""

