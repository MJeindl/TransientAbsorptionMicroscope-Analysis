from scipy.io import loadmat
from scipy.optimize import curve_fit, OptimizeWarning
import numpy as np
import ExponentialFunctions
import os
import matplotlib.pyplot as plt
import json
import sys
sys.path.insert(1, r"C:\Users\M\Documents\Books\masterprojectinformation")
import SeriesDegradation
import fileParsingMethods

class DelayScan:
    def __init__(self, pump="", probe="", filename="") -> None:
        self.pump = pump
        self.probe = probe
        self.filename = filename
        self.startParameters = []
        
    @classmethod
    @staticmethod
    def loadFromJson(jsonDict: dict):
        scan = DelayScan()
        scan.pump = jsonDict["pump"]
        scan.probe = jsonDict["probe"]
        scan.filename = jsonDict["inputFile"]
        if "wavelength" in jsonDict.keys():
            scan.wavelength = jsonDict["wavelength"]
        if "popt" in jsonDict.keys():
            scan.popt = np.array(jsonDict["popt"])
            scan.pcov = np.array(jsonDict["pcov"])
        if "startParameters" in jsonDict.keys():
            scan.startParameters = jsonDict["startParameters"]
        #check for background parameters indicating linear background subtraction
        #implementation not yet tested
        if "bgParam" in jsonDict.keys():
            scan.backgroundParameters = jsonDict["bgParam"]
        #not yet sure how to make the call for these
        if "tauExpDegrade" in jsonDict.keys():
            scan.degExpTau = jsonDict["tauExpDegrade"]
        if "degradePowerRatio" in jsonDict.keys():
            scan.degPowerRatio = jsonDict["degradePowerRatio"]


        return scan

    def readData(self, degradeBool = False):
        if degradeBool != True:
            #not sure if this path even works tbh, but I am not testing it now
            data = loadmat(self.filename)
            if "filenames" in data.keys():
                # We have a saturation type input, average overall measurements
                # Careful: Does not handle different measurement series (e.g. polarization or wavelength scans)
                filenames = data["filenames"][0]

                #to normalise values
                meas_number = len(filenames)

                for filename in filenames:
                    if os.path.isfile(os.path.dirname(self.filename)+"/" + filename[0] + ".mat") == False:
                        #correct normalisation
                        meas_number -= 1
                        continue
                        
                    currentData = loadmat(os.path.dirname(self.filename)+"/" + filename[0] + ".mat")
                    if not hasattr(self, "delay"):
                        self.delay = currentData["delay"][0, :]
                        self.T = currentData["d_vec"][0, :]
                    else:
                        # Assume we have equal lengths (should be valid for all series)
                        #self.delay += currentData["delay"][0, :]
                        self.T += currentData["d_vec"][0, :]
                
                #check if there are any measurements at all
                if meas_number > 0:
                    if meas_number < len(filenames):
                        print("Not all TA files for %s", self.filename)
                else:
                    print("No files for %s", self.filename)
                #This one is nonesense
                #self.delay /= meas_number

                self.T /= meas_number
                if hasattr(self, "backgroundParameters"):
                    self.removeBackgroundwithParameters()
                else:
                    self.removeBackground()
            else:
                self.delay = data["delay"][0, :]
                self.T = data["d_vec"][0, :]
            self.A = -1000 * np.log10(self.T)
        else:
            #get points in time
            deltaTimes = fileParsingMethods.parseTime(fileParsingMethods.getTimes(self.filename))
            self.T, self.delay, _ = fileParsingMethods.parseSummaryFileToArray(self.filename)
            if hasattr(self, "backgroundParameters"):
                self.removeBackgroundwithParameters()
            else:
                self.removeBackground()
            compensationFactors = SeriesDegradation.degradationCompensation(self.degExpTau, deltaTimes, decaysteps=len(self.delay),powerDensities=self.degPowerRatio)
            #degradation factor does not belong on T
            #this should recalculate the absorbance as well
            self.A = np.mean(self.A*compensationFactors, axis = 0)
            self.T = np.power(10,-self.A/1000)

            #self.T = np.power(10,-self.A/1000)
            
        #no need to return something I think? couldn't find references
        #return data

    def removeBackground(self, numEntries: int = 20) -> float:
        backgroundLevel = np.mean(self.T[:numEntries])
        self.T += 1 - backgroundLevel
        self.A = -1000 * np.log10(self.T)
        return backgroundLevel
    
    def removeBackgroundwithParameters(self) -> None:
        self.T += 1 - np.polyval(self.backgroundParameters, self.delay)
        self.A = -1000 * np.log10(self.T)
        return None

    def hasFits(self) -> bool:
        return hasattr(self, "popt")

    def fitDataFixed(self, position, slope, delayGuesses):
        if not hasattr(self, "A"):
            self.readData()

        if len(np.shape(self.A)) > 1:
            self.A = np.mean(self.A, axis=0)
            self.T = np.mean(self.T, axis=0)

        invalidFilter = np.isfinite(self.A)
        self.delay = self.delay[invalidFilter]
        self.A = self.A[invalidFilter]
        self.T = self.T[invalidFilter]

        startParameters = [position, slope]
        for guess in delayGuesses:
            startParameters.append(guess["amplitude"]) # A
            startParameters.append(guess["decay"]) # tau
        
        self.popt, self.pcov, infodict, mesg, self.ier = curve_fit(
            ExponentialFunctions.universalExponentialParameters,
            self.delay,
            self.A,
            startParameters,
            method="trf",
            full_output = True
        )

        self.terms = len(delayGuesses)

        return (self.popt, self.pcov)

    def fitDegradedData(self, filenames, degConst, p_ratio, dirPath = r"", backgroundLen = 15, linearBg=False):
        '''this uses Seriesdegradation to adjust single traces as well as interactive fitting to fit as per usual'''
        self.degExpTau = degConst
        self.degPowerRatio = p_ratio
        self.backgroundLen = backgroundLen
        absorbance, transmittance, correctionFactors, _, bgParameters = SeriesDegradation.autoCompensation(filenames, degConst, p_ratio, dirPath, backgroundMean = backgroundLen, linearBg=linearBg)
        #values here are already background subtracted
        self.A = np.mean(absorbance*correctionFactors, axis = 0)
        self.T = np.power(10,-self.A/1000)
        self.removeBackground()
        self.backgroundParameters = bgParameters
        popt, pcov = self.fitData()
        return popt, pcov
    
    def refitDegradedData(self, backgroundLen = 20):
        '''finally works as intended with minimal input'''
        #not sure why the autoCompensation background subtraction works so much better than what I tried with readData
        if hasattr(self, "backgroundParameters") and self.backgroundParameters[0]!=0:
            linBgBool = True
        else:
            linBgBool = False
        absorbance, transmittance, correctionFactors, _, bgParameters = SeriesDegradation.autoCompensation(self.filename, self.degExpTau, self.degPowerRatio, backgroundMean=backgroundLen, linearBg=linBgBool)
        #values here are already background subtracted
        
        self.A = np.mean(absorbance*correctionFactors, axis = 0)
        self.T = np.power(10,-self.A/1000)
        self.removeBackground(backgroundLen)
        self.backgroundParameters = bgParameters
        popt, pcov = self.fitDataFixed(self.startParameters[0]["position"], self.startParameters[0]["slope"], self.startParameters[1:])
        return popt, pcov

       

    def fitData(self, N=10):
        if not hasattr(self, "A"):
            self.readData()

        invalidFilter = np.isfinite(self.A)
        print(np.shape(invalidFilter))
        self.delay = self.delay[invalidFilter]
        self.A = self.A[invalidFilter]
        self.T = self.T[invalidFilter]

        #defaultSlope = 100
        ######---------#######
        #Instead of using a defaultslope using an interactive plot to give zero point and slope
        ######---------#######

    
        Amean = np.mean(self.A)
        #AIndMax = np.argmax(np.abs(self.A - Amean))
        #AMax = self.A[AIndMax]
        #delayMax = self.delay[AIndMax]
        delayMax, defaultSlope, decay_times, amplitudes = self.getStartParameters()
        

        invalidFilter = ~np.isnan(self.A)
        self.delay = self.delay[invalidFilter]
        self.A = self.A[invalidFilter]
        self.T = self.T[invalidFilter]

        startParameters = [delayMax, defaultSlope]

        bounds = [[-np.inf, 0], [np.inf, 5000]]
        
        for decay, amp in zip(decay_times,amplitudes):
            startParameters.append(amp) # A
            startParameters.append(decay) # tau

            bounds[0].append(-np.inf)#amplitude
            bounds[1].append(np.inf)
            bounds[0].append(0)#tau
            bounds[1].append(1e4*(np.max(self.delay) - np.min(self.delay)))
            #bounds[1].append(np.inf)



        
        self.popt, self.pcov, infodict, mesg, self.ier = curve_fit(
            ExponentialFunctions.universalExponentialParameters,
            self.delay,
            self.A,
            startParameters,
            method="trf",
            full_output = True
        )


        return (self.popt, self.pcov)
        

        # Alternate signs and start with equal amplitudes
        startParameters = [delayMax, defaultSlope]  # slope
        for index in range(len(decay_times)):
            #amplitude
            startParameters.append(amplitudes[index])
            #decay time
            startParameters.append(decay_times[index])

        bounds = [[-np.inf, 0], [np.inf, 5000]]
        for parameterIndex in range(N):
            startParameters.append(0)  # Will be set in the next step
            # AMax / N * 2 * (-1)**parameterIndex)  # total amplitude split equally
            startParameters.append(
                10 ** (2 + parameterIndex / 2)
            )  # decay: 100, 141, 1000, 1414, 10000, 14142
            bounds[0].append(-np.inf)
            bounds[1].append(np.inf)
            bounds[0].append(0)
            bounds[1].append(np.max(self.delay) - np.min(self.delay))

        # Solve linear equation using least-squares for start amplitudes
        startFitResult = np.empty(shape=(len(self.delay), N+len(decay_times)))
        for decayIndex, decay in enumerate(startParameters[3::2]):
            startFitResult[:, decayIndex] = ExponentialFunctions.singleExponential(
                self.delay, delayMax, defaultSlope, 1, decay
            )
        amplitudeApprox = np.linalg.lstsq(startFitResult, self.A, rcond=None)
        startParameters[2::2] = amplitudeApprox[0]

        # Bounds are not used at the moment, but we keep them in an array anyway
        self.popt, self.pcov, infodict, mesg, self.ier = curve_fit(
            ExponentialFunctions.universalExponentialParameters,
            self.delay,
            self.A,
            startParameters,
            method="lm",
            full_output = True
        )
        return (self.popt, self.pcov)

    def fitMaxExpTerms(self, maxN: int = 10, maxUncertainty: float = 3):
        encounteredError = True
        terms = maxN

        while encounteredError:
            if terms < 1:
                raise RuntimeError("Could not fit one exponential decay!")

            try:
                self.fitData(terms)
                uncertainty = np.sqrt(np.diag(self.pcov))
                if (
                    # Average of parameters is above uncertainty limit - high covariance, too many parameters
                    np.mean(np.abs(uncertainty / self.popt)) < maxUncertainty
                    or terms == 1   # Minimum exp terms
                    or self.ier < 1 # fit failed
                    or self.ier > 4 # fit failed
                ):
                    # Uncertainties are okay or we reached the minimum amount
                    encounteredError = False
                else:
                    terms -= 1
            except RuntimeError:
                terms -= 1

        return (self.popt, self.pcov, terms)
    
    def getJsonArray(self):
        entry = {}
        entry["pump"] = self.pump
        entry["probe"] = self.probe
        entry["startParameters"] = [] # We don't know what we want to do here yet.
        entry["inputFile"] = self.filename
        if hasattr(self, "popt"):
            entry["popt"] = self.popt.tolist()
            entry["pcov"] = self.pcov.tolist()
        if hasattr(self, "wavelength"):
            entry["wavelength"] = self.wavelength
        entry["startParameters"] = self.startParameters
        return entry
    
    def getStartParameters(self):
        if hasattr(self, 't0_start') == False or hasattr(self, 'slope_start') == False:
            fig, ax = plt.subplots(1,1)
            ax.plot(self.delay*1e-3, self.A)
            if hasattr(self, "backgroundLen"):
                ax.plot(np.ones(2)*self.delay[self.backgroundLen]*1e-3, [np.min(self.A), np.max(self.A)], "r-")
            ax.set_xlabel("t / ps")
            ax.set_ylabel("dOD / mOD")
            #ginput returns list of tuples (x,y)
            #get zero point and background level
            
            fig.suptitle("Zero point then max amplitude")
            points = plt.ginput(n=2, timeout = 0)
            slope = (points[1][1] - points[0][1]) / abs(points[1][0] - points[0][0])
            print(slope)
            slope = 1000
            self.t0_start = (points[0][0]+points[1][0])/2
            self.slope_start = abs(slope)
            user_input = input("Input decay time in fs; comma \",\" separate\n")
            decay_times = np.array(user_input.split(","),dtype=float)
            self.decay_times = decay_times
            user_input = input("Input amplitudes in mOD; comma \",\" separate\n")
            amplitudes = np.array(user_input.split(","),dtype=float)
            self.amplitudes = amplitudes

            plt.close()

        return self.t0_start, self.slope_start, self.decay_times, self.amplitudes
    
    def generateFitJSON(self):
        #have to make the "startParameters dictionary array of dictionaries first"
        numberOfParameters = len(self.popt)-2
        startParameters = [{
                    "position": self.popt[0],
                    "slope": self.popt[1]
                }]
        for parameterPairIndex in range(int(numberOfParameters/2)):
            startParameters.append({
                "amplitude": self.popt[2+parameterPairIndex*2],
                    "decay": self.popt[3+parameterPairIndex*2]
            })
        

        title = input("Measurement (Series) Title\n")
        name_type = input("Type comment (polarization)\n")
        date = input("Date of Measurement: \"YYYY-MM-DD\"\n")
        pumpWav = input("Pump wavelength \ nm\n")
        probeWav = input("Probe wavelength \ nm\n")
        #self.popt is [location(t0), slope, amplitude, decaytime, amplitude, decaytime,....]
        if hasattr(self, "degExpTau"):
            tauDegrade = self.degExpTau
        else:
            tauDegrade = None
        
        if hasattr(self, "degPowerRatio"):
            pRatio = self.degPowerRatio
        else:
            pRatio = None

        if hasattr(self, "backgroundParameters"):
            bgParam = self.backgroundParameters
            if isinstance(bgParam, np.ndarray):
                bgParam = bgParam.tolist()
        else:
            bgParam = [0,0]

        #print(startParameters)
        toBeJSON = {
            "title": title,
        "type": name_type,
        "date": date,
        "entries": [
            {
            "pump": pumpWav,
            "probe": probeWav,
            "startParameters": startParameters,
            "inputFile": os.path.split(self.filename)[1],
            "wavelength": None,
            "tauExpDegrade": tauDegrade,
            "degradePowerRatio": pRatio,
            "bgParam": bgParam
        }]
        }

        return toBeJSON
        