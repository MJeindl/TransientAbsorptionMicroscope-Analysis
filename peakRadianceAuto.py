import numpy as np
import json

def parseFitFile(filepath):
    '''takes background fitfile and returns parsed parameters with errors in pairs\\
        returns dict with number and whatever x an y fit FWHM is available (with error) as array\\
        x[0] = FWHM; x[1] = deltaFWHM'''
    #Names;µ / px;FWHM / px;Amplitude;Offset;µm/px
    #names given as "1st-x" and similar, grab first letter for number and last letter for x or y to parse
    fitnames = np.genfromtxt(filepath, delimiter = ";", dtype=str, usecols=(0), skip_header = 1)
    fitfile = np.genfromtxt(filepath, delimiter = ";", dtype=float, usecols=(2), skip_header = 1)
    

    #make a tuple array
    numberAndAxis = list()
    numbers_unique = []
    for row in range(0,int(np.shape(fitfile)[0]),2):
        number = int(fitnames[row][0]) #grab the number from first index
        xory = fitnames[row][-1]
        numberAndAxis.append((row, number, xory))
        if number not in numbers_unique:
            numbers_unique.append(number)
    
    #read names from first column and put into dictionary
    sorter = lambda array: array[:][1]
    fitdata = [dict() for i in range(len(numbers_unique))]
    numberAndAxis.sort(key=sorter)

    for tuple_loop in numberAndAxis:
        fitdata[tuple_loop[1]-1].update({"number" : tuple_loop[1], tuple_loop[2] : [fitfile[tuple_loop[0]], fitfile[tuple_loop[0]+1]]})
    return fitdata
        

#parseFitFile(r"C:\Users\M\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm\STD-680\BeamFitsPump653Probe680.txt")
#this is a terrible solution, I hate this, the absolute paths are my bane at this point, I should have known
def autoPeakRadiance(filepath, fittingPath=None):
    '''takes full filepath to summaryfile and returns peak radiances with error of any\\
    assumes power = 1 for when pump or probe power are not given\\
    ignores any errors not given'''

    f_JSON = open(filepath)
    data_JSON = json.load(f_JSON)
    pulse_power = np.zeros((2,2))
    try:
        pulse_power[1,0] = data_JSON['Pump power: W']
    except:
        pulse_power[1,0] = 1
    try:
        pulse_power[0,0] = data_JSON['Probe power: W']
    except:
        pulse_power[0,0] = 1
    try: 
        pulse_power[1,1] = data_JSON['Error Pump power: W']
    except:
        pulse_power[1,1] = 0
    try: 
        pulse_power[0,1] = data_JSON['Error Probe power: W']
    except:
        pulse_power[0,1] = 0


    umPerPx = [data_JSON["scale: Âµm/px"], 0]
    try:
        umPerPx[1] = data_JSON["Error scale: Âµm/px"]
    except:
        umPerPx[1] = 1


    try:
        if fittingPath == None:
            fittingPath = data_JSON["FittingPath"]
    except:
        Exception("no scale factor for sensor available")

    f_JSON.close()
    dataDict = parseFitFile(fittingPath)
    return peakRadiancefromInput(dataDict, pulse_power, umPerPx)

def peakRadiancefromInput(dataDict, pulse_power, umPerPx):
    #pulled it apart so I can easily parse whatever I want

    factorFWHM = 2*np.sqrt(2*np.log(2))
    peakRadiance = lambda power, sigX, sigY, umPerPx : power[0]/(sigX[0]*sigY[0]*umPerPx[0]**2*1e-12*2*np.pi)
    dPeakRadiance = lambda power, sigX, sigY, umPerPx : power[1]/(sigX[0]*sigY[0]*umPerPx[0]**2*1e-12*2*np.pi) + power[0]*(sigX[1]/sigX[0]+sigY[1]/sigY[0]+2*umPerPx[1]/umPerPx[0])/(sigX[0]*sigY[0]*umPerPx[0]**2*1e-12*2*np.pi)
    outRadiance = np.zeros(len(dataDict))
    outDeltaRadiance = np.zeros(len(dataDict))
    for index in range(len(dataDict)):
        try:
            if 'x' not in dataDict[index]:
                Exception("No xFWHM for input %d" %dataDict['number'])
            if 'y' not in dataDict[index]:
                Exception("No yFWHM for input %d" %dataDict['number'])
            
            xSig = dataDict[index]['x']/factorFWHM
            ySig = dataDict[index]['y']/factorFWHM
            outRadiance[index] = peakRadiance(pulse_power[index], xSig, ySig, umPerPx)
            outDeltaRadiance[index] = dPeakRadiance(pulse_power[index], xSig, ySig, umPerPx)

            if pulse_power[index][0] == 1:
                print("power assumed to be 1 W")
            if pulse_power[index][1] == 0:
                print("no power error given")
            if umPerPx[1] == 0:
                print("no scale error given")

            print("Peak radiance for signal %d is (%.2f +- %.1f) W/m^2" %(index, outRadiance[index], outDeltaRadiance[index]))
            

        except:
            print("No radiance possible for ...")
            continue
    return outRadiance, outDeltaRadiance


    

#outRad, outDelta = autoPeakRadiance(r"C:\Users\M\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm\STD-680\2023.12.22_PumpPowerVar\SummaryPump653Probe680_0.JSON", r"C:\Users\M\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm\STD-680\2023.12.22_PumpPowerVar\BeamFitsPump653Probe680_0.txt")