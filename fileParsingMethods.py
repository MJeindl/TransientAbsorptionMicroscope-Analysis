import numpy as np
from scipy.io import loadmat
import os
from matplotlib import pyplot as plt


#so all of these functions should take the filename without leading \ if you use a dirPath
#if they don't I am sorry and it slipped by me in one of the many revisions
#need to fix it so it gives the full filepath instead of something without directory
def parseTime(time):
    '''takes time in hh:mm:ss and returns a value in s'''
    if type(time) == str():
        times = np.zeros(len(time), dtype = int)
        hh, mm, ss = float(str(time[0]).split(":"))
    else:
        times = np.zeros(np.shape(time)[0])
        hh = time[:,0]
        mm = time[:,1]
        ss = time[:,2]

    times[0] = int(hh[0])*3600+int(mm[0])*60+int(ss[0])
    for i in range(1,np.shape(time)[0],1):
        #hh, mm, ss = str(time[i+1]).split(":")
        times[i] = int(hh[i])*3600+int(mm[i])*60+int(ss[i]) - times[0]

    times[0] = 0
    return times

#moved over from SeriesDegradation
def numberFileGenerator(start, stop, namePrefix=r"TA_fourier_", nameSuffix=r""):
    ''''generate TA_fourier_02364 style names in list form'''
    vec = np.arange(start, stop+1,1)
    filenames = []
    for number in vec:
        filenames.append(namePrefix + str(number)+ nameSuffix)
    return filenames


def parseFilenames(filenameArray, directoryPath = r""):
    if directoryPath != r"":
        pathInsert = r"\\"[0]
    else:
        pathInsert = r""
    filenames = [] 
    if type(filenameArray) != type(str()):
        #parsing for summary files
        for file in filenameArray:
            tempfile = loadmat(directoryPath + pathInsert+ file)
            if 'delay' not in tempfile.keys():
                #only works in same directory, multiple summary files
                #can swap to having a directoryPath array instead in the future
                file_current = tempfile['filenames'][0]
                for subFile in file_current:
                    #filenames.append(directoryPath + r"\\"[0] + subFile[0] + ".mat")
                    filenames.append(r"\\"[0] + subFile[0] + ".mat")
                    OutdirectoryPath = os.path.dirname(os.path.abspath(filenameArray[0]))
            else:
                #TA_fourier style file
                filenames.append(pathInsert + file + ".mat")
                OutdirectoryPath = directoryPath
        
    else:
        #if single saturation file
        summaryfile = loadmat(directoryPath + pathInsert + filenameArray)

        OutdirectoryPath = os.path.dirname(os.path.abspath(directoryPath + pathInsert + filenameArray))
        #print(summaryfile['filenames'])
        for fname in summaryfile['filenames'][0,:]:
            if isinstance(fname, np.ndarray):
                fname = fname[0]
            #print(fname)
            filenames.append(r"\\"[0] + fname)

    for i in range(len(filenames)):
        filenames[i] = OutdirectoryPath + filenames[i]
    return filenames

def removeBackground(T, numEntries: int = 20) -> float:
    '''returns T, A and Background levels\\numEntries is how many of the points are to be taken as background'''
    backgroundLevel = np.mean(T[:,:numEntries])
    T += 1 - backgroundLevel
    A = -1000 * np.log10(T)
    return T, A, backgroundLevel

def removeBackgroundLinearFit(T_vec, delay, numEntries: int = 20):
    '''returns T, A and Background levels\\numEntries is how many of the points are to be taken as background\\
    needs a delay vector as this should be related to stage movement distance, which becomes problematic for irregular delaysteps\\
    additionally subtracts a linear function from Ipumpprobe to remove an increasing pump influence within numEntries'''
    print("Linear background fitting")
    bgT_vec = np.mean(T_vec[:,:numEntries], axis=None)

    T_vec += 1 -bgT_vec#- (np.polyval(popt,delay))/bgIprobe +1
    poptTvec = np.polyfit(delay[:numEntries], np.mean(T_vec[:,:numEntries], axis = 0), deg=1)
    for meas_ind in range(np.shape(T_vec)[0]):
        T_vec[meas_ind, :] = T_vec[meas_ind,:]- np.polyval(poptTvec, delay)


    #print(T_vec)
    #print(T)
    A_vec = -1000* np.log10(T_vec)
    return T_vec, A_vec, poptTvec

        
def parseSummaryFileToArray(filenameArray, directoryPath=r"", linearBackgroundSubtract = False, backgroundLen = 15, absorbance = False):
    '''takes or filearray and returns T/A, delay, background'''

    filenames = parseFilenames(filenameArray, directoryPath)
    if linearBackgroundSubtract == True:
        #splitting it like this becaues the lower part already works and I do not want to touch it before I know how to proceed
        #I cannot have a linear subtraction of background and it be different for every measurement
        #print(len(filenames))
        for ind, file in enumerate(filenames):
            if isinstance(file, np.ndarray):
                file = file[0]
            tempFile = loadmat(file)
            delay = tempFile['delay'][0]
            if ind == 0:
                T_vec = np.zeros((len(filenames), len(delay)))

            T_vec[ind, :] = tempFile['d_vec']
            #d_vec = I1_vec/I2_vec
        datArray, aArray, backgroundParameters = removeBackgroundLinearFit(T_vec, delay, backgroundLen)
    else:
        for ind, file in enumerate(filenames):
            #if os.path.isfile(directoryPath + r'\\' + file[0]) == False:
            #    continue

            if isinstance(file, np.ndarray):
                file = file[0]
            tempFile = loadmat(file)
            delay = tempFile['delay'][0]
            if ind == 0:
                datArray = np.zeros((len(filenames), len(delay)))
                aArray = np.zeros(np.shape(datArray))
            datArray[ind, :] = tempFile['d_vec']
        datArray, aArray, backgroundParameters = removeBackground(datArray)
        backgroundParameters = np.array([0, backgroundParameters])
            

    if absorbance == False:
        return datArray, delay, backgroundParameters
    else:
        return aArray, delay, backgroundParameters

def parseSummaryFiletoRaw(filenameArray, directoryPath=r""):
    '''takes or filearray and returns Ipumpprobe, Iprobeonly, delay'''

    filenames = parseFilenames(filenameArray, directoryPath)
    #splitting it like this becaues the lower part already works and I do not want to touch it before I know how to proceed
    #I cannot have a linear subtraction of background and it be different for every measurement
    #print(len(filenames))
    for ind, file in enumerate(filenames):
        if isinstance(file, np.ndarray):
            file = file[0]
        tempFile = loadmat(file)
        delay = tempFile['delay'][0]
        if ind == 0:
            I1_array = np.zeros((len(filenames), len(delay)))
            I2_array = np.zeros((len(filenames), len(delay)))
        I1_array[ind,:] = tempFile['I1_vec']
        I2_array[ind,:] = tempFile['I2_vec']
        #d_vec = I1_vec/I2_vec

            


    return I1_array, I2_array, delay


        
def getTimes(filenameArray, dirPath=r""):
    '''moved out of degradationCompensation for general use'''
    filenames = parseFilenames(filenameArray, dirPath)
    #test if file has keys for single TA measurement
    subTimes = []
    for ind, filename in enumerate(filenames):
        tempLoad = loadmat(filename)
        if 'delay' in tempLoad.keys() and 'd_vec' in tempLoad.keys():
            #this part can easily be broken by a matlab update
            #load time
            timeString = str(tempLoad['__header__']).split("Created on: ")[1]
            timeString = timeString.split(' ')[3]

            tSplit = np.array(timeString.split(':')[:3], dtype=float)
            subTimes.append([tSplit[0], tSplit[1], tSplit[2]])

        elif 'dates' in tempLoad.keys():
            subDates = tempLoad['dates'][0]

            for ind, timeArray in enumerate(subDates):
                timeArray = np.array(timeArray[0][3:6], dtype=float)
                subTimes.append([timeArray[0], timeArray[1], timeArray[2]])
                #subTimes[ind,0] = timeArray[0]#hh
                #subTimes[ind,1] = timeArray[1]##mm
                #subTimes[ind,2] = timeArray[2]##ss with milisecond
    return np.array(subTimes, dtype = float)