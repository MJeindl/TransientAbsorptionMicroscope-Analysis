import json
from os import path
from DelayScanInteractive import DelayScan
import ExponentialFunctions

#this is the only part that needs to be changed
jsonPath = r"C:\Users\M\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm\SHG-350\TAfitPump653Probe350_DecayCorrected.JSON"
bgLen = 30
#end of part to be changed

jsonFile = open(jsonPath, 'r')
jsonDictin = json.load(jsonFile)
satFile = jsonDictin["entries"][0]["inputFile"]
#assume same directory position as the jsonPath
directoryPath = path.dirname(jsonPath)
#now redo the fitting with the same number of 



#following is adapted from fitInteractiveSingleDelay
#print(type(jsonDictin["entries"][0]))
#prescan = DelayScan()
scan = DelayScan.loadFromJson(jsonDictin["entries"][0])
scan.filename = directoryPath+r"\\"[0]+ satFile
scan.readData(degradeBool=True)

import matplotlib.pyplot as plt
import numpy as np

#parsing the different options

if hasattr(scan,"backgroundParameters") and scan.backgroundParameters[0] !=0:
    linearBgsubtraction = True
else:
    linearBgsubtraction=False

(popt, pcov) = scan.refitDegradedData(bgLen)


t_data = np.linspace(np.min(scan.delay), np.max(scan.delay), 1000)
y_data = ExponentialFunctions.universalExponentialParameters(t_data, *popt)
puncert = np.sqrt(np.diag(pcov))
fig, ax = plt.subplots(1,2)
fit_label = [
    f"$t_0$={popt[0]:.0e} $\pm$ {puncert[0]:.0e} fs",
    f"$\sigma$={popt[1]:0e} $\pm$ {puncert[1]:.0e} fs"]
for fit_index, (amplitude, decay, dAmplitude, dDecay) in enumerate(zip(popt[2::2], popt[3::2], puncert[2::2], puncert[3::2])):
    fit_label.append(f"$A_{fit_index}$={amplitude:.2e} $\pm$ {dAmplitude:.2e} $\Delta$mOD")
    fit_label.append(f"$t_{fit_index}$={decay:.1e} $\pm$ {dDecay:.0e} fs")

ax[0].plot(np.ones(2)*scan.delay[bgLen], [np.min(scan.A), np.max(scan.A)], "r-")
ax[0].plot(scan.delay, scan.A, label="experiment")
ax[1].set_ylabel("Residual")
ax[1].set_xlabel("t / fs")
ax[1].plot(scan.delay, scan.A-ExponentialFunctions.universalExponentialParameters(scan.delay, *popt))


ax[0].plot(t_data, y_data, label="\n".join(fit_label))
ax[0].legend()
ax[0].set_xlabel("$\Delta$t / fs")
ax[0].set_ylabel("$\Delta$A / $\Delta$mOD")
#ax.set_title(args.filename)
ax[0].legend()
plt.show()

#save but not destructively
correctionBool= True

fitSaveBool = input("Save Fitparameters in JSON? yes/no\n")
if fitSaveBool == "yes" or fitSaveBool == "y":
    #alter the json
    numberOfParameters = len(popt)-2
    startParameters = [{
                    "position": popt[0],
                    "slope": popt[1]
                }]
    for parameterPairIndex in range(int(numberOfParameters/2)):
        startParameters.append({
            "amplitude": popt[2+parameterPairIndex*2],
                "decay": popt[3+parameterPairIndex*2]
        })
    jsonDictin["entries"][0]["startParameters"] = startParameters
    jsonDictin["entries"][0]["bgParam"] = scan.backgroundParameters.tolist()
    #if correctionBool == True:
    #    savePath = directoryPath+"\TAfitPump" + fitJSON["entries"][0]["pump"]+"Probe" + fitJSON["entries"][0]["probe"] + "_DecayCorrected.JSON"
    #else:
    #    savePath = directoryPath+"\TAfitPump" + fitJSON["entries"][0]["pump"]+"Probe" + fitJSON["entries"][0]["probe"] + ".JSON"
    savePath = jsonPath
    if path.isfile(savePath):
        overwriteBool = input("Do you want to overwrite "+ savePath + " yes/no\n")
        if overwriteBool == "y" or overwriteBool == "yes":
            file = open(savePath, "w")
            print("Overwriting")
            json.dump(jsonDictin, file, ensure_ascii=False, indent = 4)
            file.close()
        else:
            print("aborting")
    else:
        print("Saving to " + savePath)
        file = open(savePath, "w")
        json.dump(jsonDictin, file, ensure_ascii=False, indent = 4)
        file.close()