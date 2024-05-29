from ArtrayAnalysis import ArtFit
import json
from pandas import DataFrame
import ArtrayEvaluatorFunctions
r'''
bg_map = r"C:\Users\Maximilian\Documents\CameraData_07.12.2023\Artray\ArtrayBackground_10.10.2023_exp15ms.txt"
save_path = r"C:\Users\Maximilian\Documents\CameraData_07.12.2023\Artray\JSON_ArtrayBG_27.10.2023.JSON"
ArtFit.backgroundJSON_Gen(save_path, bg_map, 780, 100, 3.5228, "27.10.2023")
'''

JSON_path = r"C:\Users\M\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm\JSON_ArtrayBG_2024.01.16.json"
pumpPath = r"C:\Users\M\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm\SHG-440\Orthogonal\Pump653forProbe440.txt"
probePath = r"C:\Users\M\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm\SHG-440\Orthogonal\Probe440.txt"
fitPath = r"C:\Users\M\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm\SHG-440\Orthogonal\BeamFitsPump653Probe440_OrthPol.txt"
summaryPath = r"C:\Users\M\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm\SHG-440\Orthogonal\SummaryPump653Probe440_OrthPol.JSON"
acquisitionDate = "2024.01.15"

### Actual settings for fitting
fitWidth = 3
pumpPower = 2.45e-6 #W


#ArtFit(probePath)
Routine = ArtFit(JSON_path, probePath, pumpPath, fitWidth)


#Set power for pump in W
Routine.setPumpPower(pumpPower)
#Set scale if not set in JSON
#currently not implemented fetch from JSON yet
#Routine.setScale(1)

#Fit data
Routine.plainFit(1)
Routine.plainFit(2)

#Plot comparison to see if all is right
Routine.plotComp()
Routine.plotGauss()

#Routine.getNumericCorrection(0.02)


#get peakRadiance
peakRadiance = Routine.getPowerDensity(2)[0]
print("%.3e W/cm^2" %(peakRadiance*1e-4))

print(Routine.getCorrectionFactor()*peakRadiance)
#save fitdata from pandas.DataFrame
Routine.readableFitData()[0].to_csv(fitPath, sep = ";", index = False, mode = 'x', encoding="utf-8", decimal = '.')



UmPerPxscale, scaleBool = Routine.getScale(True)
if scaleBool == False:
        UmPerPxscale = "NaN"

#save summary
ArtrayEvaluatorFunctions.JSON_summary(saveTo=summaryPath, pump=pumpPath, probe=probePath, JSON_path = JSON_path, fitPath = fitPath, fitWidth = fitWidth, pumpPower = pumpPower, corrFactor = Routine.getCorrectionFactor(), scale=UmPerPxscale, date = acquisitionDate, powerDensity=peakRadiance)
