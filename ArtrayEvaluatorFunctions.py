from ArtrayAnalysis import ArtFit
import json
from pandas import DataFrame

def JSON_summary(saveTo, pump, probe, JSON_path, fitPath, fitWidth, pumpPower, corrFactor, scale, date, powerDensity):
    json_template = {
    "pathPump": pump,
    "pathProbe": probe,
    "FittingPath": fitPath, 
    "BackgroundJSON": JSON_path,
    "Correction Factor: per W/m^2 pump power": corrFactor,
    "Pump power: W": pumpPower,
    "Pump power density: W/m^2": powerDensity,
    "scale: µm/px": scale,
    "FitAvgWidth / px": fitWidth,
    "Date of acquisition": date
    }
    file = open(saveTo, 'x', encoding='utf-8')
    json.dump(json_template, file, ensure_ascii=False, indent = 4)

def runRoutineJSON(JSON_summary_path):
    "only prints data and shows plots"
    summary = open(JSON_summary_path)
    summary_dat = json.load(summary)

    Routine = ArtFit(summary_dat["BackgroundJSON"], summary_dat["pathProbe"], summary_dat["pathPump"], summary_dat["FitAvgWidth / px"])

    #Set power for pump in W
    Routine.setPumpPower(summary_dat["Pump power: W"])
    #Set scale
    Routine.setScale(summary_dat["scale: µm/px"])

    #Fit data
    Routine.plainFit(1,3)
    Routine.plainFit(2,3)

    #Plot comparison to see if all is right
    Routine.plotComp()
    Routine.plotGauss()



    #get peakRadiance
    peakRadiance = Routine.getPowerDensity(2)[0]
    print("%.3e W/cm^2" %(peakRadiance*1e-4))

    #print fitdata from pandas.DataFrame
    Routine.readableFitData(True)
