#reruns artray from file

#JSON_summary = r"C:\Users\Maximilian\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm\Pump653Probe493_ArtrayCompensationControl\SummaryPump653Probe493Caution.JSON"
#ArtFit.rerunFromJSON(JSON_summary, True, True, True)


from ArtrayAnalysis import ArtFit


JSON_summary = r"C:\Users\M\Documents\phdmatlab\sqib-pmma-probe-wavelength\UV_Setup\new_parallel_pol_pump653nm\SampleAreaImagingPump653Probe680\2024.01.15\Summary653_680_2024.01.15_imaging.JSON"
#rerunFromJSON(self, JSONpath, overwrite = False, showPlots = False, bgSubtract = False, maxCorr = 3):
routine = ArtFit.rerunFromJSON(JSON_summary, True, True, True)
