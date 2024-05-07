# TransientAbsorptionMicroscope-Analysis
Contains descriptions and scripts of useful methods for analysis of data from the transient absorption microscope
The methods are listed under the headline of their file location
# fileParsingMethods
These methods are used to unify parsing of files and automate handling of different data types produced by the methods used with the TAM.
## parseFilenames:
used to easily parse TA\_fourier filenames along with saturation files into a list of full path TA\_fourier files which can be handed to any function taking arrays of filenames.
## parseSummaryFileToArray:
uses "parseFilenames" to return either transient absorbance or transient transmissivity along with a single delay vector, files must have same delay vector. Possible to remove a linear background by fitting if needed. Alternatively similar function returning the raw intensity vectors along with a delay vector is "parseSummaryFiletoRaw".
## getTimes:
used get starting times of TA\_fourier files or saturation summary files. TA\_fourier file compatibility may be broken by any matlab release changing the header of .m files.

# SeriesDegradation
This file is basis for correction of photodamage. 
## autoCompensation:
used to automatically calculate correction values to for known degradation constant and power for a filename array. Returns transient absorption and transmission corrected for background, which may be linear or static, along with an array of correction factors to be applied and a delay vector. Alternatively "degradationCompensation" calculates only the correction factors from degradation constant, power ratio and number of decay steps.
## plotTrend:
shows a map of the transient absorbances of a filename array against a delay vector.
## fitDegradation
fits degradation of the sample according to eq. \ref{eq:degradationFitting}. Allows for smoothing with a moving mean for noisy data. Needs proper t0 index of temporal pump-probe overlap to work.


# ShowDelayScan
includes methods to quickly plot raw data from either a filename array or from raw data.
# ImageAnalysis
used to show an image scan with a subtracted background. It is a direct transcription of the matlab method.
# ArtrayAnalysis - ArtFit Class
Class used to handle overlap and overlap correction factors. Intended to use raw (binary) data.
A few ArtFit classmethods that are externally accessible and useful:
## parseFitFile
returns a usable array of fitdata from fitdata file saved with ArtFit.
## calculateOverlap
calculates the LiveFitting overlap that is "one" if pump and probe are identical.
## calculateOverlapCorrection
calculates the correction factor for one dimension which multiplies onto the absorbance value.
## rerunFromJSON
runs files from summary JSON file again to allow a comparison between current fitting and past fitting. Allows overwriting of summary JSON.

