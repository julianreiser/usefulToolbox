# usefulToolbox
MATLAB toolbox with useful functions for EEG / physiological analysis.
 
# Prerequisites
### Core Requirements
1. MATLAB
 	- Ensure you have MATLAB R2008b (7.6) or later installed.
 	- While EEGLAB functions can run on Octave, full compatibility and GUI support are best with MATLAB.
2. EEGLAB (2021.1 or newer)
 	- Download from the EEGLAB website. ￼
 	- EEGLAB is essential for EEG data processing and provides a framework for various plugins and toolboxes.

### Essential EEGLAB Plugins
1. MoBILAB
  	- Facilitates importing and synchronizing multimodal data, including XDF files.
  	- Install via EEGLAB’s plugin manager or download from the MoBILAB GitHub repository. ￼
2. XDF Import Plugin
  	- Enables loading of XDF files recorded with LabRecorder.
  	- Available through EEGLAB’s plugin manager.
3. BIOSIG Toolbox
	- Supports importing various EEG data formats such as BDF and EDF.
	- Can be installed via EEGLAB’s plugin manager. ￼
4. bva-io Plugin (Optional)
	- Useful if you plan to export data to BrainVision Analyzer format. ￼
	- Install through EEGLAB’s plugin manager.
 
### Recommended MATLAB Toolboxes
While not strictly necessary, the following MATLAB toolboxes can enhance functionality and performance:
- Signal Processing Toolbox: Improves filtering and spectral analysis operations. ￼
- Statistics and Machine Learning Toolbox: Provides advanced statistical functions. ￼
- Optimization Toolbox: Useful for certain EEGLAB extensions and advanced analyses.
- Image Processing Toolbox: Required by some EEGLAB extensions like FieldTrip. ￼

# Functions:
- usefulLslAlign.m - import multiple LSL streams and/or local files (.edf, .bdf, .xdf, .xlsx) and merge them into an EEGLAB-compatible EEG struct for multimodal neuroimaging data analysis.
- usefulMergePupil.ipynb - merge pupil labs recordings from the xdf file with the downloaded pupil cloud files automatically
- usefulQuickPrepro.m - quick-and-dirty preprocessing for data measured during the workshop at the mbt conference 3.0 (see examplary data)

# Examplary data to try the code
Examplary data (recorded at the Real-time neuroergonomics workshop at the 2025 mbt conference (see presentations folder for the workshop materials) to try out the functions can be found [here](https://nextcloud.mbraintrain.com/s/GgCqTpkHmxedjGs)
