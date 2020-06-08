## Utility for processing raw Holter data
## *** NOT FOR MEDICAL DIAGNOSTIC OR OTHER CLINICAL USE ***

## Project Title
Repose: A utility for processing raw Holter monitor data

## Description
Repose is software for converting raw sensor data from a 'Holter' cardiac monitor to useable body position, heart rate and respiration rate data (saved in .xlsx format). It will be useful for groups and individuals conducting physiological research that requires these variables to be recorded from a convenient portable device. Since EDF files are not human readable, an EDF reading utility is normally used to view, edit and export signals (in CSV or EXCEL format for example).  This utility allows body position (e.g. supine, tilted, left/right side, prone etc.) and angles (tilt and rotation) to be calculated from a Bittium Faros 90, Faros 180 or Faros 360 Holter monitor, although it could be easily adapted for use with other brands.

## Technology stack
Repose is written in Python 3.X and probably works best running on a desktop IDE (such as Spyder or Pycharm).

## Dependencies
The following Python libraries should be installed prior to use:
 * pyedflib
 * datetime
 * scipy.signal
 * pandas
 * numpy
 * enum
 * os

## Usage
The input file(s) should be in EDF format with at least the following data types:

 * ECG
 * Accelerometer_X
 * Accelerometer_Y
 * Accelerometer_Z

Other data types will be ignored.  Any sampling rate (at least 12.5 Hz) is supported but we recommend 100 Hz for accelerometer and 250 Hz for ECG.

The output file (in Excel format) outputs data nominally at 1 Hz (one row per sample but this can be changed) with the following column headings:
 * Timestamp (format: YYYY-MM-DD HH-MM-SS)
 * RR (respiratory rate BPM)
 * HR (heart rate BPM)
 * Pos_Index (0 - 5)
 * Position (supine, prone, tilt, sitting, left_side, right_side)
 * Tilt (angle of the Holter body from horizontal: 0˚ - supine -> +90˚ - sitting)
 * Rotation (angle of rotation of the body: +180˚ - supine -> left side -> 0˚ - prone -> right side -> -180˚

## Configuration/settings
The following variables can be changed (see asterisks(*) in the main loop):
 - filename of directory containing EDF files.
 - Averaging window length for heart rate and resp. rate (defaults = 15s, 25s respectively)
 - Calculation frequency ('calc every') for heart rate and resp. rate (defaults = 5s, 5s respectively)
 - Output frequency (i.e. how often a row of the .xls spreadsheet is written).  Uses pandas frequency notation e.g. ('S' = every 1 second, try '0.1S', '10S' etc.)
 

## Getting Started
 * Download the .py file to a local folder.
 * In a python IDE or text editor: edit line 172 so that the 'pathname' variable points to the pathname of a folder containing EDF files.
 * Running the program will sequentially read all files in the folder (e.g. 01-01-20.EDF)
 * For each EDF file in the folder, an excel file will be written to the folder using the same base filename with a .xlsx extension (e.g. 01-01-20.xlsx)
 * Note that large EDF files may take some time to process, e.g. processing a 24 hour EDF file will take approximately 40 minutes to complete (depending on local resources).


## Known issues
Heart rate estimation uses a very lightweight algorithm that is not robust when heavy movement artifact is present on the ECG signal.  Respiratory rate is similarly susceptible to motion artifact.  Derivation of respiratory rate from acceleometer signals is not a tried and tested method in any case.  If you want robust RR, try using a respiratory belt or other monitoring device.  The accuracy nor the reliability have not been thoroughly tested or validated at any scale so treat results with caution.




