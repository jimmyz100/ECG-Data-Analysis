# Zebrafish ECG Signal Processing Algorithm
---

The enclosed program 'signal_processing_data_cleaning.m' is the main function that takes raw ECG data obtained from LabView and preprocesses, denoises, detects the PQRST waves of the ECG, and calculates ECG parameters (Wave intervals, ratios) typically used in clinical practice. 

The general workflow of the program is as follows:
1. Formats the raw txt file format taken from LabView into arrays, initially trims data based on noisy sections (fish movement, electrical noise), organizes data to allow serial processing with the program with Trimtxtfast.m
2. Preprocess and denoise data with the extractNoise.m and preprocessingWT.m programs, which segregates data into 10-second fragments, determines an amplitude threshold based on the average ECG amplitude and removes segments that exhibits amplitudes above and below the threshold. See following figure for visualization of how ECG segments are partitioned to automatically remove segments. To remove high-frequency noise, an initial Dolph-Chebyshev low-pass filter with a cutoff frequency of 40 Hz is applied, decreasing signals above 40 Hz by at least 40 Db. For low-frequency noise (baseline wander, gill motion), a wavelet denoise algorithm is applied using the Daubechies mother wavelet number 6, which is utilized for the similar characteristic shape as an ECG waveform.

![Signal Preprocessing](https://ieeexplore.ieee.org/mediastore_new/IEEE/content/media/7361/8706548/8636511/le3-2897789-large.gif)

3. Defines sampling rate, predetermined parameters for processing based on previous work: (https://www.mdpi.com/1424-8220/18/1/61)
4. Detects the R peaks, which are the easiest to visualize as the local maxima of an ECG cycle. Points that are at least 50% of the maximum signal amplitude are designated as R peaks. The RR-interval is calculated as the time interval between R peaks.
5. Detects the Q and S waves, which are the next features branching out from the R peak. They are determined by finding the minima approximately within 50 ms before and after the timepoint of the R peak, respectively.
6. Detects the P and T waves, which are the next set of features branching out. The P wave is determined by the local maximum 65% to 95% of the calculated RR-interval from the preceding R peak. The T wave is determined as the local maximum 15% to 55% of the calculated RR-interval after the R peak.
7. Output plotted ECG waveform data with the aforementioned waveforms denoted as markers.
8. Calculates parameters for subsequent analysis, including heart rate, heart rate variation, PR-interval, QRS-duration, QT-interval, ST-interval.

![ECG Acquisition Setup and Signal Processing](https://ibb.co/c3VwkKN)
(a)Figure of zebrafish ECG acquisition, processing of data from LabView, (b) Raw input data (in orange) and annotated data after utilizing aforementioned program (in blue)

![Sample ECG cycle and parameters](https://ibb.co/zZzdf2s)

For more information, see: (https://ieeexplore.ieee.org/abstract/document/8636511)
# Data Obtained from Consecutive treatments of methamphetamine promote the development of cardiac pathological symptoms in zebrafish
---

Dataset contains the tabulated data acquired from ECG recordings as well as molecular assays determining how the GPCR pathway is impacted. Zebrafish were treated with 200 uM of methamphetamine over the course of 2 weeks in order to assess the effect on heart physiology.


## Description of the Data and File Structure

The csv files titled 'Compiled ECG Data' describes the HR (heart rate), HRV (heart rate variation, calculated by RMSSD), QRS, PR, and QTc intervals, all parameters relevant for analyzing ECG. Two experimental groups were listed: control (no methamphetamine) and treated (200 uM methamphetamine), and ECG were recorded prior to any treatment (baseline), as well as subsequent days during the study, with Day 1 as the first day of the study. Both the mean and standard error are provided for the two trials conducted for this study, labeled correspondingly in the separate files. 

The supplemental figures are derived from the corresponding manuscript. Briefly, these figures include representative ECG figures, plotted ECG parameters as described above, as well as genetic expression and histological data corresponding to the consequences of methamphetamine treatment.






## Sharing/access Information

Links to publicly accessible locations of the data: https://datadryad.org/stash/dataset/doi:10.7280/D1269D
Original Manuscript: https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0294322

Was data derived from another source?
No
