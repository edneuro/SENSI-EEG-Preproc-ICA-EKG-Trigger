This folder in empty. Data will be populated here when running the file "example.m"

This repository contains example EEG data used to demonstrate the
SENSI EEG PREPROC EKG-DIN Artifact Detection Module. The data are provided
for illustration, testing, and documentation purposes only.

The example script located in the main module directory
(example.m) automatically downloads the data files from the Stanford
Digital Repository (SDR) into the local ExampleData folder when run.
Users do not need to manually download or manage the data files.

The recordings consist of real, de-identified EEG data that have been
preprocessed to support quality-control demonstrations. Specifically,
the data have been filtered and downsampled prior to release and are
not intended for primary scientific analysis or clinical use.

-----------------------
Data contents
-----------------------

Each data file contains the following variables:

- xRaw  
  EEG data matrix of size channels x samples. Originally, EGI 129-channel
  recordings with the reference and bad channels removed.

- fs  
  Sampling frequency (Hz) after downsampling.

- badCh
  List containing indices of bad channels removed.

-----------------------
Notes
-----------------------

The example data are intended to illustrate common channel-quality
issues encountered in EEG preprocessing. They are provided
to enable reproducible demonstration of the module's workflow.