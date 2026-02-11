This folder in empty. Data will be populated here when running the file "example.m"

This repository contains example EEG data used to demonstrate the
SENSI EEG PREPROC Bad Channel Detection Module. The data are provided
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

- xAll_129  
  EEG data matrix of size channels x samples. Channel 129 corresponds
  to the reference electrode and is removed prior to bad-channel
  detection, resulting in 128 scalp channels used for analysis.

- fs  
  Sampling frequency (Hz) after downsampling.

- EOG  
  Structure specifying eye-related channels, including vertical and
  horizontal EOG electrodes. These channels are retained to support
  identification and handling of ocular activity during preprocessing.

-----------------------
Notes
-----------------------

The example data are intended to illustrate common channel-quality
issues encountered in EEG preprocessing, including unstable channels,
high-amplitude artifacts, and eye-related activity. They are provided
to enable reproducible demonstration of the module's feature-based
flagging, clustering, and interactive manual review workflow.