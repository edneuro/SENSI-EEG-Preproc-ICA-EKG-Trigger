This folder is empty. Data will be populated here when running the file "example.m"

This repository contains example EEG data used to demonstrate the SENSI-EEG-Preproc-ICA-EKG-Trigger Detection Module. The data are provided for illustration, testing, and documentation purposes only.

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

## Reference (please cite)

**Preprint**
Amilcar J. Malave & Blair Kaneshiro (2026). Semi-Automated Identification of EKG and Trigger Artifacts in EEG Using ICA and Spectral Characteristics. bioRxiv. https://doi.org/10.64898/2026.04.08.717297

**GitHub repository**
Amilcar J. Malave & Blair Kaneshiro (2026). SENSI-EEG-Preproc-ICA-EKG-Trigger
(v1.0): A MATLAB framework for semi-automated identification of EKG and trigger
artifacts in EEG using ICA and spectral characteristics. Stanford University. https://doi.org/10.25740/pm515wq0875

**Dataset**
Amilcar J. Malave & Blair Kaneshiro (2026). Example Data for SENSI-EEG-Preproc-ICA-EKG-Trigger [Data set]. Stanford Digital Repository. https://doi.org/10.25740/xd818jt7842