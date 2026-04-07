# SENSI EEG Preproc — EKG & DIN Artifact Detection Module

This repository contains a MATLAB module for semi-automated detection and removal of EKG (cardiac) and DIN (digital input/trigger) artifacts from EEG data using Independent Component Analysis (ICA). The main entry point is `example.m`.

The module provides automated spectral detection of artifact components followed by interactive review UIs, allowing users to confirm or override automated flags before artifact removal. EKG components are identified via harmonic SNR analysis of the cardiac rhythm band; DIN components are identified via harmonic energy ratios in the trigger-frequency band. Artifact removal is performed by zeroing the flagged ICA components and reconstructing the signal in sensor space.

## Prerequisites

- MATLAB with the following toolboxes:
  - Signal Processing Toolbox
  - Statistics and Machine Learning Toolbox
- Pre-processed EEG data: filtered and bad channels removed (variables `xRaw`, `fs`, and `badCh` (if bad channels were removed))
- Electrode locations file (e.g., .sfp format; examples provided in `Sensor Layout/` directory)

## Example data

Example (deidentified) EEG data used by `example.m` are hosted on the Stanford Digital Repository (SDR):

**Data location:**  
https://purl.stanford.edu/xd818jt7842

The `example.m` script downloads the data files automatically into the local `ExampleData/` directory and then runs a demonstration analysis.

For full details about the dataset contents (file descriptions, variables, and notes), see the dataset README at the SDR link above.

## Quick start

1. Clone or download this repository.
2. Open MATLAB and add the repository (and subfolders) to your path.
3. Run: `example.m`

## Reference (please cite)

**Preprint**  
Amilcar J. Malave & Blair Kaneshiro (2026). "Semi-Automated Identification of EKG and Trigger Artifacts in EEG Using ICA and Spectral Characteristics". bioRxiv. [DOI: 10.64898/2026.02.04.703874](https://doi.org/10.64898/2026.02.04.703874)

**GitHub repository**  
Amilcar J. Malave & Blair Kaneshiro (2026). SENSI-EEG-Preproc-ICA-EKG-Trigger: A MATLAB framework for semi-automated identification of EKG and trigger artifacts in EEG using ICA and spectral characteristics. Stanford University. https://github.com/edneuro/SENSI-EEG-Preproc-ICA-EKG-Trigger

**Dataset**  
Amilcar J. Malave & Blair Kaneshiro (2026). Example Data for SENSI-EEG-Preproc-ICA-EKG-Trigger [Data set]. Stanford Digital Repository.  
https://doi.org/10.25740/xd818jt7842


## MIT License

Copyright (c) 2026 Amilcar J. Malave, and Blair Kaneshiro.

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
