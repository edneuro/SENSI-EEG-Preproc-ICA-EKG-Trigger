# SENSI EEG Preproc — ICA-EKG-Trigger Module

This repository contains a MATLAB module for semi-automated detection and removal of EKG (cardiac) and DIN (digital input/trigger) artifacts from EEG data using Independent Component Analysis (ICA). The main entry point is `exampleV2.m`.

The module provides automated spectral detection of artifact components followed by interactive review UIs, allowing users to confirm or override automated flags before artifact removal. EKG components are identified by two methods: harmonic SNR of the cardiac rhythm band, scored in chunks of time, and the regularity of the beat train in the time domain; DIN components are identified by a comb score over the trigger harmonics. Artifact removal is performed by zeroing the flagged ICA components and reconstructing the signal in sensor space.

## Prerequisites

- MATLAB with the following toolboxes:
  - Signal Processing Toolbox
  - Statistics and Machine Learning Toolbox
- Pre-processed EEG data: filtered and bad channels removed (variables `xRaw`, `fs`, and `badCh` (if bad channels were removed))
- Electrode locations file (e.g., .sfp format; examples provided in `Sensor Layouts/` directory)

## Example data

Example (deidentified) EEG data used by `exampleV2.m` are hosted on the Stanford Digital Repository (SDR):

**Data location:**  
https://purl.stanford.edu/xd818jt7842

The `exampleV2.m` script downloads the data files automatically into the local `ExampleData/` directory and then runs a demonstration analysis.

For full details about the dataset contents (file descriptions, variables, and notes), see the dataset README at the SDR link above.

## Quick start

1. Clone or download this repository.
2. Open MATLAB and add the repository (and subfolders) to your path.
3. Run: `exampleV2.m`, from the repository folder.

The V1 tutorial, `example.m`, has moved to `Helper Functions/V1/example.m`. It still runs the
V1 detectors and review UIs, unchanged.

**User Manual.** `User Manual.pdf` describes V1 (`Helper Functions/V1/example.m` and its functions). There is no
separate manual for V2: `exampleV2.m` follows the same steps, and its review windows are used
the same way (click a row to switch it between red, reject, and green, keep; press Done to
submit). What V2 changes is listed in the table below. V1's functions are unchanged and now
sit in `Helper Functions/V1`, so add the repository with its subfolders to the MATLAB path.

## V2 detectors and review UIs

The V1 files, now in `Helper Functions/V1`, are unchanged: `example.m`, `autoDetectEkg`,
`autoDetectDIN` and their review UIs behave exactly as before. `exampleV2.m` follows V1's
`example.m`, with the detection and review blocks switched to V2.

| v1 | V2 | what changes |
|---|---|---|
| `autoDetectEkg` | `autoDetectEkgV2` | Two arms: the released harmonic SNR, scored per 30 s chunk, and `rrReg2Score`, the regularity of the beat train in the time domain. The panel holds 6: the top 3 of each arm, then whichever components stand out most. |
| `autoDetectDIN` | `autoDetectDINV2` | `dinCombScore`'s `combC`, a weighted **sum** over the spectrum (positive on the harmonics of 1/T, negative between them), with outlying bins dropped from the peaks and from the valleys separately, instead of a harmonic-to-midpoint **ratio**. Works at any trigger period T. |
| `reviewEkgArtifactUI` / `reviewDINArtifactUI` | `reviewEkgArtifactUIV2` / `reviewDINArtifactUIV2` | Same red / green rows. Added: pages of at most 4 rows (`rowsPerPage`); a start time that is either seconds or `'strongest'` (each row at its own strongest chunk); an in-figure *Window start* box; *add IC* + *+ Add*; a status line; and an optional second window with each component's power per channel (`showICChannelsV2`, `icChannels`, on by default). |
| `reviewOtherUI` | `reviewOtherUIV2` | `reviewOtherUI`, unchanged, plus the same optional second window with each component's power per channel. |


## Reference (please cite)

**Preprint**  
Amilcar J. Malave & Blair Kaneshiro (2026). "Semi-Automated Identification of EKG and Trigger Artifacts in EEG Using ICA and Spectral Characteristics". bioRxiv. [DOI: 10.64898/2026.04.08.717297](https://doi.org/10.64898/2026.04.08.717297)

**GitHub repository**  
Amilcar J. Malave & Blair Kaneshiro (2026). SENSI-EEG-Preproc-ICA-EKG-Trigger (v1.0): A MATLAB framework for semi-automated identification of EKG and trigger artifacts in EEG using ICA and spectral characteristics. Stanford University. https://github.com/edneuro/SENSI-EEG-Preproc-ICA-EKG-Trigger. [DOI: 10.25740/pm515wq0875](https://doi.org/10.25740/pm515wq0875)

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
