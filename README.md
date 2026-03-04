# SENSI EEG Preproc — Bad Channel Detection Module

This repository contains a MATLAB module for semi-automated EEG bad-channel detection with an interactive review UI. The main entry point is `markSusChs.m`.

This introduces a MATLAB Module for bad-channel quality control that emphasizes interpretability, relational structure, and human-in-the-loop validation rather than fully automated rejection. The method operates on multi-channel EEG data and combines complementary channel-level features, including time-dependent neighbor dissimilarity and amplitude- and variance-based statistics to score and pre-label channels as good, suspicious, or bad. To expose shared artifactual structure, channels are additionally grouped using a similarity measure derived from the co-occurrence of robustly detected high-amplitude transients, allowing channels to be reviewed together. Importantly, clustering is used as an exploratory tool to reveal co-artifactual patterns rather than to impose final class labels, which are confirmed through an interactive review interface supported by summary visualizations and grouped channel displays. 

## Example data

Example (deidentified) EEG data used by `example.m` are hosted on the Stanford Digital Repository (SDR):

**Data location:**  
https://purl.stanford.edu/dg856vy8753

The `example.m` script downloads the data files from the SDR URLs into the local `ExampleData/` directory and then runs a demonstration analysis using `markSusChs.m`.

For full details about the dataset contents (file descriptions, variables, and notes), see the dataset README provided inside the `ExampleData` folder, or at the SDR link above.

## Quick start

1. Clone or download this repository.
2. Open MATLAB and add the repository (and subfolders) to your path.
3. Run: example.m


## Reference (please cite)

**Preprint**  
Amilcar J. Malave and Blair Kaneshiro (2026). “EEG Bad-Channel Detection Using Multi-Feature Thresholding and Co-Occurrence of High-Amplitude Transients”. In: bioRxiv. [DOI: 10.64898/2026.02.04.703874](https://doi.org/10.64898/2026.02.04.703874)

**GitHub repository**  
Amilcar J. Malave and Blair Kaneshiro (2026). Bad-Channel Detection Module (v1.0): A MATLAB
framework for semi-automated EEG bad-channel detection and review. Stanford University. https://github.com/edneuro/SENSI-EEG-Preproc-bad-ch

**Dataset**  
Amilcar J. Malave and Blair Kaneshiro (2025). Example EEG data for the SENSI EEG PREPROC Bad-Channel Detection Module [Data set]. Stanford Digital Repository. 
https://doi.org/10.25740/dg856vy8753


## MIT License

Copyright (c) 2025 Amilcar J. Malave, and Blair Kaneshiro.

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
