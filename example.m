% example.m
% =========================================================================
% TUTORIAL: Independent Component Analysis (ICA) for EKG and DIN Artifact Removal

% This script demonstrates a new workflow for identifying and removing
% components related to EKG (cardiac) and DIN (digital input/trigger)
% artifacts using ICA decomposition.

%%%%% PREREQUISITES AND DEPENDENCIES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% To run this script, ensure you have the following:
% (1 & 2 are provided in the sample data)
%
% 1. DATA: Requires pre-cleaned EEG/MEG data that has been FILTERED and
%    BAD CHANNELS REMOVED (variables xRaw, fs, and badCh).
% 2. ELECTRODE LOCATIONS: Requires a structure containing electrode coordinates
%    (e.g., locsEGI124.mat) for topography plotting.
% 3. MATLAB TOOLBOXES (EXTERNAL LIBRARIES):
%    - Signal Processing Toolbox
%    - Statistics and Machine Learning Toolbox
%    - Image Processing Toolbox
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% 1. SETUP AND INPUT PARAMETERS

% -------------------------------------------------------------------------
% Define paths and data identifiers. These parameters control the data
% input, where figures are saved, and the general artifact detection rules.
% -------------------------------------------------------------------------
clear; close all; clc
disp('~ * ~ * TUTORIAL: ICA Artifact Cleaning Demonstration * ~ * ~')

%%% FILE AND PATH CONFIGURATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
INFO.fileDir = 'D:\Stanford\Data\artifact_tutorial'; % Directory containing data files
INFO.figDir = './Figures';                           % Output directory for saving figures
INFO.fileName = 'example_10';                         % Base name of the input data file
saveFigs = 1;                                        % 1 to save review figures, 0 otherwise

%%% Load electrode locations %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
S = load('locsEGI124.mat','locs');

%%% ARTIFACT DETECTION PARAMETERS %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These thresholds control the automated component flagging process.

% DIN (Digital Input) Parameters
HRFOpts.fmin = 10; % Minimum frequency for detection (Hz)
HRFOpts.fmax = 50; % Maximum frequency for detection (Hz)
HRFOpts.bw = .2;   % Half-bandwidth for frequency window (Hz)

% EKG (Cardiac) Parameters
ekgOpts.fmin = .8; % Lower bound of cardiac rhythm band (Hz)
ekgOpts.fmax = 2;  % Upper bound of cardiac rhythm band (Hz)
ekgOpts.harmonics = 4; % Number of harmonics to include in spectral test

%%% DATA ASSESSMENT THRESHOLDS (Used in the final QC plot) %%%%%%%%%%%%%%%%
INFO.recUVThresh = 50;  % uV magnitude threshold for a 'bad' sample
INFO.recPctThresh = 15; % % of 'bad' samples required to flag a channel as 'recording bad'


%% 2. LOAD DATA AND PREPARE LOCATIONS

% Load the raw data file. This file must contain the variables:
load([INFO.fileDir filesep INFO.fileName]) 

% Initialize the main channel location variable (chanlocs) from the structure.
chanlocs = S.locs; 

% Remove the locations corresponding to the channels previously identified as bad.
if ~exist('badCh','var') % Ensuring badCh exists
    badCh = [];
end
if ~isempty(badCh)
    chanlocs(badCh) = []; % Removing bad channel locations
end


%% 3. ICA DECOMPOSITION (or Loading Pre-computed W Matrix)
% -------------------------------------------------------------------------
% The data must be transformed from sensor space (xRaw) into the ICA
% component space (xICA) using the unmixing matrix (W).
% For this tutorial, we will load a pre-computed W matrix for speed.
% -------------------------------------------------------------------------
INFO.DoICA = 0; % Set to 1 if you want to run doICA(xRaw), otherwise 0.

if INFO.DoICA
    W = doICA(xRaw); % Decompose the raw data
else
    load([INFO.fileDir filesep INFO.fileName '_W']); % Load precomputed W matrix
end

% Transform the raw data into ICA component space
xICA = W * xRaw;
close all;
disp('-> Data transformed to ICA space (xICA).')


%% 4. AUTOMATED DETECTION AND INTERACTIVE REVIEW

% -------------------------------------------------------------------------
% This section performs automated detection for EKG and DIN artifacts,
% followed by "Interactive Review" to confirm or correct the
% automated flags. The UIs are the primary decision points.
% -------------------------------------------------------------------------

%%% Instructions
% 1 - Run the code block.
% 2 - Two interactive figures will render: the **DIN Review UI** and the **EKG Review UI**. 
%     These figures will **not** automatically save; saving is triggered by 
%     clicking **'Done'** inside the UI, provided 'INFO.saveFigs = 1'
% 3 - Review the figures and make your final decisions for each component.

% INTERACTION (DIN UI and EKG UI):
% 4 - The UIs display suspected components. **Sources start flagged for REJECT (red)**.
% 5 - **To KEEP a source** (i.e., you do not believe it is an artifact), 
%     **click anywhere on its row** (topography, time series, or frequency plot). 
%     The row will toggle to **green**.
% 6 - Review the components:
%     - For **EKG**, look for the classic cardiac beat waveform in the time series,
%       Frequency Peaks on Harmonics of 1, and characteristic topography.
%     - For **DIN**, look for abrupt, simultaneous sharp transients in the
%       time series, and comb like harmonics in the frequency domain. 
% 7 - Click **'Done'** in each UI to submit your final selections and proceed to the next step.

% [OPTIONAL] Adjusting the Time Window:
% 8 - The time-domain plots show a 10-second segment starting at `tempStartSec` 
%     (default 5 minutes, or 5*60 seconds). In rare cases where an artifact is 
%     ambiguous (e.g., EKG-like), you may visualize a different segment of data 
%     by adjusting the value of **`tempStartSec`** in the preliminary variables 
%     section. Be sure any value you enter does not exceed the total 
%     length of the data. **In most cases, you should NOT need to adjust this 
%     parameter. Please do not spot-check different data segments by default.**


%%% PLOTTING SETUP %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tempNSources = [1:124]; % Focus on the first 124 components
tempStartSec = 5*60;   % Start time for time-series plots (e.g., 5 minutes)
tempNSec = 10;         % Duration of the time-series plots


%%%%%%%%%%%%%%%%%%%%%% DIN Detection %%%%%%%%%%%%%%%%%%%%%%%%%
% Detect DIN ICA components via harmRatioFromSignal
[dinSrc, tempDinInfo] = autoDetectDIN(W, xICA, badCh, fs, [], tempNSources, ...
    tempStartSec, tempNSec, HRFOpts);

% Interactive review of DIN candidates
dinSrc = reviewDINArtifactUI(W, xICA, fs, chanlocs, dinSrc, tempDinInfo.reviewCandidates, ...
   tempStartSec, tempNSec, [INFO.fileName '_DIN'], ...
   saveFigs, INFO.figDir);

if isempty(dinSrc); fprintf('\n\tNo DIN artifact detected\n')
else; fprintf('\n\tdinSrc = %s\n', mat2str(dinSrc));
end

%%%%%%%%%%%%%%%%%%%%%%% EKG Detection %%%%%%%%%%%%%%%%%%%%%%%%%
% 1) KG Harmonics Test (detecting beats via spectral analysis)
[tempEkgHarmScore, tempMetrics] = autoDetectEkgHarm(xICA(tempNSources,:), ...
    fs, tempNSources, ekgOpts);

% 2) ICA Topography Test (EKG like topography)
tempEkgTopoScore = autoDetectEkgTopo(W, chanlocs, tempNSources, 0);
tempEkgTopoScore = zscore(tempEkgTopoScore.');

% 3) Combining Detection Metrics
% Create a binary mask for "likely EKG" sources (high score on both metrics)
tempNoEkg = (tempEkgHarmScore > 2.5) .* (tempEkgTopoScore > 1.5);
% Simple average of Z-scores to get a composite EKG artifact score
tempEkgScore = (tempEkgHarmScore + tempEkgTopoScore)/2;
% Zero out scores for components that didn't meet the initial "likely EKG" threshold
tempEkgScore(~tempNoEkg) = 0;
% Final list of automatically flagged EKG sources (high combined score)
ekgSrc = tempNSources(tempEkgScore > 2.5);
% List of suspected EKG sources that were below the final threshold but above suspicion (for review)
tempEkgSus = tempNSources((tempEkgScore > 1.96) & ~ekgSrc);

% 4) EKG UI - Interactive review of flagged and suspected sources
ekgSrc = reviewEkgArtifactUI(W, xICA, fs, chanlocs, ekgSrc, tempEkgSus, ...
   tempStartSec, tempNSec, [INFO.fileName '_EKG'],...
   saveFigs, INFO.figDir);

if isempty(ekgSrc); fprintf('\n\tNo EKG artifact detected\n')
else; fprintf('\n\tekgSrc = %s\n', mat2str(ekgSrc));
end


%% 5. OTHER SOURCES

% -------------------------------------------------------------------------
% Use this block ONLY for components not caught by the automated EKG,
% or DIN detectors that are CLEARLY artifactual
% -------------------------------------------------------------------------

otherSrc = [];
% Identify all components NOT already flagged for removal from the first 10 ICs
tempOtherSrc = setdiff((1:10).',[ekgSrc(:); dinSrc(:)]);

% Use the updated, save-enabled UI for review. Note: This UI starts components 
% flagged for KEEP (green), as rejection here should be rare.
otherSrc = reviewOtherUI(W, xICA, fs, chanlocs, otherSrc, tempOtherSrc, ...
    [INFO.fileName '_Other'], saveFigs, INFO.figDir); % Figure 6

clear temp*


%% 6. ARTIFACT REMOVAL AND SIGNAL RECONSTRUCTION

% -------------------------------------------------------------------------
% The components flagged for removal are zeroed out in ICA space, and the
% signal is converted back to clean sensor space (xCl).
% -------------------------------------------------------------------------

% Combine all components flagged for removal
rmSrc = [ekgSrc(:); dinSrc(:); otherSrc(:)];
disp([newline 'Total Components to Remove (rmSrc) = ' mat2str(rmSrc(:)')]);

% Zero out the artifact components in the ICA data matrix
xICA(rmSrc,:) = 0;

% Convert the data back to the clean sensor space
xCl = inv(W) * xICA;

% Store component lists for documentation
INFO.W = W; INFO.EkgSrc = ekgSrc; INFO.DinSrc = dinSrc;
INFO.OtherSrc = otherSrc; INFO.RmSrc = rmSrc;


%% 7. QUALITY CONTROL: BEFORE AND AFTER PLOTS

% -------------------------------------------------------------------------
% Visualize the raw (xRaw) vs. cleaned (xCl) data to assess the impact of
% ICA removal. This is a critical quality check.
% -------------------------------------------------------------------------

% Fill bad channels with NaNs for clean plotting
tempX124_before = fillBadChRows(xRaw, badCh);
tempX124_after = fillBadChRows(xCl, badCh);

% --- Plotting: Time Series and Data Image ---
figure()
sgtitle([INFO.fileName ' : Before and After ICA Cleaning'], 'interpreter', 'none')

% Left Column: Before ICA (xRaw)
subplot(4, 2, 1); plotEEGOverlay(tempX124_before, []); title('Before ICA, Overlay');
subplot(4, 2, [3 5]); imagesc(abs(tempX124_before)); title('Before ICA, Data Image (abs)');

% Right Column: After ICA (xCl)
subplot(4, 2, 2); plotEEGOverlay(tempX124_after, []); title('After ICA, Overlay');
subplot(4, 2, [4 6]); imagesc(abs(tempX124_after)); title('After ICA, Data Image (abs)');

% --- Quality Check: Percentage of Outliers ---
% A good cleaning should significantly reduce the percentage of samples
% exceeding the high-voltage threshold (INFO.recUVThresh).
tempRecOverUVThresh = computePerChannelPctOverThresh(tempX124_after, INFO.recUVThresh);
subplot(4, 2, [7 8]);
stem(1:124, tempRecOverUVThresh); hold on; grid on
title(['Post-ICA: % of abs values exceeding ' num2str(INFO.recUVThresh) '\muV'])

% Highlight channels exceeding the recording percentage threshold (INFO.recPctThresh)
if any(tempRecOverUVThresh >= INFO.recPctThresh)
    yline(INFO.recPctThresh, 'r', 'linewidth', 1.5);
    tempChOver = find(tempRecOverUVThresh >= INFO.recPctThresh);
    stem(tempChOver, tempRecOverUVThresh(tempChOver), 'r')
    warning('One or more channels still exceed the recording bad channel threshold after ICA.')
else
    disp('- * Post-ICA QC check complete: No recording bad channels! * -')
end
clear temp*