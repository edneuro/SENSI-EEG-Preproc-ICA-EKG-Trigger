function [dinSrc, info] = autoDetectDIN(xICA, fs, nSources, Opts)
% AUTODETECTDIN  Detect DIN (trigger) artifact ICA components via harmonic
%                energy ratio and z-score ranking.
%
%   [DINSRC, INFO] = AUTODETECTDIN(XICA, FS, NSOURCES, OPTS)
%
%   Identifies ICA components that exhibit the comb-like harmonic structure
%   characteristic of periodic digital-input (DIN/trigger) artifacts. For
%   each component, a harmonic energy ratio is computed via
%   harmRatioFromSignal. Components are z-scored across nSources and flagged
%   in two stages: high-confidence DIN (|z| > thress) are flagged first,
%   then a second z-score pass (with flagged components replaced by the
%   median) identifies additional borderline candidates to present for review.
%
%   INPUTS
%     xICA      [nComp x nSamp] ICA component time series
%     fs        Scalar sampling rate (Hz)
%     nSources  Vector of component indices to scan (default: 1:30)
%     Opts      Struct with optional fields (defaults shown):
%       .T      1      DIN repetition period (sec); sets fundamental freq
%       .fmin   10     Lower bound of DIN detection band (Hz)
%       .fmax   50     Upper bound of DIN detection band (Hz)
%       .bw     0.2    Bandwidth parameter passed to harmRatioFromSignal
%
%   OUTPUTS
%     dinSrc    Row vector of component indices flagged as DIN artifacts
%               (up to MAXK=3 components, ranked by z-score)
%     info      Struct with fields:
%       .dinCheck         [nComp x 1] raw harmonic ratio per component (NaN for unscanned)
%       .z                [nComp x 1] stage-1 z-scores (NaN for unscanned)
%       .dinAdj           [nComp x 1] dinCheck with flagged components replaced by median
%       .z2               [nComp x 1] stage-2 z-scores after median replacement
%       .reviewCandidates Vector of stage-2 candidates (for interactive review)
%       .topNonDIN        Highest-ranked stage-2 candidate (or [] if none)
%
%   METHOD
%     1) harmRatioFromSignal computes the ratio of harmonic energy to
%        background energy in the DIN band for each component in nSources.
%     2) Stage-1: z-score across nSources; components with z > thress (=2.5)
%        flagged. Top MAXK (=3) by score selected as dinSrc.
%     3) Stage-2: flagged components replaced by median, z-scores recomputed.
%        Additional candidates with z2 > thress (not already flagged) are
%        collected as reviewCandidates for the interactive review UI.
%
%   DEPENDENCIES
%     harmRatioFromSignal
%
%   NOTES
%     • thress=2.5 and MAXK=3 are empirical defaults; most recordings have
%       0–1 DIN component; MAXK=3 accommodates unusual cases.
%     • Components not in nSources receive NaN in all output score fields.
%     • If no valid sources exist, dinSrc=[] and a warning is issued.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Reference (please cite):
%
% Module Citation:
% Malave, A. J., & Kaneshiro, B. (2026). SENSI-EEG-Preproc-ICA-EKG-Trigger:
% A MATLAB framework for semi-automated identification of EKG and trigger
% artifacts in EEG using ICA and spectral characteristics. Stanford University.
% https://github.com/edneuro/SENSI-EEG-Preproc-ICA-EKG-Trigger
%
% MIT License
%
% Copyright (c) 2026 Amilcar J. Malave, and Blair Kaneshiro.
%
% Permission is hereby granted, free of charge, to any person obtaining a
% copy of this software and associated documentation files (the "Software"),
% to deal in the Software without restriction, including without limitation
% the rights to use, copy, modify, merge, publish, distribute, sublicense,
% and/or sell copies of the Software, and to permit persons to whom the
% Software is furnished to do so, subject to the following conditions:
%
% The above copyright notice and this permission notice shall be included in
% all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS
% OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
% THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
% FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
% DEALINGS IN THE SOFTWARE.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% ---- Defaults
if nargin < 3 || isempty(nSources),    nSources    = 1:30; end
if nargin < 4 || isempty(Opts),        Opts     = struct('fmin',10,'fmax',50,'bw',0.2); end

% z-score threshold for flagging DIN components (empirical)
thress = 2.5;
% Maximum number of DIN components to flag across both stages
MAXK = 3;

% ---- Shapes / guards
[nComp, nSamp] = size(xICA);
nSources = unique(nSources(:)');
nSources = nSources(nSources>=1 & nSources<=nComp);
if isempty(nSources)
    dinSrc = [];
    info   = struct('dinCheck',[],'z',[],'topNonDIN',[]);
    warning('autoDetectDIN:NoValidSources','No valid sources to scan.');
    return;
end

% ---- Compute DIN metric (harmonic energy ratio per component)
dinCheck = nan(nComp,1);
for k = nSources
    out = harmRatioFromSignal(xICA(k,:), fs, Opts);
    dinCheck(k) = out.ratio_linear;
end

% ---- Stage-1: z-score over selected sources and flag high-confidence DIN
vals  = dinCheck(nSources);
zVals = zscore(vals, 0, 'omitnan');

z = nan(nComp,1);
z(nSources) = zVals;

initialMask = zVals > thress;
flag1   = nSources(initialMask);
score1  = zVals(initialMask);

% Sort stage-1 flags by score descending
if ~isempty(flag1)
    [~, ord1] = sort(score1, 'descend');
    flag1  = flag1(ord1);
end

% Take up to MAXK confirmed DIN components
dinSrc   = flag1(1:min(MAXK, numel(flag1))).';
nRemain  = MAXK - numel(dinSrc);

% ---- Stage-2: replace flagged values with median and re-score remaining
% This exposes borderline candidates that were masked by the strong DIN components.
medVal = median(dinCheck(nSources),'omitnan');
dinAdj = dinCheck;
dinAdj(flag1) = medVal;

vals2 = dinAdj(nSources);
z2    = zscore(vals2, 0, 'omitnan');

cand2     = [];
topNonDIN = [];

if nRemain > 0
    candMask  = (z2 > thress) & ~initialMask;
    cand2_all = nSources(candMask);
    score2    = z2(candMask);

    % Sort stage-2 candidates by score descending
    if ~isempty(cand2_all)
        [~, ord2] = sort(score2, 'descend');
        cand2_all = cand2_all(ord2);
        cand2     = cand2_all(1:min(nRemain, numel(cand2_all))).';
        topNonDIN = cand2(1);
    end
end

% ---- Package output info struct
info = struct('dinCheck',dinCheck,'z',z,'dinAdj',dinAdj,'z2',z2, ...
              'reviewCandidates',cand2,'topNonDIN',topNonDIN);

end
