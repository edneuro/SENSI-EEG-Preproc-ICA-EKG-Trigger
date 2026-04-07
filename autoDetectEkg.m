function [ekgSrc,ekgSus, metrics] = autoDetectEkg(xICA, fs, nSources, opts)
% AUTODETECTEKG  Score ICA components for EKG artifacts via harmonic local-SNR.
%
%   [EKGSRC, EKGSUS, METRICS] = AUTODETECTEKG(XICA, FS, NSOURCES, OPTS)
%
%   Identifies ICA components that exhibit the periodic harmonic structure
%   characteristic of cardiac (EKG) artifacts. For each component, a Welch
%   power spectrum is computed (with optional 1/f flattening), the dominant
%   frequency in the cardiac band [fmin, fmax] is located, and a local SNR
%   is measured at each harmonic. Components are ranked by summed harmonic
%   SNR and flagged as confirmed EKG (score > thrEkg) or suspected
%   (score > thrSus), capped at maxTot total components.
%
%   INPUTS
%     xICA      [nComp x nSamp] ICA component time series (DC-removed internally)
%     fs        Scalar sampling rate (Hz)
%     nSources  Vector of component indices to scan (default: 1:nComp)
%     opts      Struct with optional fields (defaults shown):
%       .window        14      Welch window length (sec)
%       .overlap       window/2  Welch overlap (sec)
%       .fmin          0.8    Lower bound of cardiac rhythm search band (Hz)
%       .fmax          2.0    Upper bound of cardiac rhythm search band (Hz)
%       .harmonics     4      Number of harmonics to score (includes h=1)
%       .peakHalfHz    0.20   Half-width of peak window around each harmonic (Hz)
%       .nbHalfHz      0.60   Half-width of neighbor baseline ring (Hz)
%       .peakAgg       'sum'  Aggregation over peak window: 'max' or 'sum'
%       .do_detrend    1      Apply log-log linear detrend for 1/f removal
%       .detrend_fmin  0.5    Lower bound for 1/f fit (Hz)
%       .detrend_fmax  30     Upper bound for 1/f fit (Hz)
%
%   OUTPUTS
%     ekgSrc    Row vector of component indices flagged as EKG artifacts
%     ekgSus    Row vector of suspected (borderline) EKG component indices
%     metrics   Struct with fields:
%       .compositeScore  [nComp x 1] robust z-score of harmSNRsum over nSources
%       .psdFreq         Frequency axis (Hz)
%       .psdPow          Welch spectra [nFreq x nComp] (after optional detrend)
%       .fPeak           Dominant frequency in [fmin,fmax] per component (Hz)
%       .harmFreq        [nComp x H] expected harmonic frequencies (Hz)
%       .harmSNR         [nComp x H] local SNR per harmonic (linear ratio)
%       .harmSNRsum      [nComp x 1] summed harmonic SNR
%
%   METHOD
%     1) Welch power spectrum computed per component; optional 1/f flattening
%        via log-log linear regression over [detrend_fmin, detrend_fmax].
%     2) Dominant frequency f1 located in [fmin, fmax].
%     3) Local SNR at each harmonic h*f1: peak window signal / median of
%        neighbor ring (±nbHalfHz excluding ±peakHalfHz).
%     4) Composite score = robust z-score of harmSNRsum over nSources.
%     5) Components with score > thrEkg (=5) flagged as ekgSrc;
%        score > thrSus (=3) flagged as ekgSus. Total capped at maxTot (=3).
%
%   DEPENDENCIES
%     fillDefault (local helper)
%
%   NOTES
%     • Components are DC-removed before spectral analysis.
%     • thrEkg=5, thrSus=3, maxTot=3 are empirical defaults tuned for
%       typical resting EEG recordings; adjust if cardiac rate is unusual.
%     • Components with no valid harmonic windows receive NaN harmSNR.

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

% ---- Defaults / checks
[nComp, ~] = size(xICA);
if nargin < 3 || isempty(nSources), nSources = 1:nComp; end
if nargin < 4, opts = struct(); end
opts = fillDefault(opts, ...
    'window',14, ...
    'fmin',0.8, 'fmax',2.0, ...
    'harmonics',4, ...
    'peakHalfHz',0.20, 'nbHalfHz',0.60, ...
    'peakAgg','sum',...
    'do_detrend', 1,...
    'detrend_fmin', 0.5, ...
    'detrend_fmax', 30);
opts.overlap = round(opts.window/2);


nSources = nSources(nSources>=1 & nSources<=nComp);

% Welch params
winN = max(8, round(opts.window*fs));
ovN  = max(0, min(winN-1, round(opts.overlap*fs)));

% ---- Center ICA components (remove DC per component)
xICA = xICA - mean(xICA, 2);

% ---- Compute Welch spectra
% Store each component spectrum in a cell array, then assemble into matrix.
PxxCell = cell(nComp,1);
f = [];
for k = 1:nComp
    x = xICA(k,:);
    [Pxx, fWelch] = pwelch(x, hann(winN), ovN, [], fs);

    if isempty(f), f = fWelch(:); end

    % Optional: spectral flattening via log-log linear detrend (1/f removal)
    if opts.do_detrend
        fitIdx = (f >= opts.detrend_fmin) & (f <= opts.detrend_fmax);
        fitIdx = fitIdx & isfinite(Pxx) & (Pxx > 0) & isfinite(f) & (f > 0);
        if nnz(fitIdx) >= 3
            lf = log10(f(fitIdx));
            lp = log10(Pxx(fitIdx) + eps);
            p = polyfit(lf, lp, 1);                 % lp ≈ p(1)*lf + p(2)
            trend = polyval(p, log10(f + eps));     % extend to all f
            Pxx = 10.^(log10(Pxx + eps) - trend);   % flattened spectrum (~baseline 1)
        end
    end

    PxxCell{k} = Pxx(:);
end

nF = numel(f);
metrics.psdFreq = f;
metrics.psdPow  = zeros(nF, nComp);
for k = 1:nComp
    metrics.psdPow(:,k) = PxxCell{k};
end

% ---- Find base frequency (fundamental) in [fmin, fmax]
fPk = nan(nComp,1);
idxBand = (f >= opts.fmin) & (f <= opts.fmax);

for k = nSources
    Pk = metrics.psdPow(:,k);
    if any(idxBand)
        [~, pi] = max(Pk(idxBand));
        fBand = f(idxBand);
        fPk(k) = fBand(pi);
    end
end
metrics.fPeak = fPk;

% ---- Harmonic SNR per component
H = opts.harmonics;
harmFreq = nan(nComp, H);
harmSNR  = nan(nComp, H);
harmSNRsum = nan(nComp,1);

useMax = strcmpi(opts.peakAgg,'max');
useSum = strcmpi(opts.peakAgg,'sum');

for k = nSources
    Pk = metrics.psdPow(:,k);
    f1 = metrics.fPeak(k);
    if isnan(f1) || f1<=0, continue; end

    % only include harmonics <= Nyquist
    maxH = min(H, floor((fs/2)/f1));
    if maxH < 1, continue; end

    snrVec = nan(1,H);

    for h = 1:maxH
        fh = h*f1;
        harmFreq(k,h) = fh;

        % Peak window
        dp = opts.peakHalfHz;
        idxPeak = (f >= fh-dp) & (f <= fh+dp);

        % Neighbor ring: ±nbHalfHz excluding the peak window
        nb = opts.nbHalfHz;
        idxNbr = ((f >= fh-nb) & (f <  fh-dp)) | ...
                 ((f >  fh+dp) & (f <= fh+nb));

        % Guard: if too close to 0 Hz or Nyquist, neighbor ring may be empty
        if ~any(idxPeak) || ~any(idxNbr)
            continue;
        end

        % Signal from peak window
        if useMax
            sigVal = max(Pk(idxPeak));
        else
            sigVal = sum(Pk(idxPeak));
        end

        % Baseline from neighbor ring
        baseVal = median(Pk(idxNbr));

        snrVec(h) = sigVal / (baseVal + eps);   % ratio SNR (not dB)
        harmSNR(k,h) = snrVec(h);
    end

    % Sum harmonic SNRs (optionally exclude the fundamental, h=1)
    harmSNRsum(k) = nansum(snrVec);

end

metrics.harmFreq    = harmFreq;
metrics.harmSNR     = harmSNR;
metrics.harmSNRsum  = harmSNRsum;

% ---- Composite score: robust z-score of harmSNRsum
v = metrics.harmSNRsum;
med  = nanmedian(v(nSources));
madv = 1.4826 * nanmedian(abs(v(nSources) - med));
compositeScore = (v - med) ./ max(madv, eps);

metrics.compositeScore  = compositeScore;


% ---- Select EKG and suspected EKG components
% thrEkg=5, thrSus=3: empirical thresholds on robust z-score of harmSNRsum.
% maxTot=3: cap on total flagged components to avoid over-rejection.
thrEkg = 5; thrSus = 3; maxTot = 3;

scores = compositeScore(nSources);
% Rank all candidates above thrSus by score (descending)
cand = nSources(scores > thrSus);
[~, ord] = sort(compositeScore(cand), 'descend');
cand = cand(ord);
% EKG = top candidates above thrEkg
ekgSrc = cand(compositeScore(cand) > thrEkg);
% Cap ekgSrc to maxTot
ekgSrc = ekgSrc(1:min(numel(ekgSrc), maxTot));
% Suspected = fill remaining slots from borderline candidates
nRemain = maxTot - numel(ekgSrc);
if nRemain > 0
    rest = setdiff(cand, ekgSrc, 'stable');
    ekgSus = rest(1:min(numel(rest), nRemain));
else
    ekgSus = [];
end

end

%% Helper to fill default struct fields
function s = fillDefault(s, varargin)
for i = 1:2:numel(varargin)
    if ~isfield(s, varargin{i}) || isempty(s.(varargin{i}))
        s.(varargin{i}) = varargin{i+1};
    end
end
end
