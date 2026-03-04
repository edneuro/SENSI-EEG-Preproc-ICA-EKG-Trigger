function [ekgSrc,ekgSus, metrics] = autoDetectEkg(xICA, fs, nSources, opts)
% AUTODETECTEKG  Score ICA components for EKG via harmonic local-SNR.
%
% Inputs:
%   xICA      : [nComp x nSamp] ICA time series
%   fs        : sampling rate (Hz)
%   nSources  : components to scan (default 1:nComp)
%   opts      : struct with fields (defaults shown):
%       .window        = 14;      % Welch window length (sec)
%       .overlap     = .window/2; % Welch overlap (sec)
%       .fmin        = 0.8;    % fundamental search band (Hz)
%       .fmax        = 2.0;    % fundamental search band (Hz)
%       .harmonics   = 4;      % number of harmonics to score (includes h=1)
%       .peakHalfHz  = 0.20;   % half-width of peak window around each harmonic (Hz)
%       .nbHalfHz    = 0.60;   % half-width of neighbor ring around each harmonic (Hz)
%       .peakAgg     = 'sum';  % 'max' or 'sum' for peak window signal
%       .do_detrend   = 1;      % spectral flattening via log-log linear detrend (1/f removal)
%       .detrend_fmin = .5;     % min freq for spectral flattening fitting
%       .detrend_fmax = 30;     % max freq for spectral flattening fitting
% 
%
% Outputs:
%   ekgSrc: sources likely to be EKG artifacts
%   ekgSus: suspecious sources
%   metrics        : struct with fields:
%       .compositeScore : [nComp x 1] robust z-score of HarmonicSNRsum across nSources
%       .psdFreq      : frequency axis (Hz)
%       .psdPow       : Welch spectra [nFreq x nComp]
%       .fPeak        : fundamental peak freq in [fmin,fmax] (Hz)
%       .harmFreq     : [nComp x H] expected harmonic freqs (Hz)
%       .harmSNR      : [nComp x H] local SNR per harmonic (ratio, not dB)
%       .harmSNRsum   : [nComp x 1] sum of harmonic SNRs (optionally excluding h=1)
%

% ---- Defaults / checks
[nComp, ~] = size(xICA);
if nargin < 3 || isempty(nSources), nSources = 1:nComp; end
if nargin < 4, opts = struct(); end
opts = fillDefault(opts, ...
    'window',14, ...
    'fmin',0.8, 'fmax',2.0, ...
    'harmonics',4, ...
    'peakHalfHz',0.20, 'nbHalfHz',0.60, ...
    'peakAgg','sum',...)
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
% We'll store each component spectrum. Note: pwelch returns one-sided by default.
Ptmp = cell(nComp,1);
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

    Ptmp{k} = Pxx(:);
end

nF = numel(f);
metrics.psdFreq = f;
metrics.psdPow  = zeros(nF, nComp);
for k = 1:nComp
    metrics.psdPow(:,k) = Ptmp{k};
end

% ---- Find base frequency (fundamental) in [fmin,fmax]
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



% Select EKG and suspected EKG (ranked by harmonic SNR, max 3 total)
thrEkg = 5; thrSus = 3; maxTot = 3;

scores = compositeScore(nSources);
% ---- Rank all candidates by score (descending)
cand = nSources(scores > thrSus);
[~, ord] = sort(compositeScore(cand), 'descend');
cand = cand(ord);
% ---- EKG = top ones above thrEkg
ekgSrc = cand(compositeScore(cand) > thrEkg);
% Cap ekgSrc to maxTot
ekgSrc = ekgSrc(1:min(numel(ekgSrc), maxTot));
% ---- Suspected = fill remaining slots from the rest
nRemain = maxTot - numel(ekgSrc);
if nRemain > 0
    rest = setdiff(cand, ekgSrc, 'stable');
    ekgSus = rest(1:min(numel(rest), nRemain));
else
    ekgSus = [];
end


% Debuggin - Plot Strong Ones
Nekg = length(ekgSrc);
figure;
for i = 1:Nekg
    subplot(Nekg,1,i)
    plot(f,Ptmp{ekgSrc(i)})
    title(sprintf('Source %d', ekgSrc(i)))
    xlabel('Freq (Hz)')
    ylabel('Amplitude')
    xlim([0,12])
    sgtitle('EKG Sources - Pwelch and 1/f flattened')
end

end

%% Helper to fill defaults
function s = fillDefault(s, varargin)
for i = 1:2:numel(varargin)
    if ~isfield(s, varargin{i}) || isempty(s.(varargin{i}))
        s.(varargin{i}) = varargin{i+1};
    end
end
end


% Legacy

% % 3) Select EKG and suspected EKG (ranked by harmonic SNR, max 3 total)
% thrEkg = 5; thrSus = 3; maxTot = 3;
% scores = tempEkgHarmSnrScore(tempNSources);
% % ---- Rank all candidates by score (descending)
% cand = tempNSources(scores > thrSus);
% [~, ord] = sort(tempEkgHarmSnrScore(cand), 'descend');
% cand = cand(ord);
% % ---- EKG = top ones above thrEkg
% ekgSrc = cand(tempEkgHarmSnrScore(cand) > thrEkg);
% % Cap ekgSrc to maxTot
% ekgSrc = ekgSrc(1:min(numel(ekgSrc), maxTot));
% % ---- Suspected = fill remaining slots from the rest
% nRemain = maxTot - numel(ekgSrc);
% if nRemain > 0
%     rest = setdiff(cand, ekgSrc, 'stable');
%     tempEkgSus = rest(1:min(numel(rest), nRemain));
% else
%     tempEkgSus = [];
% end