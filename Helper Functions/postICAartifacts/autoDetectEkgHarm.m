function [compositeScore, metrics] = autoDetectEkgHarm(xICA, fs, nSources, opts)
% AUTODETECTEKG  Flag EKG components via multi‐criterion spectral analysis.
%
% Usage:
%   [ekgSrc, metrics, compositeScore] = autoDetectEKG(xICA, fs)
%   [ekgSrc, metrics, compositeScore] = autoDetectEKG(xICA, fs, nSources, opts)
%
% Inputs:
%   xICA      : [nComp x nSamp] ICA time series
%   fs        : sampling rate in Hz
%   nSources  : vector of components to scan (default 1:nComp)
%   opts      : struct with fields:
%       .fmin         = 0.8;      % band lower bound (Hz)
%       .fmax         = 2.0;      % band upper bound (Hz)
%       .harmonics    = 3;        % number of harmonics to include

% Legacy
%       .hrRange      = [0.7 1.8];% plausible heart‐rate band (Hz)
%       .maxFlag      = 4;        % how many to return at most
%

% Outputs:
%   ekgSrc         : flagged EKG component indices (1×M, M≤maxFlag)
%   metrics        : struct with PSD, HarmonicRatio, fPeak
%   compositeScore : nComp×1 composite z‐score

% ---- Defaults and input checks
[nComp, nSamp] = size(xICA);
if nargin < 3 || isempty(nSources),   nSources = 1:nComp; end
if nargin < 4, opts = struct(); end
opts = fillDefault(opts, ...
    'fmin',0.8, 'fmax',2.0, 'harmonics',3, ...
    'hrRange',[0.7 1.8]);

nSources = nSources(nSources>=1 & nSources<=nComp);

% FFT parameters
nFFT = 2^nextpow2(nSamp);
f    = linspace(0, fs/2, nFFT/2+1)';

%% Compute PSDs
metrics.psdFreq = f;
metrics.psdPow  = zeros(numel(f), nComp);

for k = 1:nComp
    Y = abs(fft(xICA(k,:), nFFT));
    P = (Y.^2)/nSamp;
    metrics.psdPow(:,k) = P(1:numel(f));

    % if k == 9
    %     figure;
    %     subplot(2,1,1)
    %     plot(f,Y(1:numel(f)))
    %     xlim([0,15])
    %     xlabel('Frequency (Hz)');
    %     ylabel('Amplitude |Y(f)|');
    %     title('Magnitude Spectrum of Signal');
    %     subplot(2,1,2)
    %     plot(f,P(1:numel(f)))
    %     xlim([0,15])
    %     xlabel('Frequency (Hz)');
    %     ylabel('Power |Y(f)|^2');
    %     title('Power Spectrum of Signal');
    % end
end

%% (2) Harmonic ratio
H = nan(nComp,1);
for k = nSources
    Pk = metrics.psdPow(:,k);
    idx = f>=opts.fmin & f<=opts.fmax;
    [~, pi] = max(Pk(idx));
    f1 = f(idx); f1 = f1(pi);
    harmPow = 0;
    for h=1:opts.harmonics
        [~, hi] = min(abs(f - h*f1));
        harmPow = harmPow + Pk(hi);
    end
    H(k) = harmPow / max(sum(Pk),eps);
end
metrics.HarmonicRatio = H;

%% Record peak frequency in cardiac band
fPk = nan(nComp,1);
for k = nSources
    Pk = metrics.psdPow(:,k);
    idx = f>=opts.fmin & f<=opts.fmax;
    [~, pi] = max(Pk(idx));
    fBand = f(idx);
    fPk(k) = fBand(pi);
end
metrics.fPeak = fPk;

%% Composite z‐score
fields = {'HarmonicRatio'};
Z = zeros(nComp,numel(fields));
for i=1:numel(fields)
    v    = metrics.(fields{i});
    med  = nanmedian(v(nSources));
    madv = 1.4826*nanmedian(abs(v(nSources)-med));
    Z(:,i) = (v - med) ./ max(madv,eps);
end
compositeScore = sum(Z,2);

% %% Candidate selection
% % 1) top by composite
% [~, ord] = sort(compositeScore(nSources), 'descend');
% cand    = nSources(ord);
% 
% % 2) filter by peak freq in hrRange
% inRange = cand(metrics.fPeak(cand) >= opts.hrRange(1) & ...
%                metrics.fPeak(cand) <= opts.hrRange(2));
% 
% if ~isempty(inRange)
%     [~, idx] = max(compositeScore(inRange));
%     ekgSrc = inRange(idx);
% else
%     ekgSrc = cand(1);
% end
% 
% % cap to maxFlag
% ekgSrc = ekgSrc(1:min(numel(ekgSrc), opts.maxFlag));

end

%% Helper to fill defaults
function s = fillDefault(s, varargin)
for i=1:2:numel(varargin)
    if ~isfield(s,varargin{i}) || isempty(s.(varargin{i}))
        s.(varargin{i}) = varargin{i+1};
    end
end
end
