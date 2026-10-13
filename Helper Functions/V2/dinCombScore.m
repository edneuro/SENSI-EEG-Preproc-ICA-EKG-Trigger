function out = dinCombScore(xICA, fs, opts)
% DINCOMBSCORE  Comb score of every ICA component for a trigger repeating every T seconds.
%
%   OUT = DINCOMBSCORE(XICA, FS, OPTS)
%
%   A trigger artifact that repeats every T s puts a comb of narrow peaks at k/T Hz into a
%   component's magnitude spectrum. The comb score weighs that spectrum directly:
%
%       comb = w . mhat
%
%   mhat is the magnitude spectrum over [fmin, fmax] scaled to unit mean (so loudness
%   cannot win on its own); w is +raised-cosine on each harmonic (half-width 0.1*f0),
%   0 in a guard ring out to 0.2*f0, and flat negative everywhere else, with the positive
%   and the negative weights each summing to 1. A flat spectrum scores exactly 0.
%
%   Unlike a harmonic-to-midpoint RATIO, a strong off-harmonic line (on the SENSI
%   Session 1 cohort, ~23.5 Hz) can only subtract a bounded amount: it is one band
%   among hundreds in the negative weights, rather than the denominator.
%
%   combC, the default of autoDetectDINV2, removes even that bounded amount. The same
%   weights, but the bins are first split into PEAKS (within 0.1*f0 of a harmonic) and
%   VALLEYS (beyond 0.2*f0), and within EACH set the bins whose log10 magnitude lies
%   more than clipK robust z (median, 1.4826 MAD) from that set's median are dropped:
%
%       combC = weighted mean(kept peaks) - mean(kept valleys)
%
%   An interfering line between harmonics is dropped from the valleys instead of
%   lowering the score; a line that happens to sit on a harmonic is dropped from the
%   peaks instead of raising it. On the cohort below this took the comb's last missed
%   DIN component from rank 47 to rank 1 of its block-group.
%
%   INPUTS
%     xICA   [nComp x nSamp] ICA component time series
%     fs     sampling rate (Hz)
%     opts   (optional struct, defaults shown)
%       .T         1        trigger period (s); f0 = 1/T
%       .fmin      10       analysis band (Hz)
%       .fmax      50
%       .segSec    20       Welch segment (s), Hann, 50% overlap, nfft = segment;
%                           Inf = one Hann periodogram of the whole record
%       .clipK     6        combC: bins beyond K robust z (log10) dropped, per set
%       .chunkSec  []       also score each chunk of this length (for a time course)
%       .chunkScore 'combC' which score the per-chunk time course uses
%
%   OUTPUT (struct)
%     .comb    [nComp x 1]  the comb score (higher = more trigger-like; 0 = flat)
%     .combR   [nComp x 1]  harmonic-weighted / mean off-harmonic magnitude (1 = flat)
%     .combC   [nComp x 1]  comb with peaks and valleys clipped separately (0 = flat)
%     .chunk   [nComp x nChunk] chunkScore per chunk (only if opts.chunkSec), .edges
%     .df      frequency resolution achieved (Hz)
%
%   NOTES
%     • The harmonics must be resolvable: a Welch segment shorter than ~10*T leaves too
%       few bins between harmonics and the score degrades; NaN when no bin lies in the
%       negative region at all. Use segSec = Inf (whole record) for long periods.
%
%   See also autoDetectDINV2, autoDetectDIN, harmRatioFromSignal

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

if nargin < 3 || isempty(opts), opts = struct(); end
opts = fillDefault(opts, 'T',1, 'fmin',10, 'fmax',50, 'segSec',20, 'clipK',6, 'chunkSec',[], ...
    'chunkScore','combC');
[nComp, N] = size(xICA);
f0 = 1/opts.T;

[f, mag, df] = magSpectrum(xICA, fs, opts.segSec, opts.fmin, opts.fmax);
[out.comb, out.combR, out.combC] = combFromMag(f, mag, f0, opts.clipK);
out.df = df;

if ~isempty(opts.chunkSec)
    chunkN = round(opts.chunkSec*fs);
    st = 1:chunkN:max(1, N-chunkN+1);
    ed = min(st+chunkN-1, N);
    keep = (ed-st+1) >= round(10*fs); st = st(keep); ed = ed(keep);
    out.edges = [st(:) ed(:)];
    out.chunk = nan(nComp, numel(st));
    seg = min(opts.segSec, opts.chunkSec/2);           % >= 3 Welch segments per chunk
    for q = 1:numel(st)
        [fq, mq] = magSpectrum(xICA(:, st(q):ed(q)), fs, seg, opts.fmin, opts.fmax);
        [c1, c2, c3] = combFromMag(fq, mq, f0, opts.clipK);
        out.chunk(:,q) = struct('comb',c1,'combR',c2,'combC',c3).(opts.chunkScore);
    end
end
out.opts = opts;
end

% -------------------------------------------------------------------------
function [f, mag, df] = magSpectrum(x, fs, segSec, fmin, fmax)
[nComp, N] = size(x);
if isinf(segSec)
    nfft = 2^nextpow2(N); w = hann(N);
    M = floor(nfft/2)+1; fa = (0:M-1)' * (fs/nfft); b = fa >= fmin & fa <= fmax;
    mag = zeros(nComp, nnz(b));
    for c0 = 1:16:nComp
        cc = c0:min(c0+15, nComp);
        X = fft(x(cc,:).' .* w, nfft);
        mag(cc,:) = abs(X(b,:)).';
    end
    f = fa(b).'; df = fs/nfft;
else
    segN = min(round(segSec*fs), N);
    [P, fa] = pwelch(x.', hann(segN), floor(segN/2), segN, fs);
    b = fa >= fmin & fa <= fmax;
    f = fa(b).'; mag = sqrt(P(b,:)).'; df = fs/segN;
end
end

function [comb, combR, combC] = combFromMag(f, mag, f0, K)
mag = mag ./ mean(mag, 2);
d = abs(f/f0 - round(f/f0)) * f0;                  % distance to the nearest harmonic (Hz)
isPos = d <= 0.1*f0; isNeg = d > 0.2*f0;
nC = size(mag,1);
comb = nan(nC,1); combR = nan(nC,1); combC = nan(nC,1);
if ~any(isPos) || ~any(isNeg), return; end
wp = 0.5 * (1 + cos(pi * d(isPos) / (0.1*f0)));   % raised cosine on each harmonic
num = mag(:, isPos) * wp.' / sum(wp);              % weighted mean of the peaks
den = mean(mag(:, isNeg), 2);                      % mean of the valleys
comb  = num - den;
combR = num ./ max(den, eps);
for c = 1:nC
    mp = mag(c, isPos); mn = mag(c, isNeg);
    kp = keepSet(mp, K); kn = keepSet(mn, K);
    combC(c) = sum(wp(kp).*mp(kp))/sum(wp(kp)) - mean(mn(kn));
end
end

function keep = keepSet(v, K)
% Bins of one set (peaks or valleys) within K robust z of that set's median, on log10.
l = log10(max(v, realmin)); md = median(l); mad = 1.4826 * median(abs(l - md));
keep = true(size(v));
if mad > 0, keep = abs(l - md) / mad <= K; end
if ~any(keep), keep = true(size(v)); end
end

function s = fillDefault(s, varargin)
for i = 1:2:numel(varargin)
    if ~isfield(s, varargin{i}) || isempty(s.(varargin{i})), s.(varargin{i}) = varargin{i+1}; end
end
end
