function out = harmRatioFromSignal(x, fs, opts)

% harmRatioFromSignal  Quantify harmonic energy at a repetition rate vs inter-harmonic energy.
%
% This is an SNR-like metric in the frequency domain:
%   "signal"  = power near expected harmonics of f0
%   "baseline"= power between harmonics (excluding guard zones near harmonics)
%
% Usage:
%   out = harmRatioFromSignal(x, fs)
%   out = harmRatioFromSignal(x, fs, opts)
%
% Inputs
%   x    : time series (vector)
%   fs   : sampling rate (Hz)
%   opts : (struct) optional overrides:
%          .T      = [];          % din repetition period in seconds. Default: 1 sec
%          .fmin   = 10;          % analysis band low (Hz)
%          .fmax   = 50;          % analysis band high (Hz)
%
%          % Band definitions (choose either absolute Hz or fractional-of-f0):
%          .bw     = [];          % half-bandwidth around centers (Hz). Can be left empty - value calculated inside function.
%          .guard  = [];          % exclude this half-width near harmonics for "middle" bands (Hz). Can be left empty - value calculated inside function.
%          .bwFrac    = 0.15;     % bw = bwFrac*f0 (used when bw is empty)
%          .guardFrac = 0.08;     % guard = guardFrac*f0 (used when guard is empty)
%          .minBW  = 0.05;        % clamp for scaled bw (Hz)
%          .maxBW  = 0.40;        % clamp for scaled bw (Hz)
%          .minGuard = 0.02;      % clamp for scaled guard (Hz)
%          .maxGuard = 0.50;      % clamp for scaled guard (Hz)
%
%          .method = 'welch';     % 'welch' | 'fft' | 'pmtm'
%          .segLen = [];          % Welch segment length (samples). Default: 14-30 secs (calculated insde the function)
%          .overlapFrac = 0.5;    % Welch overlap fraction
%          .nfft   = [];          % NFFT (auto if empty)
%          .plot   = false;       % quick sanity plot
%
% Output (struct)
%   .f0, .T
%   .hz_list, .mid_list
%   .band_power_harm, .band_power_middle
%   .mean_harm_power, .mean_middle_power
%   .ratio_linear, .ratio_dB
%   .f, .Sxx       % frequency vector and PSD (power/Hz) used for integration
%   .bw, .guard    % actual values used (Hz)

% -------------------- defaults --------------------
if nargin < 3, opts = struct; end
defaults = struct( ...
    'T', [], ...
    'fmin',10,'fmax',50, ...
    'bw',[],'guard',[], ...
    'bwFrac',0.15,'guardFrac',0.08, ...
    'minBW',0.05,'maxBW',0.40, ...
    'minGuard',0.02,'maxGuard',0.50, ...
    'method','welch', ...
    'segLen',[],'overlapFrac',0.5,'nfft',[],'plot',false);

opts = filldefaults(opts, defaults);

x = x(:);

% -------------------- fundamental --------------------
if isempty(opts.T)
    T  = 1;      % backward-compatible
else
    T  = opts.T;
end
if ~isscalar(T) || ~isfinite(T) || T <= 0
    error('opts.T must be a positive scalar (seconds) if provided.');
end
f0 = 1/T; % fundamental harmonic

% -------------------- choose bw/guard (scale with f0 unless user provided absolute Hz) --------------------
if isempty(opts.bw)
    bw = opts.bwFrac * f0;
    bw = min(opts.maxBW, max(opts.minBW, bw));
else
    bw = opts.bw;
end

if isempty(opts.guard)
    guard = opts.guardFrac * f0;
    guard = min(opts.maxGuard, max(opts.minGuard, guard));
else
    guard = opts.guard;
end

% Guardrails so the "middle" region doesn't disappear between harmonics
guard = min(guard, 0.49*f0);      % must be < f0/2
bw    = min(bw,    0.49*f0);      % avoid overlapping adjacent harmonic centers too much (still allows wide peaks)

% -------------------- spectrum --------------------
switch lower(opts.method)
    case 'welch'
        if isempty(opts.segLen)
            segLenSec = max(14, 4*T);     % at least 14s, slightly longer for slower repetition
            segLenSec = min(segLenSec, 30);
            opts.segLen = round(segLenSec * fs);
        end
        win   = hann(opts.segLen);
        nover = round(opts.overlapFrac*opts.segLen);
        if isempty(opts.nfft)
            nfft = max(2^nextpow2(opts.segLen), opts.segLen);
        else
             nfft = max(opts.nfft, opts.segLen);
        end
        % pwelch
        [Sxx, f] = pwelch(x, win, nover, nfft, fs, 'onesided'); % power/Hz

    case 'pmtm'
        if isempty(opts.nfft), opts.nfft = []; end
        [Sxx, f] = pmtm(x, [], opts.nfft, fs, 'onesided'); % power/Hz

    case 'fft'
        N  = length(x);
        xw = x .* hann(N);
        if isempty(opts.nfft), nfft = 2^nextpow2(N); else, nfft = opts.nfft; end
        X  = fft(xw, nfft);

        M  = floor(nfft/2)+1;
        X1 = X(1:M);
        f  = (0:M-1)' * (fs/nfft);

        U  = sum(hann(N).^2)/N;        % window power normalization
        P2 = (abs(X1).^2)/(fs*N*U);    % power/Hz

        if mod(nfft,2)==0 % even
            Sxx = P2;
            Sxx(2:end-1) = 2*Sxx(2:end-1);
        else % odd
            Sxx = P2;
            Sxx(2:end)   = 2*Sxx(2:end);
        end

    otherwise
        error('Unknown method: %s', opts.method);
end

% -------------------- restrict analysis band --------------------
fmin = opts.fmin; fmax = opts.fmax;
inBand = (f >= max(0, fmin-1.5*bw)) & (f <= (fmax+1.5*bw));
f   = f(inBand);
Sxx = Sxx(inBand);

% -------------------- build harmonic and middle bands --------------------
k_list = ceil(fmin/f0):floor(fmax/f0);
hz_list = k_list * f0;

if numel(hz_list) < 2
    warning('Not enough harmonics in [fmin,fmax] for f0=%.4f Hz. Returning NaNs.', f0);
    out = struct();
    out.T = T; out.f0 = f0;
    out.hz_list = hz_list;
    out.mid_list = [];
    out.band_power_harm = nan(size(hz_list));
    out.band_power_middle = nan(0,1);
    out.mean_harm_power = NaN;
    out.mean_middle_power = NaN;
    out.ratio_linear = NaN;
    out.ratio_dB = NaN;
    out.f = f;
    out.Sxx = Sxx;
    out.bw = bw;
    out.guard = guard;
    return
end

mid_list = (hz_list(1:end-1) + hz_list(2:end))/2; % = hz_list(1:end-1) + f0/2

% Integrator
bandpow = @(fL,fH) integrate_band(f, Sxx, fL, fH);

% harmonic bands (signal)
band_power_harm = nan(size(hz_list));
for i = 1:numel(hz_list)
    h = hz_list(i);
    band_power_harm(i) = bandpow(h - bw, h + bw);
end

% middle bands (baseline), clipped away from harmonics by guard
band_power_middle = nan(size(mid_list));
for i = 1:numel(mid_list)
    m    = mid_list(i);

    % Candidate middle window: [m-bw, m+bw]
    % Then clip it to stay at least "guard" away from the adjacent harmonics:
    low  = max(m - bw, hz_list(i)   + guard);
    high = min(m + bw, hz_list(i+1) - guard);

    if high > low
        band_power_middle(i) = bandpow(low, high);
    end
end
band_power_middle = band_power_middle(~isnan(band_power_middle));

mean_harm = mean(band_power_harm,   'omitnan');
mean_mid  = mean(band_power_middle, 'omitnan');

out = struct();
out.T                 = T;
out.f0                = f0;
out.hz_list           = hz_list;
out.mid_list          = mid_list;
out.band_power_harm   = band_power_harm;
out.band_power_middle = band_power_middle;
out.mean_harm_power   = mean_harm;
out.mean_middle_power = mean_mid;
out.ratio_linear      = mean_harm / mean_mid;
out.ratio_dB          = 10*log10(out.ratio_linear);
out.f                 = f;
out.Sxx               = Sxx;
out.bw                = bw;
out.guard             = guard;

% -------------------- optional sanity plot --------------------
if opts.plot
    figure('Name','Harmonics vs Middles');
    loglog(f, Sxx); grid on; xlabel('Frequency (Hz)'); ylabel('PSD (power/Hz)');
    hold on
    yl = ylim;
    for h = hz_list
        plot([h h], yl, '--');
    end
    for m = mid_list
        plot([m m], yl, ':');
    end
    title(sprintf('f0=%.3f Hz (T=%.3fs) | Ratio = %.2f dB (harmonics / middles)', f0, T, out.ratio_dB));
end
end

% ---------- helpers ----------
function P = integrate_band(f,Sxx,fL,fH)
    idx = (f >= fL) & (f <= fH);
    if nnz(idx) < 2
        P = NaN;
    else
        P = trapz(f(idx), Sxx(idx)); % integrate PSD -> power in band
    end
end

function o = filldefaults(o, d)
    f = fieldnames(d);
    for k=1:numel(f)
        if ~isfield(o,f{k}) || isempty(o.(f{k}))
            o.(f{k}) = d.(f{k});
        end
    end
end
