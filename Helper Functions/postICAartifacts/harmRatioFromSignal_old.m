function out = harmRatioFromSignal_old(x, fs, opts)

% harmonic_ratio_from_signal  Quantify 1 Hz harmonic energy vs inter-harmonic energy.
%
% Usage:
%   out = harmRatioFromSignal(x, fs)
%   out = harmRatioFromSignal(x, fs, opts)
%
% Inputs
%   x    : time series (vector)
%   fs   : sampling rate (Hz)
%   opts : (struct) optional overrides:
%          .fmin   = 10;          % analysis band low (Hz)
%          .fmax   = 50;          % analysis band high (Hz)
%          .bw     = 0.30;        % half-bandwidth around centers (Hz)
%          .guard  = 0.15;        % exclude this half-width near harmonics for "middle" bands
%          .method = 'welch';     % 'welch' | 'fft' | 'pmtm'
%          .segLen = round(4*fs); % Welch segment length (samples)
%          .overlapFrac = 0.5;    % Welch overlap fraction
%          .nfft   = [];          % NFFT (auto if empty)
%          .detrend = 'constant'; % 'constant'|'none' for welch/fft
%          .plot   = false;       % quick sanity plot
%
% Output (struct)
%   .hz_list, .mid_list
%   .band_power_harm, .band_power_middle
%   .mean_harm_power, .mean_middle_power
%   .ratio_linear, .ratio_dB
%   .f, .Sxx       % frequency vector and PSD (power/Hz) used for integration


% -------------------- defaults --------------------
if nargin < 3, opts = struct; end
defaults = struct('fmin',10,'fmax',50,'bw',0.30,'guard',0.15,'method','welch', ...
                  'segLen',[],'overlapFrac',0.5,'nfft',[],'detrend','constant','plot',false);
opts = filldefaults(opts, defaults);

x = x(:);
% -------------------- spectrum --------------------
switch lower(opts.method)
    case 'welch'
        if isempty(opts.segLen), opts.segLen = max(round(4*fs), 256); end
        win   = hann(opts.segLen);
        nover = round(opts.overlapFrac*opts.segLen);
        if isempty(opts.nfft)
            nfft = max(2^nextpow2(opts.segLen), opts.segLen);
        else
            nfft = opts.nfft;
        end
        [Sxx, f] = pwelch(x, win, nover, nfft, fs, 'onesided');
        % Sxx is power/Hz (PSD) already

    case 'pmtm'
        % Multitaper (robust to leakage). Uses default time-halfbandwidth.
        if isempty(opts.nfft), opts.nfft = []; end
        [Sxx, f] = pmtm(x,[],opts.nfft,fs,'onesided'); % power/Hz

    case 'fft'
        % Simple periodogram with a Hann window for leakage control.
        % Convert |X|^2 to power/Hz to keep integration meaningful.
        N  = length(x);
        xw = x .* hann(N);
        if isempty(opts.nfft), nfft = 2^nextpow2(N); else, nfft = opts.nfft; end
        X  = fft(xw, nfft);
        % one-sided
        M  = floor(nfft/2)+1;
        X1 = X(1:M);
        f  = (0:M-1)' * (fs/nfft);
        % Window power correction and scaling to power/Hz:
        U  = sum(hann(N).^2)/N;        % window power normalization
        P2 = (abs(X1).^2)/(fs*N*U);    % power/Hz for one-sided except DC/Nyquist factor
        % Double non-DC/non-Nyquist bins to conserve total power:
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
fmin = opts.fmin; fmax = opts.fmax; bw = opts.bw; guard = opts.guard;
inBand = (f >= max(0,fmin-1.5*bw)) & (f <= (fmax+1.5*bw));
f   = f(inBand);
Sxx = Sxx(inBand);

% -------------------- build harmonic and middle bands --------------------
hz_list  = ceil(fmin):floor(fmax);        % integer-Hz harmonics
mid_list = (hz_list(1:end-1)+hz_list(2:end))/2;  % e.g., h+0.5

% Integrators
bandpow = @(fL,fH) integrate_band(f,Sxx,fL,fH);

% harmonic bands
band_power_harm = nan(size(hz_list));
for i=1:numel(hz_list)
    h = hz_list(i);
    band_power_harm(i) = bandpow(h-bw, h+bw);
end

% middle bands (exclude guard near each harmonic)
band_power_middle = nan(size(mid_list));
for i=1:numel(mid_list)
    m    = mid_list(i);
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
out.hz_list            = hz_list;
out.mid_list           = mid_list;
out.band_power_harm    = band_power_harm;
out.band_power_middle  = band_power_middle;
out.mean_harm_power    = mean_harm;
out.mean_middle_power  = mean_mid;
out.ratio_linear       = mean_harm / mean_mid;
out.ratio_dB           = 10*log10(out.ratio_linear);
out.f                  = f;
out.Sxx                = Sxx;

% -------------------- optional sanity plot --------------------
if opts.plot
    figure('Name','Harmonics vs Middles'); 
    loglog(f, Sxx); grid on; xlabel('Frequency (Hz)'); ylabel('PSD (power/Hz)');
    hold on
    yl = ylim;
    for h = hz_list
        plot([h h], yl, '--'); %#ok<*PLOTNC>
    end
    for m = mid_list
        plot([m m], yl, ':');
    end
    title(sprintf('Ratio = %.2f dB (harmonics / middles)', out.ratio_dB));
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
