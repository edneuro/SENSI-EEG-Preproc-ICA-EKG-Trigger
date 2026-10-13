function [dinSrc, dinCand, info] = autoDetectDINV2(xICA, fs, nSources, opts)
% AUTODETECTDINV2  Trigger (DIN) detection with the comb score, at any repetition period.
%
%   [DINSRC, DINCAND, INFO] = AUTODETECTDINV2(XICA, FS, NSOURCES, OPTS)
%
%   V2 of autoDetectDIN. Every component is scored with dinCombScore's combC -- a weighted
%   SUM over its magnitude spectrum, positive on each harmonic k/T and negative between
%   them, with outlying bins dropped from the peaks and from the valleys separately --
%   instead of autoDetectDIN's harmonic-to-midpoint RATIO.
%
%   PANEL: every component whose robust z (within nSources) exceeds minZ, highest score
%   first, at most nPanel of them -- so a block-group with one clear DIN shows few rows and
%   one with several shows up to the cap. A member starts flagged (red) when its robust z
%   exceeds redZ. z = (score - median) / (1.4826 * MAD) over the scanned components, so the
%   rule is how far a component stands out from the rest of its own recording.
%
%   Why not the ratio: on the SENSI Session 1 cohort a narrowband line at ~23.5 Hz (half
%   way between two DIN harmonics) sits in the ratio's denominator. The component with the
%   most DIN often also carries the most 23.5 Hz, so its ratio falls below 1 and its
%   z-score went as low as -10: the true DIN component ranked last of 124. A sum can only
%   lose a bounded amount to one band. Measured on 128 labelled DIN components: 98.4% in
%   the top 3 of their block-group for combC (all of them in the top 5), 92.2% for the
%   released ratio.
%
%   INPUTS
%     xICA      [nComp x nSamp] ICA component time series
%     fs        sampling rate (Hz)
%     nSources  components to scan (default 1:nComp -- autoDetectDIN defaults to 1:30)
%     opts      (optional struct, defaults shown)
%       .T         1       trigger repetition period (s)
%       .fmin      10      analysis band (Hz)
%       .fmax      50
%       .segSec    []      Welch segment (s); Inf = whole-record periodogram. [] = by T:
%                          20 s for T <= 3, 40 s for T <= 10, the whole record beyond
%       .score     'combC' 'combC' | 'comb' | 'combR' (see dinCombScore)
%       .nPanel    5       the cap: at most this many rows on the panel
%       .minZ      2       a component needs robust z > minZ to be on the panel
%       .redZ      15      a panel member starts red when its robust z > redZ
%       .stripChunkSec 60  chunk length (s) used to find each component's strongest chunk
%                          (the review UIs' 'strongest' window start)
%
%   OUTPUTS
%     dinSrc   panel members pre-flagged as DIN (start red in the review UI)
%     dinCand  the rest of the panel (start green)
%     info     struct: .score / .z / .rank [nComp x 1], .panel, .why, .strip (.t, .s:
%              the score per chunk, from which .bestStartSec comes),
%              .bestStartSec [nComp x 1], .all (all three comb variants), .opts
%
%   DEFAULTS, measured on SENSI Session 1 (128 labelled DIN components at T = 1, and a
%   semi-synthetic test: a real DIN waveform injected at T = 1-60 s, 0/1/5 ms jitter, into
%   36 non-DIN components, through the recordings' real epoch structure).
%     T = 1, combC on a 20 s Welch: recall @1/2/3/4/5 = 82.8/96.1/98.4/99.2/100% -- every
%     labelled DIN inside the top 5 of its block-group. The plain comb: 82.8/96.1/99.2/
%     99.2/99.2%, 100% only at N = 47; the released ratio 77.3/89.1/92.2/94.5/95.3%.
%     The Welch segment must hold several periods. At the real DIN's median strength
%     (top 3, no jitter) the 20 s spectrum holds to T = 5 s (100%, 53% at 7 s), 40 s to
%     T = 10 s (92%, 22% at 15 s); the whole record holds to T = 10 s and still finds
%     92% at 15 s and 81% at 20 s. At T = 1 the 20 s / 40 s spectra need half the
%     amplitude the whole record does. On the injected trains combC = comb (no
%     interferer to clip); the clip only matters on real data.
%     Panel (z > 2, cap 5) and start red (z > 15), 20 s spectrum: all 128 DIN on the
%     panel, 4.6 rows per block-group on average; 152 rows start red, 76% of them DIN,
%     90% of DIN start red. 40 s: 126/128 on the panel (2 DIN have z <= 2), whole record
%     125/128. Labelled DIN have a median robust z of 67 (20 s); non-DIN components that
%     reach a 3-row panel, 4.7; all other components, 0.0 (90th percentile 1.3).
%
%   See also dinCombScore, autoDetectDIN, reviewDINArtifactUIV2, autoDetectEkgV2

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

[nComp, ~] = size(xICA);
if nargin < 3 || isempty(nSources), nSources = 1:nComp; end
if nargin < 4 || isempty(opts), opts = struct(); end
opts = fillDefault(opts, 'T',1, 'fmin',10, 'fmax',50, 'segSec',[], 'score','combC', ...
    'nPanel',5, 'minZ',2, 'redZ',15, 'stripChunkSec',60);
if isempty(opts.segSec)
    if opts.T <= 3, opts.segSec = 20; elseif opts.T <= 10, opts.segSec = 40; else, opts.segSec = Inf; end
end
nSources = unique(nSources(nSources >= 1 & nSources <= nComp));
inS = false(nComp,1); inS(nSources) = true;

co = struct('T',opts.T,'fmin',opts.fmin,'fmax',opts.fmax,'segSec',opts.segSec, ...
    'chunkSec',opts.stripChunkSec, 'chunkScore',opts.score);
A = dinCombScore(xICA, fs, co);
s = A.(opts.score);
u = s(inS & ~isnan(s)); md = median(u); mad = 1.4826*median(abs(u - md));
z = (s - md) ./ max(mad, eps);
rank = inf(nComp,1);
[~, o] = sort(s(nSources), 'descend', 'MissingPlacement','last'); o = nSources(o);
rank(o) = 1:numel(o);
panel = o(z(o) > opts.minZ);                      % stands out from its own recording ...
panel = panel(1:min(opts.nPanel, numel(panel)));  % ... up to the cap
why = arrayfun(@(c) sprintf('%s #%d (%.2f, z %.1f)', opts.score, rank(c), s(c), z(c)), ...
    panel, 'UniformOutput', false);
red = z(panel) > opts.redZ;
dinSrc = panel(red); dinCand = panel(~red);

bestStartSec = zeros(nComp,1);
for c = 1:nComp
    [v, q] = max(A.chunk(c,:));
    if ~isempty(v) && isfinite(v), bestStartSec(c) = (A.edges(q,1)-1)/fs; end
end
info = struct('score',s,'z',z,'rank',rank,'panel',panel,'why',{why}, ...
    'strip',struct('t',mean(A.edges,2)'/fs,'s',A.chunk), 'bestStartSec',bestStartSec, ...
    'all',struct('comb',A.comb,'combR',A.combR,'combC',A.combC),'opts',opts);
end

function s = fillDefault(s, varargin)
for i = 1:2:numel(varargin)
    if ~isfield(s, varargin{i}) || isempty(s.(varargin{i})), s.(varargin{i}) = varargin{i+1}; end
end
end
