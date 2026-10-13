function [ekgSrc, ekgCand, info] = autoDetectEkgV2(xICA, fs, nSources, opts)
% AUTODETECTEKGV2  Two-arm EKG detection: released harmonic SNR + beat-train regularity.
%
%   [EKGSRC, EKGCAND, INFO] = AUTODETECTEKGV2(XICA, FS, NSOURCES, OPTS)
%
%   V2 of autoDetectEkg (which it calls, unmodified). Two independent arms rank every
%   component, and a fixed-size review panel is built from both:
%
%     released arm  autoDetectEkg's summed harmonic SNR -- over the whole record, or per
%                   chunk of relChunkSec s and aggregated (relAgg) over chunks
%     rrReg2 arm    rrReg2Score: regularity of the beat train in the time domain, per
%                   chunk, 75th percentile over chunks
%
%   Panel (nPanel slots, default 6): the top topK (3) of EACH arm by rank, repeats merged,
%   then filled to nPanel with the components that stand out most in their own arm's
%   distribution -- max(robust z of the released arm, robust z of rrReg2). On the SENSI
%   Session 1 cohort (114 labelled EKG components, found in 92 of its 168 block-groups)
%   this panel caught 113, against 96 for the released score alone given the same 6 slots.
%
%   DEFAULTS, measured on SENSI Session 1 (168 block-groups; 114 labelled EKG components,
%   113 the analyst must catch). Recall = labels inside the top N of their own block-group.
%     released arm per 30 s chunk (Welch 14 s, p75): recall @1/@3/@5 61/86/92%, against
%       46/77/82% for the whole record. In the 6-slot panel both reach 113/114; with
%       4 slots the chunked arm keeps 113, the whole-record arm 112.
%     rrReg2 per 30 s chunk, p75: 71/95/97% (15 s and 60 s chunks: <= 71/94/96%).
%     start red when the component is in the top 4 of BOTH arms: on the panel, 58% of red
%       rows are labelled EKG and 87% of labelled EKG start red. Both arms' z > 4
%       (redRule 'both'): 76% / 82%. The released rule (released z > 5): 27% / 85%.
%
%   WHY CHUNKS: a block-group is several recordings concatenated, but ICA fits one
%   stationary unmixing over all of it, so the heart can move between components over
%   time. A whole-record spectrum dilutes a component that carries the heart for only
%   part of the record; a per-chunk score does not.
%
%   INPUTS
%     xICA      [nComp x nSamp] ICA component time series
%     fs        sampling rate (Hz)
%     nSources  components to scan (default 1:nComp)
%     opts      (optional struct, defaults shown)
%       .relChunkSec  30     released arm chunk (s); Inf = whole record (as released)
%       .relWindowSec 14     released arm Welch window (s); .window works too
%       .relAgg       'p75'  released arm aggregate over chunks: 'p75' | 'max'
%       autoDetectEkg's own options can be put straight in OPTS, as in V1:
%         .fmin .fmax .harmonics .peakHalfHz .nbHalfHz .peakAgg .do_detrend
%         .detrend_fmin .detrend_fmax  (and .window, = .relWindowSec)
%       They go to the released arm. (.ekgOpts, a struct of them, still works too.)
%       .rr           struct()  passed to rrReg2Score (.chunkSec 30, ...)
%       .topK         3      per-arm top-k guaranteed a slot
%       .nPanel       6      panel size
%       .redRule      'topK' a panel member starts flagged (red) when: 'topK' -- it ranks
%                            in the top redTopK of BOTH arms | 'both' -- both arms'
%                            robust z exceed redZ | 'released' | 'rrReg2' | 'either' (z)
%       .redTopK      4
%       .redZ         4      (z rules only)
%       .stripChunkSec 30    chunk length (s) of the released arm's per-chunk scores when
%                            that arm is whole-record (relChunkSec = Inf); used only to
%                            find each component's strongest chunk
%
%   OUTPUTS
%     ekgSrc   panel members pre-flagged as EKG (start red in the review UI)
%     ekgCand  the rest of the panel (start green)
%     info     struct: .rel / .relZ / .rr / .rrZ / .best  [nComp x 1] scores and robust z,
%              .rankRel / .rankRr  ranks within nSources, .panel  panel order,
%              .why  {1 x nPanel} which arm(s) brought each member in,
%              .strip  each arm's score per chunk (.tRel, .rel, .tRr, .rr), from which
%                      .bestStartSec comes,
%              .bestStartSec [nComp x 1]  start of the chunk where it is most heart-like
%
%   See also autoDetectEkg, rrReg2Score, reviewEkgArtifactUIV2, autoDetectDINV2

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

[nComp, nSamp] = size(xICA);
if nargin < 3 || isempty(nSources), nSources = 1:nComp; end
if nargin < 4 || isempty(opts), opts = struct(); end
opts = fillDefault(opts, 'relChunkSec',30, 'relWindowSec',14, 'relAgg','p75', ...
    'ekgOpts',struct(), 'rr',struct(), 'topK',3, 'nPanel',6, 'redRule','topK', 'redTopK',4, 'redZ',4, ...
    'stripChunkSec',30);
nSources = unique(nSources(nSources >= 1 & nSources <= nComp));
nSources = nSources(:)';
inS = false(nComp,1); inS(nSources) = true;

% ---- released arm
% The released arm's options: autoDetectEkg's names given straight in OPTS (as V1 takes
% them), or nested under opts.ekgOpts; a top-level value wins.
eo = opts.ekgOpts;
for nm = {'fmin','fmax','harmonics','peakHalfHz','nbHalfHz','peakAgg','do_detrend', ...
          'detrend_fmin','detrend_fmax'}
    if isfield(opts, nm{1}) && ~isempty(opts.(nm{1})), eo.(nm{1}) = opts.(nm{1}); end
end
if isfield(opts, 'window') && ~isempty(opts.window), opts.relWindowSec = opts.window; end
eo.window = opts.relWindowSec;
if isinf(opts.relChunkSec)
    [~,~,m] = autoDetectEkg(xICA, fs, 1:nComp, eo);
    rel = m.harmSNRsum(:);
    [s0, e0] = chunks(nSamp, round(opts.stripChunkSec*fs), fs);   % only for the 'strongest' start
    relChunk = perChunkHarm(xICA, fs, eo, s0, e0);
else
    [s0, e0] = chunks(nSamp, round(opts.relChunkSec*fs), fs);
    relChunk = perChunkHarm(xICA, fs, eo, s0, e0);
    rel = aggChunks(relChunk, opts.relAgg);
end
relT = (s0 + e0)/2/fs;

% ---- rrReg2 arm
[rr, rrChunk, rrEdges] = rrReg2Score(xICA, fs, opts.rr);

% ---- robust z within the scanned sources, ranks, panel
relZ = robustZ(rel, inS); rrZ = robustZ(rr, inS);
best = max(relZ, rrZ);
rankRel = rankIn(rel, nSources, nComp); rankRr = rankIn(rr, nSources, nComp);
[~, o1] = sort(rel(nSources), 'descend', 'MissingPlacement','last'); o1 = nSources(o1);
[~, o2] = sort(rr(nSources),  'descend', 'MissingPlacement','last'); o2 = nSources(o2);
k = min(opts.topK, numel(nSources));
panel = unique([o1(1:k), o2(1:k)], 'stable');
[~, ob] = sort(best(nSources), 'descend'); ob = nSources(ob);
panel = [panel, setdiff(ob, panel, 'stable')];
panel = panel(1:min(opts.nPanel, numel(panel)));
why = cell(1, numel(panel));
for i = 1:numel(panel)
    c = panel(i); w = {};
    if rankRel(c) <= opts.topK, w{end+1} = sprintf('released #%d', rankRel(c)); end %#ok<AGROW>
    if rankRr(c)  <= opts.topK, w{end+1} = sprintf('rrReg2 #%d', rankRr(c)); end %#ok<AGROW>
    if isempty(w), w = {sprintf('stands out (best z %.1f)', best(c))}; end
    why{i} = strjoin(w, ' + ');
end
switch lower(opts.redRule)
    case 'topk',     red = rankRel(panel) <= opts.redTopK & rankRr(panel) <= opts.redTopK;
    case 'both',     red = min(relZ(panel), rrZ(panel)) > opts.redZ;
    case 'either',   red = max(relZ(panel), rrZ(panel)) > opts.redZ;
    case 'released', red = relZ(panel) > opts.redZ;
    case 'rrreg2',   red = rrZ(panel) > opts.redZ;
    otherwise, error('autoDetectEkgV2: unknown redRule ''%s''.', opts.redRule);
end
ekgSrc = panel(red); ekgCand = panel(~red);

% ---- where to look: the chunk in which each component is most heart-like
bestStartSec = zeros(nComp,1);
for c = 1:nComp
    [v, q] = max(rrChunk(c,:));
    if ~isempty(v) && isfinite(v)
        bestStartSec(c) = (rrEdges(q,1)-1)/fs;
    else
        [v, q] = max(relChunk(c,:));
        if ~isempty(v) && isfinite(v), bestStartSec(c) = (s0(q)-1)/fs; end
    end
end

info = struct('rel',rel,'relZ',relZ,'rr',rr,'rrZ',rrZ,'best',best,'rankRel',rankRel, ...
    'rankRr',rankRr,'panel',panel,'why',{why}, ...
    'strip',struct('tRel',relT(:)','rel',relChunk,'tRr',mean(rrEdges,2)'/fs,'rr',rrChunk), ...
    'bestStartSec',bestStartSec,'opts',opts);
end

% -------------------------------------------------------------------------
function [st, ed] = chunks(nSamp, chunkN, fs)
st = 1:chunkN:max(1, nSamp-chunkN+1);
ed = min(st+chunkN-1, nSamp);
keep = (ed-st) >= round(10*fs); st = st(keep); ed = ed(keep);
end

function H = perChunkHarm(xICA, fs, eo, st, ed)
% autoDetectEkg, unmodified, on each chunk; a chunk it cannot score stays NaN.
H = nan(size(xICA,1), numel(st));
for q = 1:numel(st)
    try
        [~,~,m] = autoDetectEkg(xICA(:, st(q):ed(q)), fs, 1:size(xICA,1), eo);
        H(:,q) = m.harmSNRsum(:);
    catch
    end
end
end

function s = aggChunks(v, agg)
switch lower(agg)
    case 'p75', s = prctile(v, 75, 2);
    case 'max', s = max(v, [], 2);
    otherwise, error('autoDetectEkgV2: unknown relAgg ''%s''.', agg);
end
end

function z = robustZ(v, inS)
% (v - median) / (1.4826 MAD) over the scanned sources; NaN -> -Inf (never stands out).
u = v(inS & ~isnan(v));
z = -inf(size(v));
if isempty(u), return; end
md = median(u); mad = 1.4826 * median(abs(u - md));
z = (v - md) ./ max(mad, eps); z(isnan(v)) = -Inf;
end

function r = rankIn(v, nSources, nComp)
r = inf(nComp,1);
[~, o] = sort(v(nSources), 'descend', 'MissingPlacement','last');
r(nSources(o)) = 1:numel(o);
end

function s = fillDefault(s, varargin)
for i = 1:2:numel(varargin)
    if ~isfield(s, varargin{i}) || isempty(s.(varargin{i})), s.(varargin{i}) = varargin{i+1}; end
end
end
