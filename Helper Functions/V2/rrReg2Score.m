function [score, reg, edges] = rrReg2Score(xICA, fs, opts)
% RRREG2SCORE  Heartbeat-train regularity of every ICA component, per chunk of time.
%
%   [SCORE, REG, EDGES] = RRREG2SCORE(XICA, FS, OPTS)
%
%   EKG in an ICA component is a train of beats at a steady rate. rrReg2 finds that train
%   directly in the time domain, chunk by chunk, so a component that carries the heart for
%   only part of the recording still scores in the chunks where it does:
%
%     1) A = |x - movmean(x, maWinSec)|   (removes slow activity, keeps sharp deflections)
%     2) per chunk: peaks above nSigma robust sd (median(A)/0.6745), >= refracSec apart
%     3) a chunk scores only if it has >= minPeaks peaks at a median rate inside bpmRange
%     4) chunk score = 1 - CV(inter-peak intervals), floored at 0
%     5) component score = the agg (default 75th percentile) over its chunks
%
%   INPUTS
%     xICA   [nComp x nSamp] ICA component time series
%     fs     sampling rate (Hz)
%     opts   (optional struct, defaults shown)
%       .chunkSec   30     chunk length (s); chunks are non-overlapping, runts < 10 s dropped
%       .maWinSec   0.3    moving-average window (s)
%       .nSigma     2.5    peak threshold in robust sd
%       .refracSec  0.3    minimum peak spacing (s)
%       .minPeaks   5      peaks a chunk needs to be scored
%       .bpmRange   [40 120]  accepted median rate
%       .agg        'p75'  'p75' | 'max' | 'median' over chunks
%
%   OUTPUTS
%     score  [nComp x 1]  higher = more heart-like; NaN if no chunk passed the gates
%     reg    [nComp x nChunk] per-chunk regularity (NaN = chunk not scored)
%     edges  [nChunk x 2] first/last sample of every chunk
%
%   NOTES
%     • The rate gate is wide, so a 1 Hz (60 bpm) trigger train passes it; a DIN
%       component can therefore score high here. The DIN detector is what sorts that out.
%     • Peaks are taken on a moving-average residual, not a band-passed signal: on the
%       SENSI Session 1 cohort that change (with nSigma 2.5) is what lifted recall at
%       the top 10 from 95% to 99% (114 labelled EKG components).
%
%   See also autoDetectEkgV2, autoDetectEkg

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
opts = fillDefault(opts, 'chunkSec',30, 'maWinSec',0.3, 'nSigma',2.5, 'refracSec',0.3, ...
    'minPeaks',5, 'bpmRange',[40 120], 'agg','p75');

x = xICA - mean(xICA, 2);
A = abs(x - movmean(x, round(opts.maWinSec*fs), 2));
[nComp, nSamp] = size(A);
chunkN = round(opts.chunkSec*fs);
st = 1:chunkN:max(1, nSamp-chunkN+1);
ed = min(st+chunkN-1, nSamp);
keep = (ed-st+1) >= round(10*fs); st = st(keep); ed = ed(keep);
edges = [st(:) ed(:)];
refrac = max(1, round(opts.refracSec*fs));
reg = nan(nComp, numel(st));
for c = 1:nComp
    for q = 1:numel(st)
        a = A(c, st(q):ed(q));
        sd = median(a)/0.6745;
        if sd <= 0, continue; end
        [~, pk] = findpeaks(a, 'MinPeakHeight', opts.nSigma*sd, 'MinPeakDistance', refrac);
        if numel(pk) < opts.minPeaks, continue; end
        rr = diff(pk)/fs;
        bpm = 60/median(rr);
        if bpm < opts.bpmRange(1) || bpm > opts.bpmRange(2), continue; end
        reg(c,q) = max(0, 1 - std(rr)/max(mean(rr),eps));
    end
end
score = aggChunks(reg, opts.agg);
end

function s = aggChunks(v, agg)
switch lower(agg)
    case 'p75',    s = prctile(v, 75, 2);
    case 'max',    s = max(v, [], 2);
    case 'median', s = median(v, 2, 'omitnan');
    otherwise, error('rrReg2Score: unknown agg ''%s''.', agg);
end
end

function s = fillDefault(s, varargin)
for i = 1:2:numel(varargin)
    if ~isfield(s, varargin{i}) || isempty(s.(varargin{i})), s.(varargin{i}) = varargin{i+1}; end
end
end
