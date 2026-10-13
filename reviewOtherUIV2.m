function otherSrcFinal = reviewOtherUIV2(W, xICA, fs, chanlocs, otherInit, candidates, figName, saveFigs, saveFolder, icChannels)
% REVIEWOTHERUIV2  reviewOtherUI, unchanged, with the optional IC power-per-channel window.
%
%   OTHERSRCFINAL = REVIEWOTHERUIV2(W, XICA, FS, CHANLOCS, OTHERINIT, CANDIDATES,
%                                   FIGNAME, SAVEFIGS, SAVEFOLDER, ICCHANNELS)
%
%   Calls the released reviewOtherUI with the same arguments -- same panel, same
%   behaviour, same output. When ICCHANNELS is true (default) it first opens
%   showICChannelsV2 on the same components (OTHERINIT then CANDIDATES), the second,
%   non-blocking window the V2 EKG / DIN reviews also open, and closes it when the
%   review returns.
%
%   INPUTS   as reviewOtherUI, plus
%     icChannels  true (default) open the IC power-per-channel window beside the review
%
%   See also reviewOtherUI, showICChannelsV2, reviewEkgArtifactUIV2, reviewDINArtifactUIV2

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

if nargin < 7, figName = []; end
if nargin < 8, saveFigs = []; end
if nargin < 9, saveFolder = []; end
if nargin < 10 || isempty(icChannels), icChannels = true; end

comps = unique([otherInit(:); candidates(:)].', 'stable');
if icChannels && ~isempty(comps)
    [~, icFig] = showICChannelsV2(W, xICA, fs, chanlocs, comps, struct('label', 'Other'));
    closeIC = onCleanup(@() closeIfOpen(icFig)); %#ok<NASGU>   % also on an error or Ctrl+C
end
otherSrcFinal = reviewOtherUI(W, xICA, fs, chanlocs, otherInit, candidates, figName, saveFigs, saveFolder);
end

function closeIfOpen(f)
if ~isempty(f) && ishghandle(f), delete(f); end
end
