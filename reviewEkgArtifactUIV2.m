function [ekgSrcFinal, added] = reviewEkgArtifactUIV2(W, xICA, fs, chanlocs, ekgSrc, candidates, info, StartSec, NSec, figName, saveFigs, saveFolder, rowsPerPage, icChannels)
% REVIEWEKGARTIFACTUIV2  V2 interactive review of EKG (cardiac) ICA components.
%
%   [EKGSRCFINAL, ADDED] = REVIEWEKGARTIFACTUIV2(W, XICA, FS, CHANLOCS, EKGSRC, CANDIDATES,
%                              INFO, STARTSEC, NSEC, FIGNAME, SAVEFIGS, SAVEFOLDER,
%                              ROWSPERPAGE, ICCHANNELS)
%
%   The released reviewEkgArtifactUI, plus what autoDetectEkgV2 knows about each component:
%   one row per panel component, at most ROWSPERPAGE rows per page ('< Prev' / 'Next >'):
%   scalp map | NSEC s of it | |FFT|, as released -- red = reject, green = keep, click a
%   row to toggle. The window starts at STARTSEC (seconds, or 'strongest': each row at its
%   own strongest chunk); 'add IC' + '+ Add' brings any other component onto the panel (it
%   starts green); a second window shows each component's power per channel. See
%   reviewArtifactUIV2 for the details.
%
%   INPUTS   as reviewEkgArtifactUI, plus
%     info       third output of autoDetectEkgV2
%     StartSec   'strongest' (default) = each row at its own strongest chunk; or a start
%                time in seconds, the same for every row
%     rowsPerPage  4 (default) rows per page; more components = more pages
%     icChannels   true (default) open the IC power-per-channel window (showICChannelsV2)
%   OUTPUTS
%     ekgSrcFinal  components left red: confirmed EKG artifacts
%     added      components added with '+ Add' (whatever colour they ended)
%
%   See also autoDetectEkgV2, reviewArtifactUIV2, reviewEkgArtifactUI

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

if nargin < 8 || isempty(StartSec), StartSec = 'strongest'; end
if nargin < 9 || isempty(NSec), NSec = 10; end
if nargin < 10 || isempty(figName), figName = 'Review EKG Artifacts (V2)'; end
if nargin < 11 || isempty(saveFigs), saveFigs = 0; end
if nargin < 12 || isempty(saveFolder), saveFolder = pwd; end
if nargin < 13 || isempty(rowsPerPage), rowsPerPage = 4; end
if nargin < 14 || isempty(icChannels), icChannels = true; end
[ekgSrcFinal, added] = reviewArtifactUIV2('ekg', W, xICA, fs, chanlocs, ekgSrc, candidates, info, ...
    struct('figName', figName, 'saveFigs', saveFigs, 'saveFolder', saveFolder, ...
           'nSec', NSec, 'startSec', StartSec, 'rowsPerPage', rowsPerPage, 'icChannels', icChannels));
end
