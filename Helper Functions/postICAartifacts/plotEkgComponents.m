function [figDef, figSus] = plotEkgComponents(W, xICA, fs, chanlocs, ekgSrc, susEkg, ekgStartSec, ekgNSec)
% plotEkgComponents  Plot definite and suspect EKG ICA components
%
% Syntax:
%   [figDef,figSus] = plotEkgComponents(W, xICA, fs, chanlocs, ekgSrc, susEkg)
%   [figDef,figSus] = plotEkgComponents(W, xICA, fs, chanlocs, ekgSrc, susEkg, ekgStartSec, ekgNSec)
%
% Inputs:
%   W            - ICA unmixing matrix (components × channels)
%   xICA         - ICA component time series (components × samples)
%   fs           - Sampling frequency (Hz)
%   chanlocs     - Channel location structs for topoplotStandalone
%   ekgSrc       - Vector of definite EKG component indices
%   susEkg       - Vector of suspect EKG component indices
%   ekgStartSec  - (optional) start time for time-domain snippet, default = 5*60
%   ekgNSec      - (optional) duration in seconds for snippet, default = 10
%
% Outputs:
%   figDef       - Handle to the figure of definite EKG components
%   figSus       - Handle to the figure of suspect EKG components

if nargin < 7 || isempty(ekgStartSec), ekgStartSec = 5*60; end
if nargin < 8 || isempty(ekgNSec),     ekgNSec     = 10;   end

Wi = inv(W);
[nComp, nSamp] = size(xICA);
maxFreq = 50;

startIdx = round(ekgStartSec * fs) + 1;
endIdx   = min(nSamp, startIdx + round(ekgNSec * fs) - 1);
idxPlot  = startIdx:endIdx;

% Initialize output handles
figDef = [];
figSus = [];

    function hFig = plotRows(sources, figName)
        if isempty(sources)
            hFig = [];
            return;
        end
        nRow = numel(sources);
        hFig = figure('Name', figName, 'Color', 'w', 'Units', 'pixels');
        hFig.Position(3) = 900;
        hFig.Position(4) = nRow*150 + 50;
        T = tiledlayout(nRow, 7, 'TileSpacing','compact','Padding','compact');
        for r = 1:nRow
            src = sources(r);
            % Topography
            nexttile(T, (r-1)*7 + 1);
            topoplotStandalone(Wi(:,src), chanlocs);
            title(sprintf('IC %d', src), 'Interpreter','none');
            % Time Domain snippet
            nexttile(T, (r-1)*7 + 2, [1 3]);
            plot(idxPlot/fs, xICA(src, idxPlot), 'LineWidth',1); grid on;
            xlabel('Time (s)'); ylabel('\muV');
            title('Time Domain','Interpreter','none');
            % Frequency Domain
            nexttile(T, (r-1)*7 + 5, [1 3]);
            data = xICA(src, :);
            fAx  = computeFFTFrequencyAxis(numel(data), fs);
            plot(fAx, abs(fft(data)), 'LineWidth',1); grid on;
            xlim([0 maxFreq]); xlabel('Hz');
            title('Frequency Domain','Interpreter','none');
        end
        sgtitle(T, figName, 'Interpreter','none');
    end

warning('off','all');
figDef = plotRows(ekgSrc, sprintf('Definite EKG Components [%s]', mat2str(ekgSrc)));
figSus = plotRows(susEkg, sprintf('Suspect EKG Components [%s]', mat2str(susEkg)));
warning('on','all');

end