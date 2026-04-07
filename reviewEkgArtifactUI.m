function ekgSrcFinal = reviewEkgArtifactUI(W, xICA, fs, chanlocs, ekgSrc, candidates, StartSec, NSec, figName, saveFigs, saveFolder)
% REVIEWEKGARTIFACTUI  Interactive review UI for EKG artifact ICA components.
%
%   EKGSRCFINAL = REVIEWEKGARTIFACTUI(W, XICA, FS, CHANLOCS, EKGSRC,
%                   CANDIDATES, STARTSEC, NSEC, FIGNAME, SAVEFIGS, SAVEFOLDER)
%
%   Opens a MATLAB figure displaying each candidate component as a row of 7
%   tiles: topoplot | time-series window | full-recording |FFT|. Components
%   in EKGSRC start flagged red (REJECT); components in CANDIDATES start
%   green (KEEP). Click anywhere in a row to toggle its status. Click "Done"
%   to finalize and return the list of rejected components.
%
%   INPUTS
%     W           [nCh x nComp] ICA unmixing matrix
%     xICA        [nComp x nSamp] ICA component time series
%     fs          Scalar sampling rate (Hz)
%     chanlocs    Electrode location struct (from chanlocsFromFile)
%     ekgSrc      Vector of component indices pre-flagged for REJECT (red)
%     candidates  Vector of component indices to review (green)
%     StartSec    Start time (sec) for the time-series plot window
%     NSec        Duration (sec) of the time-series plot window
%     figName     String used as figure title and output filename base
%     saveFigs    (Optional) 1 = save figure on Done, 0 = no save (default 0)
%     saveFolder  (Optional) Output folder for saved figure (default: pwd)
%
%   OUTPUTS
%     ekgSrcFinal  Vector of component indices confirmed as EKG artifacts
%
%   NOTES
%     • If both ekgSrc and candidates are empty, returns [] immediately.
%     • Figure is saved as a PNG at 300 DPI when saveFigs=1 and Done is clicked.
%     • The frequency plot uses the full IC time series (not just the window)
%       and displays 0–10 Hz to highlight cardiac harmonic structure.

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

% -------- Input Argument Handling for New Inputs --------
if nargin < 11 || isempty(saveFolder), saveFolder = pwd; end
if nargin < 10 || isempty(saveFigs), saveFigs = 0; end
if nargin < 9 || isempty(figName), figName = 'Review EKG Artifacts (tiled)'; end

if isempty(StartSec), StartSec = 0; end
if isempty(NSec) || NSec <= 0, NSec = max(1, min(10, floor(size(xICA,2)/fs))); end
% --- Prep data ---
Wi = inv(W);                              % mixing (columns are topographies)
allList = unique([ekgSrc(:); candidates(:)].','stable');
if isempty(allList)
    ekgSrcFinal = [];
    return;
end
isEKG = ismember(allList, ekgSrc(:).');
N = size(xICA,2);
i0 = max(1, floor(StartSec*fs) + 1);
i1 = min(N, i0 + max(1, floor(NSec*fs)) - 1);
if i1 < i0, i1 = i0; end
idxPlot = i0:i1;
t = (0:numel(idxPlot)-1) / fs;            % 0-based window axis
maxFreq = min(60, fs/2);
% --- Colors & styles ---
cEKG   = [0.75 0.10 0.10];   % red
cCAND  = [0.00 0.55 0.00];   % green
lw     = 1.0;
% --- Figure sizing (no pagination) ---
nRows = numel(allList);
ctrlH = 54;
rowH_est = 170;
figW = 1100;
figH = ctrlH + (nRows*rowH_est + 60);
fig = figure('Name', figName, 'NumberTitle','off', 'Color','w', 'Units','pixels');
fig.Position(3:4) = [figW, figH];
% Initial center + clamp to current screen
centerAndClamp(fig);
% Re-clamp if MATLAB resizes the window
set(fig, 'SizeChangedFcn', @(~,~)centerAndClamp(fig));
% --- Control bar ---
ctrl = uipanel(fig,'Units','pixels','Position',[0 figH-ctrlH figW ctrlH], ...
    'BackgroundColor',[.97 .97 .97],'BorderType','none');
uicontrol(ctrl,'Style','text','Units','normalized','BackgroundColor',[.97 .97 .97], ...
    'Position',[.02 .12 .72 .76], 'FontSize',12, 'FontWeight','bold', ...
    'HorizontalAlignment','left', ...
    'String', sprintf('%s',figName));
uicontrol(ctrl,'Style','pushbutton','String','Done','Units','normalized', ...
    'Position',[.82 .12 .16 .76],'FontWeight','bold','Callback',@(~,~)onDone);
% --- Tiled layout content (7 cols per row: 1 topo | 3 time | 3 freq) ---
contentH = figH - ctrlH;
content = uipanel(fig,'Units','pixels','Position',[0 0 figW contentH], ...
    'BackgroundColor','w','BorderType','none');
T = tiledlayout(content, nRows, 7, 'Padding','compact', 'TileSpacing','compact');
sgtitle(T, 'EKG Review (red = EKG, green = not EKG)', 'Interpreter','none');
% Keep handles per row index for recoloring on toggle
hTime = gobjects(nRows,1);
hFreq = gobjects(nRows,1);
hBadge= gobjects(nRows,1);
% Suppress specific tiledlayout position warning
warning('off','all');
for r = 1:nRows
    comp = allList(r);
    col = tern(isEKG(r), cEKG, cCAND);
    % --- (1) Topoplot ---
    ax1 = nexttile(T, (r-1)*7 + 1);
    topoplot_new(Wi(:,comp), chanlocs);
    title(ax1, sprintf('IC %d', comp), 'Interpreter','none');
    axis(ax1, 'off');
    hold(ax1,'on');
    if isEKG(r)
        hBadge(r) = text(ax1, 0.95, 0.95, '✓', 'Units','normalized', ...
            'HorizontalAlignment','right','VerticalAlignment','top', ...
            'FontWeight','bold','Color',col,'FontSize',12,'HitTest','off');
    else
        hBadge(r) = text(ax1, 0.95, 0.95, '•', 'Units','normalized', ...
            'HorizontalAlignment','right','VerticalAlignment','top', ...
            'FontWeight','bold','Color',col,'FontSize',12,'HitTest','off');
    end
    hold(ax1,'off');
    % --- (2) Time (windowed) ---
    ax2 = nexttile(T, (r-1)*7 + 2, [1 3]);
    seg = xICA(comp, idxPlot);
    hTime(r) = plot(ax2, t, seg, 'LineWidth', lw, 'Color', col);
    grid(ax2,'on');
    xlim(ax2, [0, max(eps, (numel(idxPlot)-1)/fs)]);
    xlabel(ax2,'Time (s)'); ylabel(ax2,'Amp');
    title(ax2, 'Time (window)', 'Interpreter','none');
    % --- (3) Frequency: |FFT| (full time series) ---
    ax3 = nexttile(T, (r-1)*7 + 5, [1 3]);
    seg = xICA(comp, :); % NOTE: This uses the FULL series, unlike the time plot!
    L = numel(seg);
    if L < 2
        hFreq(r) = plot(ax3, 0, abs(seg), 'o', 'Color', col, 'LineWidth', lw);
        xlim(ax3,[0, 10]);
    else
        nfft = 2^nextpow2(L);
        S = fft(double(seg), nfft);
        Mag = abs(S(1:nfft/2+1));
        F   = (0:nfft/2) * (fs/nfft);
        hFreq(r) = plot(ax3, F, Mag, 'LineWidth', lw, 'Color', col);
        xlim(ax3, [0, 10]);
    end
    grid(ax3,'on');
    xlabel(ax3,'Hz'); ylabel('|X(f)|');
    title(ax3, '|FFT| (full series)', 'Interpreter','none'); % Updated title for clarity
    % Row-wide toggle
    set(ax1, 'ButtonDownFcn', @(~,~)toggleRow(r));
    set(ax2, 'ButtonDownFcn', @(~,~)toggleRow(r));
    set(ax3, 'ButtonDownFcn', @(~,~)toggleRow(r));
    % Let clicks pass through line objects
    set(hTime(r),'HitTest','off','PickableParts','none');
    set(hFreq(r),'HitTest','off','PickableParts','none');
end
uiwait(fig);
ekgSrcFinal = allList(isEKG);
if ishghandle(fig), close(fig); end
% % Re-enable warnings
warning('on','all');
return;
% ================= nested callbacks =================
    function toggleRow(r)
        isEKG(r) = ~isEKG(r);
        if isEKG(r)
            newCol = cEKG; badgeChar = '✓';
        else
            newCol = cCAND; badgeChar = '•';
        end
        if ishghandle(hTime(r)), set(hTime(r), 'Color', newCol); end
        if ishghandle(hFreq(r)), set(hFreq(r), 'Color', newCol); end
        if ishghandle(hBadge(r))
            set(hBadge(r), 'String', badgeChar, 'Color', newCol);
        end
    end

    % *** UPDATED onDone function with save logic ***
    function onDone(~,~)
        if saveFigs && ishghandle(fig)
            saveCurrentFigurePage(fig, saveFolder, figName);
        end
        uiresume(fig);
    end
end
% ================= nested save function (renamed for clarity) =================
function saveCurrentFigurePage(figh, folder, name)
    % Saves the current (and only) figure page.
    if ~isfolder(folder), mkdir(folder); end
    try
        % Create safe filename
        baseName = strrep(name, ' ', '_');
        fileName = sprintf('%s.png', baseName); % No page number since there is only one page
        savePath = fullfile(folder, fileName);

        exportgraphics(figh, savePath, 'Resolution', 300);
        fprintf('Saved figure to: %s\n', savePath);
    catch ME
        warning('reviewEkgArtifactUI:SaveError', 'Could not save figure "%s" to path "%s". Error: %s', name, folder, ME.message);
    end
end
% -------- utilities --------
function centerAndClamp(fig)
% CENTERANDCLAMP  Center the figure on the current monitor and clamp it
%   fully on-screen. On multi-monitor setups, targets the monitor with the
%   greatest overlap with the current figure position.
    if ~ishghandle(fig), return; end
    oldUnits = get(fig,'Units');
    set(fig,'Units','pixels');
    pos = get(fig,'Position');     % [x y w h]
    mon = get(0,'MonitorPositions');  % Nx4: [x y w h] in pixels
    % Choose the monitor with max overlap (so multi-monitor setups behave)
    [~, idx] = max(arrayfun(@(i) overlapArea(pos, mon(i,:)), 1:size(mon,1)));
    m = mon(idx,:);
    % Desired center on that monitor
    desiredX = m(1) + (m(3) - pos(3))/2;
    desiredY = m(2) + (m(4) - pos(4))/2;
    % Clamp so the whole window stays visible
    newX = min(max(desiredX, m(1)), m(1) + m(3) - pos(3));
    newY = min(max(desiredY, m(2)), m(2) + m(4) - pos(4));
    set(fig,'Position',[newX, newY, pos(3), pos(4)]);
    set(fig,'Units', oldUnits);
end
function a = overlapArea(aPos, bPos)
% OVERLAPAREA  Return the pixel overlap area between two rectangles.
%   aPos, bPos : [x y w h] in pixels.
    ax1=aPos(1); ay1=aPos(2); ax2=ax1+aPos(3); ay2=ay1+aPos(4);
    bx1=bPos(1); by1=bPos(2); bx2=bx1+bPos(3); by2=by1+bPos(4);
    iw = max(0, min(ax2,bx2) - max(ax1,bx1));
    ih = max(0, min(ay2,by2) - max(ay1,by1));
    a = iw*ih;
end
function out = tern(cond, a, b)
% TERN  Inline ternary: returns A if COND is true, B otherwise.
if cond, out = a; else, out = b; end
end