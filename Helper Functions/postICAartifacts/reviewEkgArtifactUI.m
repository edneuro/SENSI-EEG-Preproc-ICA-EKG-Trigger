function ekgSrcFinal = reviewEkgArtifactUI(W, xICA, fs, chanlocs, ekgSrc, candidates, StartSec, NSec, figName, saveFigs, saveFolder)
% reviewEkgArtifactUI
% Interactive EKG artifact reviewer with TILEDLAYOUT and color toggling.
% **NOW INCLUDES OPTIONAL FIGURE SAVING**
% Layout per row (7 tiles): [1] Topoplot | [2-4] Time (windowed) | [5-7] |FFT| (full time series)
% Click anywhere in a row to toggle membership: EKG (red) <-> candidate (green).
% Auto-centers and clamps figure position to stay on-screen.
%
% Inputs (New Inputs Bolded):
%   W, xICA, fs, chanlocs : as in your pipeline
%   ekgSrc      : components pre-flagged for REJECT (start red)
%   candidates  : components to review (start green)
%   StartSec    : start time (s) for the time window
%   NSec        : window length (s) for time plot
%   figName     : title in control bar
%   saveFigs    : (Optional) Boolean flag (0 or 1). If 1, figure is saved. Defaults to 0.
%   saveFolder  : (Optional) Path to the saving folder. Defaults to pwd.

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
    topoplotStandalone(Wi(:,comp), chanlocs);
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

        saveas(figh, savePath, 'png');
        fprintf('Saved figure to: %s\n', savePath);
    catch ME
        warning('reviewEkgArtifactUI:SaveError', 'Could not save figure "%s" to path "%s". Error: %s', name, folder, ME.message);
    end
end
% -------- utilities --------
function centerAndClamp(fig)
% Center the figure on the *current* monitor and clamp it fully on-screen.
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
% Return overlap area between rect A and B (both [x y w h]).
    ax1=aPos(1); ay1=aPos(2); ax2=ax1+aPos(3); ay2=ay1+aPos(4);
    bx1=bPos(1); by1=bPos(2); bx2=bx1+bPos(3); by2=by1+bPos(4);
    iw = max(0, min(ax2,bx2) - max(ax1,bx1));
    ih = max(0, min(ay2,by2) - max(ay1,by1));
    a = iw*ih;
end
function out = tern(cond, a, b)
if cond, out = a; else, out = b; end
end