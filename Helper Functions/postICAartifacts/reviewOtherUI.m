function otherSrcFinal = reviewOtherUI(W, xICA, fs, chanlocs, otherInit, candidates, figName, saveFigs, saveFolder)
% REVIEWOTHERUI (Paged, ≤4 rows/page)
% Interactive reviewer for "OTHER" sources with **automatic figure saving**.
% Each row (4 tiles): [1] Topoplot | [2-4] Time (FULL time series)
% Click a row to toggle: Reject (red ✗) <-> Keep (green •).
% Figures are automatically saved upon clicking 'Done' if requested.
%
% Inputs:
%   W, xICA, fs, chanlocs : as in your pipeline
%   otherInit   : components pre-flagged for REJECT (start red)
%   candidates  : components to review (start green)
%   figName     : title in control bar (default 'Review OTHER Sources (fixed)')
%   saveFigs    : (Optional) Boolean flag (0 or 1). If 1, figures are saved. Defaults to 0.
%   saveFolder  : (Optional) Path to the saving folder. Defaults to pwd.
%
% Output:
%   otherSrcFinal : row vector of final components to REJECT
%
% --------------------------------------------------------------------------

% -------- Defaults / guards for NEW saving inputs --------
if nargin < 9 || isempty(saveFolder), saveFolder = pwd; end
if nargin < 8 || isempty(saveFigs),   saveFigs = 0; end
if nargin < 7 || isempty(figName),    figName = 'Review OTHER Sources (fixed)'; end
if nargin < 6 || isempty(candidates), candidates = []; end
if nargin < 5 || isempty(otherInit),  otherInit  = []; end

Wi = inv(W);

allList = unique([otherInit(:); candidates(:)].','stable');

if isempty(allList)
    otherSrcFinal = [];
    return;
end

% Initial states: otherInit = reject (red), candidates = keep (green)
isOther = ismember(allList, otherInit(:).');

% Time series parameters (FULL SERIES PLOT)
N = size(xICA, 2);
fullIdxPlot = 1:N;
tFull = (0:N-1) / fs; % Full time vector
timeRange = max(eps, tFull(end)); % Total length in seconds

% Colors / style
cREJ  = [0.75 0.10 0.10];   % red
cKEEP = [0.00 0.55 0.00];   % green
lw    = 1.0;

% Paging
nRowsTotal = numel(allList);
pageCap    = 4;
nPages     = ceil(nRowsTotal / pageCap);
curPage    = 1;

% Figure sizing (per page height)
ctrlH   = 54;
rowH    = 170;
figW    = 1100;
basePad = 60;   % a little bottom padding below tiledlayout
figH    = 0;    % Will be set by applyPageSize

% Build figure + control bar
fig = figure('Name', figName, 'NumberTitle','off', 'Color','w', 'Units','pixels');
applyPageSize();
ctrl = uipanel(fig,'Units','pixels','BackgroundColor',[.97 .97 .97],'BorderType','none');
content = uipanel(fig,'Units','pixels','BackgroundColor','w','BorderType','none');

repositionPanels();

% Control bar UI
uicontrol(ctrl,'Style','pushbutton','String','< Back','Units','normalized', ...
    'Position',[.01 .12 .10 .76],'Callback',@(~,~)changePage(-1));
uicontrol(ctrl,'Style','pushbutton','String','Next >','Units','normalized', ...
    'Position',[.12 .12 .10 .76],'Callback',@(~,~)changePage(+1));
titleTxt = uicontrol(ctrl,'Style','text','Units','normalized','BackgroundColor',[.97 .97 .97], ...
    'Position',[.24 .12 .52 .76],'FontSize',12,'FontWeight','bold', ...
    'HorizontalAlignment','center','String','');
uicontrol(ctrl,'Style','pushbutton','String','Done','Units','normalized', ...
    'Position',[.82 .12 .16 .76],'FontWeight','bold','Callback',@(~,~)onDone);

centerAndClamp(fig);
set(fig, 'SizeChangedFcn', @(~,~)centerAndClamp(fig));

% Suppress specific tiledlayout position warning
warning('off','all');

% First page render
drawPage();
uiwait(fig);
otherSrcFinal = allList(isOther);
if ishghandle(fig), close(fig); end

% Re-enable warnings
warning('on','all');

return;


% ================ nested helpers ================
    function drawPage()
        % Which rows on this page
        i0 = (curPage-1)*pageCap + 1;
        i1 = min(curPage*pageCap, nRowsTotal);
        rows = allList(i0:i1);
        nRows = numel(rows);
        
        % Resize to fit this page
        applyPageSize();
        repositionPanels();
        
        % Clear previous content
        delete(get(content,'Children'));
        
        % Tiled layout: nRows x 4 (1 topo | 3 time)
        T = tiledlayout(content, nRows, 4, 'Padding','compact', 'TileSpacing','compact');
        sgtitle(T, 'OTHER Review (red = Reject, green = Keep)', 'Interpreter','none');
        set(titleTxt,'String',sprintf('%s — Page %d/%d', figName, curPage, nPages));
        
        % Keep handles to recolor when toggling
        hTime  = gobjects(nRows,1);
        hBadge = gobjects(nRows,1);

        for r = 1:nRows
            comp = rows(r);
            isR  = isOther(allList == comp);
            col  = tern(isR, cREJ, cKEEP);
            
            % (1) Topoplot (Tile 1)
            ax1 = nexttile(T, (r-1)*4 + 1); % Index: 1
            topoplot_new(Wi(:,comp), chanlocs);
            title(ax1, sprintf('IC %d', comp), 'Interpreter','none');
            axis(ax1, 'off');
            hold(ax1,'on');
            hBadge(r) = text(ax1, 0.95, 0.95, tern(isR,'✗','•'), 'Units','normalized', ...
                'HorizontalAlignment','right','VerticalAlignment','top', ...
                'FontWeight','bold','Color',col,'FontSize',12,'HitTest','off');
            hold(ax1,'off');
            
            % (2) Time (FULL time series) (Tiles 2-4)
            ax2 = nexttile(T, (r-1)*4 + 2, [1 3]); % Index: 2, spans 3 columns
            
            seg = xICA(comp, fullIdxPlot);
            hTime(r) = plot(ax2, tFull, seg, 'LineWidth', lw, 'Color', col);
            grid(ax2,'on');
            xlim(ax2, [0, timeRange]);
            xlabel(ax2,'Time (s)'); ylabel(ax2,'Amp');
            
            % Row toggle (click anywhere in the row)
            set(ax1, 'ButtonDownFcn', @(~,~)toggleRow(comp, hTime(r), hBadge(r)));
            set(ax2, 'ButtonDownFcn', @(~,~)toggleRow(comp, hTime(r), hBadge(r)));
            set(hTime(r),'HitTest','off','PickableParts','none');
        end
    end
    function toggleRow(comp, hT, hB)
        idx = find(allList == comp, 1, 'first');
        if isempty(idx), return; end
        isOther(idx) = ~isOther(idx);
        if isOther(idx)
            newCol = cREJ; badgeChar = '✗';
        else
            newCol = cKEEP; badgeChar = '•';
        end
        if ishghandle(hT), set(hT, 'Color', newCol); end
        if ishghandle(hB), set(hB, 'String', badgeChar, 'Color', newCol); end
    end
    function changePage(delta)
        curPage = max(1, min(nPages, curPage + delta));
        drawPage();
    end
    
    % *** NEW onDone function with save logic ***
    function onDone(~,~)
        if saveFigs && ishghandle(fig)
            saveAllPages(fig, saveFolder, figName);
        end
        uiresume(fig);
    end
    
    % *** NEW function to iterate and save all pages ***
    function saveAllPages(figh, folder, name)
        % Save all pages by iterating through them and redrawing
        
        % Ensure folder exists
        if ~isfolder(folder), mkdir(folder); end
        
        % Save current page index
        originalPage = curPage;
        
        for p = 1:nPages
            % Change page and force redraw
            curPage = p;
            drawPage();
            
            try
                % Create safe filename
                baseName = strrep(name, ' ', '_');
                fileName = sprintf('%s_Page%d.png', baseName, p);
                savePath = fullfile(folder, fileName);
                
                % Save figure (saving as PNG is robust for UI)
                exportgraphics(figh, savePath, 'Resolution', 300);
            catch ME
                warning('ReviewOTHERUI:SaveError', 'Could not save figure page %d: %s', p, ME.message);
            end
        end
        
        % Restore the original page view
        curPage = originalPage;
        drawPage();
    end
    
    function applyPageSize()
        % Height based on rows on current page
        nRows = min(pageCap, nRowsTotal - (curPage-1)*pageCap);
        figH = ctrlH + (nRows * rowH) + basePad;
        fig.Position(3:4) = [figW, figH];
    end
    function repositionPanels()
        pos = get(fig,'Position');
        figHcurr = pos(4);
        set(ctrl,   'Position', [0, figHcurr-ctrlH, figW, ctrlH]);
        set(content,'Position', [0, 0, figW, figHcurr-ctrlH]);
    end
end

% -------- utilities --------
function centerAndClamp(fig)
    if ~ishghandle(fig), return; end
    oldUnits = get(fig,'Units'); set(fig,'Units','pixels');
    pos = get(fig,'Position');     % [x y w h]
    mon = get(0,'MonitorPositions');
    [~, idx] = max(arrayfun(@(i) overlapArea(pos, mon(i,:)), 1:size(mon,1)));
    m = mon(idx,:);
    desiredX = m(1) + (m(3) - pos(3))/2;
    desiredY = m(2) + (m(4) - pos(4))/2;
    newX = min(max(desiredX, m(1)), m(1) + m(3) - pos(3));
    newY = min(max(desiredY, m(2)), m(2) + m(4) - pos(4));
    set(fig,'Position',[newX, newY, pos(3), pos(4)]);
    set(fig,'Units', oldUnits);
end
function a = overlapArea(aPos, bPos)
    ax1=aPos(1); ay1=aPos(2); ax2=ax1+aPos(3); ay2=ay1+aPos(4);
    bx1=bPos(1); by1=bPos(2); bx2=bx1+bPos(3); by2=by1+bPos(4);
    iw = max(0, min(ax2,bx2) - max(ax1,bx1));
    ih = max(0, min(ay2,by2) - max(ay1,by1));
    a = iw*ih;
end
function out = tern(cond, a, b)
    if cond, out = a; else, out = b; end
end