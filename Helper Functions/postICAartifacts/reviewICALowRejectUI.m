function lowRejectFinal = reviewICALowRejectUI(W, xICA, corrV, corrH, badCh, lowReject, figName, saveFigs, saveFolder)

% REVIEWICALOWREJECTUI Interactive review and rejection of ICA components.
%
%   lowRejectFinal = REVIEWICALOWREJECTUI(W, xICA, corrV, corrH, badCh, lowReject, figName, saveFigs, saveFolder)
%
%   Opens a highly customized, interactive graphical user interface (GUI)
%   to visually review and manually toggle the acceptance or rejection status
%   of pre-flagged Independent Components (ICs). The UI displays the IC
%   topography, time course, and vertical/horizontal EOG correlations for
%   each component.
%
% Inputs:
%   W           (matrix): The ICA unmixing matrix (channels x components). Used
%                       to calculate the inverse (mixing) matrix for topoplots.
%   xICA        (matrix): The ICA component time series (components x time points).
%                       Used to plot the time course of each component.
%   corrV       (vector): Vector of vertical EOG correlation values for all ICs.
%                       Indices must correspond to component numbers.
%   corrH       (vector): Vector of horizontal EOG correlation values for all ICs.
%                       Indices must correspond to component numbers.
%   badCh       (vector): Indices of bad channels to be excluded from the
%                       topoplot visualizations. Pass [] if no bad channels.
%   lowReject   (vector): The **initial list of component indices** that have been
%                       flagged for rejection (the 'low-reject' set). These
%                       components are initialized in the GUI as 'Suspicious' (Yellow).
%   figName     (char):   A string title for the figure window and the saved file.
%   saveFigs    (scalar): Boolean flag (0 or 1). If 1, the final state of the
%                       figure is saved upon clicking 'Done'. Defaults to 0.
%   saveFolder  (char):   Path to the folder where the figure should be saved
%                       if 'saveFigs' is 1. Defaults to the current working
%                       directory (`pwd`).
%
% Output:
%   lowRejectFinal (vector): A vector of component indices that the user
%                           ultimately chose to **reject**. This is the
%                           final list of components to remove from the data.
%
% Usage:
%   % Assuming ICA results are available
%   finalRejects = reviewICALowRejectUI(W, S_ica, V_corr, H_corr, bad_channels, initial_rejects, 'ICA Low-Reject Review');
%
% This function requires 'locsEGI124.mat' and an internal 'topoplotStandalone'
% function (not included in the source snippet) to run correctly.
%
% The UI features pagination, keyboard shortcuts (Next/Back), and direct
% clicking on component rows to toggle the status between REJECT (Red ✗)
% and KEEP (Green ✓), updating the line color and badge glyph immediately.
%
% All components in the input `lowReject` list are initially set to **Suspicious**.
%
% ---------- Early exit ----------
if isempty(lowReject)
    fprintf('\nNo low-reject components detected.\n');
    lowRejectFinal = [];
    return;
end
% --------------------------------

% ---------------- Layout knobs ----------------
pageCap     = 6;
figW        = 830;    % fixed width (px)
figH        = 620;    % initial height (px)
ctrlHpx     = 56;     % control bar height (px)

mLpx        = 28;     % left margin
mRpx        = 18;     % right margin
mTpx_min    = 12;     % min top margin
mBpx        = 22;     % bottom margin
vGapPx      = 10;     % vertical gap

% Color definitions (used in plot and legend)
GREEN_KEEP   = [0.00 0.45 0.00]; % 0
YELLOW_SUSP  = [0.80 0.62 0.10]; % 1 
RED_REJECT   = [0.70 0.25 0.25]; % 2

% Column ratios: [text, topo, time] with total 1+0.5+4.5 = 6
colRat      = [1.0, 0.5, 4.5] / 6;

maxRowHpx   = 105;    % cap row height so few-row pages aren’t huge
txtFS       = 13;     % correlation text font size (larger)
badgeFS     = 16;     % badge font size (larger, more visible)
lineW       = 1.0;    % time-course line width
% ------------------------------------------------
% --- Input Argument Handling for Optional Save Variables ---
if nargin < 9
    saveFolder = pwd; % Default to current directory
end
if nargin < 8
    saveFigs = 0; % Default to NOT saving figures
end
% ------------------------------------------------------------

C        = lowReject(:);
nC       = numel(C);
P        = max(1, ceil(nC / pageCap));

% CRITICAL CHANGE: State is now 0=Keep, 1=Suspicious, 2=Reject
% All components in C start as Suspicious (1).
state = ones(size(C)); 
curPage  = 1;
rowHandles = cell(pageCap, 2);

% Precompute Wi + locs once
Wi = inv(W);
S  = load('locsEGI124.mat','locs'); locs = S.locs;
if ~isempty(badCh), locsPlot = locs; locsPlot(badCh) = []; else, locsPlot = locs; end

% Suppress specific tiledlayout position warning
warning('off','all');

% Figure + panels
fig = figure('Name',figName,'NumberTitle','off','Color','w', ...
    'Units','pixels','Position',[100 100 figW figH], ...
    'WindowKeyPressFcn',@onKey);
ctrl    = uipanel(fig,'Units','pixels','Position',[0 figH-ctrlHpx figW ctrlHpx], ...
    'BackgroundColor',[.97 .97 .97],'BorderType','none');
content = uipanel(fig,'Units','pixels','Position',[0 0 figW figH-ctrlHpx], ...
    'BackgroundColor','w','BorderType','none');
uicontrol(ctrl,'Style','pushbutton','String','< Back','Units','normalized', ...
    'Position',[.01 .12 .10 .76],'Callback',@(~,~)changePage(-1));
uicontrol(ctrl,'Style','pushbutton','String','Next >','Units','normalized', ...
    'Position',[.12 .12 .10 .76],'Callback',@(~,~)changePage(+1));
titleTxt = uicontrol(ctrl,'Style','text','Units','normalized','BackgroundColor',[.97 .97 .97], ...
    'Position',[.24 .12 .52 .76],'FontSize',12,'FontWeight','bold', ...
    'HorizontalAlignment','center','String','');
uicontrol(ctrl,'Style','pushbutton','String','Done (return)','Units','normalized', ...
    'Position',[.80 .12 .19 .76],'FontWeight','bold','Callback',@(~,~)onDone);

% First paint
drawPage();
uiwait(fig);

% lowRejectFinal is now components where state == 2 (Reject)
lowRejectFinal = C(state == 2); 
if ishghandle(fig), close(fig); end

% % Re-enable warnings
warning('on','all');

return;

% =================== Nested helpers ===================
    function drawPage()
        % Determine components on this page
        i0 = (curPage-1)*pageCap + 1;
        i1 = min(curPage*pageCap, nC);
        comps  = C(i0:i1);
        % CRITICAL CHANGE: Get state for components on this page
        pageStates = state(i0:i1); 
        nRows  = numel(comps);
        
        % Clear handles for the new page
        rowHandles = cell(pageCap, 2);
        % 1. Define space for the instruction text banner
        instrHpx = 25; 
        % Adjust figure height to rows (keep width fixed)
        contentHpx_nom = figH - ctrlHpx;
        totalVGap      = max(0,nRows-1) * vGapPx;
        
        % 2. CRITICAL CHANGE: Subtract instruction height from available space
        rowHpx_nom     = (contentHpx_nom - mTpx_min - mBpx - totalVGap - instrHpx) / max(1,nRows); 
        
        rowHpx         = min(rowHpx_nom, maxRowHpx);
        
        % CRITICAL CHANGE: Add instruction height back to used height
        usedH          = nRows*rowHpx + totalVGap + mTpx_min + mBpx + instrHpx; 
        
        figH_new       = ctrlHpx + max(usedH, 200);
        % Resize figure + panels if needed
        pos = get(fig,'Position');
        if pos(4) ~= figH_new
            set(fig,'Position',[pos(1) pos(2) figW figH_new]);
        end
        set(ctrl,   'Position',[0 figH_new-ctrlHpx figW ctrlHpx]);
        set(content,'Position',[0 0 figW figH_new-ctrlHpx]);
        % Recompute content metrics
        contentHpx = figH_new - ctrlHpx;
        usedRowsH  = nRows*rowHpx + totalVGap + mBpx + instrHpx; 
        
        % mTpx is now the space *above* the instructions
        mTpx       = max(mTpx_min, contentHpx - usedRowsH); 
        % Pixel columns from ratios
        Wcontent   = figW - (mLpx + mRpx);
        x_text     = mLpx;
        w_text     = colRat(1) * Wcontent;     % 1/6
        x_topo     = x_text + w_text;
        w_topo     = colRat(2) * Wcontent;     % 0.5/6
        x_time     = x_topo + w_topo;
        w_time     = colRat(3) * Wcontent;     % 4.5/6
        % Clear content and update title
        delete(get(content,'Children'));
        set(titleTxt,'String',sprintf('%s — Page %d/%d', figName, curPage, P));
        % =============================================
        % 💡 SIMPLIFIED IN-FIGURE INSTRUCTIONS/LEGEND 💡
        % This is placed in the newly reserved space (instrHpx).
        % =============================================
        
        instrX = mLpx; 
        instrW = figW * 0.45; 
        instrY = contentHpx - mTpx - instrHpx; 
        uicontrol(content, 'Style', 'text', 'Units', 'pixels', ...
            'Position', [instrX, instrY, instrW, instrHpx], ... 
            'BackgroundColor', 'w', 'HorizontalAlignment', 'left', ...
            'String', sprintf('Toggle status (Low Reject=Red ✗ / Keep=Green ✓)'), ...
            'FontSize', 11, 'FontWeight', 'bold', 'ForegroundColor', [0 0 0]);
        
        % =============================================
        
        for r = 1:nRows
            % Top-aligned rows: r=1 is top row
            % 4. CRITICAL CHANGE: Push component rows down by the instruction height
            yPos = contentHpx - mTpx - instrHpx - r*rowHpx - (r-1)*vGapPx; 
            
            % Get initial state properties
            current_state = pageStates(r);
            switch current_state
                case 0, col=GREEN_KEEP; glyph='✓'; % Keep
                case 1, col=YELLOW_SUSP; glyph='?'; % Suspicious <-- NEW LOGIC
                case 2, col=RED_REJECT; glyph='✗'; % Reject
            end
            % [1] Correlation text (bold, larger)
            ax1 = axes('Parent',content,'Units','pixels', ...
                       'Position',[x_text, yPos, w_text, rowHpx]);
            axis(ax1,'off');
            txtLines = {sprintf('Source %d', comps(r)), ...
                        sprintf('corrV = %.4f', corrV(comps(r))), ...
                        sprintf('corrH = %.4f', corrH(comps(r)))};
            text(ax1, 0.02, 0.62, txtLines, 'Units','normalized', ...
                 'FontSize', txtFS, 'Interpreter','none', 'FontWeight','bold');
            % Place badge and SAVE HANDLE
            hold(ax1,'on');
            h_badge = text(ax1, 0.92, 0.90, glyph, 'Units','normalized', 'FontSize', badgeFS, ...
                 'FontWeight','bold', 'Color', col, 'HorizontalAlignment','center', ...
                 'VerticalAlignment','top', 'HitTest','off');
            hold(ax1,'off');
            % [2] Topoplot (smaller column)
            ax2 = axes('Parent',content,'Units','pixels', ...
                       'Position',[x_topo, yPos, w_topo, rowHpx]);
            cla(ax2);
            topo = Wi(:, comps(r));
            topoplotStandalone(topo, locsPlot);
            % [3] Time course (wide)
            ax3 = axes('Parent',content,'Units','pixels', ...
                       'Position',[x_time, yPos, w_time, rowHpx]);
            
            % Plot line and SAVE HANDLE
            h_line = plot(ax3, xICA(comps(r),:), 'LineWidth', lineW, 'Color', col);
            grid(ax3,'on'); xlim(ax3,'tight');
            % Hide X tick labels except on the bottom row (avoid overlap)
            if r < nRows
                ax3.XTickLabel = [];
            end
            
            % SAVE HANDLES FOR DIRECT UPDATE
            rowHandles{r, 1} = h_badge;
            rowHandles{r, 2} = h_line;
            % Click anywhere in the row toggles (Pass the page row index)
            set(ax1,'ButtonDownFcn', @(~,~)toggle(comps(r), r)); 
            set(ax2,'ButtonDownFcn', @(~,~)toggle(comps(r), r));
            set(ax3,'ButtonDownFcn', @(~,~)toggle(comps(r), r));
        end
        drawnow; 
    end % END for drawPage
    function toggle(src, row_idx)
        % OPTIMIZED: Update properties directly (no full page redraw)
        
        % 1. Update the main state variable
        i = find(C==src,1,'first'); 
        if isempty(i)
            return; 
        end
        
        % CRITICAL CHANGE: Cycle between Keep (0) and Reject (2). 
        % If the state is 1 (Suspicious/Yellow), the first click moves it to Keep (0).
        if state(i) == 0 % Current state is Keep (Green) -> next is Reject (Red)
            state(i) = 2; 
        else % Current state is Suspicious (Yellow, 1) or Reject (Red, 2) -> next is Keep (Green)
            state(i) = 0;
        end
        
        % 2. Calculate new color and glyph
        current_state = state(i);
        switch current_state
            case 0, col=GREEN_KEEP; glyph='✓'; % Keep
            case 1, col=YELLOW_SUSP; glyph='?'; % This should only happen at initial load, not during a toggle
            case 2, col=RED_REJECT; glyph='✗'; % Reject
        end
        
        % 3. Retrieve saved handles for the clicked row
        h_badge = rowHandles{row_idx, 1};
        h_line  = rowHandles{row_idx, 2};
        
        % 4. Update graphics object properties directly
        if ishghandle(h_badge)
            set(h_badge, 'String', glyph, 'Color', col);
        end
        if ishghandle(h_line)
            set(h_line, 'Color', col);
        end
        
        % 5. Force immediate GUI update (Crucial for modern MATLAB graphics)
        drawnow; 
        refresh(fig); 
    end % END for toggle
    function changePage(delta)
        curPage = max(1, min(P, curPage + delta));
        drawPage();
    end % END for changePage
    function onKey(~,ev)
        switch lower(ev.Key)
            case {'rightarrow','n'}, changePage(+1);
            case {'leftarrow','p'},  changePage(-1);
            case {'space','return'}  % click rows to toggle
        end
    end % END for onKey
    function onDone(~,~)
        % Figure saving logic
        if saveFigs && ishghandle(fig)
            try
 
                % Saving as PNG as requested
                savePath = fullfile(saveFolder, [figName, '.png']);
                saveas(fig, savePath, 'png');
            catch ME
                warning('reviewICARejectUI:SaveError', 'Could not save figure: %s\nError: %s', figName, ME.message);
            end
        end
        
        uiresume(fig);
    end 
end