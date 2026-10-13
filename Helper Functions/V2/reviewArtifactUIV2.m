function [finalSrc, added] = reviewArtifactUIV2(kind, W, xICA, fs, chanlocs, src, cand, info, opts)
% REVIEWARTIFACTUIV2  Shared review UI behind reviewEkgArtifactUIV2 / reviewDINArtifactUIV2.
%
%   [FINALSRC, ADDED] = REVIEWARTIFACTUIV2(KIND, W, XICA, FS, CHANLOCS, SRC, CAND, INFO, OPTS)
%
%   KIND 'ekg' or 'din'. One row per panel component, in panel order, as in the released
%   UIs:
%
%     scalp map | 10 s of the component | |FFT| of the whole component
%
%   at most OPTS.rowsPerPage rows per page; '< Prev' / 'Next >' move between pages.
%   Components in SRC start red (REJECT), those in CAND green (KEEP), exactly as in the
%   released UIs; click anywhere on a row to toggle it. What V2 adds:
%
%     * WHY each row is there, and its rank under every arm, in the row's title.
%     * The time window starts at OPTS.startSec: a number (s, the same for every row),
%       or the code string 'strongest' -- each row at THAT component's strongest chunk
%       (INFO.bestStartSec). The figure's 'Window start' box takes either, + 'Show'.
%     * A second, non-blocking window (OPTS.icChannels) with each component's power per
%       channel -- scalp map of its share, and its top channels with the component's part
%       overlaid (showICChannelsV2).
%     * 'add IC' + '+ Add' appends any component as a new row (on the last page),
%       starting green, so a component you decide to examine because of what another
%       row showed can be reviewed and rejected here. Its number is returned in ADDED.
%     * A status line in the figure; every action is also printed to the Command
%       Window, so a diary records the review.
%
%   OPTS (struct, all optional)
%     .figName     figure name and output file base     .saveFigs   0
%     .saveFolder  pwd                                    .nSec       10
%     .rowsPerPage 4   rows per page; more components = more pages
%     .startSec    'strongest'  window start: a number (s) or 'strongest'
%     .icChannels  true  open the IC power-per-channel window beside the review
%     .maxFreq     |FFT| x-limit (Hz): ekg 10, din 60 (as released)
%
%   Clicking Done saves each page that holds a component flagged for REJECT (red)
%   (when saveFigs: <figName>.png for one page, <figName>_p1.png, _p2.png ... for
%   more, numbered by page; pages with none are not saved) and returns. Closing the window with its X
%   also returns, with the colours as they stand, and saves nothing.
%
%   See also reviewEkgArtifactUIV2, reviewDINArtifactUIV2, autoDetectEkgV2, autoDetectDINV2

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

if nargin < 9 || isempty(opts), opts = struct(); end
isEkg = strcmpi(kind, 'ekg');
opts = fillDefault(opts, 'figName', ['Review ' upper(kind) ' Artifacts (V2)'], 'saveFigs', 0, ...
    'saveFolder', pwd, 'nSec', 10, 'startSec', 'strongest', 'maxFreq', tern(isEkg, 10, 60), ...
    'rowsPerPage', 4, 'icChannels', true);
startSec = parseStart(opts.startSec);
opts.maxFreq = min(opts.maxFreq, fs/2);
rpp = max(1, round(opts.rowsPerPage));

[nComp, N] = size(xICA);
Wi = inv(W);
list = unique([info.panel(ismember(info.panel, [src(:); cand(:)])), src(:)', cand(:)'], 'stable');
isRed = ismember(list, src(:)');
added = zeros(1,0);
page = 1; icFig = [];
cRed = [0.75 0.10 0.10]; cGrn = [0.00 0.55 0.00];
busy = false; pending = false;
hLines = {}; hTtl = gobjects(0); hLayout = [];

% ---- figure + two-row control bar
ctrlH = 72; rowH = 165; figW = 1300;
fig = figure('Name', opts.figName, 'NumberTitle','off', 'Color','w', 'Units','pixels', ...
    'CloseRequestFcn', @onClose);
fig.Position(3:4) = [figW, ctrlH + rpp*rowH + 40];
centerAndClamp(fig);
ctrl = uipanel(fig, 'Units','pixels', 'BackgroundColor',[.97 .97 .97], 'BorderType','none');
bg = [.97 .97 .97];
uicontrol(ctrl, 'Style','text', 'Units','normalized', 'Position',[.01 .55 .45 .40], ...
    'BackgroundColor',bg, 'FontSize',12, 'FontWeight','bold', 'HorizontalAlignment','left', ...
    'String', opts.figName);
hStatus = uicontrol(ctrl, 'Style','text', 'Units','normalized', 'Position',[.47 .55 .52 .40], ...
    'BackgroundColor',bg, 'FontSize',10, 'HorizontalAlignment','left', 'ForegroundColor',[.2 .2 .6], ...
    'String', 'Click a row to toggle red (reject) / green (keep).');
uicontrol(ctrl, 'Style','text', 'Units','normalized', 'Position',[.01 .08 .10 .36], ...
    'BackgroundColor',bg, 'HorizontalAlignment','right', 'String','Window start:');
hStart = uicontrol(ctrl, 'Style','edit', 'Units','normalized', 'Position',[.115 .08 .07 .40], ...
    'String', startString(startSec), 'TooltipString', 'seconds, or ''strongest''');
uicontrol(ctrl, 'Style','pushbutton', 'Units','normalized', 'Position',[.19 .08 .06 .40], ...
    'String','Show', 'Callback', @onShow);
uicontrol(ctrl, 'Style','text', 'Units','normalized', 'Position',[.27 .08 .05 .36], ...
    'BackgroundColor',bg, 'HorizontalAlignment','right', 'String','add IC:');
hAdd = uicontrol(ctrl, 'Style','edit', 'Units','normalized', 'Position',[.325 .08 .06 .40], 'String','');
uicontrol(ctrl, 'Style','pushbutton', 'Units','normalized', 'Position',[.39 .08 .06 .40], ...
    'String','+ Add', 'Callback', @onAdd);
hPrev = uicontrol(ctrl, 'Style','pushbutton', 'Units','normalized', 'Position',[.53 .08 .07 .40], ...
    'String','< Prev', 'Callback', @onPrev);
hPage = uicontrol(ctrl, 'Style','text', 'Units','normalized', 'Position',[.605 .08 .08 .36], ...
    'BackgroundColor',bg, 'HorizontalAlignment','center', 'String','');
hNext = uicontrol(ctrl, 'Style','pushbutton', 'Units','normalized', 'Position',[.69 .08 .07 .40], ...
    'String','Next >', 'Callback', @onNext);
uicontrol(ctrl, 'Style','pushbutton', 'Units','normalized', 'Position',[.86 .08 .12 .40], ...
    'String','Done', 'FontWeight','bold', 'Callback', @onDone);
content = uipanel(fig, 'Units','pixels', 'BackgroundColor','w', 'BorderType','none');
set(fig, 'SizeChangedFcn', @(~,~) layout());
layout();
redraw();
openICChannels(1);
if isempty(list), status('No component on the panel: use add IC to review one, or Done.'); end
uiwait(fig);
finalSrc = list(isRed);
if ishghandle(fig), delete(fig); end
return;

% ================= nested functions =================
    function n = nPages()
        n = max(1, ceil(numel(list) / rpp));
    end

    function layout()
        if ~ishghandle(fig), return; end
        p = fig.Position;
        ctrl.Position = [0 p(4)-ctrlH p(3) ctrlH];
        content.Position = [0 0 p(3) p(4)-ctrlH];
    end

    function redraw()
        % Re-entrancy guard: a callback arriving mid-render is coalesced into one more pass.
        if busy, pending = true; return; end
        busy = true;
        try
            delete(content.Children);
            page = min(max(1, page), nPages());
            set(hPage, 'String', sprintf('page %d of %d', page, nPages()));
            set(hPrev, 'Enable', tern(page > 1, 'on', 'off'));
            set(hNext, 'Enable', tern(page < nPages(), 'on', 'off'));
            hLines = cell(1, numel(list)); hTtl = gobjects(1, numel(list));
            rows = (page-1)*rpp + 1 : min(page*rpp, numel(list));
            if ~isempty(rows)
                T = tiledlayout(content, rpp, 7, 'Padding','compact', 'TileSpacing','compact'); hLayout = T;
                title(T, sprintf('%s review (red = artifact, removed; green = kept)   --   page %d of %d', ...
                    upper(kind), page, nPages()), 'Interpreter','none');
                for k = 1:numel(rows), drawRow(T, k, rows(k)); end
            end
        catch ME
            busy = false; pending = false;
            rethrow(ME);
        end
        busy = false;
        if pending, pending = false; redraw(); end
    end

    function drawRow(T, k, i)
        % k = row on this page, i = index into list
        c = list(i); col = tern(isRed(i), cRed, cGrn);
        s0 = rowStart(c);
        i0 = max(1, floor(s0*fs) + 1); i1 = min(N, i0 + round(opts.nSec*fs) - 1);
        % (1) scalp map -- topoplot_new draws into gca, so pin figure and axes first
        ax1 = nexttile(T, (k-1)*7 + 1);
        set(0, 'CurrentFigure', fig); set(fig, 'CurrentAxes', ax1);
        % topoplot_new resizes its axes, which a tiled layout refuses with a warning on
        % every map (the map still draws). Silence warnings for this call only and put
        % the caller's warning state back -- not warning('off','all') for the whole UI.
        ws = warning('off', 'all');
        try, topoplot_new(Wi(:,c), chanlocs); catch, end
        warning(ws);
        ax1 = gca; axis(ax1, 'off');
        hTtl(i) = title(ax1, sprintf('IC %d', c), 'Color', col, 'FontWeight','bold');
        % (2) the component, nSec s
        ax2 = nexttile(T, (k-1)*7 + 2, [1 3]);
        h2 = plot(ax2, (i0:i1)/fs, xICA(c, i0:i1), 'Color', col, 'LineWidth', 0.9);
        grid(ax2, 'on'); xlim(ax2, [i0 max(i1, i0+1)]/fs); xlabel(ax2, 's');
        title(ax2, rowTitle(c), 'Interpreter','none', 'FontSize', 9);
        % (3) |FFT| of the whole component
        ax3 = nexttile(T, (k-1)*7 + 5, [1 3]);
        x = double(xICA(c,:)) - mean(xICA(c,:)); nf = 2^nextpow2(numel(x));
        X = abs(fft(x, nf)); F = (0:nf/2)*fs/nf; X = X(1:nf/2+1);
        kf = F <= opts.maxFreq;
        h3 = plot(ax3, F(kf), X(kf), 'Color', col, 'LineWidth', 0.7); grid(ax3, 'on');
        xlim(ax3, [0 opts.maxFreq]); xlabel(ax3, 'Hz'); title(ax3, '|FFT| (full series)', 'FontSize', 9);
        hLines{i} = [h2 h3];
        for a = [ax1 ax2 ax3]
            set(a, 'ButtonDownFcn', @(~,~) toggleRow(i));
            set(findall(a, '-not', 'Type','axes'), 'HitTest','off', 'PickableParts','none');
        end
    end

    function openICChannels(startIdx)
        % Reopened (not updated) when the list changes: it is a reference window.
        if ~opts.icChannels || isempty(list), return; end
        if ~isempty(icFig) && ishghandle(icFig), delete(icFig); end
        [~, icFig] = showICChannelsV2(W, xICA, fs, chanlocs, list, ...
            struct('label', upper(kind), 'start', startIdx));
        if ishghandle(fig), figure(fig); end          % the review stays in front
    end

    function closeICChannels()
        if ~isempty(icFig) && ishghandle(icFig), delete(icFig); end
    end

    function s = rowStart(c)
        if ischar(startSec), s = info.bestStartSec(c); else, s = startSec; end
        s = min(max(0, s), max(0, N/fs - opts.nSec));
    end

    function t = rowTitle(c)
        if isEkg
            t = sprintf('released #%s (z %.1f)  |  rrReg2 #%s (z %.1f)', rk(info.rankRel(c)), ...
                info.relZ(c), rk(info.rankRr(c)), info.rrZ(c));
        else
            t = sprintf('%s #%s (%.2f, z %.1f)', info.opts.score, rk(info.rank(c)), info.score(c), info.z(c));
        end
        j = find(info.panel == c, 1);
        if ~isempty(j)
            if ~isEkg, t = info.why{j};                                          % DIN's why is its rank
            elseif startsWith(info.why{j}, 'stands out'), t = ['stands out   --   ' t];
            end                                                                  % else the ranks say it
        elseif ismember(c, added)
            t = ['added   --   ' t];
        end
        t = sprintf('%s   [%.0f-%.0f s]', t, rowStart(c), rowStart(c) + opts.nSec);
    end

    function toggleRow(i)
        isRed(i) = ~isRed(i);
        col = tern(isRed(i), cRed, cGrn);
        if i <= numel(hLines) && ~isempty(hLines{i})
            set(hLines{i}(isgraphics(hLines{i})), 'Color', col);
        end
        if i <= numel(hTtl) && isgraphics(hTtl(i)), set(hTtl(i), 'Color', col); end
        status(sprintf('IC %d -> %s', list(i), tern(isRed(i), 'REJECT (red)', 'KEEP (green)')));
    end

    % Nested, not anonymous: an anonymous @() goPage(page - 1) would freeze PAGE at the
    % value it had when the button was made.
    function onPrev(~,~), goPage(page - 1); end
    function onNext(~,~), goPage(page + 1); end

    function goPage(p)
        if p < 1 || p > nPages(), return; end
        page = p; redraw();
        status(sprintf('Page %d of %d.', page, nPages()));
    end

    function onShow(~,~)
        v = strtrim(get(hStart, 'String'));
        try
            startSec = parseStart(v);
        catch
            status(sprintf('"%s" is neither a start time in seconds nor ''strongest''.', v)); return;
        end
        if ischar(startSec), status('Each row shows its own strongest chunk.');
        else, status(sprintf('All rows show %g-%g s.', startSec, startSec + opts.nSec)); end
        set(hStart, 'String', startString(startSec));
        redraw();
    end

    function onAdd(~,~)
        v = strtrim(get(hAdd, 'String')); c = str2double(v);
        if ~isfinite(c) || c ~= round(c) || c < 1 || c > nComp
            status(sprintf('"%s" is not a component number (1-%d).', v, nComp)); return;
        end
        if ismember(c, list)
            status(sprintf('IC %d is already on the panel.', c)); return;
        end
        list(end+1) = c; isRed(end+1) = false; added(end+1) = c;
        set(hAdd, 'String', '');
        page = nPages(); redraw(); openICChannels(numel(list));
        status(sprintf('IC %d added on page %d (green -- click its row to reject it).', c, page));
    end

    function onDone(~,~)
        status(sprintf('Done: %s = %s', upper(kind), mat2str(list(isRed))));
        if opts.saveFigs && ishghandle(fig), savePages(); end
        closeICChannels(); uiresume(fig);
    end

    function onClose(~,~)
        fprintf('[%s V2 review] window closed: %s = %s (figure not saved)\n', upper(kind), ...
            upper(kind), mat2str(list(isRed)));
        closeICChannels(); uiresume(fig); delete(fig);
    end

    function savePages()
        if isempty(list), return; end                 % nothing drawn, nothing to save
        if ~isfolder(opts.saveFolder), mkdir(opts.saveFolder); end
        base = strrep(opts.figName, ' ', '_'); p0 = page; nSaved = 0;
        for p = 1:nPages()
            % only pages holding a component flagged for REJECT (red)
            if ~any(isRed((p-1)*rpp + 1 : min(p*rpp, numel(list)))), continue; end
            nSaved = nSaved + 1;
            page = p; redraw(); drawnow;
            if nPages() == 1, fn = fullfile(opts.saveFolder, [base '.png']);
            else, fn = fullfile(opts.saveFolder, sprintf('%s_p%d.png', base, p)); end
            try
                % The page's layout, not the figure: exportgraphics warns that the
                % figure's control bar (UI components) cannot be included.
                exportgraphics(hLayout, fn, 'Resolution', 200);
                fprintf('Saved figure to: %s\n', fn);
            catch ME
                warning('reviewArtifactUIV2:SaveError', 'Could not save "%s": %s', fn, ME.message);
            end
        end
        if nSaved == 0, fprintf('No component flagged for REJECT: no page saved.\n'); end
        page = p0; redraw();
    end

    function status(msg)
        if ishghandle(hStatus), set(hStatus, 'String', msg); end
        fprintf('[%s V2 review] %s\n', upper(kind), msg);
    end
end

% -------- utilities --------
function s = rk(r)
if isfinite(r), s = sprintf('%d', r); else, s = '-'; end
end

function v = parseStart(v)
% A start time in seconds, or the code string 'strongest' (also: empty).
if isstring(v), v = char(v); end
if isempty(v) || (ischar(v) && strcmpi(strtrim(v), 'strongest')), v = 'strongest'; return; end
if ischar(v), v = str2double(v); end
if ~isscalar(v) || ~isfinite(v) || v < 0, error('reviewArtifactUIV2:start', 'bad start'); end
end

function s = startString(v)
if ischar(v), s = v; else, s = num2str(v); end
end

function centerAndClamp(fig)
% Centre on the monitor holding most of the figure and keep it fully on screen.
pos = fig.Position; mon = get(0, 'MonitorPositions');
ov = arrayfun(@(i) overlap(pos, mon(i,:)), 1:size(mon,1)); [~, i] = max(ov); m = mon(i,:);
pos(3) = min(pos(3), m(3)); pos(4) = min(pos(4), m(4) - 60);
pos(1) = m(1) + (m(3) - pos(3))/2; pos(2) = m(2) + max(0, (m(4) - pos(4))/2);
fig.Position = pos;
end

function a = overlap(p, q)
a = max(0, min(p(1)+p(3), q(1)+q(3)) - max(p(1), q(1))) * max(0, min(p(2)+p(4), q(2)+q(4)) - max(p(2), q(2)));
end

function out = tern(cond, a, b)
if cond, out = a; else, out = b; end
end

function s = fillDefault(s, varargin)
for i = 1:2:numel(varargin)
    if ~isfield(s, varargin{i}) || (isempty(s.(varargin{i})) && ~strcmp(varargin{i}, 'startSec'))
        s.(varargin{i}) = varargin{i+1};
    end
end
end
