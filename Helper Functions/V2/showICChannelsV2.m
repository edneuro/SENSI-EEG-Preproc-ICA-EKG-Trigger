function [R, fig] = showICChannelsV2(W, xICA, fs, chanlocs, comps, opts)
% SHOWICCHANNELSV2  Second window: where each component's power lands, channel by channel.
%
%   [R, FIG] = SHOWICCHANNELSV2(W, XICA, FS, CHANLOCS, COMPS, OPTS)
%
%   For each component of COMPS (one at a time, '< Prev IC' / 'Next IC >' or the arrow
%   keys), the left half maps the component's share of power on every channel -- the
%   sensor map (sensorMapAxes) for the EGI HydroCel GSN net (129, or its cropped 128 /
%   124), the regular topoplot_new on CHANLOCS for any other montage -- and the right half
%   shows its top channels: the real channel trace in grey
%   with the component's own part of it, A(c,k)*s_k(t), overlaid in colour. If the two
%   coincide on an artifact, removing the component removes it from that channel; a
%   feature in the grey the colour does not follow will NOT be removed.
%
%   Two numbers per channel, because they answer different questions:
%     % of IC   100 * A(c,k)^2 / sum_c A(c,k)^2   where this component lives (sums to 100)
%     % of ch   100 * var(A(c,k)*s_k) / var(x_c)  how much of this channel it explains
%   The channels shown are ranked by % of IC.
%
%   READ ONLY and NON-BLOCKING: it never calls uiwait, so it can stay open beside a
%   review UI (which does). The channel data are rebuilt as inv(W)*xICA. Data rows are
%   matched to electrodes through CHANLOCS labels ('E12' -> electrode 12), so removed
%   bad channels simply stay dim on the sensor map.
%
%   INPUTS
%     W         [nComp x nCh] unmixing matrix (square)
%     xICA      [nComp x nSamp] component time series
%     fs        sampling rate (Hz)
%     chanlocs  channel locations, one per row of the channel data
%     comps     components to show, in this order
%     opts      (optional struct)
%       .nCh     4      channels per component (2 across)
%       .label   ''     text in the window title (e.g. 'DIN', 'EKG')
%       .visible true   false = compute R only, open no window
%       .start   1      which entry of COMPS to show first
%
%   OUTPUTS
%     R    table: comp, ch, label, pctOfIC, pctOfCh, corrWithCh -- every component x
%          shown-channel pair
%     fig  the window's handle ([] when not visible)
%
%   See also reviewArtifactUIV2, reviewEkgArtifactUIV2, reviewDINArtifactUIV2

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
%
% Ported from the SENSI CS pipeline's am_showICChannels, with its sensor map
% (sensorMapAxes) for EGI nets and topoplot_new for any other montage.

if nargin < 6 || isempty(opts), opts = struct(); end
opts = fillDefault(opts, 'nCh', 4, 'label', '', 'visible', true, 'start', 1);
comps = reshape(comps, 1, []); comps = comps(comps >= 1 & comps <= size(xICA,1));
fig = [];
A = inv(W);
x = A * double(xICA);
[nCh, nSamp] = size(x);
labs = channelLabels(chanlocs, nCh);
elecOf = double(extractAfter(labs, 1));          % 'E12' -> 12; NaN if not an EGI label
% The sensor map is built for the EGI HydroCel GSN 128 layout only. Use it when the channels
% are that net -- every label E1..E129 (or Cz), at most 129 channels: the full 129, the
% cropped 128 / 124, or any of them minus bad channels -- and its montage file is on the
% path. Otherwise the left panel is the regular topoplot on the user's chanlocs.
isEgi = all(~isnan(elecOf) & elecOf >= 1 & elecOf <= 129 | labs == "Cz") && nCh <= 129;
useSensorMap = isEgi && ~isempty(which('Hydrocel GSN 128 1.0.sfp'));
varCh = var(x, 0, 2);
t = (0:nSamp-1) / fs;

comp = []; ch = []; label = strings(0,1); pctOfIC = []; pctOfCh = []; corrWithCh = [];
top = cell(1, numel(comps));
for ii = 1:numel(comps)
    k = comps(ii); a = A(:,k); p = a.^2 / sum(a.^2);
    [~, ord] = sort(p, 'descend'); top{ii} = ord(1:min(opts.nCh, nCh))';
    for c = top{ii}
        bp = a(c) * double(xICA(k,:)); cc = corrcoef(bp, x(c,:));
        comp(end+1,1) = k; ch(end+1,1) = c; label(end+1,1) = labs(c);            %#ok<AGROW>
        pctOfIC(end+1,1) = 100*p(c); pctOfCh(end+1,1) = 100*var(bp)/max(varCh(c), eps); %#ok<AGROW>
        corrWithCh(end+1,1) = cc(1,2);                                           %#ok<AGROW>
    end
end
R = table(comp, ch, label, pctOfIC, pctOfCh, corrWithCh);
if ~opts.visible || isempty(comps), return; end

tag = 'IC power per channel'; if ~isempty(opts.label), tag = [tag ' - ' opts.label]; end
idx = min(max(1, opts.start), numel(comps));
fig = figure('Name', tag, 'NumberTitle','off', 'Color','w', 'Units','normalized', ...
    'Position', [.05 .06 .88 .82], 'WindowKeyPressFcn', @onKey);
bg = [.97 .97 .97];
ctrl = uipanel(fig, 'Units','normalized', 'Position', [0 .95 1 .05], 'BackgroundColor', bg, 'BorderType','none');
uicontrol(ctrl, 'Style','pushbutton', 'String','< Prev IC', 'Units','normalized', ...
    'Position', [.01 .1 .09 .8], 'Callback', @onPrev);
uicontrol(ctrl, 'Style','pushbutton', 'String','Next IC >', 'Units','normalized', ...
    'Position', [.11 .1 .09 .8], 'Callback', @onNext);
lbl = uicontrol(ctrl, 'Style','text', 'Units','normalized', 'Position', [.22 .1 .77 .8], ...
    'BackgroundColor', bg, 'HorizontalAlignment','left', 'FontSize', 10, 'String', '');
mapPanel = uipanel(fig, 'Units','normalized', 'Position', [0 0 .30 .95], 'BackgroundColor','w', 'BorderType','none');
content  = uipanel(fig, 'Units','normalized', 'Position', [.30 0 .70 .95], 'BackgroundColor','w', 'BorderType','none');
draw();

    function draw()
        if ~ishghandle(fig), return; end
        k = comps(idx); sel = top{idx}; a = A(:,k); p = 100 * a.^2 / sum(a.^2);
        delete(mapPanel.Children); delete(content.Children);
        ax = axes(mapPanel, 'Position', [.02 .04 .80 .86]);
        if useSensorMap
            % The sensor map, coloured by this component's share of EVERY electrode. Values
            % are indexed by ELECTRODE number; removed channels (and 125-128) stay dim. The
            % colour scale is clipped at the 97th percentile, since a focal component puts
            % ~40% on one electrode and would turn every other marker the same colour; the
            % true maximum goes in the title.
            try
                ok = ~isnan(elecOf) & elecOf >= 1 & elecOf <= 128;
                vMap = nan(1, 128); vMap(elecOf(ok)) = p(ok);
                vFin = vMap(~isnan(vMap)); vTop = max(prctile(vFin, 97), eps);
                sensorMapAxes(ax, 'values', vMap, 'valueLabel', '% of IC', 'valueClim', [0 vTop], ...
                    'figSizeForFont', [4.4 8.0], 'axFracForFont', 0.80);
                [mx, im] = max(vMap);
                title(ax, {sprintf('IC %d: share per electrode', k), ...
                    sprintf('max %.1f%% at E%d  (scale clipped at %.1f%%)', mx, im, vTop)}, ...
                    'Interpreter','none', 'FontSize', 9);
            catch ME
                % A map failure must not stop the traces from being shown.
                axis(ax, 'off'); text(ax, .5, .5, sprintf('sensor map unavailable:\n%s', ME.message), ...
                    'HorizontalAlignment','center', 'Units','normalized', 'FontSize', 8, 'Color', [.6 .2 .2]);
            end
        else
            % Not the EGI HydroCel 128 layout the sensor map is built for: the regular
            % topoplot of the same share, on the user's own chanlocs. topoplot_new draws
            % into gca, so pin this figure and axes first.
            set(0, 'CurrentFigure', fig); set(fig, 'CurrentAxes', ax);
            try
                topoplot_new(p, chanlocs); ax = gca; axis(ax, 'off');
                clim(ax, [0 max(max(p), eps)]);       % a share is never negative; topoplot_new centres on 0
                cb = colorbar(ax, 'southoutside'); cb.Label.String = '% of IC';
            catch ME
                axis(ax, 'off'); text(ax, .5, .5, sprintf('scalp map unavailable:\n%s', ME.message), ...
                    'HorizontalAlignment','center', 'Units','normalized', 'FontSize', 8, 'Color', [.6 .2 .2]);
            end
            [mx, im] = max(p);
            title(ax, {sprintf('IC %d: share per channel', k), sprintf('max %.1f%% at %s', mx, labs(im))}, ...
                'Interpreter','none', 'FontSize', 9);
        end
        TL = tiledlayout(content, ceil(numel(sel)/2) + 1, 2, 'Padding','compact', 'TileSpacing','compact');
        title(TL, sprintf('IC %d   --   grey = the real channel, colour = this IC''s part of it', k), ...
            'Interpreter','none', 'FontWeight','bold');
        ax0 = nexttile(TL, 1, [1 2]);
        plot(ax0, t, xICA(k,:), 'Color', [.2 .3 .7], 'LineWidth', .6); grid(ax0, 'on');
        xlim(ax0, t([1 end])); ylabel(ax0, 'IC'); title(ax0, sprintf('IC %d', k), 'Interpreter','none');
        for r = 1:numel(sel)
            c = sel(r); axc = nexttile(TL);
            plot(axc, t, x(c,:), 'Color', [.6 .6 .62], 'LineWidth', .5); hold(axc, 'on');
            plot(axc, t, a(c) * double(xICA(k,:)), 'Color', [.85 .33 .1], 'LineWidth', .6); hold(axc, 'off');
            grid(axc, 'on'); xlim(axc, t([1 end]));
            q = R(R.comp == k & R.ch == c, :);
            title(axc, sprintf('%s   %.1f%% of IC   %.1f%% of channel   r = %.2f', labs(c), ...
                q.pctOfIC(1), q.pctOfCh(1), q.corrWithCh(1)), 'Interpreter','none', 'FontSize', 9);
            if r > numel(sel) - 2, xlabel(axc, 's'); end
        end
        set(lbl, 'String', sprintf('  IC %d   (%d of %d)   top channels: %s   [arrow keys page]', ...
            k, idx, numel(comps), strjoin(labs(sel), ', ')));
    end

    function onPrev(~,~), idx = max(1, idx - 1); draw(); end
    function onNext(~,~), idx = min(numel(comps), idx + 1); draw(); end
    function onKey(~, ev)
        switch lower(ev.Key)
            case {'rightarrow','n'}, onNext();
            case {'leftarrow','p'},  onPrev();
        end
    end
end

% -------------------------------------------------------------------------
function labs = channelLabels(chanlocs, nCh)
labs = "ch" + string(1:nCh);
if isstruct(chanlocs) && isfield(chanlocs, 'labels') && numel(chanlocs) == nCh
    for i = 1:nCh
        if ~isempty(chanlocs(i).labels), labs(i) = string(chanlocs(i).labels); end
    end
end
end

function s = fillDefault(s, varargin)
for i = 1:2:numel(varargin)
    if ~isfield(s, varargin{i}) || isempty(s.(varargin{i})), s.(varargin{i}) = varargin{i+1}; end
end
end
