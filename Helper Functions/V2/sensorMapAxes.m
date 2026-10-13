function M = sensorMapAxes(ax, opts)
% M = sensorMapAxes(ax, ...)
% -----------------------------
% Draw the labelled sensor map into an existing axes and return per-channel
% handles, so a caller can recolour markers or attach click callbacks.
%
% This is the drawing core. In this module showICChannelsV2 uses it to colour the
% markers by each component's share of power per electrode (opts.values).
%
% ############################ 128 CHANNELS ONLY ############################
% Hard-coded to the EGI HydroCel GSN 128 montage. It reads
% 'Hydrocel GSN 128 1.0.sfp' and relies on that file's specific layout:
%   entries 1..3 are fiducials (FidNz, FidT9, FidT10)
%   entries 4..131 are E1..E128
%   entry 132 is Cz
% so electrode k is at index k+3. A different net means a different file, a
% different offset, and different packing constants. NOT portable as written --
% asserted at runtime rather than left to fail quietly.
% ###########################################################################
%
% Inputs
% - ax: target axes.
% - opts.flagged / bad / eog : channel lists drawn in the suspicious / bad / EOG colours
% - opts.layout  : 'exact' | 'packed'   (default 'packed')
% - opts.compact : 0..1                 (default 0.20)
% - opts.values  : [optional] 1xN vector, one value per ELECTRODE NUMBER (not
%                  per data row), NaN where there is no data. Switches the
%                  markers to a continuous colormap with a colorbar instead of
%                  the good/suspicious/bad palette. Display only -- setState
%                  still writes categorical colours, so do not mix the modes.
% - opts.valueLabel : colorbar label.
% - opts.valueClim  : [lo hi] colour limits. Without it the scale spans the
%                  full data range, which one dominant electrode can swamp.
% - opts.dimmed  : channels to draw faded and NOT clickable (e.g. channels the
%                  review UI never pages through). Default [].
%
% Output struct M
% - M.h(ch)    : rectangle handle for channel ch
% - M.txt(ch)  : text handle for channel ch
% - M.px, M.py : final marker centres (after packing/compaction)
% - M.r        : marker radius
% - M.col      : the palette, so callers recolour consistently
% - M.setState(ch, state) : recolour one marker, state 0=good 1=susp 2=bad
% - M.setRing(ch, tf)     : draw/remove the "current page" highlight
% - M.setRef(ch)          : mark ch as the panels' reference-overlay channel
%                           ([] or 0 clears it)
%
% See also sensorMapPack, showICChannelsV2

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
% Ported from the SENSI CS pipeline (untracked/am_sensorMapAxes.m), renamed for this module.

arguments
    ax
    opts.flagged double = []
    opts.bad     double = []
    opts.eog     double = []
    opts.dimmed  double = []
    opts.layout  (1,:) char {mustBeMember(opts.layout, {'exact','packed'})} = 'packed'
    opts.compact (1,1) double {mustBeInRange(opts.compact, 0, 1)} = 0.20
    opts.figSizeForFont (1,2) double = [10 9.4]
    opts.axFracForFont  (1,1) double = 0.96
    opts.values     double = []
    opts.valueLabel (1,:) char = ''
    opts.valueClim  double = []
end

N_CH = 128;

%% Electrode coordinates (128 only -- see the banner above)

locFile = 'Hydrocel GSN 128 1.0.sfp';
shift   = 3;

% C must be declared before the evalc: this function contains nested functions,
% which gives it a STATIC workspace, and eval/evalc cannot introduce a new
% variable into one ("Attempt to add C to a static workspace").
C = [];
evalc('C = chanlocsFromFile(locFile);');   % suppress the readlocs banner

assert(numel(C) == N_CH + 4, ...
    ['%s has %d entries; expected %d (3 fiducials + %d electrodes + Cz). ' ...
     'This function is written for the GSN 128 montage only.'], ...
    locFile, numel(C), N_CH + 4, N_CH);
assert(strcmp(char(C(1 + shift).labels), 'E1') && ...
       strcmp(char(C(N_CH + shift).labels), sprintf('E%d', N_CH)), ...
    'Unexpected layout in %s: index %d is "%s", not "E1".', ...
    locFile, 1 + shift, char(C(1 + shift).labels));

theta  = [C((1:N_CH) + shift).theta];
radius = [C((1:N_CH) + shift).radius];

% theta = 0 is the nasion, so this puts the nose at the top. Left hemisphere
% ends up on the LEFT (FidT9 projects to negative x), as in EEGLAB topoplot.
px = radius .* sind(theta);
py = radius .* cosd(theta);

%% Packing

[px, py, r] = sensorMapPack(px, py, opts.layout, opts.compact);

%% Axes setup

hold(ax, 'on'); axis(ax, 'equal'); axis(ax, 'off');

pad = r + 0.030;
xlo = min(px) - pad; xhi = max(px) + pad;
ylo = min(py) - pad; yhi = max(py) + pad + 0.030;
xlim(ax, [xlo xhi]); ylim(ax, [ylo yhi]);

pointsPerUnit = 72 * opts.axFracForFont * opts.figSizeForFont(1) / (xhi - xlo);

text(ax, 0, yhi - 0.004, 'front', 'FontSize', 9, 'FontAngle', 'italic', ...
    'Color', [0.35 0.35 0.40], 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', 'HitTest', 'off', 'PickableParts', 'none');
text(ax, xlo + 0.004, 0, 'L', 'FontSize', 11, 'FontWeight', 'bold', ...
    'Color', [0.55 0.55 0.60], 'HorizontalAlignment', 'left', ...
    'VerticalAlignment', 'middle', 'HitTest', 'off', 'PickableParts', 'none');
text(ax, xhi - 0.004, 0, 'R', 'FontSize', 11, 'FontWeight', 'bold', ...
    'Color', [0.55 0.55 0.60], 'HorizontalAlignment', 'right', ...
    'VerticalAlignment', 'middle', 'HitTest', 'off', 'PickableParts', 'none');

%% Palette -- muted fills so the channel NUMBER stays the readable thing

col.plain  = [0.90 0.90 0.92];
col.flag   = [0.99 0.94 0.78];
col.bad    = [0.97 0.84 0.81];
col.dim    = [0.965 0.965 0.972];
col.eog    = [0.15 0.40 0.80];
col.edge   = [0.55 0.55 0.58];
col.edgeF  = [0.62 0.50 0.20];
col.edgeB  = [0.55 0.25 0.22];
col.edgeDim= [0.85 0.85 0.87];
col.page   = [0.10 0.55 0.25];   % "showing this page" highlight
col.ref    = [0.30 0.30 0.34];   % reference-overlay channel indicator

%%% AM: set up continuous-value mode before the marker loop.
useValues = ~isempty(opts.values) && any(~isnan(opts.values));
vals = nan(1, N_CH); vLo = 0; vHi = 1; cmap = [];
if useValues
    v = reshape(opts.values, 1, []);
    assert(numel(v) <= N_CH, ...
        'opts.values has %d entries; the montage has %d electrodes.', numel(v), N_CH);
    vals(1:numel(v)) = v;
    if ~isempty(opts.valueClim)
        %%% Caller-supplied limits. Needed because a single dominant electrode
        %%% otherwise sets the whole scale: on a focal component 120 of 123
        %%% markers came out the same dark blue and all spatial structure was
        %%% lost. Clipping to a percentile restores the mid-range; values above
        %%% the top limit saturate, so state the true maximum in the title.
        vLo = opts.valueClim(1);
        vHi = opts.valueClim(2);
    else
        vLo = min(vals(~isnan(vals)));
        vHi = max(vals(~isnan(vals)));
    end
    if vHi <= vLo, vHi = vLo + eps; end
    cmap = parula(256);
end

flagged = reshape(intersect(opts.flagged, 1:N_CH), 1, []);
bad     = reshape(intersect(opts.bad,     1:N_CH), 1, []);
eog     = reshape(intersect(opts.eog,     1:N_CH), 1, []);
dimmed  = reshape(intersect(opts.dimmed,  1:N_CH), 1, []);
flagged = reshape(setdiff(flagged, bad), 1, []);

%% Markers

h   = gobjects(1, N_CH);
txt = gobjects(1, N_CH);

for k = 1:N_CH
    st = 0;
    if ismember(k, bad),          st = 2;
    elseif ismember(k, flagged),  st = 1;
    end
    isDim = ismember(k, dimmed);

    [fc, ec, lw] = faceEdge(st, isDim, ismember(k, eog), col);

    %%% AM: continuous-value mode. When opts.values is supplied, the face colour
    %%% comes from a colormap scaled across the supplied values instead of the
    %%% good/suspicious/bad palette -- used to show how much of an ICA component
    %%% lands on each electrode. NaN means "no data for this electrode", which is
    %%% the normal case here: Step 4 works on 124 channels minus the bad ones, so
    %%% 125-128 and every removed electrode have nothing to show and stay dim.
    %%% DISPLAY ONLY -- setState still writes categorical colours, so do not mix
    %%% the two modes on one map.
    if useValues
        if isnan(vals(k))
            fc = col.dim; ec = col.edgeDim; lw = 0.4;
        else
            fc = valueColour(vals(k), vLo, vHi, cmap);
            ec = [0.35 0.35 0.38]; lw = 0.6;
        end
    end

    h(k) = rectangle(ax, 'Position', [px(k)-r, py(k)-r, 2*r, 2*r], ...
        'Curvature', [1 1], 'FaceColor', fc, 'EdgeColor', ec, 'LineWidth', lw);

    % Label. Font is computed from the data-to-point scale so "fits inside"
    % does not depend on window size.
    lbl = num2str(k);
    availPts = 1.75 * r * pointsPerUnit;
    f = min(max(availPts / (0.60 * numel(lbl)), 5.0), 8.5);
    if isDim, tcol = [0.62 0.62 0.66]; else, tcol = [0.16 0.16 0.20]; end
    % PickableParts 'none', not just HitTest 'off'. HitTest only decides what
    % happens once an object has CAPTURED a click; with the default
    % PickableParts the label captures it and passes it to the AXES, never to
    % the marker underneath. Since the label sits dead centre on the marker,
    % that silently swallowed every click on a channel number.
    txt(k) = text(ax, px(k), py(k), lbl, 'FontSize', f, ...
        'FontWeight', ternary(st > 0 && ~isDim, 'bold', 'normal'), ...
        'Color', tcol, 'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', 'HitTest', 'off', ...
        'PickableParts', 'none', 'Clipping', 'off');
end

%% Output

M = struct();
%%% AM: colorbar for continuous-value mode. Placed after the markers so it
%%% does not shift the axes the labels were sized against.
if useValues
    colormap(ax, cmap);
    caxis(ax, [vLo vHi]);
    cb = colorbar(ax, 'eastoutside');
    cb.FontSize = 8;
    if ~isempty(opts.valueLabel)
        cb.Label.String = opts.valueLabel;
        cb.Label.FontSize = 8;
    end
    M.cb = cb;
end

M.h = h; M.txt = txt;
M.px = px; M.py = py; M.r = r;
M.col = col;
M.nCh = N_CH;
M.dimmed = dimmed;
M.ax = ax;
M.ring = gobjects(1, N_CH);
M.refRing = gobjects(1);

M.setState = @setStateFcn;
M.setRing  = @setRingFcn;
M.setRef   = @setRefFcn;

    function setStateFcn(ch, state)
        if ch < 1 || ch > N_CH || ~isvalid(h(ch)), return; end
        isDim = ismember(ch, dimmed);
        [fc, ec, lw] = faceEdge(state, isDim, ismember(ch, eog), col);
        set(h(ch), 'FaceColor', fc, 'EdgeColor', ec, 'LineWidth', lw);
        if isvalid(txt(ch))
            set(txt(ch), 'FontWeight', ternary(state > 0 && ~isDim, 'bold', 'normal'));
        end
    end

    function setRingFcn(ch, tf)
        % Green ring marking "this channel is on the page currently shown".
        if ch < 1 || ch > N_CH, return; end
        if tf
            if ~isgraphics(M.ring(ch))
                rr = r * 1.30;
                M.ring(ch) = rectangle(ax, ...
                    'Position', [px(ch)-rr, py(ch)-rr, 2*rr, 2*rr], ...
                    'Curvature', [1 1], 'EdgeColor', col.page, ...
                    'LineWidth', 2.0, 'FaceColor', 'none', 'HitTest', 'off', ...
                    'PickableParts', 'none');
            end
        elseif isgraphics(M.ring(ch))
            delete(M.ring(ch));
            M.ring(ch) = gobjects(1);
        end
    end

    function setRefFcn(ch)
        % Dashed ring marking the channel drawn as the grey reference overlay
        % in the review panels. Pass [] or 0 to clear.
        %
        % Radius 1.52r reaches slightly into a neighbouring marker, which is
        % why the EOG cue is an edge rather than a ring (several EOG rings
        % overlapped each other). Here it is safe: there is only ever ONE
        % reference channel, so at most one pair of shapes can touch, and the
        % dashed line plus the 1.30r solid green page ring stay tellable apart.
        if isgraphics(M.refRing)
            delete(M.refRing);
            M.refRing = gobjects(1);
        end
        if isempty(ch) || ch < 1 || ch > N_CH, return; end
        rr = r * 1.52;
        M.refRing = rectangle(ax, ...
            'Position', [px(ch)-rr, py(ch)-rr, 2*rr, 2*rr], ...
            'Curvature', [1 1], 'EdgeColor', col.ref, 'LineStyle', '--', ...
            'LineWidth', 1.5, 'FaceColor', 'none', 'HitTest', 'off', ...
            'PickableParts', 'none');
    end

end

% ------------------------------------------------------------------------

function [fc, ec, lw] = faceEdge(state, isDim, isEog, col)
switch state
    case 2, fc = col.bad;   ec = col.edgeB; lw = 0.9;
    case 1, fc = col.flag;  ec = col.edgeF; lw = 0.9;
    otherwise, fc = col.plain; ec = col.edge; lw = 0.5;
end
if isDim
    fc = col.dim; ec = col.edgeDim; lw = 0.4;
end
if isEog
    % EOG as a heavy edge on the marker, not a ring around it: a ring at 1.55x
    % the radius reached into neighbouring markers and the rings overlapped
    % each other (128/32 and 1/125 are adjacent).
    ec = col.eog; lw = 2.2;
end
end

function c = valueColour(v, lo, hi, cmap)
% Map one value onto the colormap. Clamped, so a value at either end of the
% range cannot index outside the table.
f = (v - lo) / (hi - lo);
f = min(max(f, 0), 1);
c = cmap(1 + round(f * (size(cmap,1) - 1)), :);
end

function out = ternary(cond, a, b)
if cond, out = a; else, out = b; end
end
