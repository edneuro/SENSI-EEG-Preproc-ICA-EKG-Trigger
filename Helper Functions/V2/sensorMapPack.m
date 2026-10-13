function [x, y, r] = sensorMapPack(x0, y0, layout, compact)
% [x, y, r] = sensorMapPack(x0, y0, layout, compact)
% -----------------------------------------------------
% Choose one uniform marker radius and, if asked, nudge marker positions so the
% markers can be larger and the sparse rim of the montage is less gappy.
% Used by sensorMapAxes.
%
% WHY THIS IS NOT TRIVIAL
% Nearest-neighbour spacing in the GSN 128 montage ranges from 0.0586 (ch15 to
% ch16) to 0.1788 -- a 3.1x spread. A single marker size safe for the tightest
% pair is therefore far too small to hold a three-digit number, and most of the
% figure ends up as white space.
%
%   layout 'exact'   r = 0.45 * min(spacing). True positions, small markers.
%   layout 'packed'  r = 0.46 * median(spacing) (1.58x larger), then:
%                      1. separate to remove the overlaps that creates,
%                      2. optionally compact the sparse rim (see `compact`),
%                      3. separate again,
%                    shrinking r and repeating until zero overlaps is actually
%                    achieved. Convergence is verified, not assumed: an earlier
%                    version left 35 overlapping pairs because the spring
%                    holding markers near their true positions never let the
%                    relaxation settle.
%
% compact (0..1), 'packed' only. Measured on this montage:
%   0.00  extent 1.48 x 1.31, mean move 0.01 marker diameters, 99% of 6-NN kept
%   0.20  extent 1.07 x 1.03, mean move 0.92,                  88%   <- default
%   1.00  extent 0.98 x 0.93, mean move 1.63,                  81%   blob
% Nearly all the gain arrives by 0.20; beyond that distortion rises for almost
% no extra tightness.
%
% After compaction, positions are RELATIVE, not anatomical: a channel keeps its
% neighbours and its side of the head but is not at its true coordinate.
%
% See also sensorMapAxes, showICChannelsV2

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
% Ported from the SENSI CS pipeline (untracked/am_sensorMapPack.m), renamed for this module.

arguments
    x0 (1,:) double
    y0 (1,:) double
    layout (1,:) char {mustBeMember(layout, {'exact','packed'})} = 'packed'
    compact (1,1) double {mustBeInRange(compact, 0, 1)} = 0.20
end

% MEMOISED. The layout is a pure function of (coords, layout, compact) and the
% montage never changes, but the relaxation costs ~24 s for compact > 0. Step 2
% draws the map twice per file, so without this the pipeline would spend ~47 s
% per file -- over two hours across 170 files -- recomputing an identical answer.
persistent CACHE
if isempty(CACHE), CACHE = containers.Map('KeyType','char','ValueType','any'); end
key = sprintf('%s|%.6f|%d|%.10g|%.10g', layout, compact, numel(x0), ...
    sum(x0) + 7*sum(y0), sum(x0.^2) + 7*sum(y0.^2));
if isKey(CACHE, key)
    v = CACHE(key);
    x = v{1}; y = v{2}; r = v{3};
    return
end

n = numel(x0);
D = hypot(x0' - x0, y0' - y0);
D(1:n+1:end) = inf;
nnDist = min(D, [], 2)';

x = x0; y = y0;

switch layout
    case 'exact'
        r = 0.45 * min(nnDist);

    case 'packed'
        r = 0.46 * median(nnDist);
        ok = false;
        for attempt = 1:14
            [tx, ty] = separate(x0, y0, r, 1.0);
            if compact > 0
                [tx, ty] = compactRim(tx, ty, r, compact);
            end
            [tx, ty] = separate(tx, ty, r, 1e6);   % pure separation, no spring
            if minCentreDist(tx, ty) >= 2*r - 1e-9
                x = tx; y = ty; ok = true;
                break
            end
            r = r * 0.96;
        end
        assert(ok, 'Could not find a non-overlapping layout (last r = %.5f).', r);
end

assert(minCentreDist(x, y) >= 2*r - 1e-9, ...
    'Post-condition failed: markers overlap.');

CACHE(key) = {x, y, r};

end

% ------------------------------------------------------------------------

function [x, y] = separate(x0, y0, r, strength, nIter)
% Push equal-radius circles apart until they stop overlapping, with a spring
% holding each near its original position. The spring DECAYS to zero: a constant
% spring fights the separation forever and never reaches zero overlaps.
%
% nIter is capped low when called from inside compactRim's loop: full
% convergence on every one of 250 compaction steps is wasted work, and the final
% separation plus the assert enforce the invariant anyway. Uncapped, this cost
% ~24 s per layout.
x = x0(:); y = y0(:);
n = numel(x);
minD = 2 * r;
springK0 = 0.08 / max(strength, 1e-6);
if nargin < 5 || isempty(nIter), nIter = 600; end

for iter = 1:nIter
    springK = springK0 * max(0, 1 - iter/(0.75*nIter));

    dx = x - x'; dy = y - y';
    d = hypot(dx, dy);
    d(1:n+1:end) = inf;

    overlap = minD - d;
    overlap(overlap < 0) = 0;
    if ~any(overlap(:)) && springK == 0
        break
    end

    coincident = d < 1e-9;
    d(coincident) = 1e-9;
    ux = dx ./ d; uy = dy ./ d;
    if any(coincident(:))
        % Deterministic direction for exactly coincident points. No RNG, so
        % repeated runs give identical layouts.
        [ii, jj] = find(coincident);
        ang = 2*pi * mod((ii + 3*jj), n) / n;
        ux(coincident) = cos(ang);
        uy(coincident) = sin(ang);
    end

    push = 0.5 * overlap;
    x = x + sum(push .* ux, 2);
    y = y + sum(push .* uy, 2);

    x = x + springK * (x0(:) - x);
    y = y + springK * (y0(:) - y);
end

x = x'; y = y';
end

function [x, y] = compactRim(x0, y0, r, strength)
% Pull markers inward wherever there is room, leaving dense areas alone.
%
% Slack is blended between "distance to nearest neighbour" (conservative, keeps
% the silhouette) and "mean distance to the 6 nearest" (frees isolated PAIRS
% like 48/43 and 119/120, which touch each other and so read as having no room
% under the conservative measure). The blend is what `compact` controls.
x = x0(:); y = y0(:);
n = numel(x);
K = 6;
minD = 2 * r;

for iter = 1:250
    D = hypot(x - x', y - y');
    D(1:n+1:end) = inf;

    Ds = sort(D, 2);
    slackNear = Ds(:, 1) - minD;
    slackWide = mean(Ds(:, 1:K), 2) - minD;
    slack = (1-strength)*slackNear + strength*slackWide;
    slack(slack < 0) = 0;

    if max(slack) < 1e-4
        break
    end

    nb = zeros(n, K);
    [~, ord] = sort(D, 2);
    nb(:) = ord(:, 1:K);
    cx = mean(reshape(x(nb), n, K), 2);
    cy = mean(reshape(y(nb), n, K), 2);
    gx = mean(x); gy = mean(y);

    % Global pull scales with strength: light settings tighten locally and keep
    % the shape; heavy ones drag peripheral groups into the mass.
    wg = 0.55 * strength;
    vx = (1-wg)*(cx - x) + wg*(gx - x);
    vy = (1-wg)*(cy - y) + wg*(gy - y);
    L = hypot(vx, vy); L(L < 1e-12) = 1;

    step = min(0.35 * strength * slack, 0.5 * slack);
    x = x + step .* vx ./ L;
    y = y + step .* vy ./ L;

    [xr, yr] = separate(x', y', r, 1e6, 60);   % capped; see separate()
    x = xr(:); y = yr(:);
end

x = x'; y = y';
end

function d = minCentreDist(x, y)
D = hypot(x(:) - x(:)', y(:) - y(:)');
D(1:numel(x)+1:end) = inf;
d = min(D(:));
end
