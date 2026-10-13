function scores = autoDetectEkgTopo(W, chanlocs, nSources, doPlot)
% autoDetectEkgTopo - Automatically detect EKG ICA component via topography
%
% Syntax:
%   scores = autoDetectEkgTopo(W, chanlocs, nSources)
%   scores = autoDetectEkgTopo(W, chanlocs, nSources, doPlot)
%
% Inputs:
%   W         - ICA unmixing matrix (components x channels)
%   chanlocs  - vector of channel location structs for topoplot_new
%   nSources  - vector of component indices to evaluate (e.g., 1:20)
%   doPlot    - (optional) logical flag to generate diagnostic plots (default: false)
%
% Output:
%   scores    - vector of composite length×smoothness scores for each component in nSources
%
if nargin < 4
    doPlot = false;
end

Wi = inv(W);  % mixing matrix (channels x components)

filterSize = 25;  % Gaussian kernel size (odd)
sigma      = 5;   % Gaussian standard deviation
h          = fspecial('gaussian', filterSize, sigma);
numComp    = numel(nSources);

scores        = zeros(1, numComp);
dipoleLength = zeros(1, numComp);
smoothness   = zeros(1, numComp);

f = figure('visible', 'off');

for idx = 1:numComp
    c = nSources(idx);

    % Interpolate scalp map for component c
    [~, Zi, ~, Xi, Yi] = topoplot_new(Wi(:,c), chanlocs, ...
                                      'noplot','on','gridscale',100);
    mask = ~isnan(Zi);

    % Smooth to find poles
    Zf = Zi;
    Zf(~mask) = griddata(Xi(mask), Yi(mask), Zi(mask), Xi(~mask), Yi(~mask), 'nearest');
    Zs = imfilter(Zf, h, 'replicate');
    Zs(~mask) = NaN;

    % Locate positive and negative poles
    [~, pIdx] = max(Zs(:));
    [~, nIdx] = min(Zs(:));
    [rP, cP]  = ind2sub(size(Zs), pIdx);
    [rN, cN]  = ind2sub(size(Zs), nIdx);
    xP = Xi(rP,cP); yP = Yi(rP,cP);
    xN = Xi(rN,cN); yN = Yi(rN,cN);

    % Compute pole-to-pole distance
    L = hypot(xP - xN, yP - yN);
    dipoleLength(idx) = L;

    % Normalize raw map to max absolute amplitude
    maxAbs = max(abs(Zi(mask)));
    Zraw = NaN(size(Zi));
    Zraw(mask) = Zi(mask)/maxAbs;

    % Sample normalized voltages along dipole axis
    N = max(round(L)*2 + 1, 50);
    t = linspace(0,1,N)';
    xLine = xN + t*(xP - xN);
    yLine = yN + t*(yP - yN);
    Zline = interp2(Xi, Yi, Zraw, xLine, yLine, 'linear');

    % Fill remaining NaNs by nearest neighbor interpolation
    nanM = isnan(Zline);
    if any(nanM)
        ok = find(~nanM);
        Zline(nanM) = interp1(ok, Zline(ok), find(nanM), 'nearest');
    end

    % Compute physical step lengths and voltage drops
    dx        = diff(xLine);
    dy        = diff(yLine);
    distStep  = hypot(dx, dy);
    distStep(distStep==0) = eps;
    dZ        = abs(diff(Zline));
    dropUnit  = dZ ./ distStep;

    % Mean and std of drop per unit distance (omit NaNs)
    meanDrop = mean(dropUnit, 'omitnan');
    stdDrop  = std(dropUnit,  'omitnan');

    % Smoothness metrics
    S_mean = 1/(1 + meanDrop);
    S_std  = 1/(1 + stdDrop);
    smoothness(idx) = S_mean * S_std;

    % Composite score
    scores(idx) = L * smoothness(idx);

    % Optional diagnostic plots
    if doPlot
        figure;
        subplot(1,2,1);
        imagesc(Zs); axis off tight; colormap(jet);
        hold on; plot([cP,cN],[rP,rN],'w--','LineWidth',1); hold off;
        title(sprintf('Comp %d: L=%.1f', c, L));
        subplot(1,2,2);
        distances = [0; cumsum(distStep)];
        plot(distances, Zline, '-o','MarkerSize',4);
        xlabel('Distance'); ylabel('Norm voltage');
        title(sprintf('\muΔ=%.3f \sigmaΔ=%.3f', meanDrop, stdDrop));
    end
end

close(f);

end