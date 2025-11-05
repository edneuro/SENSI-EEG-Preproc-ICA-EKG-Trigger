function [dinSrc, info] = autoDetectDIN(W, xICA, badCh, fs, INFO, nSources, ekgStartSec, ekgNSec, HRFOpts)
% AUTODETECTDIN
% Detect DIN ICA components via harmRatioFromSignal + |z|>thress on the selected sources.
% Plot DIN sources together (topo + wide time + wide spectrum), and a separate
% figure for the single strongest non-DIN candidate for review.
%
% Usage:
%   [dinSrc, info] = autoDetectDIN(W, xICA, badCh, fs, INFO, nSources, ekgStartSec, ekgNSec, HRFOpts)

% ---- Defaults
if nargin < 6 || isempty(nSources),    nSources    = 1:30; end
if nargin < 7 || isempty(ekgStartSec), ekgStartSec = 5*60; end
if nargin < 8 || isempty(ekgNSec),     ekgNSec     = 10;   end
if nargin < 9 || isempty(HRFOpts),     HRFOpts     = struct('fmin',10,'fmax',50,'bw',0.2); end

thress = 2.5;

% ---- Shapes / guards
[nComp, nSamp] = size(xICA);
nSources = unique(nSources(:)');
nSources = nSources(nSources>=1 & nSources<=nComp);
if isempty(nSources)
    dinSrc = [];
    info   = struct('dinCheck',[],'z',[],'topNonDIN',[]);
    warning('autoDetectDIN:NoValidSources','No valid sources to scan.');
    return;
end

% ---- Compute DIN metric
dinCheck = nan(nComp,1);
for k = nSources
    out = harmRatioFromSignal(xICA(k,:), fs, HRFOpts);
    dinCheck(k) = out.ratio_linear;
end

% ---- Z-score over selected sources
vals  = dinCheck(nSources);
zVals = zscore(vals, 0, 'omitnan');  

z = nan(nComp,1);
z(nSources) = zVals;

% ---- Flag DIN components
dinSrc = nSources(zVals > thress);
dinSrc = dinSrc(:).';

% ---- Top non-DIN candidate via two-stage z-score
% Stage 1: initial flags
initialMask = zVals > thress;
dinSrc = nSources(initialMask)';
% Stage 2: remap flagged to median and recompute
medVal = median(dinCheck(nSources),'omitnan');
dinAdj = dinCheck;
dinAdj(dinSrc) = medVal;
vals2 = dinAdj(nSources);
z2    = zscore(vals2, 0, 'omitnan');
% review candidates: z2>thress & not initially flagged

cand2 = nSources(z2>thress & ~initialMask);

if ~isempty(cand2)
    [~, idxMax] = max(z2(ismember(nSources, cand2)));
    topNonDIN = cand2(idxMax);
else
    topNonDIN = [];
end
% capture info
info = struct('dinCheck',dinCheck,'z',z,'dinAdj',dinAdj,'z2',z2,'reviewCandidates',cand2,'topNonDIN',topNonDIN);


% % ---- Prepare plotting parameters
% plotIdx = round(ekgStartSec*fs) + (1:round(ekgNSec*fs));
% plotIdx = plotIdx(plotIdx>=1 & plotIdx<=nSamp);
% maxFreq = 50;
% 
% % ---- Load and adjust locations
% load('locsEGI124.mat','locs');
% goodBad = badCh(badCh>=1 & badCh<=numel(locs));
% locsAdj = locs;
% locsAdj(goodBad) = [];
% Wi = inv(W);
% 
% % ---- Plot DIN sources
% % Suppress layout position warnings under tiledlayout
% warning('off','all');
% if ~isempty(dinSrc)
%     nRow = numel(dinSrc);
%     fig = figure('Name','DIN flagged','Color','w','Units','pixels');
%     fig.Position(3) = 900;
%     fig.Position(4) = nRow*125 + 50;
%     T = tiledlayout(nRow, 7, 'TileSpacing','compact','Padding','compact');
%     for r = 1:nRow
%         src = dinSrc(r);
%         % Topo (col 1)
%         nexttile(T, (r-1)*7 + 1, [1 1]);
%         topoplotStandalone(Wi(:,src), locsAdj);
%         title(sprintf('Source %d', src), 'Interpreter','none');
%         % Time (cols 2–4)
%         nexttile(T, (r-1)*7 + 2, [1 3]);
%         plot(plotIdx/fs, xICA(src,plotIdx), 'LineWidth',1); grid on;
%         xlabel('Time (s)'); ylabel('Amplitude');
%         title('Time Domain','Interpreter','none');
%         % Freq (cols 5–7)
%         nexttile(T, (r-1)*7 + 5, [1 3]);
%         data = xICA(src,:);
%         fAx  = computeFFTFrequencyAxis(length(data), fs);
%         plot(fAx, abs(fft(data)), 'LineWidth',1); grid on;
%         xlim([0 maxFreq]); xlabel('Freq (Hz)');
%         title('Frequency Domain','Interpreter','none');
%     end
%     fn   = INFO.preproc4_MIR_fNameLoaded;
%     base = fn(1:end-4);
%     sgtitle(T, sprintf('%s: DIN flagged [%s]', base, mat2str(dinSrc')), 'Interpreter','none');
% end
% 
% % ---- Plot review candidate
% % Suppress warnings again for review figure
% warning('off','MATLAB:layout:UnableToSetPosition');
% if ~isempty(topNonDIN)
%     fig = figure('Name','Review (not DIN)','Color','w','Units','pixels');
%     fig.Position(3) = 900;
%     fig.Position(4) = 125 + 50;
%     T2 = tiledlayout(1, 7, 'TileSpacing','compact','Padding','compact');
%     src = topNonDIN;
%     nexttile(T2, 1, [1 1]);
%     topoplotStandalone(Wi(:,src), locsAdj);
%     title(sprintf('Source %d', src), 'Interpreter','none');
%     nexttile(T2, 2, [1 3]);
%     plot(plotIdx/fs, xICA(src,plotIdx), 'LineWidth',1); grid on;
%     xlabel('Time (s)'); ylabel('Amplitude');
%     title('Time Domain','Interpreter','none');
%     nexttile(T2, 5, [1 3]);
%     data = xICA(src,:);
%     fAx  = computeFFTFrequencyAxis(length(data), fs);
%     plot(fAx, abs(fft(data)), 'LineWidth',1); grid on;
%     xlim([0 maxFreq]); xlabel('Freq (Hz)');
%     title('Frequency Domain','Interpreter','none');
%     sgtitle(T2, sprintf('%s: Review candidate (not DIN)', base), 'Interpreter','none');
% end
% 
% % Re-enable warnings
% warning('on','all');

end
