function reviewICAHighRejectSAVE(W, xICA, corrV, corrH, badCh, highReject, figName, saveFigs, saveFolder)
% REVIEWICAHIGHREJECTSAVE Plotting and automatic saving of high-reject ICA components.
%
%   REVIEWICAHIGHREJECTSAVE(W, xICA, corrV, corrH, badCh, highReject, figName, saveFigs, saveFolder)
%
%   Generates one or more non-interactive figures to display topographies,
%   time courses (in RED), and EOG correlations for components flagged for
%   high-frequency rejection. Figures are split every 6 components.
%
% Inputs:
%   W           (matrix): The ICA unmixing matrix.
%   xICA        (matrix): The ICA component time series.
%   corrV       (vector): Vertical EOG correlation values.
%   corrH       (vector): Horizontal EOG correlation values.
%   badCh       (vector): Indices of bad channels (use [] if none).
%   highReject  (vector): The list of component indices flagged for rejection.
%   figName     (char):   A string title for the figure window(s) and file prefix.
%   saveFigs    (scalar): Boolean flag (0 or 1). If 1, figures are saved. Defaults to 0.
%   saveFolder  (char):   Path to the saving folder. Defaults to pwd.
%
% Outputs:
%   (None) - This function is purely for visualization and saving.
%
% --------------------------------------------------------------------------

% --- Input Argument Handling (For optional save parameters) ---
if nargin < 9
    saveFolder = pwd;
end
if nargin < 8
    saveFigs = 0;
end

% --- Early exit ---
if isempty(highReject)
    fprintf('\nNo high-reject components detected.\n');
    return;
end

% ---------------- Recursive Logic for Figure Splitting -------------------

nSources = length(highReject);
pageCap = 6; % Max components per figure

if nSources > pageCap
    nFigs = ceil(nSources / pageCap);
    
    for i = 1:(nFigs - 1)
       % Recursive call for all pages except the last one
       reviewICAHighRejectSAVE(W, xICA, corrV, corrH, badCh,...
           highReject((1:pageCap) + (i-1)*pageCap),...
           [figName ' - Page ' num2str(i) ' of ' num2str(nFigs)],...
           saveFigs, saveFolder);
    end
    
    % Recursive call for the final page
    reviewICAHighRejectSAVE(W, xICA, corrV, corrH, badCh,...
        highReject(((nFigs-1)*pageCap + 1):end),...
        [figName ' - Page ' num2str(nFigs) ' of ' num2str(nFigs)],...
        saveFigs, saveFolder);
    
    return
end

% ----------------- Plotting Logic (Executed for each figure) -----------------

% Load the pre-saved locs file (Required for topoplotStandalone)
load('locsEGI124.mat', 'locs');

% Remove coordinates for bad channels
if length(badCh) > 0
    locs(badCh) = [];
end

% Calculate inverse matrix for topoplots
Wi = inv(W);

% --- Plotting Parameters ---
nrow = length(highReject);
ncol = 10; % Columns: 1 (Text), 1 (Topo), 8 (Timecourse)
corrFontSize = 18;
RED_REJECT = [0.70 0.25 0.25]; % Color for time course

% --- Create Figure ---
fig = figure('name', figName, 'numbertitle', 'off');

% NOTE: Requires 'tight_subplot' to be in the MATLAB path.
% Using the original outer margin '0' is crucial for the text positioning to work.
hh = tight_subplot(nrow, ncol, 0.0005, 0); 

% --- Loop to Draw Components ---
for i = 1:nrow
    sourceNum = highReject(i);
    
    % The column index for the start of the current row (i-1)*ncol + 1
    rowStartIdx = (i-1) * ncol + 1;

    % 1. Correlation text (Column 1)
    % Use the tight_subplot handle and the original negative X position
    textAxes = hh(rowStartIdx + 0); 
    axes(textAxes);
    set(gca, 'Visible', 'off');
    
    % Correct text placement: -0.5 puts the text in the figure's left margin
    text(-0.5, 0.5, ...
        sprintf('Source %d \ncorrV = %0.4f \ncorrH = %0.4f', ...
        sourceNum, corrV(sourceNum), corrH(sourceNum)), ...
        'fontsize', corrFontSize);
    
    % 2. Topoplot (Column 2)
    topoAxes = hh(rowStartIdx + 1);
    axes(topoAxes)
    tempTopo = Wi(:, sourceNum);
    % NOTE: Assumes topoplotStandalone is available
    topoplotStandalone(tempTopo, locs);
    
    % 3. Timecourse of the source (Columns 3 through 10)
    % Use the standard subplot function to cleanly span the required columns.
    subplot(nrow, ncol, (3:ncol) + (i-1) * ncol)
    
    % Plot in RED, as requested
    plot(xICA(sourceNum, :), 'Color', RED_REJECT); 
    grid on; 
    xlim('tight');
    
    % Hide X tick labels on all but the bottom row
    if i < nrow
        set(gca, 'XTickLabel', []);
    end
end

% Set the overall figure title
sgtitle(figName, 'interpreter', 'none', 'fontsize', corrFontSize);

% ------------------------- Saving Logic -------------------------
if saveFigs
    if ishghandle(fig)
        try
            % Ensure the file name is safe and includes the page number from figName
            safeName = strrep(figName, ' ', '_'); 
            safeName = strrep(safeName, '-', '');
            safeName = strrep(safeName, '/', '_');
            
            saveFileName = sprintf('%s.png', safeName);
            savePath = fullfile(saveFolder, saveFileName);
            
            saveas(fig, savePath, 'png');
        catch ME
            warning('reviewICAHighRejectSAVE:SaveError', 'Could not save figure "%s" to path "%s". Error: %s', figName, saveFolder, ME.message);
        end
    end
end
% ----------------------------------------------------------------

end % END of function