function getExampleData_ekgdin(target_folder)
% GETEXAMPLEDATA_EKGDIN Downloads hardcoded files into a target folder,
% automatically extracting ZIP archives.
%
% This function iterates through a hardcoded structure of URLs and 
% corresponding filenames. For each entry, it checks if the file already 
% exists in the target folder. If it does not exist, the file is downloaded 
% using websave. If the downloaded file is a .zip archive, it is
% automatically extracted and the archive is deleted.
%
% INPUTS:
%   target_folder - A string specifying the path to the directory where
%                   the files should be saved. This folder MUST exist, 
%                   otherwise the function will throw an error.
%
% Example Usage:
%
%   % 1. Define the download location (and ensure it exists)
%   download_path = 'MyStanfordData'; 
%
%   % 2. Run the function
%   getExampleData_ekgdin(download_path);

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

      % 1. Define the data structure (Hardcoded)
      % NOTE: info_struct(1) is set to a hypothetical ZIP file to demonstrate extraction.
      info_struct(1).url = 'https://stacks.stanford.edu/file/xd818jt7842/data1.mat';
      info_struct(1).fn = 'data1.mat';

      info_struct(2).url = 'https://stacks.stanford.edu/file/xd818jt7842/data1_W.mat';
      info_struct(2).fn = 'data1_W.mat';
      
      info_struct(3).url = 'https://stacks.stanford.edu/file/xd818jt7842/data2.mat';
      info_struct(3).fn = 'data2.mat';

      info_struct(4).url = 'https://stacks.stanford.edu/file/xd818jt7842/data2_W.mat';
      info_struct(4).fn = 'data2_W.mat';

      info_struct(5).url = 'https://stacks.stanford.edu/file/xd818jt7842/data3.mat';
      info_struct(5).fn = 'data3.mat';

      info_struct(6).url = 'https://stacks.stanford.edu/file/xd818jt7842/data3_W.mat';
      info_struct(6).fn = 'data3_W.mat';
      
    % --- 1. Validate Target Folder Existence ---
    if ~exist(target_folder, 'dir')
        % Pause execution altogether with a warning/error
        error('getExampleData:FolderNotFound', ...
              'Target folder "%s" does not exist. Please create it manually before running the function.', target_folder);
    else
        fprintf('Data directory: %s\n', target_folder);
    end
    num_files = length(info_struct);
    fprintf('\nProcessing %d files...\n', num_files);
    
    % --- 2. Loop Through Files and Download if Missing ---
    for i = 1:num_files
        current_info = info_struct(i);
        url = current_info.url;
        filename = current_info.fn;
        
        % Construct the full path where the file should be saved
        full_filepath = fullfile(target_folder, filename);
        
        fprintf('  File %d/%d: %s\n', i, num_files, filename);
        
        % Check if the file already exists locally
        % NOTE: If this is a ZIP file and it's deleted after extraction (as requested), 
        % this check will fail on subsequent runs, and the file will be re-downloaded.
        if exist(full_filepath, 'file') == 2
            fprintf('    --> File already exists. Skipping download.\n');
            continue; % Skip to the next file
        end
        
        % File does not exist, proceed with download
        fprintf('    --> File missing. Downloading from: %s\n', url);
        try
            % Use websave to download the file
            websave(full_filepath, url);
            fprintf('    --> Successfully downloaded: %s\n', filename);

            % --- Extraction Logic (Only for ZIP) ---
            
            if endsWith(filename, '.zip', 'IgnoreCase', true)
                
                fprintf('    --> Detected ZIP file. Starting extraction...\n');
                % Unzip the file into the target folder
                unzip(full_filepath, target_folder);
                
                fprintf('    --> Extraction complete. Deleting ZIP file.\n');
                % Delete the original ZIP file after successful extraction
                delete(full_filepath);
                
            else
                
                % Standard file (e.g., .mat, .txt)
                fprintf('    --> File saved successfully as a standard data file.\n');
            end
            % --- END: Extraction Logic ---
            
        catch ME
            % Handle potential errors (e.g., network issue, bad URL, 404)
            warning('GETEXAMPLEDATA:DownloadFailed', ...
                    'Failed to download file %s. Error: %s', ...
                    filename, ME.message);
        end
    end
    fprintf('\nAll file checks and necessary downloads complete.\n');
end



