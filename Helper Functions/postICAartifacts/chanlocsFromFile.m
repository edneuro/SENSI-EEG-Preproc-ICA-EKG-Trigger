% Stand-alone version of EEGLab's popular readlocs function.  Based on
% EEGLab version 6.03b, with some mods by Jason Ki ....
% 
% Adapted from topoplot_new from the Parra Lab. https://github.com/dmochow/SRC/blob/master/topoplot_new.m

% Amilcar Malave, SENSI Stanford. February 28, 2026

% READLOCS - read electrode location coordinates and other information from a file. 
%              Several standard file formats are supported. Users may also specify 
%              a custom column format. Defined format examples are given below 
%              (see File Formats).
% Usage:
%   >>  eloc = readlocs( filename );

% Inputs:
%   filename   - Name of the file containing the electrode locations
%                {default: 2-D polar coordinates} (see >> help topoplot )

% Outputs:
%   eloc        - structure containing the channel names and locations (if present).
%                 It has three fields: 'eloc.labels', 'eloc.theta' and 'eloc.radius' 
%                 identical in meaning to the EEGLAB struct 'EEG.chanlocs'.


% Copyright (C) Colin Humphries & Scott Makeig, CNL / Salk Institute, Aug, 1996
%                                          
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


function [handle,Zi,grid,Xi,Yi] = chanlocsFromFile(loc_file)

% -------------------------------------------------------------------------
% Special mode: return chanlocs struct (so caller can edit before plotting)
% Usage:
%   chanlocs = topoplot_new('getchanlocs', 'my_chan.locs');
%   chanlocs = topoplot_new('getchanlocs', EEG.chanlocs); % also allowed
% -------------------------------------------------------------------------

% accept both "string" and 'char'
if isstring(loc_file), loc_file = char(loc_file); end

if ischar(loc_file)
    % filename -> struct
    chanlocs = readlocs(loc_file);   % uses the nested readlocs in this file
elseif isstruct(loc_file)
    % already a struct -> just return it (optionally you could normalize it)
    chanlocs = loc_file;
else
    error('getchanlocs: loc_file must be a filename (char/string) or a chanlocs struct.');
end

% Return chanlocs as first output and exit early
handle = chanlocs;
Zi = [];
grid = [];
Xi = [];
Yi = [];
return


%=========================================================================
%                       readlocs fundtion
%==========================================================================


function [eloc, labels, theta, radius, indices] = readlocs( filename, varargin ); 

if nargin < 1
	help readlocs;
	return;
end;

% NOTE: To add a new channel format:
% ----------------------------------
% 1) Add a new element to the structure 'chanformat' (see 'ADD NEW FORMATS HERE' below):
% 2)  Enter a format 'type' for the new file format, 
% 3)  Enter a (short) 'typestring' description of the format
% 4)  Enter a longer format 'description' (possibly multiline, see ex. (1) below)
% 5)  Enter format file column labels in the 'importformat' field (see ex. (2) below)
% 6)  Enter the number of header lines to skip (if any) in the 'skipline' field
% 7)  Document the new channel format in the help message above.
% 8)  After testing, please send the new version of readloca.m to us
%       at eeglab@sccn.ucsd.edu with a sample locs file.
% The 'chanformat' structure is also used (automatically) by the writelocs() 
% and pop_readlocs() functions. You do not need to edit these functions.

chanformat(1).type         = 'polhemus';
chanformat(1).typestring   = 'Polhemus native .elp file';
chanformat(1).description  = [ 'Polhemus native coordinate file containing scanned electrode positions. ' ...
                               'User must select the direction ' ...
                               'for the nose after importing the data file.' ];
chanformat(1).importformat = 'readelp() function';
% ---------------------------------------------------------------------------------------------------
chanformat(2).type         = 'besa';
chanformat(2).typestring   = 'BESA spherical .elp file';
chanformat(2).description  = [ 'BESA spherical coordinate file. Note that BESA spherical coordinates ' ...
                               'are different from Matlab spherical coordinates' ];
chanformat(2).skipline     = 0; % some BESA files do not have headers
chanformat(2).importformat = { 'type' 'labels' 'sph_theta_besa' 'sph_phi_besa' 'sph_radius' };
% ---------------------------------------------------------------------------------------------------
chanformat(3).type         = 'xyz';
chanformat(3).typestring   = 'Matlab .xyz file';
chanformat(3).description  = [ 'Standard 3-D cartesian coordinate files with electrode labels in ' ...
                               'the first column and X, Y, and Z coordinates in columns 2, 3, and 4' ];
chanformat(3).importformat = { 'channum' '-Y' 'X' 'Z' 'labels'};
% ---------------------------------------------------------------------------------------------------
chanformat(4).type         = 'sfp';
chanformat(4).typestring   = 'BESA or EGI 3-D cartesian .sfp file';
chanformat(4).description  = [ 'Standard BESA 3-D cartesian coordinate files with electrode labels in ' ...
                               'the first column and X, Y, and Z coordinates in columns 2, 3, and 4.' ...
                               'Coordinates are re-oriented to fit the EEGLAB standard of having the ' ...
                               'nose along the +X axis.' ];
chanformat(4).importformat = { 'labels' '-Y' 'X' 'Z' };
chanformat(4).skipline     = 0;
% ---------------------------------------------------------------------------------------------------
chanformat(5).type         = 'loc';
chanformat(5).typestring   = 'EEGLAB polar .loc file';
chanformat(5).description  = [ 'EEGLAB polar .loc file' ];
chanformat(5).importformat = { 'channum' 'theta' 'radius' 'labels' };
% ---------------------------------------------------------------------------------------------------
chanformat(6).type         = 'sph';
chanformat(6).typestring   = 'Matlab .sph spherical file';
chanformat(6).description  = [ 'Standard 3-D spherical coordinate files in Matlab format' ];
chanformat(6).importformat = { 'channum' 'sph_theta' 'sph_phi' 'labels' };
% ---------------------------------------------------------------------------------------------------
chanformat(7).type         = 'asc';
chanformat(7).typestring   = 'Neuroscan polar .asc file';
chanformat(7).description  = [ 'Neuroscan polar .asc file, automatically recentered to fit EEGLAB standard' ...
                               'of having ''Cz'' at (0,0).' ];
chanformat(7).importformat = 'readneurolocs';
% ---------------------------------------------------------------------------------------------------
chanformat(8).type         = 'dat';
chanformat(8).typestring   = 'Neuroscan 3-D .dat file';
chanformat(8).description  = [ 'Neuroscan 3-D cartesian .dat file. Coordinates are re-oriented to fit ' ...
                               'the EEGLAB standard of having the nose along the +X axis.' ];
chanformat(8).importformat = 'readneurolocs';
% ---------------------------------------------------------------------------------------------------
chanformat(9).type         = 'elc';
chanformat(9).typestring   = 'ASA .elc 3-D file';
chanformat(9).description  = [ 'ASA .elc 3-D coordinate file containing scanned electrode positions. ' ...
                               'User must select the direction ' ...
                               'for the nose after importing the data file.' ];
chanformat(9).importformat = 'readeetraklocs';
% ---------------------------------------------------------------------------------------------------
chanformat(10).type         = 'chanedit';
chanformat(10).typestring   = 'EEGLAB complete 3-D file';
chanformat(10).description  = [ 'EEGLAB file containing polar, cartesian 3-D, and spherical 3-D ' ...
                               'electrode locations.' ];
chanformat(10).importformat = { 'channum' 'labels'  'theta' 'radius' 'X' 'Y' 'Z' 'sph_theta' 'sph_phi' ...
                               'sph_radius' 'type' };
chanformat(10).skipline     = 1;
% ---------------------------------------------------------------------------------------------------
chanformat(11).type         = 'custom';
chanformat(11).typestring   = 'Custom file format';
chanformat(11).description  = 'Custom ASCII file format where user can define content for each file columns.';
chanformat(11).importformat = '';
% ---------------------------------------------------------------------------------------------------
% ----- ADD MORE FORMATS HERE -----------------------------------------------------------------------
% ---------------------------------------------------------------------------------------------------

listcolformat = { 'labels' 'channum' 'theta' 'radius' 'sph_theta' 'sph_phi' ...
      'sph_radius' 'sph_theta_besa' 'sph_phi_besa' 'gain' 'calib' 'type' ...
      'X' 'Y' 'Z' '-X' '-Y' '-Z' 'custom1' 'custom2' 'custom3' 'custom4' 'ignore' 'not def' };

% ----------------------------------
% special mode for getting the info
% ----------------------------------
if isstr(filename) & strcmp(filename, 'getinfos')
   eloc = chanformat;
   labels = listcolformat;
   return;
end;

g = finputcheck( varargin, ...
   { 'filetype'	   'string'  {}                 '';
     'importmode'  'string'  { 'eeglab' 'native' } 'eeglab';
     'defaultelp'  'string'  { 'besa'   'polhemus' } 'polhemus';
     'skiplines'   'integer' [0 Inf] 			[];
     'elecind'     'integer' [1 Inf]	    	[];
     'format'	   'cell'	 []					{} }, 'readlocs');
if isstr(g), error(g); end;  

if isstr(filename)
   
   % format auto detection
	% --------------------
   if strcmpi(g.filetype, 'autodetect'), g.filetype = ''; end;
   g.filetype = strtok(g.filetype);
   periods = find(filename == '.');
   fileextension = filename(periods(end)+1:end);
   g.filetype = lower(g.filetype);
   if isempty(g.filetype)
       switch lower(fileextension),
        case {'loc' 'locs' }, g.filetype = 'loc';
        case 'xyz', g.filetype = 'xyz'; 
          fprintf( [ 'WARNING: Matlab Cartesian coord. file extension (".xyz") detected.\n' ... 
                  'If importing EGI Cartesian coords, force type "sfp" instead.\n'] );
        case 'sph', g.filetype = 'sph';
        case 'ced', g.filetype = 'chanedit';
        case 'elp', g.filetype = g.defaultelp;
        case 'asc', g.filetype = 'asc';
        case 'dat', g.filetype = 'dat';
        case 'elc', g.filetype = 'elc';
        case 'eps', g.filetype = 'besa';
        case 'sfp', g.filetype = 'sfp';
        otherwise, g.filetype =  ''; 
       end;
       fprintf('readlocs(): ''%s'' format assumed from file extension\n', g.filetype); 
   else 
       if strcmpi(g.filetype, 'locs'),  g.filetype = 'loc'; end
       if strcmpi(g.filetype, 'eloc'),  g.filetype = 'loc'; end
   end;
   
   % assign format from filetype
   % ---------------------------
   if ~isempty(g.filetype) & ~strcmpi(g.filetype, 'custom') ...
           & ~strcmpi(g.filetype, 'asc') & ~strcmpi(g.filetype, 'elc') & ~strcmpi(g.filetype, 'dat')
      indexformat = strmatch(lower(g.filetype), { chanformat.type }, 'exact');
      g.format = chanformat(indexformat).importformat;
      if isempty(g.skiplines)
         g.skiplines = chanformat(indexformat).skipline;
      end;
      if isempty(g.filetype) 
         error( ['readlocs() error: The filetype cannot be detected from the \n' ...
                 '                  file extension, and custom format not specified']);
      end;
   end;
   
   % import file
   % -----------
   if strcmp(g.filetype, 'asc') | strcmp(g.filetype, 'dat')
       eloc = readneurolocs( filename );
       eloc = rmfield(eloc, 'sph_theta'); % for the conversion below
       eloc = rmfield(eloc, 'sph_theta_besa'); % for the conversion below
       if isfield(eloc, 'type')
           for index = 1:length(eloc)
               type = eloc(index).type;
               if type == 69,     eloc(index).type = 'EEG';
               elseif type == 88, eloc(index).type = 'REF';
               elseif type >= 76 & type <= 82, eloc(index).type = 'FID';
               else eloc(index).type = num2str(eloc(index).type);
               end;
           end;
       end;
   elseif strcmp(g.filetype, 'elc')
       eloc = readeetraklocs( filename );
       %eloc = read_asa_elc( filename ); % from fieldtrip
       %eloc = struct('labels', eloc.label, 'X', mattocell(eloc.pnt(:,1)'), 'Y', ...
       %                        mattocell(eloc.pnt(:,2)'), 'Z', mattocell(eloc.pnt(:,3)'));
       eloc = convertlocs(eloc, 'cart2all');
       eloc = rmfield(eloc, 'sph_theta'); % for the conversion below
       eloc = rmfield(eloc, 'sph_theta_besa'); % for the conversion below
   elseif strcmp(lower(g.filetype(1:end-1)), 'polhemus') | ...
           strcmp(g.filetype, 'polhemus')
       try, 
           [eloc labels X Y Z]= readelp( filename );
           if strcmp(g.filetype, 'polhemusy')
               tmp = X; X = Y; Y = tmp;
           end;
           for index = 1:length( eloc )
               eloc(index).X = X(index);
               eloc(index).Y = Y(index);	
               eloc(index).Z = Z(index);	
           end;
       catch, 
           disp('readlocs(): Could not read Polhemus coords. Trying to read BESA .elp file.');
           [eloc, labels, theta, radius, indices] = readlocs( filename, 'defaultelp', 'besa', varargin{:} );
       end;
   else      
       % importing file
       % --------------
       if isempty(g.skiplines), g.skiplines = 0; end;
       array = load_file_or_array( filename, g.skiplines);
       if size(array,2) < length(g.format)
           fprintf(['readlocs() warning: Fewer columns in the input than expected.\n' ...
                    '                    See >> help readlocs\n']);
       elseif size(array,2) > length(g.format)
           fprintf(['readlocs() warning: More columns in the input than expected.\n' ...
                    '                    See >> help readlocs\n']);
       end;
       
       % removing lines BESA
       % -------------------
       if isempty(array{1,2})
           disp('BESA header detected, skipping three lines...');
           array = load_file_or_array( filename, g.skiplines-1);
           if isempty(array{1,2})
               array = load_file_or_array( filename, g.skiplines-1);
           end;
       end;
       
       % removing comments and empty lines
       % ---------------------------------
       indexbeg = 1;
       while isempty(array{indexbeg,1}) | ...
               (isstr(array{indexbeg,1}) & array{indexbeg,1}(1) == '%' )
           indexbeg = indexbeg+1;
       end;
       array = array(indexbeg:end,:);
       
       % converting file
       % ---------------
       for indexcol = 1:min(size(array,2), length(g.format))
           [str mult] = checkformat(g.format{indexcol});
           for indexrow = 1:size( array, 1)
               if mult ~= 1
                   eval ( [ 'eloc(indexrow).'  str '= -array{indexrow, indexcol};' ]);
               else
                   eval ( [ 'eloc(indexrow).'  str '= array{indexrow, indexcol};' ]);
               end;
           end;
       end;
   end;
   
   % handling BESA coordinates
   % -------------------------
   if isfield(eloc, 'sph_theta_besa')
       if isfield(eloc, 'type')
           if isnumeric(eloc(1).type)
               disp('BESA format detected ( Theta | Phi )');
               for index = 1:length(eloc)
                   eloc(index).sph_phi_besa   = eloc(index).labels;
                   eloc(index).sph_theta_besa = eloc(index).type;
                   eloc(index).labels         = '';
                   eloc(index).type           = '';
               end;
               eloc = rmfield(eloc, 'labels');
           end;
       end;
       if isfield(eloc, 'labels')       
           if isnumeric(eloc(1).labels)
               disp('BESA format detected ( Elec | Theta | Phi )');
               for index = 1:length(eloc)
                   eloc(index).sph_phi_besa   = eloc(index).sph_theta_besa;
                   eloc(index).sph_theta_besa = eloc(index).labels;
                   eloc(index).labels         = eloc(index).type;
                   eloc(index).type           = '';
                   eloc(index).radius         = 1;
               end;           
           end;
       end;
       
       try
           eloc = convertlocs(eloc, 'sphbesa2all');
           eloc = convertlocs(eloc, 'topo2all'); % problem with some EGI files (not BESA files)
       catch, disp('Warning: coordinate conversion failed'); end;
       fprintf('Readlocs: BESA spherical coords. converted, now deleting BESA fields\n');   
       fprintf('          to avoid confusion (these fields can be exported, though)\n');   
       eloc = rmfield(eloc, 'sph_phi_besa');
       eloc = rmfield(eloc, 'sph_theta_besa');

       % converting XYZ coordinates to polar
       % -----------------------------------
   elseif isfield(eloc, 'sph_theta')
       try
           eloc = convertlocs(eloc, 'sph2all');  
       catch, disp('Warning: coordinate conversion failed'); end;
   elseif isfield(eloc, 'X')
       try
           eloc = convertlocs(eloc, 'cart2all');  
       catch, disp('Warning: coordinate conversion failed'); end;
   else 
       try
           eloc = convertlocs(eloc, 'topo2all');  
       catch, disp('Warning: coordinate conversion failed'); end;
   end;
   
   % inserting labels if no labels
   % -----------------------------
   if ~isfield(eloc, 'labels')
       fprintf('readlocs(): Inserting electrode labels automatically.\n');
       for index = 1:length(eloc)
           eloc(index).labels = [ 'E' int2str(index) ];
       end;
   else 
       % remove trailing '.'
       for index = 1:length(eloc)
           if isstr(eloc(index).labels)
               tmpdots = find( eloc(index).labels == '.' );
               eloc(index).labels(tmpdots) = [];
           end;
       end;
   end;
   
   % resorting electrodes if number not-sorted
   % -----------------------------------------
   if isfield(eloc, 'channum')
       if ~isnumeric(eloc(1).channum)
           error('Channel numbers must be numeric');
       end;
       allchannum = [ eloc.channum ];
       if any( sort(allchannum) ~= allchannum )
           fprintf('readlocs(): Re-sorting channel numbers based on ''channum'' column indices\n');
           [tmp newindices] = sort(allchannum);
           eloc = eloc(newindices);
       end;
       eloc = rmfield(eloc, 'channum');      
   end;
else
    if isstruct(filename)
        eloc = filename;
    else
        disp('readlocs(): input variable must be a string or a structure');
    end;        
end;
if ~isempty(g.elecind)
	eloc = eloc(g.elecind);
end;
if nargout > 2
    tmptheta          = { eloc.theta }; % check which channels have (polar) coordinates set
    indices           = find(~cellfun('isempty', tmptheta));
    tmpx              = { eloc.X }; % check which channels have (polar) coordinates set
    indices           = intersect(find(~cellfun('isempty', tmpx)), indices);
    indices           = sort(indices);
    
    indbad            = setdiff(1:length(eloc), indices);
    tmptheta(indbad)  = { NaN };
    theta             = [ tmptheta{:} ];
end;
if nargout > 3
    tmprad            = { eloc.radius };
    tmprad(indbad)    = { NaN };
    radius            = [ tmprad{:} ];
end;
%tmpnum = find(~cellfun('isclass', { eloc.labels }, 'char'));
%disp('Converting channel labels to string');
for index = 1:length(eloc)
    if ~isstr(eloc(index).labels)
        eloc(index).labels = int2str(eloc(index).labels);
    end;
end;
labels = { eloc.labels };
if isfield(eloc, 'ignore')
    eloc = rmfield(eloc, 'ignore');
end;

% process fiducials if any
% ------------------------
fidnames = { 'nz' 'lpa' 'rpa' };
for index = 1:length(fidnames)
    ind = strmatch(fidnames{index}, lower(labels), 'exact');
    if ~isempty(ind), eloc(ind).type = 'FID'; end;
end;

return;

% interpret the variable name
% ---------------------------

function array = load_file_or_array( varname, skiplines );
	 if isempty(skiplines),
       skiplines = 0;
    end;
    if exist( varname ) == 2
        array = loadtxt(varname,'verbose','off','skipline',skiplines);
    else % variable in the global workspace
         % --------------------------
         try, array = evalin('base', varname);
	     catch, error('readlocs(): cannot find the named file or variable, check syntax');
		 end;
    end;     
return;


function array = loadtxt( filename, varargin );

if nargin < 1
	help loadtxt;
	return;
end;	
if ~isempty(varargin)
   try, g = struct(varargin{:});
   catch, disp('Wrong syntax in function arguments'); return; end;
else
    g = [];
end;

g = finputcheck( varargin, { 'convert'   'string'   { 'on' 'off' 'force' }   'on';
                             'skipline'  'integer'  [0 Inf]          0;
                             'verbose'   'string'   { 'on' 'off' }   'on';
                             'delim'     { 'integer' 'string' } []               [9 32];
                             'nlines'    'integer'  []               Inf });
if isstr(g), error(g); end;
g.convert = lower(g.convert);
g.verbose = lower(g.verbose);
g.delim = char(g.delim);

% open the file
% -------------
if exist(filename) ~=2, error( ['file ' filename ' not found'] ); end;  
fid=fopen(filename,'r','ieee-le');
if fid<0, error( ['file ' filename ' found but error while opening file'] ); end;  

index = 0;
while index < abs(g.skipline)
    tmpline = fgetl(fid); 
    if g.skipline > 0 | ~isempty(tmpline)
        index = index + 1;
    end;    
end; % skip lines ---------

inputline = fgetl(fid);
linenb = 1;
if strcmp(g.verbose, 'on'), fprintf('Reading file (lines): '); end;
while isempty(inputline) | inputline~=-1
     colnb = 1;
     if ~isempty(inputline)
	     switch g.convert
	        case 'off',
			     while ~isempty(deblank(inputline))
			         % 07/29/04 Petr Janata added following line to
			         % mitigate problem of strtok ignoring leading
			         % delimiters and deblanking residue in the event
			         % of only space existing between delimiters
			         inputline = strrep(inputline,[g.delim g.delim],[g.delim ' ' g.delim]);
                     
			         [array{linenb, colnb} inputline] = strtok(inputline, g.delim);
			         colnb = colnb+1;
			     end;
	        case 'on',
			     while ~isempty(deblank(inputline))
			         [tmp inputline] = mystrtok(inputline, g.delim);
			         if ~isempty(tmp) & tmp(1) > 43 & tmp(1) < 59, tmp2 = str2num(tmp);
                     else tmp2 = []; end;
			         if isempty( tmp2 )  , array{linenb, colnb} = tmp;
			         else                  array{linenb, colnb} = tmp2;
			         end;
			         colnb = colnb+1;
			     end;
	        case 'force',
			     while ~isempty(deblank(inputline))
			         [tmp inputline] = mystrtok(inputline, g.delim);
			         array{linenb, colnb} = str2double( tmp );
			         colnb = colnb+1;
			     end;
	        otherwise, error('Unrecognized conversion option');
	     end;   
	     linenb = linenb +1;
     end;
     inputline = fgetl(fid);
     if linenb > g.nlines
         inputline = -1;
     end;
     if ~mod(linenb,10) & strcmp(g.verbose, 'on'), fprintf('%d ', linenb); end;
end;        
if strcmp(g.verbose, 'on'),  fprintf('%d\n', linenb-1); end;
if strcmp(g.convert, 'force'), array = [ array{:} ]; end;
fclose(fid); 

% problem strtok do not consider tabulation
% -----------------------------------------
function [str, strout] = mystrtok(strin, delim);
    if delim == 9 % tab
        if length(strin) > 1 & strin(1) == 9 & strin(2) == 9 
            str = '';
            strout = strin(2:end);
        else
            [str, strout] = strtok(strin, delim);
        end;
    else
        [str, strout] = strtok(strin, delim);
    end;

% check field format
% ------------------
function [str, mult] = checkformat(str)
	mult = 1;
	if strcmpi(str, 'labels'),         str = lower(str); return; end;
	if strcmpi(str, 'channum'),        str = lower(str); return; end;
	if strcmpi(str, 'theta'),          str = lower(str); return; end;
	if strcmpi(str, 'radius'),         str = lower(str); return; end;
	if strcmpi(str, 'ignore'),         str = lower(str); return; end;
	if strcmpi(str, 'sph_theta'),      str = lower(str); return; end;
	if strcmpi(str, 'sph_phi'),        str = lower(str); return; end;
	if strcmpi(str, 'sph_radius'),     str = lower(str); return; end;
	if strcmpi(str, 'sph_theta_besa'), str = lower(str); return; end;
	if strcmpi(str, 'sph_phi_besa'),   str = lower(str); return; end;
	if strcmpi(str, 'gain'),           str = lower(str); return; end;
	if strcmpi(str, 'calib'),          str = lower(str); return; end;
	if strcmpi(str, 'type') ,          str = lower(str); return; end;
	if strcmpi(str, 'X'),              str = upper(str); return; end;
	if strcmpi(str, 'Y'),              str = upper(str); return; end;
	if strcmpi(str, 'Z'),              str = upper(str); return; end;
	if strcmpi(str, '-X'),             str = upper(str(2:end)); mult = -1; return; end;
	if strcmpi(str, '-Y'),             str = upper(str(2:end)); mult = -1; return; end;
	if strcmpi(str, '-Z'),             str = upper(str(2:end)); mult = -1; return; end;
	if strcmpi(str, 'custom1'), return; end;
	if strcmpi(str, 'custom2'), return; end;
	if strcmpi(str, 'custom3'), return; end;
	if strcmpi(str, 'custom4'), return; end;
    error(['readlocs(): undefined field ''' str '''']);
   


function chans = convertlocs(chans, command, varargin);

if nargin < 1
   help convertlocs;
   return;
end;

if nargin < 2
   command = 'auto';
end;
if nargin == 4 & strcmpi(varargin{2}, 'on')
    verbose = 1;
else
    verbose = 0; % off
end;

% test if value exists for default
% --------------------------------
if strcmp(command, 'auto')
    if isfield(chans, 'X') & ~isempty(chans(1).X)
        command = 'cart2all';
        if verbose
            disp('Make all coordinate frames uniform using Cartesian coords');
        end;
    else
        if isfield(chans, 'sph_theta') & ~isempty(chans(1).sph_theta)
            command = 'sph2all';
            if verbose
                disp('Make all coordinate frames uniform using spherical coords');
            end;
        else
            if isfield(chans, 'sph_theta_besa') & ~isempty(chans(1).sph_theta_besa)
                command = 'sphbesa2all';
                if verbose
                    disp('Make all coordinate frames uniform using BESA spherical coords');
                end;
            else
                command = 'topo2all';
                if verbose
                    disp('Make all coordinate frames uniform using polar coords');
                end;
            end;
        end;
    end;
end;

% convert
% -------         
switch command
 case 'topo2sph',
   theta  = {chans.theta};
   radius = {chans.radius};
   indices = find(~cellfun('isempty', theta));
   [sph_phi sph_theta] = topo2sph( [ [ theta{indices} ]' [ radius{indices}]' ] );
   if verbose
       disp('Warning: electrodes forced to lie on a sphere for polar to 3-D conversion');
   end;
   for index = 1:length(indices)
      chans(indices(index)).sph_theta  = sph_theta(index);
      chans(indices(index)).sph_phi    = sph_phi  (index);
   end;
   if isfield(chans, 'sph_radius'),
       meanrad = mean([ chans(indices).sph_radius ]);
       if isempty(meanrad), meanrad = 1; end;
   else
       meanrad = 1;
   end;
   sph_radius(1:length(indices)) = {meanrad};
case 'topo2sphbesa',
   chans = convertlocs(chans, 'topo2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
case 'topo2cart'
   chans = convertlocs(chans, 'topo2sph', varargin{:}); % search for spherical coords
   if verbose
       disp('Warning: spherical coordinates automatically updated');
   end;
   chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords
case 'topo2all',
   chans = convertlocs(chans, 'topo2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords
case 'sph2cart',
   sph_theta  = {chans.sph_theta};
   sph_phi    = {chans.sph_phi};
   indices = find(~cellfun('isempty', sph_theta));
   if ~isfield(chans, 'sph_radius'), sph_radius(1:length(indices)) = {1};
   else                              sph_radius = {chans.sph_radius};
   end;
   inde = find(cellfun('isempty', sph_radius));
   if ~isempty(inde)
       meanrad = mean( [ sph_radius{:} ]);
       sph_radius(inde) = { meanrad };
   end;
   [x y z] = sph2cart([ sph_theta{indices} ]'/180*pi, [ sph_phi{indices} ]'/180*pi, [ sph_radius{indices} ]');
   for index = 1:length(indices)
      chans(indices(index)).X = x(index);
      chans(indices(index)).Y = y(index);
      chans(indices(index)).Z = z(index);
   end;
case 'sph2topo',
 if verbose
     % disp('Warning: all radii constrained to one for spherical to topo transformation');
 end;
 sph_theta  = {chans.sph_theta};
 sph_phi    = {chans.sph_phi};
 indices = find(~cellfun('isempty', sph_theta));
 [chan_num,angle,radius] = sph2topo([ ones(length(indices),1)  [ sph_phi{indices} ]' [ sph_theta{indices} ]' ], 1, 2); % using method 2
 for index = 1:length(indices)
     chans(indices(index)).theta  = angle(index);
     chans(indices(index)).radius = radius(index);
     if ~isfield(chans, 'sph_radius') | isempty(chans(indices(index)).sph_radius)
         chans(indices(index)).sph_radius = 1;
     end;
 end;
case 'sph2sphbesa',
   % using polar coordinates
   sph_theta  = {chans.sph_theta};
   sph_phi    = {chans.sph_phi};
   indices = find(~cellfun('isempty', sph_theta));
   [chan_num,angle,radius] = sph2topo([ones(length(indices),1)  [ sph_phi{indices} ]' [ sph_theta{indices} ]' ], 1, 2);
   [sph_theta_besa sph_phi_besa] = topo2sph([angle radius], 1, 1);
   for index = 1:length(indices)
      chans(indices(index)).sph_theta_besa  = sph_theta_besa(index);
      chans(indices(index)).sph_phi_besa    = sph_phi_besa(index);
   end;   
case 'sph2all',
   chans = convertlocs(chans, 'sph2topo', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords
case 'sphbesa2sph',
   % using polar coordinates
   sph_theta_besa  = {chans.sph_theta_besa};
   sph_phi_besa    = {chans.sph_phi_besa};
   indices = find(~cellfun('isempty', sph_theta_besa));
   [chan_num,angle,radius] = sph2topo([ones(length(indices),1)  [ sph_theta_besa{indices} ]' [ sph_phi_besa{indices} ]' ], 1, 1);
   %for index = 1:length(chans)
   %   chans(indices(index)).theta  = angle(index);
   %   chans(indices(index)).radius = radius(index);
   %   chans(indices(index)).labels = int2str(index);
   %end;   
   %figure; topoplot([],chans, 'style', 'blank', 'electrodes', 'labelpoint');
   
   [sph_phi sph_theta] = topo2sph([angle radius], 2);
   for index = 1:length(indices)
      chans(indices(index)).sph_theta  = sph_theta(index);
      chans(indices(index)).sph_phi    = sph_phi  (index);      
   end;
case 'sphbesa2topo',
   chans = convertlocs(chans, 'sphbesa2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2topo', varargin{:}); % search for spherical coords
case 'sphbesa2cart',
   chans = convertlocs(chans, 'sphbesa2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2cart', varargin{:}); % search for spherical coords   
case 'sphbesa2all',
   chans = convertlocs(chans, 'sphbesa2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2all', varargin{:}); % search for spherical coords
case 'cart2topo',
   chans = convertlocs(chans, 'cart2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2topo', varargin{:}); % search for spherical coords
case 'cart2sphbesa',
   chans = convertlocs(chans, 'cart2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2sphbesa', varargin{:}); % search for spherical coords
case 'cart2sph',
    if verbose
        disp('WARNING: If XYZ center has not been optimized, optimize it using Edit > Channel Locations');
	end;
    X  = {chans.X};
    Y  = {chans.Y};
    Z  = {chans.Z};
    indices = find(~cellfun('isempty', X));
    [th phi radius] = cart2sph( [ X{indices} ], [ Y{indices} ], [ Z{indices} ]);
	for index = 1:length(indices)
		 chans(indices(index)).sph_theta     = th(index)/pi*180;
		 chans(indices(index)).sph_phi       = phi(index)/pi*180;
		 chans(indices(index)).sph_radius    = radius(index);
	end;
case 'cart2all',
   chans = convertlocs(chans, 'cart2sph', varargin{:}); % search for spherical coords
   chans = convertlocs(chans, 'sph2all', varargin{:}); % search for spherical coords
end;



function g = fieldtest( fieldname, fieldtype, fieldval, tmpval, callfunc );
	NAME = 1;
	TYPE = 2;
	VALS = 3;
	DEF  = 4;
	SIZE = 5;
    g = [];
    
    switch fieldtype
     case { 'integer' 'real' 'boolean' 'float' }, 
      if ~isnumeric(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be numeric' ]; return;
      end;
      if strcmpi(fieldtype, 'boolean')
          if tmpval ~=0 & tmpval ~= 1
              g = [ callfunc 'error: argument ''' fieldname ''' must be 0 or 1' ]; return;
          end;  
      else 
          if strcmpi(fieldtype, 'integer')
              if ~isempty(fieldval)
                  if (isnan(tmpval) & ~any(isnan(fieldval))) ...
                          & (~ismember(tmpval, fieldval))
                      g = [ callfunc 'error: wrong value for argument ''' fieldname '''' ]; return;
                  end;
              end;
          else % real or float
              if ~isempty(fieldval)
                  if tmpval < fieldval(1) | tmpval > fieldval(2)
                      g = [ callfunc 'error: value out of range for argument ''' fieldname '''' ]; return;
                  end;
              end;
          end;
      end;  
      
      
     case 'string'
      if ~isstr(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be a string' ]; return;
      end;
      if ~isempty(fieldval)
          if isempty(strmatch(lower(tmpval), lower(fieldval), 'exact'))
              g = [ callfunc 'error: wrong value for argument ''' fieldname '''' ]; return;
          end;
      end;

      
     case 'cell'
      if ~iscell(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be a cell array' ]; return;
      end;
      
      
     case 'struct'
      if ~isstruct(tmpval)
          g = [ callfunc 'error: argument ''' fieldname ''' must be a structure' ]; return;
      end;
      
      
     case '';
     otherwise, error([ 'finputcheck error: unrecognized type ''' fieldname '''' ]);
    end;

% remove duplicates in the list of parameters
% -------------------------------------------
function cella = removedup(cella)
% make sure if all the values passed to unique() are strings, if not, exist
%try
    [tmp indices] = unique(cella(1:2:end));
    if length(tmp) ~= length(cella)/2
        fprintf('Note: duplicate ''key'', ''val'' parameter(s), keeping the last one(s)\n');
    end;
    cella = cella(sort(union(indices*2-1, indices*2)));

    function [c, h] = topo2sph(eloc_locs,eloc_angles, method, unshrink)

MAXCHANS = 1024;

if nargin < 1
    help topo2sph;
    return;
end;
if nargin > 1 & ~isstr(eloc_angles)
	if nargin > 2
		unshrink = method;
	end;
	method = eloc_angles;
else
	method = 2;
end;

if isstr(eloc_locs)
	fid = fopen(eloc_locs);
	if fid<1,
	    fprintf('topo2sph()^G: cannot open eloc_loc file (%s)\n',eloc_locs)
	    return
	end
	E = fscanf(fid,'%d %f %f  %s',[7 MAXCHANS]);
	E = E';
	fclose(fid);
else
    E = eloc_locs;
    E = [ ones(size(E,1),1) E ];
end;
    
if nargin > 1 & isstr(eloc_angles)
	if exist(eloc_angles)==2,
	   fprintf('topo2sph: eloc_angles file (%s) already exists and will be erased.\n',eloc_angles);
	end

	fid = fopen(eloc_angles,'a');
	if fid<1,
	    fprintf('topo2sph()^G: cannot open eloc_angles file (%s)\n',eloc_angles)
	    return
	end
end;

if method == 2
	t = E(:,2); % theta
	r = E(:,3); % radius
	h = -t;  % horizontal rotation
	c = (0.5-r)*180;
else
	for e=1:size(E,1)
		% (t,r) -> (c,h)
		
		t = E(e,2); % theta
		r = E(e,3); % radius
		r = r*unshrink;
		if t>=0
			h(e) = 90-t; % horizontal rotation
		else
			h(e) = -(90+t);
		end
		if t~=0
			c(e) = sign(t)*180*r; % coronal rotation
		else
			c(e) = 180*r;
		end
	end;
	t = t';
	r = r';
end;

for e=1:size(E,1)
   if nargin > 1 & isstr(eloc_angles)
        chan = E(e,4:7);
        fprintf('%d	%g	%g	%s\n',E(e,1),c(e),h(e),chan);
        fprintf(fid,'%d	%g	%g	%s\n',E(e,1),c(e),h(e),chan);
   end;     
end


function [g, varargnew] = finputcheck( vararg, fieldlist, callfunc, mode )

	if nargin < 2
		help finputcheck;
		return;
	end;
	if nargin < 3
		callfunc = '';
	else 
		callfunc = [callfunc ' ' ];
	end;
    if nargin < 4
        mode = 'do not ignore';
    end;
	NAME = 1;
	TYPE = 2;
	VALS = 3;
	DEF  = 4;
	SIZE = 5;
	
	varargnew = {};
	% create structure
	% ----------------
	if ~isempty(vararg)
		for index=1:length(vararg)
			if iscell(vararg{index})
				vararg{index} = {vararg{index}};
			end;
		end;
		try
			g = struct(vararg{:});
		catch
            vararg = removedup(vararg);
            try,
                g = struct(vararg{:});
            catch
                g = [ callfunc 'error: bad ''key'', ''val'' sequence' ]; return;
            end;
		end;
	else 
		g = [];
	end;
	
	for index = 1:size(fieldlist,NAME)
		% check if present
		% ----------------
		if ~isfield(g, fieldlist{index, NAME})
			g = setfield( g, fieldlist{index, NAME}, fieldlist{index, DEF});
		end;
		tmpval = getfield( g, {1}, fieldlist{index, NAME});
		
		% check type
		% ----------
        if ~iscell( fieldlist{index, TYPE} )
            res = fieldtest( fieldlist{index, NAME},  fieldlist{index, TYPE}, ...
                           fieldlist{index, VALS}, tmpval, callfunc );
            if isstr(res), g = res; return; end;
        else 
            testres = 0;
            tmplist = fieldlist;
            for it = 1:length( fieldlist{index, TYPE} )
                if ~iscell(fieldlist{index, VALS})
                     res{it} = fieldtest(  fieldlist{index, NAME},  fieldlist{index, TYPE}{it}, ...
                                           fieldlist{index, VALS}, tmpval, callfunc );
                else res{it} = fieldtest(  fieldlist{index, NAME},  fieldlist{index, TYPE}{it}, ...
                                           fieldlist{index, VALS}{it}, tmpval, callfunc );
                end;
                if ~isstr(res{it}), testres = 1; end;
            end;
            if testres == 0,
                g = res{1};
                for tmpi = 2:length(res)
                    g = [ g 10 'or ' res{tmpi} ];
                end;
                return; 
            end;
        end;
	end;
    
    % check if fields are defined
	% ---------------------------
	allfields = fieldnames(g);
	for index=1:length(allfields)
		if isempty(strmatch(allfields{index}, fieldlist(:, 1)', 'exact'))
			if ~strcmpi(mode, 'ignore')
				g = [ callfunc 'error: undefined argument ''' allfields{index} '''']; return;
			end;
			varargnew{end+1} = allfields{index};
			varargnew{end+1} = getfield(g, {1}, allfields{index});
		end;
	end;



function [channo,angle,radius] = sph2topo(input,factor, method)

chans = size(input,1);
angle = zeros(chans,1);
radius = zeros(chans,1);

if nargin < 1
   help sph2topo
   return
end
   
if nargin< 2
  factor = 0;
end
if factor==0
  factor = 1;
end
if factor < 1
  help sph2topo
  return
end

if size(input,2) ~= 3
   help sph2topo
   return
end

channo = input(:,1);
az = input(:,2);
horiz = input(:,3);

if exist('method')== 1 & method == 1
  radius = abs(az/180)/factor;
  i = find(az>=0);
  angle(i) = 90-horiz(i);
  i = find(az<0);
  angle(i) = -90-horiz(i);
else
  angle  = -horiz;
  radius = 0.5 - az/180;
end;