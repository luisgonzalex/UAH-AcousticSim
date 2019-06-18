%function setupmatfile = EDB2readsetup(EDsetupfile,cadgeofiletoopen);
function setupmatfile = EDB2readsetup(cadgeofiletoopen,api,infos);
% EDB2readsetup - Runs the EDsetup m-file and saves all the settings in a mat-file.
% If the sources and receivers should be taken from the
% CAD-file, they will be extracted and saved in the setupmatfile.
%
% Input parameters:
%   EDsetupfile         The name of the EDsetupfile (an m-file) that
%                       should be read.
%   cadgeofiletoopen    (optional) If the CAD file has already been read
%                       and stored as a cadgeofile, the cadgeofile could
%                       be used instead of the CAD-file.
%   SHOWTEXT (global)   Determines whether there will be a pause or just
%                       printed information when another setting is
%                       recommended.
%
% Output parameters:
%   setupmatfile        The name of the output file. The name is created
%                       from the Filestem in the setupfile, with
%                       "_setup.mat" added.
% 
% ----------------------------------------------------------------------------------------------
%   This file is part of the Edge Diffraction Toolbox by Peter Svensson.                       
%                                                                                              
%   The Edge Diffraction Toolbox is free software: you can redistribute it and/or modify       
%   it under the terms of the GNU General Public License as published by the Free Software     
%   Foundation, either version 3 of the License, or (at your option) any later version.        
%                                                                                              
%   The Edge Diffraction Toolbox is distributed in the hope that it will be useful,       
%   but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS  
%   FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.             
%                                                                                              
%   You should have received a copy of the GNU General Public License along with the           
%   Edge Diffraction Toolbox. If not, see <http://www.gnu.org/licenses/>.                 
% ----------------------------------------------------------------------------------------------
% Peter Svensson 23 April 2015 (svensson@iet.ntnu.no)
%
% setupmatfile = EDB2readsetup(EDsetupfile,cadgeofile);

% 16 July 2010 Version published in the EDB2toolbox
% 23 April 2015 Cleaned up the check of "calcirs" and "calctfs"
%               Also, added the filesep to the end of Filepath if filesep
%               wasn't included by the user.

global SHOWTEXT
if api == 1
   directsound=infos.directsound;
   elemsize=infos.elemsize;
   nedgesubs=infos.nedgesubs;
   calcpaths=infos.calcpaths;
   calcirs=infos.calcirs;
   calctfs=infos.calctfs;
   skipcorners=infos.skipcorners;
   Rstart=infos.Rstart;
   CADfile=infos.CADfile;
   Filestem=infos.Filestem;
   Filepath=infos.Filepath;
   int_or_ext_model=infos.int_or_ext_model;
   open_or_closed_model=infos.open_or_closed_model;
   EDcalcmethod=infos.EDcalcmethod;
   specorder=infos.specorder;
   difforder=infos.difforder;
   sources=infos.sources;
   receivers=infos.receivers;
   FSAMP=infos.FSAMP;
   CAIR=infos.CAIR;
   RHOAIR=infos.RHOAIR;
   SHOWTEXT=infos.SHOWTEXT;
else
    load('info.mat');
    load('sources_receivers.mat');
end

% if nargin == 1,
% 	cadgeofiletoopen = [];
% end

%[EDsetupfilepath,EDsetupfile] = EDB2strppath(EDsetupfile);
% [EDsetupfilepath,EDsetupfile,fileext] = fileparts(EDsetupfile); 
% EDsetupfilepath = [EDsetupfilepath,filesep];
%%%EDsetupfile = [filestem,fileext];

%eval(['cd ''',EDsetupfilepath,''''])

% comptype = computer;
% if length(comptype) == 3,
%     comptype = [comptype,'P'];
% end
% 
% if comptype(1:4) == 'MACI',
%     [temp1,tempfilestem,temp2] = fileparts(EDsetupfile);
%     eval([tempfilestem])
% else    
%     if comptype(1:3) == 'MAC',  % | comptype(1:3) == 'PCW',
%     %	eval([EDB2strpext(EDsetupfile)])
%         eval([EDsetupfile])
%     else
%     %	eval(['run ''',EDB2strpext(EDsetupfile),''''])
%         eval(['run ''',EDsetupfile,''''])
%     end
% end

if exist('Filestem') ~= 1,
	error(['ERROR: The setup file does not specify the Filestem'])
else
	Varlist = [' Filestem'];
end

%-----------------------------------------------------------------------------------------
% Check if the EDsetupfile specified the needed parameters

if exist('calcirs') ~= 1
	error(['ERROR: calcirs was not specified in ',EDsetupfile,'.m',' Valid values are 0 and 1'])       
else
    if calcirs == 1 & isempty(FSAMP)
        error(['ERROR: You can not set calcirs = 1 without specifying FSAMP'])
    end
end

if exist('calctfs') ~= 1
	error(['ERROR: calctfs was not specified in ',EDsetupfile,'.m',' Valid values are 0 and 1'])       
else
    if calctfs == 1 & isempty(frequencies)
        error(['ERROR: You can not set calctfs = 1 without specifying frequencies'])
    end
end

% if isempty(RHOAIR),
% 	error(['ERROR: RHOAIR was not specified in ',EDsetupfile,'.m'])
% else
% 	Varlist = [Varlist,' RHOAIR'];
% end
% 
% if isempty(CAIR),
% 	error(['ERROR: CAIR was not specified in ',EDsetupfile,'.m'])
% else
% 	Varlist = [Varlist,' CAIR'];
% end

if exist('Filepath') ~= 1,
	error(['ERROR: Filepath was not specified in ',EDsetupfile,'.m'])
else
	Varlist = [Varlist,' Filepath'];
    if Filepath(end) ~= filesep
       disp('WARNING: Filepath did not end with the fileseparator. It was added accordingly');
       Filepath = [Filepath,filesep];
    end
    
end

if exist('Filestem') ~= 1,
	error(['ERROR: Filestem was not specified in ',EDsetupfile,'.m'])
else
	Varlist = [Varlist,' Filestem'];
end

if (exist('CADfile') ~= 1) & isempty(cadgeofiletoopen),
	error(['ERROR: CADfile or cadgeofile was not specified in ',EDsetupfile,'.m'])
else
    if exist('CADfile') == 1
    	Varlist = [Varlist,' CADfile'];	
    else
        if exist('cadgeofiletoopen')
            cadgeofile = cadgeofiletoopen;
        	Varlist = [Varlist,' cadgeofile'];
        end
    end
end

if exist('open_or_closed_model') ~= 1,
	error(['ERROR: open_or_closed_model was not specified in ',EDsetupfile,'.m',' Valid values are ''open'' and ''closed'''])    
else
	Varlist = [Varlist,' open_or_closed_model'];	    
end

if exist('int_or_ext_model') ~= 1,
	error(['ERROR: int_or_ext_model was not specified in ',EDsetupfile,'.m',' Valid values are ''int'' and ''ext'''])    
else
	Varlist = [Varlist,' int_or_ext_model'];	    
end

if lower(open_or_closed_model(1)) == 'o' & lower(int_or_ext_model(1)) == 'i',
    error(['ERROR: When a model is open, it must be an external model. int_or_ext_model was set to ''int'''])
end

if exist('EDcalcmethod') ~= 1,
	error(['ERROR: EDcalcmethod was not specified in ',EDsetupfile,'.m',' Valid values are n,k,v'])
else
	Varlist = [Varlist,' EDcalcmethod'];
end

if exist('directsound') ~= 1,
	error(['ERROR: directsound was not specified in ',EDsetupfile,'.m',' Valid values are 0 and 1'])
else
	Varlist = [Varlist,' directsound'];
end

if exist('specorder') ~= 1,
	error(['ERROR: specorder was not specified in ',EDsetupfile,'.m',' Valid values are integers >= 0'])
else
	Varlist = [Varlist,' specorder'];
end

if exist('difforder') ~= 1,
	error(['ERROR: difforder was not specified in ',EDsetupfile,'.m',' Valid values are integers >= 0'])
else
	Varlist = [Varlist,' difforder'];
end

if specorder < difforder,
    error(['ERROR: specorder must be the same as, or larger than, difforder in this version of EDtoolbox!'])    
end

if difforder >=2 & calcirs == 1
	if exist('elemsize') ~= 1,
		error(['ERROR: elemsize was not specified in ',EDsetupfile,'.m'])
	else
		if min(elemsize) <= 0,
			disp(['WARNING: elemsize was set smaller than zero'])
		else
			if length(elemsize) < difforder,
				error(['ERROR: elemsize must have one value for each diffraction order'])
			else
				Varlist = [Varlist,' elemsize']; 
			end
		end
	end
end

if difforder >2 & calctfs == 1
    error(['ERROR: Maximum diffraction order = 2, for the transfer function calculation'])
end

if exist('sources') ~= 1,
	error(['ERROR: sources were not specified in ',EDsetupfile,'.m'])
else
	Varlist = [Varlist,' sources'];
	[nsou,slask] = size(sources);
end

if exist('receivers') ~= 1,
	error(['ERROR: receivers were not specified in ',EDsetupfile,'.m'])
else
	Varlist = [Varlist,' receivers'];
end

%---------------------------------------------------------------------
% Optional parameters that should get a default value

if isempty(SHOWTEXT),
	disp(['   The parameter SHOWTEXT was not set. It will get the default value 0.'])
	SHOWTEXT = 0;
end
Varlist = [Varlist,' SHOWTEXT'];

if exist('SUPPRESSFILES') ~= 1
	disp(['   The parameter SUPPRESSFILES was not set. It will get the default value 0.'])
	SUPPRESSFILES = 0;
end
Varlist = [Varlist,' SUPPRESSFILES'];


if exist('nedgesubs') ~= 1,
    if SHOWTEXT > 0,
    	disp(['   The parameter nedgesubs was not set. It will get the default value 2.'])
    end
    nedgesubs = 2;
else
    if nedgesubs == 1,
	    nedgesubs = 2;
        if SHOWTEXT > 0,
	        disp('WARNING! nedgesubs was set to 1 in the setup file but it is forced to 2')
        end    
     end
end
Varlist = [Varlist,' nedgesubs'];

if exist('saveindividualdiffirs') ~= 1,
    if SHOWTEXT > 0,
    	disp(['   The parameter saveindividualdiffirs was not set. It will get the default value [0 0].'])
    end
    saveindividualdiffirs = [0 0];
else
    if length(saveindividualdiffirs) == 1
        if SHOWTEXT > 0,
        	disp(['   The parameter saveindividualdiffirs had just one value. It was expanded with a second default value 0.'])
        end        
        saveindividualdiffirs = [saveindividualdiffirs 0];
    end
end
Varlist = [Varlist,' saveindividualdiffirs'];

if exist('firstskipcorner') ~= 1,
    if SHOWTEXT > 0,
    	disp(['   The parameter firstskipcorner was not set. It will get the default value 1000000.'])
    end
    firstskipcorner = 1000000;
end
Varlist = [Varlist,' firstskipcorner'];
	
if exist('doaddsources') ~= 1,
	if nsou > 1,
		disp(['   The parameter doaddsources was not set. It will get the default value 0.'])
	end
	doaddsources = 0;
end
Varlist = [Varlist,' doaddsources'];

if exist('dobackscatter') ~= 1,
    dobackscatter = 0;
end
Varlist = [Varlist,' dobackscatter'];

if dobackscatter*doaddsources == 1
    error(['ERROR: You have set both dobackscatter and doaddsources to 1. This does not make sense!'])
end

if exist('Rstart') ~= 1,
    if SHOWTEXT > 0,
    	disp(['   The parameter Rstart was not set. It will get the default value 0.'])
    end
    Rstart = 0;
end
Varlist = [Varlist,' Rstart'];

if exist('planeseesplanestrategy') ~= 1,
    if SHOWTEXT > 0,
    	disp(['   The parameter planeseesplanestrategy was not set. It will get the default value 0.'])
    end
    planeseesplanestrategy = 0;
end
Varlist = [Varlist,' planeseesplanestrategy'];

%---------------------------------------------------------------------
% Optional parameters 

%if exist('edgeofile')
%	Varlist = [Varlist,' edgeofile'];
%end
if exist('eddatafile')
	Varlist = [Varlist,' eddatafile'];
end
if exist('srdatafile')
	Varlist = [Varlist,' srdatafile'];
end
% % % if exist('edpathsfile')
% % % 	Varlist = [Varlist,' edpathsfile'];
% % % end
if exist('calcpaths')
	Varlist = [Varlist,' calcpaths'];
end
if exist('calcirs')
	Varlist = [Varlist,' calcirs'];
end
if exist('calctfs')
	Varlist = [Varlist,' calctfs'];
end
if exist('calcinteq')
	Varlist = [Varlist,' calcinteq inteq_solmethod inteq_niter inteq_ngauss'];
end

%---------------------------------------------------------------------------------------------
% Find source and receiver coordinates either in the setup-file or in the CAD-file

[nsou,nsoucols] = size(sources);
if nsoucols < 3,
	if isempty(cadgeofiletoopen),
		error(['ERROR: sources should be taken from the CAD file, but no cadgeofile was specified as input parameter to EDB2readsetup'])
	end
	eval(['load ',cadgeofiletoopen])
	[n1,n2] = size(Snames);
	Snames = reshape(Snames.',1,n1*n2);
	strpos = zeros(nsou,1);
	for ii = 1:nsou,
      strpos_source = findstr(sources(ii,:),Snames);
      if ~isempty(strpos_source),
	      strpos(ii) = strpos_source;
      else
			error(['ERROR: One of the sources could not be found in the CAD-file: ',sources(ii,:)])
      end
      
   end
	if prod(strpos) == 0,
		error('ERROR: One of the sources could not be found in the CAD-file')
	end
	sounumbers = floor(strpos/n2)+1;
	sources = Scoords(sounumbers,:);
end

[nrec,nreccols] = size(receivers);
if nreccols ~= 3,
	if isempty(cadgeofiletoopen),
		error(['ERROR: receivers should be taken from the CAD file, but no cadgeofile was specified as input parameter to EDB2readsetup'])
	end
	eval(['load ',cadgeofiletoopen])
	strpos = zeros(nrec,1);
	for ii = 1:nrec,
		strpos_rec = find(Rnumbers==receivers(ii));
		if ~isempty(strpos_rec),
			strpos(ii) = strpos_rec;
		else
			error(['ERROR: One of the receivers could not be found in the CAD-file: ',int2str(receivers(ii))])
		end
	end
	if prod(strpos) == 0,
		error('ERROR: One of the receivers could not be found in the CAD-file')
	end
	receivers = Rcoords(strpos,:);
else
    if nrec == 1,
        if abs(round(receivers(1))) == receivers(1) & abs(round(receivers(2)))==receivers(2) & abs(round(receivers(3))) == receivers(3),
            if SHOWTEXT > 0,
                disp(['WARNING: You have specified the receivers as a row of three integer values, which is ambiguous'])
                disp(['         It is interpreted as one receiver, with given coordinates'])
            end
        end    
    end
end

%---------------------------------------------------------------------------------------------
% Check the optional parameters nSbatches etc

nsources = size(sources,1);
nreceivers = size(receivers,1);

if exist('soulist_for_findpaths') ~= 1,
    soulist_for_findpaths = [1:nsources];
else
    if soulist_for_findpaths(1) < 1,
        error('ERROR: soulist_for_findpaths must start with an integer > 0')    
    end
    if soulist_for_findpaths(end) > nsources,
        error('ERROR: soulist_for_findpaths must end with an integer <= the number of sources')    
    end    
end
if exist('reclist_for_findpaths') ~= 1,
    reclist_for_findpaths = [1:nreceivers];
else
    if reclist_for_findpaths(1) < 1,
        error('ERROR: reclist_for_findpaths must start with an integer > 0')    
    end
    if reclist_for_findpaths(end) > nreceivers,
        error('ERROR: reclist_for_findpaths must end with an integer <= the number of receivers')    
    end    
end
if exist('soulist_for_calcirs') ~= 1,
    soulist_for_calcirs = [1:nsources];
else
    if soulist_for_calcirs(1) < 1,
        error('ERROR: soulist_for_calcirs must start with an integer > 0')    
    end
    if soulist_for_calcirs(end) > nsources,
        error('ERROR: soulist_for_calcirs must end with an integer <= the number of sources')    
    end    
end
if exist('reclist_for_calcirs') ~= 1,
    reclist_for_calcirs = [1:nreceivers];
else
    if reclist_for_calcirs(1) < 1,
        error('ERROR: reclist_for_calcirs must start with an integer > 0')    
    end
    if reclist_for_calcirs(end) > nreceivers,
        error('ERROR: reclist_for_calcirs must end with an integer <= the number of receivers')    
    end    
end
if exist('soulist_for_calctfs') ~= 1,
    soulist_for_calctfs = [1:nsources];
else
    if soulist_for_calctfs(1) < 1,
        error('ERROR: soulist_for_calctfs must start with an integer > 0')    
    end
    if soulist_for_calctfs(end) > nsources,
        error('ERROR: soulist_for_calctfs must end with an integer <= the number of sources')    
    end    
end
if exist('reclist_for_calctfs') ~= 1,
    reclist_for_calctfs = [1:nreceivers];
else
    if reclist_for_calctfs(1) < 1,
        error('ERROR: reclist_for_calctfs must start with an integer > 0')    
    end
    if reclist_for_calctfs(end) > nreceivers,
        error('ERROR: reclist_for_calctfs must end with an integer <= the number of receivers')    
    end    
end

Varlist = [Varlist,' soulist_for_findpaths reclist_for_findpaths soulist_for_calcirs reclist_for_calcirs soulist_for_calctfs reclist_for_calctfs'];

%---------------------------------------------------------------------
% Optional parameters for batch processing

eval(['load ',cadgeofiletoopen])
batchsize = 6000;

suggestednSbatches = ceil(size(planecorners,1)*size(sources,1)/batchsize);
suggestednRbatches = ceil(size(planecorners,1)*size(receivers,1)/batchsize);

if exist('nRbatches') ~= 1,
    if suggestednRbatches > 1
        disp(['WARNING!!! The number of planes and receivers is high enough that it is recommended'])
        disp(['           to set nRbatches to ',int2str(suggestednRbatches)])
        disp(['           Set SHOWTEXT to 0 or 1 to skip the pause'])
        if SHOWTEXT > 1
           disp(['Take the chance to stop the processing and change the value of nRbatches'])  
           pause
        end
        
    else
        nRbatches = 1;
    end
else
    if nRbatches < 0,        
        error('ERROR: The parameter nRbatches must be a positive integer')
    else
        if suggestednRbatches > nRbatches
            disp(['WARNING!!! The number of planes and receivers is high enough that it is recommended'])
            disp(['           to set nRbatches to ',int2str(suggestednRbatches)])
            disp(['           Set SHOWTEXT to 0 or 1 to skip the pause'])
            if SHOWTEXT > 1
                disp(['Take the chance to stop the processing and change the value of nRbatches'])  
                pause
            end        
        else
            nRbatches = round(nRbatches);    
        end
    end
end
if exist('nSbatches') ~= 1,
    if suggestednSbatches > 1
        disp(['WARNING!!! The number of planes and sources is high enough that it is recommended'])
        disp(['           to set nSbatches to ',int2str(suggestednSbatches)])        
        disp(['           Set SHOWTEXT to 0 or 1 to skip the pause'])
        if SHOWTEXT > 1
           disp(['Take the chance to stop the processing and change the value of nSbatches'])  
           pause
        end
    else
        nSbatches = 1;
    end
else
    if nSbatches < 0,        
        error('ERROR: The parameter nSbatches must be a positive integer')
    else
        if suggestednSbatches > nSbatches
            disp(['WARNING!!! The number of planes and sources is high enough that it is recommended'])
            disp(['           to set nSbatches to ',int2str(suggestednSbatches)])
            disp(['           Set SHOWTEXT to 0 or 1 to skip the pause'])
            if SHOWTEXT > 1
                disp(['Take the chance to stop the processing and change the value of nSbatches'])  
                pause
            end        
        else
            nSbatches = round(nSbatches);    
        end
    end
end

Varlist = [Varlist,' nSbatches nRbatches'];

%---------------------------------------------------------------------
% Optional parameters for path editing

if exist('symmetricedgepairs') ~= 1,
    symmetricedgepairs = [];
else
    if size(symmetricedgepairs,2) ~= 2,        
        error('ERROR: The matrix symmetricedgepairs must have two columns')
    end
end

Varlist = [Varlist,' symmetricedgepairs'];

%---------------------------------------------------------------------
% Save in the output file

setupmatfile = [Filepath,Filestem,'_setup.mat'];

eval(['save ',setupmatfile,Varlist])
