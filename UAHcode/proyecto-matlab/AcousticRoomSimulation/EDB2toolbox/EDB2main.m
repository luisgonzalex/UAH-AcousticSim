function EDB2main(EDsetupfile)
% EDB2main - Calculates the specular and edge diffraction IRs and saves them in a file.
% Calculates the specular and edge diffraction IRs and saves them in a file or
% a series of files, if there are many sources and receivers.The
% calculation goes through three stages:
%   1. Geometrical precalculations: define edges and determine visibility
%      between edges, planes, sources and receivers.
%   2. Find the valid sound paths and save a list of them
%   3. Go through the list and construct IRs
%
% This version, EDB2main, includes the option of dividing the
% calculations into batches. When there are many receivers, the geometrical
% precalculations might overload the available memory since all receivers
% are treated simulataneously for each source.
%
% Input parameters:
%	EDsetupfile 	(optional)	This file contains all parameter settings for the calculations.
%					If no filename is given, a file open window will appear and the
%					desired EDsetupfile can be chosen.
%                   The file could either be a .m file (the "EDsetupfile")
%                   or a .mat file (the "setupmatfile").
% Input parameters in the EDsetupfile:
%	FSAMP, CAIR, RHOAIR, SHOWTEXT
%						(global)	The sampling frequency, speed of sound, density of air
%									and the control of how much text should
%									be plotted (SHOWTEXT could be 0,1,2,3,4).
%   SUPPRESSFILES (global, optional) If this variable is introduced in the
%                       setupfile as a global variable, with the value 1, then the files
%                       'xxx_ISEStree.mat' will not be generated. Instead values are passed
%                       directly.
%	Filepath			The path to the folder where all output data files will be stored.
%	Filestem			The first part of the output data files will be given this name
%	CADfile			    The name of the .CAD file where the geometry is given, including path
%   open_or_closed_model      Specifies whether the geometrical model is open or closed:
%                       'o'     Open model
%                       'c'     Closed model
%   int_or_ext_model    Specifies whether the geometrical model is interior or exterior:
%                       'i'     Interior model (NB! An interior model must be closed, 
%                               otherwise it should be called an exterior model)
%                       'e'     Exterior model
%	EDcalcmethod	    Specifies the edge diffraction method. Must have one of the following values:
%						'v'		Vanderkooy
%						'n'		Svensson et al
%						'k'		Kirchhoff approximation NB! The Kirchhoff
%						        approximation gives only first-order diffraction.
%	directsound		0 or 1
%	specorder		The maximum order of specular reflections. (To get all possible combinations,
%                   specorder and difforder should be given the same value).
%	difforder		The maximum order of diffraction components. (difforder can never be higher
%                   than specorder, i.e., if you want difforder = 2, you must set specorder to 2 or higher)
%   elemsize		A vector of size [1,difforder] which specifies the size of the edge subdivision,
%                   for edge diffraction calculations. A value of 0.2 means that the edge elements
%                   are 5 times (1/0.2 =5) the wavelength at FSAMP.
%				    A higher value gives more accurate edge diffraction but is also more time consuming. 
%                   Larger and larger elements could be employed for higher-order diffraction.
%                   Suggested values for high accuracy: elemsize = [1 0.5 0.25 0.125 ....].
%   nedgesubs       All edges are subdivided very coarsely into edge segments, for visibility checks.
%                   The value of nedgesubs states how many such segments each edge will be divided into.
%                   The minimum value is 2. NB! All edges will be divided into equally many segments
%                   in order to employ efficient Matlab vector-based processing.
%	sources			A vector of the source coordinates, size [nsources,3], or source names/numbers,
%                   size [nsources,1] or [nsources,2] (if they should be taken from the CADfile.)
%	receivers		A vector of the receiver coordinates or receiver numbers (if they should be taken
%					from the CADfile. Same syntax as for sources.
%   nSbatches,nRbatches
%                   (optional) An integer >= 1 which defines how many chunks the number of sources/receivers
%                   will be divided into. Default = 1. nSbatches/nRbatches = 1 means that all sources/receivers are 
%                   treated simultaneously. nSbatches/nRbatches = 2 means that the list of sources/receivers
%                   is divided into two chunks etc. It is fastest with nSbatches/nRbatches = 1 but it also
%                   requires the largest amount of memory. (For large problems, there might even be no time
%                   penalty in choosing nRbatches = the number of receivers!)
%   soulist_for_findpaths/reclist_for_findpaths
%                   (optional) A list of integers that specifies which of the sources/receivers that
%                   should be run for the "find paths" stage. Default = [1:nsources]/[1:nreceivers].
%                   This is useful for huge problems where one might want to stop (with CTRL-C) and restart calculations.
%   soulist_for_calcirs/reclist_for_calcirs
%                   (optional) A list of integers that specifies which of the sources/receivers that
%                   should be run for the "calc irs" stage. Default = [1:nsources]/[1:nreceivers].
%                   This is useful for huge problems where one might want to stop (with CTRL-C) and restart calculations.
%	firstskipcorner (optional) All edges including at least one corner with this number or
%             		higher will be excluded from the calculations. Default: 1e6
%	Rstart 			(optional) The distance that should correspond to the time zero for
%             		all IRs. Default = 0. NB! With Rstart = 0, all IRs will include the initial time
%                   delay that is caused by the propagation time from the source to the receiver. For
%             		long distances and high FSAMP, this can lead to many initial samples that are zero.
%                   If you know that all source-receiver combinations have a certain minimum
%             		distance, this minimum distance (or a value slightly below) can be given as Rstart.
%                   If Rstart is given a value which is smaller than the shortest distance encountered, 
%                   an error will occur.
%
% Optional input parameters in the EDsetupfile:
%	cadgeofile	    If this parameter is specified, this file will be read rather than the CADfile.
%                   Saves some time for big geometries.
%	eddatafile	    If this parameter is specified, the function EDB2edgeo is not called. Saves
%					some time.
%	srdatafile	    If this parameter is specified, the function EDB2srgeo is not called. Saves
%					some time.
%   restartnumber   If this parameter is specified, the calculations will start at receiver number
%                   restartnumber, rather than number 1. This is useful if a calculation for many
%                   receivers can not be completed in one run.
%   stopnumber      If this parameter is specified, the calculations will stop for receiver number
%                   stopnumber.
%   saveindividualdiffirs   A vector with two values, 0 or 1.
%					If the first value is given the value 1, the diffraction
%                   IR will have one column for each order. Default is that
%                   all orders are summed into a single irdiff.
%					If the second value is given the value 1, all individual diffraction irs
%					will be saved on individual lines in a large matrix called 'alldiffirs'.
%				 	This matrix has only single precision.
%					Default value: [0 0]
%   doaddsources    0 or 1, (default value 0) indicating whether:
%                       0: every individual S-to-R IR should be saved
%                       1: a single IR is saved for each R, which is
%                       constructed by first computing all individual
%                       S-to-R IRs, saving them in a temp file, and then
%                       add IRs from all sources together. This way a
%                       loudspeaker element can be simulated by a
%                       distribution of point sources.
%   dobackscatter   0 or 1, (default value 0), indicating whether:
%                       0: every individual S-to-R IR should be saved
%                       1: only S1 to R1, S2 to R2, S3 to R3 combos
%                       are calculated.
%
% Output parameters:
% 	The following IRs are calculated and saved in the output file(s). Any of them can be
% 	excluded from the calculations by setting the parameters in the setup file.
%
%	irtot    		the total impulse response
%	irdirect 		the direct sound
%  	irgeom   		the purely geometrical acoustics components, i.e. specular reflections.
%  	irdiff   		the diffracted components. This will innclude all combinations of specular and  
%                   diffraction, and first-order as well as higher-order combinations.
%
% Note that these variables, together with those in the EDsetupfile are saved in a results file
%	which gets the name in Filestem (see below) + '_SS_RR_ir.mat' (if addsources = 0) or
%	'_RR_ir.mat' (if addsources = 1).
%
% Note also that the impulse responses are scaled such that the direct sound gets an amplitude of
% 1/distance. This means that the impulse response represents a system with the sound pressure as
% output signal and an input signal = RHOAIR*Volumeacc(t)/(4*pi) where Volumeacc(t) is
% the volume acceleration of the monopole source.
%
% A logfile containing the computation time for the three stages is
% created, with the name Filestem + '_log.mat'
%
% Uses functions EDB2readcad EDB2readsetup EDB2edgeo
% EDB2Sgeo EDB2Rgeo EDB2findISEStree EDB2findpathsISES EDB2makeirs
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
% Peter Svensson (svensson@iet.ntnu.no) 23 April 2015
%
% EDB2main(EDsetupfile);

% 16 August 2013 Released version in the EDB2toolbox
% 23 April 2015  Bug found and fixed: EDB2main erroneously called
%                EDB2makeirsISESx, which was replaced by a call to EDB2makeirs

global FSAMP CAIR RHOAIR SHOWTEXT JJ JJnumbofchars SUPPRESSFILES
global POTENTIALISES ISCOORDS IVNDIFFMATRIX
global IVNSPECMATRIX ORIGINSFROM ISESVISIBILITY REFLORDER

global IRDIFFVEC

FSAMP = []; CAIR = []; RHOAIR = [];  SHOWTEXT = [];

%--------------------------------------------------------------
% First run the setup m file, just to get access to the CAD-file name
% so that the CAD-file can be converted into a cadgeofile.
% (This conversion must be done before the setup m-file is converted
% to a setup mat-file since the source and receiver coordinates should
% be specified in the setupmatfile, and the sources and receivers might
% have to be read from the CAD-file).
%
% First we make a little check of a common mistake: check if the paths
% seem to be for the correct computer type.
% The EDsetupfile-path does not need to be checked if a file window is
% presented. The CAD-file path always needs to be checked.

% We check if the type of setupfile:
%   .m    setupfiletype = 1;
%   .mat  setupfiletype = 2;

setupfiletype = 1;

if nargin == 0,
	[EDsetupfile,EDsetupfilepath] = uigetfile('*.m','Please select the EDsetupfile');
    [temp1,EDsetupfile,temp2] = fileparts([EDsetupfilepath,EDsetupfile]);

    EDsetupfile = [EDsetupfilepath,EDsetupfile];
    
	if ~isstr(EDsetupfile) | isempty(EDsetupfile),
		return
	end
else
    [EDsetupfilepath,tempfilestem,fileext] = fileparts(EDsetupfile);
    EDsetupfilepath = [EDsetupfilepath,filesep];
    if length(fileext) == 4
        if fileext == '.mat'
            setupfiletype = 2; 
        end
    end
end

eval(['cd ''',EDsetupfilepath,''''])

if setupfiletype == 1
    comptype = computer;
    if length(comptype) == 3,
        comptype = [comptype,'P'];
    end
    if comptype(1:4) == 'MACI', 
        [temp1,tempfilestem,temp2] = fileparts(EDsetupfile);
        eval([tempfilestem])
    else    
        eval(['run ''',EDsetupfile,''''])
    end
else
    eval(['load ''',EDsetupfile,''''])
end

if SHOWTEXT >= 1,
	disp(' ')
	disp('####################################################################')
	disp('#')
	disp('#  EDB2main. Version B2, 16 August 2013')
	disp('#')
end

if SHOWTEXT >= 2,
	disp('#')
	disp('# This version does the following:')
	disp('# *  Caclulates impulse responses that include specular reflections up to any order and') 
    disp('#    multiple diffraction up to order 6. Combinations are also taken into account, such as')
    disp('#    specular-diffraction-specular but higher-order diffractions must be')
    disp('#    in a single sequence, possibly with specular reflections before and after.')
    disp('#    There will be memory problems for higher-order specular reflections due to a')
    disp('#    memory-hungry vectorized version of the image source method.')
    disp('# *  Combinations with two diffractions that have specular reflections inbetween')
    disp('#    are handled, but not with any specular reflections after the last diffraction.')
	disp('# *  Partly visible edges are handled in first-order diffraction')
	disp('#    calculations. The accuracy depends on how high value of')
	disp('#    nedgesubs that is chosen. Edges that are split into two or more subedges are handled properly.')
    disp('# *  Specular reflections in the IR are handled by splitting a pulse between two consecutive time slots.')
    disp('#    This gives quite a severe low-pass filter roll-off, but zero phase-error.')
    disp('#    Choose a high over-sampling ratio to minimize this magnitude error.')
    disp('#    In addition, the low-pass filter effect can be avoided altogether for the direct sound by setting')
    disp('#    Rstart to exactly the source-receiver distance. Thereby the direct sound pulse will occur exactly in the')
    disp('#    first time sample. For more tweaking, the three values Rstart and CAIR and FSAMP should be possible')
    disp('#    to adjust so that one, two or three geometrical pulses arrive exactly at sample times.')
    disp('# *  An accurate integration technique is used for first-order diffraction so that receiver positions close')
    disp('#    to zone boundaries are handled without problem. Exactly at the zone boundary the specular reflection')
    disp('#    or direct sound gets half the amplitude (See below though).')
    disp('# *  Frequency-domain diffraction is implemented for up to second-order diffraction components.')
    disp('#    Matlab''s quadgk integration routine is used. Receivers near zone boundaries are not handled well in this version. ')    
	disp('#')
	disp('# Limitations with this version:')
	disp('# 1. Some, but not all, possible combinations where two or more diffracted')
	disp('#    reflections have specular reflections in-between are calculated. ')
    disp('#    This requires a bit of work to complete.')
	disp('# 2. In-plane visibility of edges is only partly checked, i.e., if a plane has indents')
	disp('#    so that an edge can not see all other edges of its own plane.')
	disp('# 3. Part-visibility of higher-order diffraction is not checked. These components are treated')
    disp('#    as fully visible or totally invisible. However, if the first or last edge in the')
    disp('#    sequence is partly visible, this part-visibility will be used.')
	disp('# 4. Non-rigid specular reflections or edge diffraction for non-rigid planes are not implemented. ')
	disp('#    Pressure-release surfaces would be quite easy to implement.')
    disp('#    Other surface impedances could be implemented since the hitpoints for all specular')
    disp('#    reflections are stored in the edpaths output file.')
    disp('# 5. Other receivers than omnidirectional microphones are not implemented. HRTFs or directional microphones')
    disp('#    could be implemented since the hitpoints for all specular reflections are stored in the edpaths')
    disp('#    so the incidence angles are easy to calculate. Reflections that end with a diffraction need some')
    disp('#    more consideration.')
    disp('# 6. The special cases where the receiver is exactly at a zone boundary is handled correctly only')
    disp('#    for first-order diffraction, and only for IRs. This takes a bit of work to implement for higher-order diffraction.')
    disp('# 7. Sources and receivers can not be place directly at a surface. They must be raised by 0.1 mm or so.')
	disp('#')
	disp('# This version could be made more efficient by:')
	disp('# 1. Employing ray-tracing or beam-tracing for finding the visibility of')
	disp('#    plane to plane etc. Now, a simple plane-to-plane visibility check')
	disp('#    is performed: are planes in front of each other? The final visibility')
	disp('#    check is done for each image source to each receiver. An improved method would be more efficient')
	disp('#    for geometries with lots of indents and obscuring planes, e.g., as in a city')
	disp('#    geometry.')
	disp('# 2. Developing an algorithm that calculates the diffraction IR in an iterative process')
	disp('#    instead of restarting for each new order of diffraction. This will be implemented in a coming')
	disp('#    version of the toolbox which is based on the integral equation formulation.')
    disp('# 3. The direct sound visibility test could be done earlier (in EDB2SRgeo).')
    disp('# 4. Calculation of Image Receiver coordinates could probably be made more efficiently.')
    disp('#')
end

if SHOWTEXT >= 1,
    disp(' ')
    disp(['   Using the setupfile: ',EDsetupfile])
    disp(' ')
end

% #########################################################################
% #########################################################################
%   
%      Geometry pre-calculations
%   
% #########################################################################

%---------------------------------------------------------------------------------------------
% Convert the CAD-file into a mat file - or read the cadgeofile if it is specified

t00 = clock;

if SHOWTEXT >= 1,
	disp(' ')
	disp('####################################################################')
	disp('#')
	disp('#  Geometry pre-calculations')
	disp(' ')
end

if exist('cadgeofile') ~= 1,
    desiredname = [Filepath,Filestem,'_cadgeo'];
    if SHOWTEXT >= 2,
		disp(['   Reading CAD-file: ',CADfile,' and converting it to a cadgeofile: ',desiredname])
    end
    cadgeofile = EDB2readcad(CADfile,desiredname)
else
	if SHOWTEXT >= 2,
		disp(['   Using existing cadgeofile: ',cadgeofile])
	end
end

%--------------------------------------------------------------
% Convert the EDsetupfile into a setupmatfile. If sources and receivers should
% be read in the CADfile, add these to the setupmatfile.
%
% This is needed only if the setupfile was of type 1, that is, a .m file.

if setupfiletype == 1
    setupmatfile = EDB2readsetup(EDsetupfile,cadgeofile);

    eval(['load ',setupmatfile])
end

tcadgeo = etime(clock,t00);

%--------------------------------------------------------------
% Derive the edge parameters.

t00 = clock;

if calcpaths == 1,
    
    if exist('eddatafile') ~= 1,
	    desiredfile = [Filepath,Filestem,'_eddata'];
		if SHOWTEXT >= 2,
			disp(['   Calculating edges and creating an eddatafile: ',desiredfile])
		end
		eddatafile = EDB2edgeox(cadgeofile,desiredfile,lower(open_or_closed_model(1)),lower(int_or_ext_model(1)),...
            specorder,difforder,nedgesubs,firstskipcorner,planeseesplanestrategy);
	else
		if SHOWTEXT >= 2,
			disp(['   Using existing eddatafile: ',eddatafile])
		end
	end
else
	if exist('eddatafile') ~= 1,
	    desiredfile = [Filepath,Filestem,'_eddata'];
    	eddatafile = desiredfile;
	else
		if SHOWTEXT >= 2,
			disp(['   Using existing eddatafile: ',eddatafile])
		end
	end
end

tedgeo = etime(clock,t00);

%--------------------------------------------------------------
% If batches should be used, divide the list of receivers into batches.

nsources = size(sources,1);
nreceivers = size(receivers,1);

if nSbatches == 1,
    doSbatches = 0;    
elseif nSbatches > 1,
    if nSbatches > nsources,
        error('ERROR: nSbatches can not be larger than the number of sources')    
    end
    doSbatches = 1;
    Sbatchlist = zeros(nSbatches,2);
    nperbatch = ceil(nsources/nSbatches);
    Sbatchlist(:,2) = [1:nSbatches].'*nperbatch;
    Sbatchlist(:,1) = [0:nSbatches-1].'*nperbatch+1;
    Sbatchlist(nSbatches,2) = nsources;
    reftoSbatchlist = [1:nsources].';
    reftoSbatchlist = reftoSbatchlist(:,ones(1,nSbatches));
    compvec = Sbatchlist(:,1).';
    reftoSbatchlist = reftoSbatchlist >= compvec(ones(nsources,1),:);
    reftoSbatchlist = sum(reftoSbatchlist.').';
end

if nRbatches == 1,
    doRbatches = 0;    
elseif nRbatches > 1,
    if nRbatches > nreceivers,
        error('ERROR: nRbatches can not be larger than the number of receivers')    
    end
    doRbatches = 1;
    Rbatchlist = zeros(nRbatches,2);
    nperbatch = ceil(nreceivers/nRbatches);
    Rbatchlist(:,2) = [1:nRbatches].'*nperbatch;
    Rbatchlist(:,1) = [0:nRbatches-1].'*nperbatch+1;
    Rbatchlist(nRbatches,2) = nreceivers;
    reftoRbatchlist = [1:nreceivers].';
    reftoRbatchlist = reftoRbatchlist(:,ones(1,nRbatches));
    compvec = Rbatchlist(:,1).';
    reftoRbatchlist = reftoRbatchlist >= compvec(ones(nreceivers,1),:);
    reftoRbatchlist = sum(reftoRbatchlist.').';
end

%---------------------------------------------------------------------------------------------
% Find out the S and R related parameters, make obstruction checks etc by calling 
% EDB2Sgeo/EDB2Rgeo or load existing sdata/rdata files.

t00 = clock;

if calcpaths == 1,
	if exist('sdatafile') ~= 1,
        if doSbatches == 0,
		    desiredname = [Filepath,Filestem,'_sdata'];
			if SHOWTEXT >= 2,
				disp(['   Calculating S parameters and creating an sdatafile: ',desiredname])
			end

            sdatafile = EDB2SorRgeo(eddatafile,desiredname,sources,'S',difforder,nedgesubs);
        else
            error(['ERROR: S-batches not implemented yet'])
            for ii = 1:nSbatches,
			    desiredname = [Filepath,Filestem,'_sdata'];
				if SHOWTEXT >= 2,
					disp(['   Calculating S parameters and creating an sdatafile: ',desiredname,'_',int2str(ii)])
				end
        		sdatafile = EDB2SorRgeo(eddatafile,[desiredname,'_',int2str(ii)],sources(Sbatchlist(ii,1):Sbatchlist(ii,2),:),'S',difforder,nedgesubs);
                sdatafile = desiredname;   % To trim off the _1 or _2 or...    
            end
        end
	else
		if SHOWTEXT >= 2,
			disp(['   Using existing sdatafile: ',sdatafile])
		end
	end
else
	if exist('sdatafile') ~= 1,
	    desiredname = [Filepath,Filestem,'_sdata'];
    	sdatafile = desiredname;
	else
		if SHOWTEXT >= 2,
			disp(['   Using existing sdatafile: ',sdatafile])
		end
	end
end

tsgeo = etime(clock,t00);

%------------------------------------
% We check that there is not a small geometrical mistake which makes the
% source be behind all planes - but only if the user has set SHOWTEXT.

if SHOWTEXT >= 2,
    eval(['load ',sdatafile])
    if size(visplanesfroms,2) == 1,
        if ~any(visplanesfroms),
            disp(' ')
            disp('WARNING!!! The source(s) can not see a single plane!')
            disp('           Check if this is correct! ')
            disp('           A common cause is that the source is placed very close to a plane,')
            disp('           but behind it. (CR => continue calculations)')
            disp(' ')
            pause
        end                
    else
        iv = find(sum(visplanesfroms) == 0);
        if ~isempty(iv),
            disp(' ')
            disp('WARNING!!! Some of the source(s) can not see a single plane!')
            disp('           These are the sources number:')
           disp(['           ',int2str(iv(:).')])
            disp('           Check if this is correct! ')
            disp('           A common cause is that the source is placed very close to a plane,')
            disp('           but behind it. (CR => continue calculations)')
            disp(' ')
        end
    end
end

t00 = clock;

if calcpaths == 1,
	if exist('rdatafile') ~= 1,
        if doRbatches == 0,
		    desiredname = [Filepath,Filestem,'_rdata'];
			if SHOWTEXT >= 2,
				disp(['   Calculating R parameters and creating an rdatafile: ',desiredname])
			end
    		rdatafile = EDB2SorRgeo(eddatafile,desiredname,receivers,'R',difforder,nedgesubs);
        else
            for ii = 1:nRbatches,
			    desiredname = [Filepath,Filestem,'_rdata'];
				if SHOWTEXT >= 2,
					disp(['   Calculating R parameters and creating an rdatafile: ',desiredname,'_',int2str(ii)])
				end
                
                rdatafile = EDB2SorRgeo(eddatafile,[desiredname,'_',int2str(ii)],receivers(Rbatchlist(ii,1):Rbatchlist(ii,2),:),'R',difforder,nedgesubs);
                rdatafile = desiredname;   % To trim off the _1 or _2 or...    
            end
        end
	else
		if SHOWTEXT >= 2,
			disp(['   Using existing rdatafile: ',rdatafile])
		end
	end
else
	if exist('rdatafile') ~= 1,
	    desiredname = [Filepath,Filestem,'_rdata'];
    	rdatafile = desiredname;
	else
		if SHOWTEXT >= 2,
			disp(['   Using existing rdatafile: ',rdatafile])
		end
	end
end

trgeo = etime(clock,t00);

%------------------------------------
% We check that there is not a small geometrical mistake which makes the
% receiver(s) be behind all planes - but only if the user has set SHOWTEXT.

if SHOWTEXT >= 2,
    eval(['load ',rdatafile])
    if size(visplanesfromr,2) == 1,
        if ~any(visplanesfromr),
            disp(' ')
            disp('WARNING!!! The receiver can not see a single plane!')
            disp('           Check if this is correct! ')
            disp('           A common cause is that the receiver is placed very close to a plane,')
            disp('           but behind it. (CR => continue calculations)')
            disp(' ')
            pause
        end                
    else
        iv = find(sum(visplanesfromr) == 0);
        if ~isempty(iv),
            disp(' ')
            disp('WARNING!!! Some of the receiver(s) can not see a single plane!')
            disp('           These are the receiver numbers:')
           disp(['           ',int2str(iv(:).')])
            disp('           Check if this is correct! ')
            disp('           A common cause is that a receiver is placed very close to a plane,')
            disp('           but behind it. (CR => continue calculations)')
            disp(' ')
        end
    end
end

%---------------------------------------------------------------------------------------------
% The edgeseesedge test is run separately, for the cases where difforder >= 2. 
% The data is stored in the eddata file again.

if difforder >= 2 & calcpaths == 1,
% if difforder >= 1 & calcpaths == 1,
    if exist('eddata2file') ~= 1,
	    desiredfile = [Filepath,Filestem,'_eddata'];
    	if SHOWTEXT >= 2,
		    disp(['   Adding edge-to-edge visibility to the eddatafile: ',eddatafile])
    	end
        if exist('ndiff2batches') ~= 1,
            ndiff2batches = 1;    
        end
    	eddatafile = EDB2ed2geox(eddatafile,sdatafile,rdatafile,doSbatches,doRbatches,specorder,difforder,nedgesubs,ndiff2batches);    

	else
		if SHOWTEXT >= 2,
			disp(['   Using existing eddata2file: ',eddata2file])
		end
	end
    
end

%---------------------------------------------------------------------------------------------
% Make a big loop running through the source and receiver coordinates.

t00 = clock;

eval(['load ',eddatafile])
clear cornerinfrontofplane planeseesplane edgeseesplane

[nplanes,slask]      = size(planecorners);
[nedges,slask]      = size(edgecorners);

if exist('calcpaths')~=1,
	calcpaths = 1;
end
if exist('calcirs')~=1,
	calcirs = 1;
end

% Display receiver numbers with an interval that depends on the 
% value of SHOWTEXT.

if SHOWTEXT == 1,
    idisp = ceil(nreceivers/20);
elseif SHOWTEXT > 1,
    idisp = 1;
end

% We construct the string version of counters and let them be global
% variables so we don't have to call int2str all the time

ncountersneeded = max(max([specorder difforder nplanes nedges]));
if ncountersneeded < 10,
    JJ = setstr(32*ones(ncountersneeded,1));
else
    JJ = setstr(32*ones(ncountersneeded,2));
end
for jj=1:ncountersneeded,
   jjstr = int2str(jj);
   JJ(jj,1:length(jjstr)) = jjstr;
end
[n1,n2] = size(JJ);
JJnumbofchars = ones(n1,1);
if n1>9,
    JJnumbofchars(10:n1) = JJnumbofchars(10:n1)+1;
    if n1 > 99,
        JJnumbofchars(100:n1) = JJnumbofchars(100:n1)+1;    
        if n1 > 999,
            JJnumbofchars(1000:n1) = JJnumbofchars(1000:n1)+1;    
        end
    end
end

NSOU = int2str(nsources);
NREC = int2str(nreceivers);

tsetup = etime(clock,t00);

% #########################################################################
% #########################################################################
%   
%      Find the paths
%   
% #########################################################################

t00 = clock;
tISEStree = zeros(nsources,1);

if calcpaths == 1,

    if SHOWTEXT >= 1,
        disp(' ')
        disp('####################################################################')
        disp('#')
        disp('#  Finding the paths')
        disp(' ')
    end

    % Some variables that are needed for empty combinations

	pathtypevec = [];
	reflpaths = [];
	pathlengthvec = [];
	specextradata = [];
	edgeextradata = [];
	mainlistguide = [];
    Sinsideplanenumber = [];
    Rinsideplanenumber = [];
    mainlistguide = [];
    mainlistguidepattern = [];
    directsoundrow = [];
    allspecrows = [];
    firstdiffrow = [];
	Varlist = [' pathtypevec reflpaths specextradata edgeextradata S R mainlistguide'];		
    Varlist = [Varlist,' Sinsideplanenumber Rinsideplanenumber mainlistguide mainlistguidepattern directsoundrow allspecrows firstdiffrow'];

	if doSbatches == 0,
        eval(['load ',sdatafile])   
	else
        latestsdatafile = reftoSbatchlist(1);
        [sdatafilepath,filestem,temp1] = fileparts(sdatafile);
        eval(['load ',[sdatafilepath,filesep,filestem],'_',int2str(latestsdatafile),'.mat'])     
	end
% % % 	% We create a souinsideplane matrix from the visplanesfroms matrix    
% % % 	% for use in EDB2makeirs
	
	if doRbatches == 0,
        eval(['load ',rdatafile])
	else
        latestrdatafile = reftoRbatchlist(1);
        [rdatafilepath,filestem,temp1] = fileparts(rdatafile);
        eval(['load ',[rdatafilepath,filesep,filestem],'_',int2str(latestrdatafile),'.mat'])     
	end
% % % 	% We create a recinsideplane matrix from the visplanesfromr matrix    
% % % 	% for use in EDB2makeirs
    
    for isou = soulist_for_findpaths,
    
        t00_sou = clock;
        ISOU = int2str(isou);
        if SHOWTEXT >= 1,
			disp(['Calculating for source ',ISOU,' (of ',NSOU,') '])
		end
        if doSbatches == 1,
            if latestsdatafile ~= reftoSbatchlist(isou),
                latestsdatafile = reftoSbatchlist(isou);
                [sdatafilepath,filestem,temp1] = fileparts(sdatafile);
                eval(['load ',[sdatafilepath,filesep,filestem],'_',int2str(latestsdatafile),'.mat'])     
            end
            Scolnumber = isou-Sbatchlist(reftoSbatchlist(isou),1)+1
        else
            Scolnumber = isou;
        end
        S = sources(Scolnumber,:);
        %%%Sdirection = [1 0 0];
	
        %   ###########################################################
        %   #
        %   #   Calculate the ISES tree for each source
        %   #   IS = image sources, for specular reflections
        %   #   ES = edge sources, for edge diffraction
        %   #
        %   ###########################################################

	    if exist('ISEStreefile') ~= 1,
            if SHOWTEXT >= 2,
                disp(['   Calculating ISES tree'])    
            end
            if difforder >= 1,
                usedISEStreefile = EDB2findISEStree(eddatafile,S,isou,specorder,difforder,visplanesfroms(:,isou),vispartedgesfroms(:,isou),nedgesubs);
            else
                usedISEStreefile = EDB2findISEStree(eddatafile,S,isou,specorder,difforder,visplanesfroms(:,isou),[],nedgesubs);
            end
            
        else
            [ISEStreefilepath,filestem,temp1] = fileparts(ISEStreefile);
            usedISEStreefile = [[ISEStreefilepath,filesep,filestem],'_',ISOU,'_ISEStree.mat;'];   
            if SHOWTEXT >= 2,
                disp(['   Using existing ISEStreefile: ',ISEStreefile,'_',ISOU,'_ISEStree.mat'])    
            end
        end
        if SUPPRESSFILES == 0
            eval(['load ',usedISEStreefile])
        else
            POTENTIALISES = usedISEStreefile.POTENTIALISES;
            ORIGINSFROM = usedISEStreefile.ORIGINSFROM;
            ISCOORDS = usedISEStreefile.ISCOORDS;
            ISESVISIBILITY = usedISEStreefile.ISESVISIBILITY;
            IVNSPECMATRIX = usedISEStreefile. IVNSPECMATRIX;
            lengthNspecmatrix = usedISEStreefile.lengthNspecmatrix;
            IVNDIFFMATRIX = usedISEStreefile.IVNDIFFMATRIX;
            lengthNdiffmatrix = usedISEStreefile.lengthNdiffmatrix;
            singlediffcol = usedISEStreefile.singlediffcol;
            REFLORDER = usedISEStreefile.REFLORDER;
            startindicessinglediff = usedISEStreefile.startindicessinglediff;
            endindicessinglediff = usedISEStreefile.endindicessinglediff;
            ndecimaldivider = usedISEStreefile.ndecimaldivider;
            PointertoIRcombs = usedISEStreefile.PointertoIRcombs;
            IRoriginsfrom = usedISEStreefile.IRoriginsfrom;
        end
        tISEStree(isou) = etime(clock,t00_sou);

        if dobackscatter ~= 1, 
            reclist_in_forloop = reclist_for_findpaths;
        else
            reclist_in_forloop = isou;            
        end

        for irec = reclist_in_forloop,

            %   ###########################################################
            %   #
            %   #   Find the valid paths for each source-receiver combination
            %   #   _edpaths files are created
            %   #
            %   ###########################################################

            IREC = int2str(irec);
        	if SHOWTEXT >= 1,
                if round(irec/idisp)*idisp==irec,  
			        disp(['...receiver no. ',IREC,' (of ',NREC,')'])
                end    
            end
            if doRbatches == 1,
                if latestrdatafile ~= reftoRbatchlist(irec),
                    latestrdatafile = reftoRbatchlist(irec);
                    [rdatafilepath,filestem,temp1] = fileparts(rdatafile);
                    eval(['load ',[rdatafilepath,filesep,filestem],'_',int2str(latestrdatafile),'.mat'])     
                end
                Rcolnumber = irec-Rbatchlist(reftoRbatchlist(irec),1)+1;
            else
                Rcolnumber = irec;
            end
            R = receivers(Rcolnumber,:);

            if exist('edpathsfile') ~= 1,
			    desirededpathsfile = [Filepath,Filestem,'_',ISOU,'_',IREC,'_edpaths.mat'];
                if soutsidemodel(Scolnumber) == 0 & routsidemodel(Rcolnumber) == 0,        
                    if difforder > 0,
                        usededpathsfile = EDB2findpathsISESx(eddatafile,lengthNspecmatrix,lengthNdiffmatrix,singlediffcol,startindicessinglediff,...
                            endindicessinglediff,ndecimaldivider,PointertoIRcombs,IRoriginsfrom,...
                            S,R,Scolnumber,Rcolnumber,directsound,specorder,difforder,...
                            nedgesubs,visplanesfroms(:,Scolnumber),visplanesfromr(:,Rcolnumber),...
                            vispartedgesfroms(:,Scolnumber),vispartedgesfromr(:,Rcolnumber),desirededpathsfile);
                    else

                        % Should it really be Scolnumber/Rcolnumber
                        % below??? They will always have the value 1 if
                        % nRbatches = nreceivers??
                        
                        usededpathsfile = EDB2findpathsISESx(eddatafile,lengthNspecmatrix,[],[],[],...
                            [],[],[],[],S,R,Scolnumber,Rcolnumber,directsound,specorder,difforder,...
                            nedgesubs,visplanesfroms(:,Scolnumber),visplanesfromr(:,Rcolnumber),[],[],desirededpathsfile);
                    end
                else
                    usededpathsfile = desirededpathsfile;
                    eval(['save ',usededpathsfile,Varlist])
                end
			else
				if SHOWTEXT >= 2,
					disp(['   Using existing edpathsfile: ',edpathsfile,'_',ISOU,'_',IREC,'_edpaths.mat'])
                    desirededpathsfile = [edpathsfile,'_',ISOU,'_',IREC,'_edpaths.mat'];
				end
			end
            
	
		end		% of the for isou = ,
	
	end		% of the for irec = ,
end

tfindpaths = etime(clock,t00);

% #########################################################################
% #########################################################################
%   
%      Edit paths
%   
% #########################################################################

t00 = clock;

if ~isempty(symmetricedgepairs) & exist('edpathsfile') == 1,

    disp(['WARNING: You specified symmetric edge pairs, but path pruning is not'])
    disp(['         allowed when you have specified an edpaths file'])
    
end
if ~isempty(symmetricedgepairs) & exist('edpathsfile') ~= 1,

    if SHOWTEXT >= 1,
        disp(' ')
        disp('####################################################################')
        disp('#')
        disp('#  Pruning paths: removing redundant (symmetric) diffraction combinations')
        disp(' ')
    end
        
	for isou = soulist_for_calcirs,
        ISOU = int2str(isou);
        if SHOWTEXT >= 1,
			disp(['Pruning paths for source ',ISOU,' (of ',NSOU,') '])
		end

        if dobackscatter ~= 1, 
            reclist_in_forloop = reclist_for_calcirs;
        else
            reclist_in_forloop = isou;            
        end

        for irec = reclist_in_forloop,
            IREC = int2str(irec);
        	if SHOWTEXT >= 1,
                if round(irec/idisp)*idisp==irec,  
			        disp(['...receiver ',IREC,' (of ',NREC,')'])
                end    
            end

			edpathsfiletoedit = [Filepath,Filestem,'_',ISOU,'_',IREC,'_edpaths.mat'];

            EDB2editpathlist(edpathsfiletoedit,[],symmetricedgepairs,[]);    
    
        end
    end
end

tprunepaths = etime(clock,t00);

% #########################################################################
% #########################################################################
%   
%      Construct IRs
%   
% #########################################################################

t00 = clock;

approxplanemidpoints = 0.5*(maxvals+minvals);

if calcirs == 1,

    if SHOWTEXT >= 1,
        disp(' ')
        disp('####################################################################')
        disp('#')
        disp('#  Constructing IRs')
        disp(' ')
    end
    
% % %     % We create a souinsideplane matrix from the visplanesfroms matrix
% % %     
% % %     souinsideplane = (visplanesfroms==4);
% % %     recinsideplane = (visplanesfromr==4);

    % If the parameter guiderowstouse wasn't specified in the setup file,
    % we give it the default value.
    
    if exist('guiderowstouse') ~= 1,
        guiderowstouse = [];    
    end

    % Some variables that are needed for empty combinations
    irgeom = [];
    irtot = [];
    irdirect = [];
    irdiff = [];
    Varlist = [' irdirect irgeom irdiff irtot FSAMP Rstart CAIR elemsize nedgesubs'];
    Varlistaccum = [' irdirectaccum irgeomaccum irdiffaccum irtotaccum FSAMP Rstart CAIR elemsize nedgesubs'];

    reclist_in_forloop = reclist_for_calcirs;

	for irec = reclist_in_forloop,
		IREC = int2str(irec);
		if SHOWTEXT >= 1,
			disp(['...receiver ',IREC,' (of ',NREC,')'])
		end
		
		irdirectaccum = [];
        irgeomaccum = [];
        irtotaccum = [];
        irdiffaccum = [];

		if dobackscatter == 1, 
			soulist_for_calcirs = irec;
		end
		
		for isou = soulist_for_calcirs,
    	    ISOU = int2str(isou);
        	if SHOWTEXT >= 1,
                if round(isou/idisp)*idisp==isou,  
					disp(['Constructing IRs for source ',ISOU,' (of ',NSOU,') '])
				end
			end
			
			if doaddsources == 1
           		desiredirfile = [Filepath,Filestem,'_',IREC,'_ir.mat'];
           	else
           		desiredirfile = [Filepath,Filestem,'_',ISOU,'_',IREC,'_ir.mat'];
			end
			
            if exist('edpathsfile') ~= 1,
			    usededpathsfile = [Filepath,Filestem,'_',ISOU,'_',IREC,'_edpaths.mat'];
			else
				if SHOWTEXT >= 2,
                    if round(irec/idisp)*idisp==irec,  
					    disp(['   Using existing edpathsfile: ',edpathsfile,'_',ISOU,'_',IREC,'_edpaths.mat'])
                    end
				end
                usededpathsfile = [edpathsfile,'_',ISOU,'_',IREC,'_edpaths.mat'];
			end
           
% % %                sourceonplane = sum(double(souinsideplane(:,isou)).*reflfactors) > 0;
% % %                receiveronplane = sum(double(recinsideplane(:,irec)).*reflfactors) > 0;
% % %                souspecboost = sourceonplane+1;
% % %                recspecboost = receiveronplane+1;
                
           if difforder >= 2,
                
      		    edirfile = EDB2makeirs(usededpathsfile,specorder,...
                    Rstart,EDcalcmethod,edgestartcoords,edgeendcoords,edgenvecs,...
                    edgelengthvec,planeeqs,approxplanemidpoints,reflfactors,closwedangvec,planesatedge,elemsize,...
                    reftoshortlistE,re1sho,re2sho,thetae1sho,thetae2sho,ze1sho,ze2sho,edgeseespartialedge,edgeplaneperptoplane1,desiredirfile,guiderowstouse,directsound,saveindividualdiffirs);
            else

      		    edirfile = EDB2makeirs(usededpathsfile,specorder,...
                    Rstart,EDcalcmethod,edgestartcoords,edgeendcoords,edgenvecs,...
                    edgelengthvec,planeeqs,[],reflfactors,closwedangvec,planesatedge,elemsize,[],[],[],[],[],[],[],[],[],desiredirfile,guiderowstouse,directsound,saveindividualdiffirs);
                
            end
 
	        if doaddsources == 1
				eval(['load ',edirfile])
			
				nirnew = length(irtot);
				nirold = length(irtotaccum);
				if nirnew > nirold,
					irtotaccum = [irtotaccum;zeros(nirnew-nirold,1)]; 
					irdirectaccum = [irdirectaccum;zeros(nirnew-nirold,1)]; 
					irgeomaccum = [irgeomaccum;zeros(nirnew-nirold,1)]; 
					irdiffaccum = [irdiffaccum;zeros(nirnew-nirold,size(irdiff,2))];             
				end
				irtotaccum(1:nirnew) = irtotaccum(1:nirnew) + irtot;
				irdirectaccum(1:nirnew) = irdirectaccum(1:nirnew) + irdirect;
				irgeomaccum(1:nirnew) = irgeomaccum(1:nirnew) + irgeom;
				irdiffaccum(1:nirnew,:) = irdiffaccum(1:nirnew,:) + irdiff;
			end
			
        end		% of the for isou = 1:nsou,
        
        if doaddsources == 1
	        eval(['save ',edirfile,Varlistaccum])
    	end    
	
	end		% of the for irec = 1:nrec,
end

tmakeirs = etime(clock,t00);

% #########################################################################
% #########################################################################
%   
%      Calculate frequency-domain results (TFs) 
%   
% #########################################################################

t00 = clock;

if calctfs == 1,

    if SHOWTEXT >= 1,
        disp(' ')
        disp('####################################################################')
        disp('#')
        disp('#  Calculating frequency-domain results (TFs)')
        disp(' ')
    end
        
    % If the parameter guiderowstouse wasn't specified in the setup file,
    % we give it the default value.
    
    if exist('guiderowstouse') ~= 1,
        guiderowstouse = [];    
    end

    % Some variables that are needed for empty combinations
    tfgeom = [];
    tftot = [];
    tfdirect = [];
    tfdiff = [];
    Varlist = [' tfdirect tfgeom tfdiff tftot frequencies Rstart CAIR elemsize nedgesubs'];
    Varlistaccum = [' tfdirectaccum tfgeomaccum tfdiffaccum tftotaccum frequencies Rstart CAIR elemsize nedgesubs'];

    reclist_in_forloop = reclist_for_calctfs;

	for irec = reclist_in_forloop,
		IREC = int2str(irec);
		if SHOWTEXT >= 1,
			disp(['...receiver ',IREC,' (of ',NREC,')'])
		end
		
		tfdirectaccum = [];
        tfgeomaccum = [];
        tftotaccum = [];
        tfdiffaccum = [];

		if dobackscatter == 1, 
			soulist_for_calctfs = irec;
		end
		
		for isou = soulist_for_calctfs,
    	    ISOU = int2str(isou);
        	if SHOWTEXT >= 1,
                if round(isou/idisp)*idisp==isou,  
					disp(['Calculating TFs for source ',ISOU,' (of ',NSOU,') '])
				end
			end
			
			if doaddsources == 1
           		desiredtffile = [Filepath,Filestem,'_',IREC,'_tf.mat'];
           	else
           		desiredtffile = [Filepath,Filestem,'_',ISOU,'_',IREC,'_tf.mat'];
			end
			
            if exist('edpathsfile') ~= 1,
			    usededpathsfile = [Filepath,Filestem,'_',ISOU,'_',IREC,'_edpaths.mat'];
			else
				if SHOWTEXT >= 2,
                    if round(irec/idisp)*idisp==irec,  
					    disp(['   Using existing edpathsfile: ',edpathsfile,'_',ISOU,'_',IREC,'_edpaths.mat'])
                    end
				end
                usededpathsfile = [edpathsfile,'_',ISOU,'_',IREC,'_edpaths.mat'];
            end
                        
            if difforder >= 2,
                
                    edtffile = EDB2maketfs(usededpathsfile,specorder,frequencies,...
                    Rstart,EDcalcmethod,edgestartcoords,edgeendcoords,edgenvecs,...
                    edgelengthvec,planeeqs,approxplanemidpoints,reflfactors,closwedangvec,planesatedge,elemsize,...
                    reftoshortlistE,re1sho,re2sho,thetae1sho,thetae2sho,ze1sho,ze2sho,edgeseespartialedge,edgeplaneperptoplane1,desiredtffile,guiderowstouse,directsound,saveindividualdiffirs);
            else
                
      		    edtffile = EDB2maketfs(usededpathsfile,specorder,...
                    frequencies,Rstart,EDcalcmethod,edgestartcoords,edgeendcoords,edgenvecs,...
                    edgelengthvec,planeeqs,[],reflfactors,closwedangvec,planesatedge,elemsize,[],[],[],[],[],[],[],[],[],desiredtffile,guiderowstouse,directsound,saveindividualdiffirs);
                
            end
 
	        if doaddsources == 1
				eval(['load ',edtffile])
			
% 				nirnew = length(tftot);
% 				nirold = length(tftotaccum);
% 				if nirnew > nirold,
% 					irtotaccum = [irtotaccum;zeros(nirnew-nirold,1)]; 
% 					irdirectaccum = [irdirectaccum;zeros(nirnew-nirold,1)]; 
% 					irgeomaccum = [irgeomaccum;zeros(nirnew-nirold,1)]; 
% 					irdiffaccum = [irdiffaccum;zeros(nirnew-nirold,size(irdiff,2))];             
% 				end
				tftotaccum = tftotaccum + tftot;
				tfdirectaccum = tfdirectaccum + tfdirect;
				tfgeomaccum = tfgeomaccum + tfgeom;
				tfdiffaccum = tfdiffaccum + tfdiff;
			end
			
        end		% of the for isou = 1:nsou,
        
        if doaddsources == 1
	        eval(['save ',edtffile,Varlistaccum])
    	end    
	
	end		% of the for irec = 1:nrec,
end

tmaketfs = etime(clock,t00);


% #########################################################################
% #########################################################################
%   
%      Save the timing data
%   
% #########################################################################

logfilename = [Filepath,Filestem,'_log.mat'];
Varlist = [' tcadgeo tedgeo tsgeo trgeo tsetup tISEStree tfindpaths tmakeirs tmaketfs tprunepaths nsources nreceivers'];
eval(['save ',logfilename,Varlist])
