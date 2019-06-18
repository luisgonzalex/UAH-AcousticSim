function [outputfile] = EDB2readac(CADfile,outputfile,planecornerstype,checkgeom)
% EDB2readac - Reads a file of type .AC (made by e.g. Invis AC3D) and saves all the geometry data in a mat-file. 
%
% Input parameters:
%  CADfile     	(optional) The input file, with or without the .AC extension.
%					If this file is not specified, a file opening window will appear.
%  outputfile      (optional) The output file. If not specified, it will be given the
%                   same name as the CAD file with an '_cadgeo.mat' ending instead of the '.AC'.
%  planecornerstype (optional) Could have the value 'zero' or 'circ'.Default: 'circ'.
%               	Affects the matrix planecorners, see below.
%  checkgeom		(optional) If this parameter is given the value 'check', then a few checks
%                   of the geometry consistency will be done: a check for duplicate corners
%                   for redundant corners and for corners that are connected to only one plane.
%					Only warnings are given. As default, no check is done.
%  SHOWTEXT (global)	If this global parameter has the value 3 or higher, informative text will
%                   be printed on the screen.
%
% Output parameters:
%	outputfile		The name of the outputfile, generated either automatically or specified as 
%					an input parameter.
%
% Output parameters in the outputfile:
% 	corners         Matrix [ncorners,3] with the corner coordinates
%   cornernumbers   Vector [ncorners,1] with the corner numbers used in the CAD-file. That is
%                   in the outputfile, the corners will effectively be renumbered from 1 to ncorners
%                   in the order they appeared in the CAD-file. 
%   planecorners    Matrix [planes,nmaxcornersperplane] with the corner numbers that make up each plane.
%                   Since planes can have different numbers of corners, the number of columns is the 
%                   maximum number of corners that any plane has. NB! The corner numbers are in the
%                   renumbered system, not in the CAD-file system.
%   planenames      Matrix [nplanes,nmaxcharacters1] (sparse), with the planenames in the CAD file.
%   planeabstypes   Matrix [nplanes,nmaxcharacters2] (sparse), with the absorber names in the CAD file. 
% 	planenumbers    Vector [nplanes,1] with the plane numbers used in the CAD file.
%   cornercrossref  Vector [nmaxnumberinCADfile,1] (sparse), which gives the matlab-file corner
%                   numbers for each CAD-file corner number. For instance, if a CAD-file corner
%                   number 1200 is translated to corner number 72 in the matlab-file, then
%                   cornercrossref(1000) = 72.
%   Snumbers        Vector [nsources,1] of the source numbers in the CAD-file. NB! Only CAD-files
%                   v6 use source numbers; later versions use source names instead of numbers. 
%                   For CAD-file v7 and v8, Snumbers is empty.
%   Snames          Matrix [nsources,2] of the source names in the CAD-file. NB! Only CAD-files v7
%                   or later use source names (such as 'A0' etc). Older versions use source numbers, 
%                   in which case Snames is empty.
%   Sdirectivitynames   Matrix [nsources,nmaxcharacters3] (sparse) with the source directivity names 
%                   in the CAD file. 
%   Sdelays         Vector [nsources,1] of the source delays in the CAD-file, in milliseconds.
% 	Scoords         Matrix [nsources,3] of the source coordinates.
%   Sdirections     Matrix [nsources,3] of the aim point coordinates for the sources.
%   Srotations      Vector [nsources,1] of the source rotations, in degrees.
%   Slevels         Matrix [nsources,6/8] of the source levels for six (CAD-files v6) or eight
%                   octave bands (CAD-files v7 or later).
%   Rnumbers        Vector [nreceivers,1] of the receiver numbers in the CAD-file.
%   Rcoords         Matrix [nreceivers,3] of the receiver coordinates.
%   CATTversionnumber   6,7 or 8
%	planeeqs        Matrix [nplanes,4] of the plane equations as derived  from the plane definitions. 
%                   Each row has the values [A B C D] of the plane equation on the form Ax + By + Cz = D
%   planenvecs      Matrix [nplanes,3] of the normalized normal vectors of the planes.
%   ncornersperplanevec     Vector [nplanes,1] which gives the number of corners for each plane.
%   planesatcorners     Matrix [ncorners,4] which gives the plane numbers that each corner is connected to.
%                   NB! Here it is assumed that each corner is connected to maximum four planes.
%   nplanespercorners   Vector [ncorners,1] which gives the number of planes that each corner is connected to.
% 	minvals         Matrix [nplanes,3] which for each plane gives the smallest coordinate in the x-, y- and z-direction.
% 	maxvals         Matrix [nplanes,3] which for each plane gives the largest coordinate in the x-, y- and z-direction.
%   planehasindents Vector [nplanes,1] which for each plane gives the
%                   number of indeting corners
%   indentingcorners  Matrix [nplanes,max(ncornersperplanevec)] which for each plane gives the number to the first corner
%                   in a corner triplet which identifies an indenting corner. The number of the corner is the
%                   order given in planecorners for that plane.
%
% Uses the functions EDB2fistr, EDB2extrnums, EDB2cross
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
% Peter Svensson (svensson@iet.ntnu.no) 18Jul09
%
% ----------------------------------------------------------------------------------------------
%   This File is created based on EDB2eadcad.m to support the Edge Diffraction Toolbox
% ----------------------------------------------------------------------------------------------
% Alexander Pohl (alexander.pohl@hcu-hamburg.de) 25May2011
%
% [outputfile] = EDB2readac(CADfile,outputfile,planecornerstype,checkgeom);

% Modified by Peter Svensson on 3Dec2011

global SHOWTEXT

if nargin == 0
	CADfile = '';
	outputfile = '';
	planecornerstype = '';
	checkgeom = '';
elseif nargin == 1
	outputfile = '';
	planecornerstype = '';
	checkgeom = '';
elseif nargin == 2
	planecornerstype = '';
	checkgeom = '';
elseif nargin == 3
	checkgeom = '';
end

geomacc = 1e-10;

%---------------------------------------------------------------
% If no AC-file was specified, present a file opening window

if isempty(CADfile)
	[CADfile,filepath] = uigetfile('*.ac','Please select the AC3Dfile');
    
    [filepath,temp1] = fileparts(filepath);
	
	if ~isstr(CADfile)
		return
	end

    [temp1,filestem] = fileparts(CADfile);


    CADfile = [[filepath,filesep],filestem,'.ac'];
else

    [filepath,filestem,fileext] = fileparts(CADfile);
    CADfile = [[filepath,filesep],filestem,fileext];
end


if exist(CADfile) ~= 2
	error(['ERROR: AC3D-file: ',CADfile,' can not be found'])
end

%---------------------------------------------------------------
% If no output file was specified, construct an automatic file name

if isempty(outputfile)
	outputfile = [[filepath,filesep],filestem,'_cadgeo.mat'];
end

%---------------------------------------------------------------
% Read in the entire file into the string B
% Find the line starts and ends

fid = fopen([CADfile],'r');
if fid == 0
	error(['ERROR: The AC3Dfile ',CADfile,' could not be opened.'])
end

%---------------------------------------------------------------
% read header "AC3Db"
tline = fgetl(fid);
if (~strcmp(tline,'AC3Db'))
    disp('AC3Db header missing');
end

% read materials
p_NumberOfMaterials=0;
while(true)
    tline = fgetl(fid);
    tokens = readLineToToken(tline);
    
    if(strcmp(tokens{1},'MATERIAL'))
        p_MaterialList{p_NumberOfMaterials+1}=tokens{2};
        p_NumberOfMaterials=p_NumberOfMaterials+1;
        
    else
        break;
    end  
end

% read world object (all objects are kids of this)
m_LocalMatrix = [1.0 0.0 0.0 0.0; 
                 0.0 1.0 0.0 0.0; 
                 0.0 0.0 1.0 0.0; 
                 0.0 0.0 0.0 1.0];

% Read OBJECT world
tokens = readLineToToken(tline);
if(~strcmp(tokens{1},'OBJECT'))    
    disp('Section does not start with "object"');
end
if(~strcmp(tokens{2},'world'))    
    disp('Section does not start with "world"');
end

% Read kids
tline = fgetl(fid);
tokens = readLineToToken(tline);
if(~strcmp(tokens{1},'kids'))    
    disp('root object section does not start with kids');
end
p_Kids = str2num(tokens{2});

for i = 1:p_Kids
    p_Objects{i}=read_ac3d_object(fid, m_LocalMatrix);
end

fclose(fid);

 %% Export
eval(['save a.mat p_Objects']);
% Combine Vertices of EVERY Object to on corner list (incl. Counting PlanesPerPoly)
for i = 1:p_Kids
    if(i==1)
        corners=p_Objects{i}{1}';
        p_CornersPerPoly = size(p_Objects{i}{1},2);
        p_PlanesPerPoly = size(p_Objects{i}{2},2);
    else
        corners = [corners; p_Objects{i}{1}'];
        p_PlanesPerPoly = [p_PlanesPerPoly size(p_Objects{i}{2},2)];
        p_CornersPerPoly = [p_CornersPerPoly size(p_Objects{i}{1},2)];
    end
end

% Enumerate each corner from 1:end
ncorners = size(corners,1);
cornernumbers=(1:ncorners)';

%Number of Planes is sum over all PlanesPerPoly
nplanes=sum(p_PlanesPerPoly);

ncornersperplanevec=zeros(nplanes,1);
currentPlaneIndex=1;
for i = 1:p_Kids
    for j=1:p_PlanesPerPoly(i);
        ncornersperplanevec(currentPlaneIndex)=length(p_Objects{i}{2}{j});
        currentPlaneIndex=currentPlaneIndex+1;
    end
end

% We need to compute the maximum number of vertices for any plane 
% (others are filled up with 0)
p_MaximumNumberOfVerticesOfPlane = max(ncornersperplanevec);

% Compute indices of corners for each plane
% Simple for first Object, as here numbers are identical
% For every further object, the number of previous used corners has be be
% taken as offset
planecorners=zeros(nplanes,p_MaximumNumberOfVerticesOfPlane);
currentPlaneIndex=1;
p_Offset=0;
for i = 1:p_Kids
    for j=1:p_PlanesPerPoly(i);
        for k=1:length(p_Objects{i}{2}{j})
            planecorners(currentPlaneIndex,k)=p_Objects{i}{2}{j}(k)+p_Offset;
        end
        currentPlaneIndex=currentPlaneIndex+1;
    end
    p_Offset=p_Offset+p_CornersPerPoly(i);
end

planenames=zeros(nplanes,30);
currentPlaneIndex=1;
for i = 1:p_Kids
    for j=1:p_PlanesPerPoly(i);
        p_Text = sparse(double(p_MaterialList{p_Objects{i}{4}{j}+1}));
        
        if length(p_Text) > size(planenames,2)
           planenames = [planenames zeros(planenumbers(currentPlaneIndex),length(p_Text)-size(planenames,2))]; 
        end
        planenames(currentPlaneIndex,1:length(p_Text)) = p_Text;
        currentPlaneIndex=currentPlaneIndex+1;
    end
end

planeabstypes = sparse(planenames(:,1:6));

% Addition on 2.12.2011 by Peter Svensson
% The absorber names had an initial " which must be removed.
% We remove columns until A-Z,a-z

if size(planeabstypes,1) > 1
    columnhasalphabet = 1 - sign( prod(sign(64-full(planeabstypes))+1));
else
    columnhasalphabet = 1 - sign( (sign(64-full(planeabstypes))+1));    
end
firstOKcolumn = find(columnhasalphabet);
if isempty(firstOKcolumn)
    error('ERROR: All plane absorber types have names without any letters. They must start with a letter')
end
firstOKcolumn = firstOKcolumn(1);
planeabstypes = planeabstypes(:,firstOKcolumn:end);


planenumbers=1:nplanes;

eval(['save ',outputfile,' corners planecorners ncornersperplanevec'])

%% Begining from here, it is code based of EDB2readcad.m without changes - except SOURCE - RECIEVER Values, as not in AC3Dfile
%---------------------------------------------------------------
% Change the numbering in the CAD file to a contiguous one

% Corner numbers
% First we make an inelegant crossreference vector giving the
% resulting corner number in the position given by the CAD-file
% corner number.

% PS 20120414: CHECK Corners with number zero are not handled correctly.
% Insert code from EDB2readcad here.

maxnumb = max(cornernumbers);
cornercrossref = sparse(zeros(1,maxnumb));
cornercrossref(cornernumbers) = [1:ncorners];

[nplanes,ncornersperplane] = size(planecorners);
iv = find(planecorners~=0);
planecorners(iv) = full(cornercrossref(planecorners(iv)));

%---------------------------------------------------------------
% Go through all planes. If there is a plane definition including
% zeros, and planecornerstype == 'circ', expand it repeating the
% same corner order again.

if isempty(planecornerstype)
	planecornerstype = 'circ';
else
	planecornerstype = setstr(lower(planecornerstype(1)));
	if planecornerstype(1) == 'z'
		planecornerstype = 'zero';
	else
		planecornerstype = 'circ';
	end
end

if planecornerstype == 'circ'
	for ii = 1:nplanes
		iv = find( planecorners(ii,:) ~= 0);
		ncornersatplane = length(iv);
		if ncornersatplane ~= ncornersperplane
			pattern = planecorners(ii,iv);
			nrepeatings = ceil(ncornersperplane/ncornersatplane);
			for jj = 1:nrepeatings-1
				pattern = [pattern planecorners(ii,iv)];
			end
			planecorners(ii,:) = pattern(1:ncornersperplane);
		end
	end
end


%---------------------------------------------------------------
% Find the normal vectors to the planes using the cross products
%
% The normal vector is basically the cross product between two vectors
% connecting three consecutive points. If there are indents though
% they will cause reversed normal vectors, so one must go through all
% combinations and choose the majority normal vector.
%
% 26mar09  Use the fact described above for detecting indention corners

planenvecs = zeros(nplanes,3);
planehasindents = zeros(nplanes,1);
indentingcorners = sparse(zeros(nplanes,max(ncornersperplanevec)));

for ii = 1:nplanes
    
	iv = find(planecorners(ii,:)~=0);
	cornerlist = planecorners(ii,iv);
	iv = find(cornerlist == cornerlist(1));
	if length(iv) > 1
		cornerlist = cornerlist(1:iv(2)-1);
	end
	ncorners = length( cornerlist );
	cornerlist = [cornerlist cornerlist(1) cornerlist(2)];

	nvectorlist = zeros(ncorners,3);
	nveclen = zeros(ncorners,1);	
	
	for jj = 1:ncorners
		co1numb = cornerlist(jj);
		co2numb = cornerlist(jj+1);
		co3numb = cornerlist(jj+2);
		vec1 = corners(co1numb,:) - corners(co2numb,:);
		vec2 = corners(co3numb,:) - corners(co2numb,:);

		nvec = EDB2cross(vec1.',vec2.').';
		nveclen(jj) = norm(nvec);
		if nveclen(jj) > 0
			nvectorlist(jj,:) = nvec./nveclen(jj);
		end
	end
	
	iv = find(nveclen < max(nveclen)*0.001);
	nvectorlist(iv,:) = [];
	nvecref = nvectorlist(1,:);

    [n1,slask] = size(nvectorlist);
	nvecsigns = round(sum(     (nvectorlist.').*(nvecref(ones(n1,1),:).')      ));
    
    if sum(nvecsigns) == 0
        disp(' ')        
       error(['ERROR: Plane ',int2str(planenumbers(ii)),' (plane numbering as in the CAD file) seems to be twisted.'])        
    end
    
    
    if abs(sum(nvecsigns)) ~= n1
       disp(['Plane ',int2str(ii),' has indents!']) 
       nindents = (n1 - abs(sum(nvecsigns)))/2;
       planehasindents(ii) = nindents;
       ivindent = find(nvecsigns == -1);
       if length(ivindent) == nindents
            indentingcorners(ii,ivindent) = 1;
       else
           if length(ivindent) == ncorners - nindents
               ivindent = find(nvecsigns == 1);
                indentingcorners(ii,ivindent) = 1;           
           else
              error(['ERROR: An unexpected problem. Please report to the developer']) 
           end
       end

    end
    
    nvecdiff = [nvectorlist(2:n1,1).*nvecsigns(2:n1).' nvectorlist(2:n1,2).*nvecsigns(2:n1).' nvectorlist(2:n1,3).*nvecsigns(2:n1).'] - nvecref(ones(n1-1,1),:);
        
    if n1 > 2
        nvecdiff = sum(nvecdiff.'.^2).';    
    else
        nvecdiff = norm(nvecdiff);
    end

    %APO Change
    if any(nvecdiff>1e-3),        
        nvecdiff
        error(['ERROR: Normal vectors for plane ',int2str(planenumbers(ii)),' (in the CAD file, = ',int2str(ii),' in the EDB2 file), get different normal vectors for consecutive corner triplets. Check the geometry in the CAD-file'])
    elseif any(nvecdiff>1e-8)
        nvecdiff
        disp(['WARNING: Normal vectors for plane ',int2str(planenumbers(ii)),' (in the CAD file, = ',int2str(ii),' in the EDB2 file), get somewhat different normal vectors for consecutive corner triplets. Check the geometry in the CAD-file'])
    end
    
	if ncorners > 5 & abs(sum(nvecsigns)) <= 1
		disp(['WARNING for plane number ',int2str(planenumbers(ii)),' in the CAD-file'])
		disp(['   with the name ',strtrim(setstr(full(planenames(ii,:))))])
		disp('   The normal vector can not be determined for this plane because there are')
		disp('   the same number of inwards and outwards corners')
		disp('   This is a list of the possible normal vectors:')
		[nv1,slask] = size(nvectorlist);
		for kk = 1:nv1
			vecstr = ['   ',int2str(kk),'. ',num2str(-nvectorlist(kk,1)),' ',num2str(-nvectorlist(kk,2)),' ',num2str(-nvectorlist(kk,3))];
			disp(vecstr)
		end
		disp(' ')
	
      preferredsign = input('   Give the number of a correct normal vector for this plane please ');
      switchsign = nvecsigns.'./nvecsigns(preferredsign);
		nvectorlist = nvectorlist.*switchsign(:,ones(1,3));	
	else
		mostcommonsign = sign(sum(nvecsigns));

		switchsign = nvecsigns.'./mostcommonsign;
		nvectorlist = nvectorlist.*switchsign(:,ones(1,3));
    end

	planenvecs(ii,:) = mean(nvectorlist);
end

planenvecs = -planenvecs;

%---------------------------------------------------------------
% Plane equations, as Ax + By + Cz = D for each plane

planeeqs = zeros(nplanes,4);
planeeqs(:,1:3) = planenvecs;
planeeqs(:,4) = sum( (planenvecs.').*(corners(planecorners(:,1),:).')  ).';

%---------------------------------------------------------------
% Useful data: planesatcorners, minvals and maxvals

[ncorners,slask] = size(corners);
planesatcornerhits = zeros(ncorners,nplanes);

for ii = 1:nplanes
	cornerlist = planecorners(ii,1:ncornersperplanevec(ii));
	planesatcornerhits(cornerlist,ii) = planesatcornerhits(cornerlist,ii) + 1;
end

maxplanespercorner = 0;
for ii = 1:ncorners
	nplanes = length(find(planesatcornerhits(ii,:) ~= 0));
	if nplanes > maxplanespercorner
		maxplanespercorner = nplanes;
	end	
end

planesatcorners = zeros(ncorners,maxplanespercorner);
nplanespercorners = zeros(ncorners,1);
for ii = 1:ncorners
	iv = find(planesatcornerhits(ii,:)~=0);
	planesatcorners(ii,1:length(iv)) = iv;
	nplanespercorners(ii) = length(iv);
end

% find cubic boxes inside which the planes are placed

[nplanes,slask] = size(planeeqs);

minvals = zeros(nplanes,3);
maxvals = zeros(nplanes,3);

for ii = 1:nplanes
	cornerlist = planecorners(ii,:);
	cornerlist = cornerlist( find(cornerlist~= 0) );
	cornercoord = corners(cornerlist,:);
	minvals(ii,:) = min(cornercoord);
	maxvals(ii,:) = max(cornercoord);
end

minvals = minvals - geomacc;
maxvals = maxvals + geomacc;

%---------------------------------------------------------------
% Check the geometry a bit

if isempty(checkgeom)
	checkgeom = 0;
else
	if setstr(checkgeom(1)) == 'c'
		checkgeom = 1;
	else
		checkgeom = 0;
	end
end

if checkgeom

	badproblem = 0;
	disp(' ')
	disp('   Checking if there are corners with identical coordinates')
	for ii = 1:ncorners-1
		othercorners = [ii+1:ncorners];
		iv = find((corners(othercorners,1)==corners(ii,1)) & (corners(othercorners,2)==corners(ii,2)) & (corners(othercorners,3)==corners(ii,3)));
		if ~isempty(iv)
			disp(['      WARNING: corner ',int2str(ii),' and ',int2str(othercorners(iv(1))),' have the same coordinates'])
			disp(['      This should be fixed in the CAD-file or in the CATT-model'])
			badproblem = 1;
		end
	end

	disp('   Checking if there are unused corners')
	iv = find(planesatcorners(:,1)==0);
	if ~isempty(iv)
		disp('      WARNING: The following corners in the CAD-file are not used:')
		printvec = int2str(cornernumbers(iv(1)));
		for iiprint = 2:length(iv)
			printvec = [printvec,' ',int2str(cornernumbers(iv(iiprint)))];
		end
		disp(['            ',printvec])
		disp(['      This is not a big problem'])
		planesatcorners(iv,:) = [];
	end

	disp('   Checking if any corners seem redundant')
	iv = find(planesatcorners(:,2)==0);
	if ~isempty(iv)
		disp('      WARNING: The following corners in the CAD-file belong to only one plane:')
		printvec = int2str(cornernumbers(iv(1)));
		for iiprint = 2:length(iv)
			printvec = [printvec,' ',int2str(cornernumbers(iv(iiprint)))];
		end
		disp(['            ',printvec])
		disp(['      This should be fixed in the CAD-file or in the CATT-model'])
		badproblem = 1;
	end

	if badproblem == 1
		disp(' ');
		error(['ERROR: Problems with the geometry defined in the CAD-file: ',CADfile])
	end

end

%---------------------------------------------------------------
% Save the relevant variables in the output file

% NB! Sparse matrices can not get a non-double format

if ncorners < 256
    planecorners = uint8(planecorners);
elseif ncorners < 65536
    planecorners = uint16(planecorners);    
end   

if max(ncornersperplanevec) <= 255
    ncornersperplanevec = uint8(ncornersperplanevec);
    planehasindents = uint8(planehasindents);
else
    ncornersperplanevec = uint16(ncornersperplanevec);
    planehasindents = uint16(planehasindents);
end

Varlist = [' corners cornernumbers planecorners planenames planeabstypes '];
Varlist = [Varlist,' planenumbers cornercrossref '];
Varlist = [Varlist,' planeeqs planenvecs ncornersperplanevec planesatcorners nplanespercorners'];
Varlist = [Varlist,' minvals maxvals planehasindents indentingcorners'];
%Varlist = [Varlist,' Scoords Sdirections Srotations Slevels Rnumbers Rcoords CATTversionnumber'];

planenames = sparse(planenames+1-1);
planeabstypes = sparse(planeabstypes+1-1);
%Sdirectivitynames = sparse(Sdirectivitynames+1-1);

eval(['save ',outputfile,Varlist])

end

function [p_Plane, p_Norm, p_Mat] = read_ac3d_surface(fid, i_Vertices)
    tline = fgetl(fid);
    tokens = readLineToToken(tline);
    if(~strcmp(tokens{1},'SURF'))    
        disp(['surface section does not start with SURF',tokens{1}]);
    end
    if(~strcmp(tokens{2},'0x10'))    
        disp(['only polygons are accepted TODO?', tokens{2}]);
    end
    
    tline = fgetl(fid);
    tokens = readLineToToken(tline);
    if(~strcmp(tokens{1},'mat'))    
        disp('missing expected token "mat"');
    end
    p_Mat=str2num(tokens{2});
    
    tline = fgetl(fid);
    tokens = readLineToToken(tline);
    
    if(~strcmp(tokens{1},'refs'))    
        disp('surface section missing refs');
    end
    
    p_NumberOfRefs = str2num(tokens{2});
    
    if(tokens{2}<3)
        disp(['too less points in surface:',tokens{2}]);
        return;
    end
    
    p_Plane=zeros(p_NumberOfRefs,1);
    
    for r=1:p_NumberOfRefs
        tline = fgetl(fid);
        tokens = readLineToToken(tline);
        p_Plane(r)=str2num(tokens{1})+1; 
    end
    
    v1 = i_Vertices(:,p_Plane(2))-i_Vertices(:,p_Plane(1));
    v2 = i_Vertices(:,p_Plane(3))-i_Vertices(:,p_Plane(1));
    p_Norm = cross(v1,v2);
    p_Norm = p_Norm./norm(p_Norm);
end

function [token] = readLineToToken(tline)
	[token{1}, remain] = strtok(tline);
    i=1;
    while(remain~=0)
        i=i+1;
        tline = remain;
        [token{i}, remain] = strtok(tline);
    end
end

function [p_Object] = read_ac3d_object(fid, i_LocalMatrix)

    tline = fgetl(fid);
    tokens = readLineToToken(tline);
    if(~strcmp(tokens{1},'OBJECT'))    
        disp('Section does not start with "object"');
    end
    if(strcmp(tokens{2},'poly'))    
        p_Object = read_ac3d_poly(fid, i_LocalMatrix);
    elseif(strcmp(tokens{2},'group'))    
        disp('GROUP not implemented yet');
    elseif(strcmp(tokens{2},'light'))    
        disp('LIGHT not implemented yet');
    else        
        disp('Section does not start with valid key');
    end
end
function [p_Polygon] = read_ac3d_poly(fid, i_LocalMatrix)
    
    p_LocalMatrix = [1.0 0.0 0.0 0.0; 
                     0.0 1.0 0.0 0.0; 
                     0.0 0.0 1.0 0.0; 
                     0.0 0.0 0.0 1.0];
    
    tline = fgetl(fid);
    
    p_Vertices=0;
    while(ischar(tline))
        tokens = readLineToToken(tline);

        if(strcmp(tokens{1},'numvert'))
            
            p_ResultMatrix = i_LocalMatrix * p_LocalMatrix;
            
            p_NumVertices=str2num(tokens{2});
            p_Vertices=zeros(3,p_NumVertices);
            for l=1:p_NumVertices
               tline = fgetl(fid);
               tokens = readLineToToken(tline);
               if(~ischar(tline))
                   disp(['unexpected EOF while reading vertex (', num2str(l),'/',num2str(p_NumVertices)]);
               end
               if(size(tokens,2)~=3)
                   disp(['error while reading vertex coords (', num2str(l),'/',num2str(p_NumVertices)]);
               end
               
               p_Vec = p_ResultMatrix * [str2num(tokens{1}); str2num(tokens{2}); str2num(tokens{3}); 1];
               p_Vertices(:,l) = p_Vec(1:3);
            end
            
        elseif(strcmp(tokens{1},'numsurf'))
            p_NumSurfaces=str2num(tokens{2});
            for s=1:p_NumSurfaces
                [p_Planes{s}, p_Norms{s}, p_Mats{s}] = read_ac3d_surface(fid, p_Vertices);
            end
        elseif(strcmp(tokens{1},'loc'))
            p_LocalMatrix(1,4)=str2num(tokens{2});
            p_LocalMatrix(2,4)=str2num(tokens{3});
            p_LocalMatrix(3,4)=str2num(tokens{4});
        elseif(strcmp(tokens{1},'rot'))
            p_LocalMatrix(1,1)=str2num(tokens{2});
            p_LocalMatrix(1,2)=str2num(tokens{3});
            p_LocalMatrix(1,3)=str2num(tokens{4});
            p_LocalMatrix(2,1)=str2num(tokens{5});
            p_LocalMatrix(2,2)=str2num(tokens{6});
            p_LocalMatrix(2,3)=str2num(tokens{7});
            p_LocalMatrix(3,1)=str2num(tokens{8});
            p_LocalMatrix(3,2)=str2num(tokens{9});
            p_LocalMatrix(3,3)=str2num(tokens{10});
        elseif(strcmp(tokens{1},'kids'))
            if(str2num(tokens{2})>2)
                disp(['Warning: Kids of Polygon not handled yet:',tokens{2}])
                p_Polygon=0;
                return;
            end
            p_Polygon{1}=p_Vertices;
            p_Polygon{2}=p_Planes;
            p_Polygon{3}=p_Norms;
            p_Polygon{4}=p_Mats;
            return;
        else
            if(strcmp(tokens{1},'name'))
                % Only name of poly, irrelevant -->SKIP
            elseif(strcmp(tokens{1},'crease'))
            else
                disp(['Warning: Ignored Entry in AC3D File: ', tokens{1}]);
            end
        end
        tline = fgetl(fid);
    end
end
