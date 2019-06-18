function outputfile = EDB2SorRgeo(eddatafile,outputfile,pointcoords,typeofcoords,difforder,nedgesubs)
% EDB2SorRgeo - Calculates some source- or receiver-related geometrical parameters.
% Calculates some source- or receiver-related geometrical parameters
% based on corners,edges and planes in an eddata-file
% and on a list of source/receiver coordinates.
% The output is saved in a .mat-file.
%
% Input parameters:
%	eddatafile (optional)	The .mat file with the corner,edge and plane data.			
%							If it is not specified, a file opening window is presented.
%	outputfile (optional)	The .mat file where the output data will be stored.
%							If it is not specified, an automatic file name is constructed
%							using the same filename stem as the eddatafile, but ending with '_srdata'.
%	pointcoords				Matrix, [nrec,3], of source or receiver coordinates.
%   typeofcoords            'S' or 'R' - specifying if the point coordinates are sources
%                           or receivers. This determines what the output data in the output
%                           file will be called.
%   difforder               The highest order of diffraction that will be calculated.
%	nedgesubs (optional)	The number of subdivisions that each edge will be
%							subdivided into for visibility/obstruction tests. Default: 2
%							NB! When nedgesubs = 2, the two end points will be checked.
%	SHOWTEXT (global)		An integer value; if it is 4 or higher, then messages will be printed on the
%							screen during calculations.
% Output parameters:
%	outputfile				The name of the outputfile, generated either automatically or specified as 
%							an input parameter.
%
% Data in the outputfile:
%   sources/receivers                       (renamed) copy of the input parameter 'pointcoords'
%   visplanesfroms/visplanesfromr           Matrix, [nplanes,nsources]/[nplanes/nreceivers]
%                                           with integer values 0-5:
%                           0 means S/R behind a plane which is reflective or totabs
%                           1 means S/R is aligned with a plane, but outside it
%                           2 means S/R is in front of a plane which is reflective
%                           3 means S/R is in front of a plane which is totabs
%                           4 means S/R is inside a plane which is reflective
%                           5 means S/R is inside a plane which is totabs
%   vispartedgesfroms/vispartedgesfromr     Matrix, [nedges,nsources]/[nedges/nreceivers]
%                                           with integer values 0-2^nedgesubs-1 indicating
%                                           part-visibility of the edge.
%   soutsidemodel/routsidemodel             List, [nsources,1]/[nreceivers,1], with values
%                                           1 or 0 indicating whether the S/R is outside (1)
%                                           the model or not (0).
%                                           NB! Only simple tests are done, so an S/R could still
%                                           be outside the model even if the value here is 0.
%
% NB! The text on the screeen, and in the code refers to 'R' or 'receivers' but it should be S or R.
%
% Uses the functions 	EDB2strpend, EDB2infrontofplane, EDB2poinpla, 
%                       EDB2getedgepoints, EDB2checkobstr_pointtoedge
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
% Peter Svensson (svensson@iet.ntnu.no) 060601
%
% outputfile = EDB2SorRgeo(eddatafile,outputfile,pointcoords,typeofcoords,difforder,nedgesubs);

global SHOWTEXT

geomacc = 1e-10;

if nargin < 6
	nedgesubs = 2;
	if nargin < 4
		error('ERROR: The input parameter pointcoords must be specified')
	end
else
	if isempty(nedgesubs), nedgesubs = 2; end
end

typeofcoords = lower(typeofcoords(1));
if typeofcoords~='r' & typeofcoords~='s'
    error(['ERROR: The input parameter typeofcoords must have the value S or R'])    
end

%---------------------------------------------------------------
% If no eddatafile was specified, present a file opening window

if exist('eddatafile')~=1
	eddatafile = '';	
end
if isempty(eddatafile)
	[eddatafile,filepath] = uigetfile('*.mat','Please select the eddatafile');
    [filepath,temp1,temp2] = fileparts(filepath);
	if ~isstr(eddatafile)
		return
	end
	[temp1,filestem,temp2] = fileparts(eddatafile);
    eddatafile = [[filepath,filesep],filestem,'.mat'];
else
	[filepath,filestem,fileext] = fileparts(eddatafile);    
    eddatafile = [[filepath,filesep],filestem,'.mat'];
end


%---------------------------------------------------------------
% If no output file was specified, construct an automatic file name

if exist('outputfile')~=1
	outputfile = '';	
end
if isempty(outputfile)
	filestem = EDB2strpend(filestem,'_eddata');
	outputfile = [[filepath,filesep],filestem,'_rdata.mat'];
end

%---------------------------------------------------------------

eval(['load ',eddatafile])
clear cornerinfrontofplane

[ncorners,slask] = size(corners);
[nplanes,slask] = size(planecorners);
maxncornersperplane = double(max(ncornersperplanevec));
[nedges,slask] = size(edgecorners);
[nreceivers,slask] = size(pointcoords);

[n1,n2] = size(canplaneobstruct);
if n2>1
    canplaneobstruct = canplaneobstruct.';
end

totalmodelmin = min(minvals);
totalmodelmax = max(maxvals);

zerosvecE1 = zeros(nedges,1);
zerosvec1P = zeros(1,nplanes);
onesvec1R = uint8(ones(1,nreceivers));
onesvec1ES = ones(1,nedgesubs);
onesvecP1 = ones(nplanes,1);
onesvec1Max = ones(1,maxncornersperplane);

%##################################################################
%##################################################################
%##################################################################
%
%		PLANE RELATED PARAMETERS
%
%##################################################################

activeplanelist = find(reflfactors ~= 0);
totabsplanelist = find(reflfactors == 0);
ntotabsplanes = length(totabsplanelist);
nthinplanes = length(find(planeisthin));

if SHOWTEXT >= 3
    if lower(typeofcoords(1)) == 'r'
    	disp('         Checking visible planes from R')
    else
    	disp('         Checking visible planes from S')        
    end
end

%--------------------------------------------------------------
%
%		visplanesfromr      [nplanes,nrec]  uint8
%
%           0 means S/R behind a plane which is reflective or totabs
%           1 means S/R is aligned with a plane
%           2 means S/R is in front of a plane which is reflective
%           3 means S/R is in front of a plane which is totabs
%           4 means S/R is inside a plane which is reflective
%           5 means S/R is inside a plane which is totabs
%
% We can call EDB2infrontofplane with a single call by extracting receiver
% numbers and plane numbers for the [nplanes,nsources] matrix.
%
% EDB2infrontofplane returns -1 for behind, 0 for in-plane with and 1 for in
% front of so if we add the value 1 we get close to the final value we want
% to have. Another modification is to give totabs planes the value 3
% instead of 2. A last modification is to find which S/R are inside the
% plane.

iv = [1:nplanes*nreceivers].';                
if nreceivers < 256
    colnumb = uint8(ceil(iv/nplanes));             % This is the receiver number
elseif nreceivers < 65536
    colnumb = uint16(ceil(iv/nplanes));             % This is the receiver number
else
    colnumb = uint32(ceil(iv/nplanes));             % This is the receiver number    
end
if nplanes < 256
    rownumb = uint8(iv - (double(colnumb)-1)*nplanes);     % This is the plane number
elseif nplanes < 65536
    rownumb = uint16(iv - (double(colnumb)-1)*nplanes);     % This is the plane number
else
    rownumb = uint32(iv - (double(colnumb)-1)*nplanes);     % This is the plane number
end
clear iv 
visplanesfromr = EDB2infrontofplane(pointcoords(colnumb,:),planenvecs(rownumb,:),corners(planecorners(rownumb,1),:),corners(planecorners(rownumb,2),:)) + 1;
clear rownumb colnumb

if ntotabsplanes > 0        
    colvec = [0:nreceivers-1]*nplanes;
    iv = uint32(totabsplanelist(:,onesvec1R) + colvec(ones(ntotabsplanes,1),:));    
    visplanesfromr(iv) = uint8(visplanesfromr(iv) + (visplanesfromr(iv)==2));
    clear iv colvec
    visplanesfromr = uint8(visplanesfromr);
end

% For all the planes that the S/R are aligned with
% check if the S/R is inside or outside the plane.
%
% We can call EDB2poinpla with a single call by extracting S/R coordinates
% (=colnumb) and planenumbers (=rownumb).

ivec_visR1 = find(visplanesfromr==1);
if ~isempty(ivec_visR1)
	if nreceivers < 256
        colnumb = uint8(ceil(ivec_visR1/nplanes));             % This is the receiver number
	elseif nreceivers < 65536
        colnumb = uint16(ceil(ivec_visR1/nplanes));             % This is the receiver number
	else
        colnumb = uint32(ceil(ivec_visR1/nplanes));             % This is the receiver number    
	end
	if nplanes < 256
        rownumb = uint8(ivec_visR1 - (double(colnumb)-1)*nplanes);     % This is the plane number
	elseif nplanes < 65536
        rownumb = uint16(ivec_visR1 - (double(colnumb)-1)*nplanes);     % This is the plane number
	else
        rownumb = uint32(ivec_visR1 - (double(colnumb)-1)*nplanes);     % This is the plane number
	end
    hitvec = find(EDB2poinpla(pointcoords(colnumb,:),rownumb,minvals,maxvals,planecorners,corners,ncornersperplanevec,planenvecs));
    if ~isempty(hitvec)
        visplanesfromr(ivec_visR1(hitvec)) = uint8(ones(size(hitvec))*4);
        if ntotabsplanes > 0
            insidetotabsplane = (reflfactors(rownumb(hitvec))==0);
            ivec2 = find(insidetotabsplane);
            if ~isempty(ivec2)
                visplanesfromr(ivec_visR1(hitvec(ivec2))) = uint8(ones(size(hitvec(ivec2)))*5);
            end
        end
        ivec_visR1(hitvec) = [];
        if nthinplanes > 0
            insidethinplane = planeisthin(rownumb(hitvec));
            ivec2 = find(insidethinplane);
            if ~isempty(ivec2)
                recinthinplane = colnumb(hitvec(ivec2));
                thinplanenumber = rownumb(hitvec(ivec2));
        	    error(['ERROR: R number ',int2str(double(recinthinplane(:).')),' has been placed exactly on the thin plane ',int2str(double(thinplanenumber(:).')),...
				    ' so it is undefined which side of the plane the R is. Move the R a short distance away'])
            end
        end
    end
end

visplanesfromr = reshape(visplanesfromr,nplanes,nreceivers);

%--------------------------------------------------------------
%
%       routsidemodel       (0 or 1, size [1,nrec])
%
% We make a simple check: if the S/R is outside the big cube defined by all
% the planes, and we have an interior problem, then S/R must be outside the
% model.

routsidemodel = uint8(zeros(1,nreceivers));
if int_or_ext_model == 'i'
    ivec = find(pointcoords(:,1)<totalmodelmin(1) | pointcoords(:,2)<totalmodelmin(2) | pointcoords(:,3)<totalmodelmin(3) | ...
          pointcoords(:,1)>totalmodelmax(1) | pointcoords(:,2)>totalmodelmax(2) | pointcoords(:,3)>totalmodelmax(3));
	if ~isempty(ivec)
        routsidemodel(ivec) = ones(size(ivec));
	end
	routsidemodel = uint8(routsidemodel);
    
    % If we have a convex model, then the obstruction check is turned off
    % so we must check here if the source/receiver is inside.
    if sum(closwedangvec<pi) == 0

        nplanesvisible = sum(sign(visplanesfromr));
        ivec = find(nplanesvisible < nplanes);
        if ~isempty(ivec)
            routsidemodel(ivec) = ones(size(ivec));        
        end
    end
end

%##################################################################
%##################################################################
%##################################################################
%
%		EDGE RELATED PARAMETERS
%
%##################################################################

if difforder == 0
    if typeofcoords=='r'
        vispartedgesfromr = [];
        receivers = pointcoords;
        Varlist = ['  visplanesfromr  receivers '];
        Varlist = [Varlist,'  vispartedgesfromr  routsidemodel '];
        eval(['save ',outputfile,Varlist])
        return
    else
        vispartedgesfroms = [];
        sources = pointcoords;
        visplanesfroms = visplanesfromr;
        soutsidemodel = routsidemodel;
        Varlist = ['  visplanesfroms  sources '];
        Varlist = [Varlist,'  vispartedgesfroms  soutsidemodel '];
        eval(['save ',outputfile,Varlist])
        return
    end
end

closwedlargerthanpi = closwedangvec>pi;
closwedsmallerthanpi = closwedangvec<pi;

if SHOWTEXT >= 3
    if lower(typeofcoords(1)) == 'r'
    	disp('         Checking which edges are seen from R')
    else
    	disp('         Checking which edges are seen from S')        
    end    
end

nbigcombs = nplanes*nreceivers;

%--------------------------------------------------------------
%
%		visedgesfromr      [nedges,nrec]   uint8
%       
%       6    Edge belongs to a plane which is aligned with R and thin and
%            rigid (but R is not inside the plane)
%       5    Edge belongs to a plane which is aligned with R and not thin 
%
% These can be derived from the cases where visplanesfromr = 1

visedgesfromr = uint8(ones(nedges,nreceivers));

if ~isempty(ivec_visR1)
	if nreceivers < 256
        colnumb = uint8(ceil(ivec_visR1/nplanes));             % This is the receiver number
	elseif nreceivers < 65536
        colnumb = uint16(ceil(ivec_visR1/nplanes));             % This is the receiver number
	else
        colnumb = uint32(ceil(ivec_visR1/nplanes));             % This is the receiver number    
	end
	if nplanes < 256
        rownumb = uint8(ivec_visR1 - (double(colnumb)-1)*nplanes);     % This is the plane number
	elseif nplanes < 65536
        rownumb = uint16(ivec_visR1 - (double(colnumb)-1)*nplanes);     % This is the plane number
	else
        rownumb = uint32(ivec_visR1 - (double(colnumb)-1)*nplanes);     % This is the plane number
    end
    
    % Divide these lists into two categories: thin planes and non-thin
    % planes.
    
    iv2 = find(planeisthin(rownumb));
    if ~isempty(iv2)
        colnumb_thin = colnumb(iv2);            
        rownumb_thin = rownumb(iv2);
        colnumb(iv2) = [];
        rownumb(iv2) = [];
    else
        colnumb_thin = [];
        rownumb_thin = [];
    end
    
    if ~isempty(colnumb)
        % Select all the edges that are connected to these planes
	
        if nedges < 256
            edgeselection = uint8(edgesatplane(rownumb,1:maxncornersperplane));
        elseif nedges < 65536
            edgeselection = uint16(edgesatplane(rownumb,1:maxncornersperplane));
        else
            edgeselection = uint32(edgesatplane(rownumb,1:maxncornersperplane));
        end
        if length(rownumb) > 1
            maxcols = sum(sign(sum(double(edgeselection))));
        else
            maxcols = sum(sign(double(edgeselection)));    
        end
        edgeselection = edgeselection(:,1:maxcols);
        colnumb = colnumb(:,ones(1,maxcols));
        
        indexvec = uint32(double(edgeselection) + (double(colnumb)-1)*nedges);
        indexvec(find(edgeselection==0)) = [];
        
        visedgesfromr(indexvec) = uint8(5*ones(size(indexvec)));

    end
        
    if ~isempty(colnumb_thin)
        % Select all the edges that are connected to these planes
	
        if nedges < 256
            edgeselection = uint8(edgesatplane(rownumb_thin,1:maxncornersperplane));
        elseif nedges < 65536
            edgeselection = uint16(edgesatplane(rownumb_thin,1:maxncornersperplane));
        else
            edgeselection = uint32(edgesatplane(rownumb_thin,1:maxncornersperplane));
        end
        if length(rownumb_thin) > 1
            maxcols = sum(sign(sum(edgeselection)));
        else
            maxcols = sum(sign(edgeselection));    
        end
        edgeselection = edgeselection(:,1:maxcols);
        colnumb_thin = colnumb_thin(:,ones(1,maxcols));
 
        indexvec = uint32(double(edgeselection) + (double(colnumb_thin)-1)*nedges);
    	indexvec(find(edgeselection==0)) = [];
        
        visedgesfromr(indexvec) = uint8(6*ones(size(indexvec)));
    end
end

%       4    Edge belongs to a plane which the S/R is inside
%            and the plane is totabs
%       3    Edge belongs to a plane which the S/R is inside
%            and the plane has indents
%       2    Edge belongs to a plane which the S/R is inside
%            and the plane has no indents
%
% These can be derived from the cases where visplanesfromr = 4 and 5

ivec = find(visplanesfromr==5);
if ~isempty(ivec)
    colnumb = ceil(ivec/nplanes);               % This is the receiver number
	if nreceivers < 256
        colnumb = uint8(ceil(ivec/nplanes));             % This is the receiver number
	elseif nreceivers < 65536
        colnumb = uint16(ceil(ivec/nplanes));             % This is the receiver number
	else
        colnumb = uint32(ceil(ivec/nplanes));             % This is the receiver number    
	end
    rownumb = ivec - (double(colnumb)-1)*nplanes;       % This is the plane number
	if nplanes < 256
        rownumb = uint8(ivec - (double(colnumb)-1)*nplanes);     % This is the plane number
	elseif nplanes < 65536
        rownumb = uint16(ivec - (double(colnumb)-1)*nplanes);     % This is the plane number
	else
        rownumb = uint32(ivec - (double(colnumb)-1)*nplanes);     % This is the plane number
	end
    clear ivec
    
    % Select all the edges that are connected to these planes
	
    if nedges < 256
        edgeselection = uint8(edgesatplane(rownumb,1:maxncornersperplane));
    elseif nedges < 65536
        edgeselection = uint16(edgesatplane(rownumb,1:maxncornersperplane));
    else
        edgeselection = uint32(edgesatplane(rownumb,1:maxncornersperplane));
    end
    if length(rownumb) > 1
        maxcols = sum(sign(sum(edgeselection)));
    else
        maxcols = sum(sign(edgeselection));    
    end
    edgeselection = edgeselection(:,1:maxcols);
    colnumb = colnumb(:,ones(1,maxcols));
    
    indexvec = uint32(double(edgeselection) + (double(colnumb)-1)*nedges);
    indexvec(find(edgeselection==0)) = [];
    
    visedgesfromr(indexvec) = uint8(4*ones(size(indexvec)));
end

ivec = find(visplanesfromr==4);    
if ~isempty(ivec)
	if nreceivers < 256
        colnumb = uint8(ceil(ivec/nplanes));             % This is the receiver number
	elseif nreceivers < 65536
        colnumb = uint16(ceil(ivec/nplanes));             % This is the receiver number
	else
        colnumb = uint32(ceil(ivec/nplanes));             % This is the receiver number    
	end
	if nplanes < 256
        rownumb = uint8(ivec - (double(colnumb)-1)*nplanes);     % This is the plane number
	elseif nplanes < 65536
        rownumb = uint16(ivec - (double(colnumb)-1)*nplanes);     % This is the plane number
	else
        rownumb = uint32(ivec - (double(colnumb)-1)*nplanes);     % This is the plane number
	end
    clear ivec
    
    % Divide these lists into two categories: planes with indents and 
    % planes without
    
    iv2 = find(planehasindents(rownumb));
    if ~isempty(iv2)
        colnumb_ind = colnumb(iv2);            
        rownumb_ind = rownumb(iv2);
        colnumb(iv2) = [];
        rownumb(iv2) = [];
    else
        colnumb_ind = [];
        rownumb_ind = [];
    end
    if ~isempty(colnumb)
        % Select all the edges that are connected to these planes
	
    %    edgeselection = edgesatplane(rownumb,1:maxncornersperplane);
        if nedges < 256
            edgeselection = uint8(edgesatplane(rownumb,1:maxncornersperplane));
        elseif nedges < 65536
            edgeselection = uint16(edgesatplane(rownumb,1:maxncornersperplane));
        else
            edgeselection = uint32(edgesatplane(rownumb,1:maxncornersperplane));
        end
        if length(rownumb) > 1
            maxcols = sum(sign(sum(double(edgeselection))));
        else
            maxcols = sum(sign(double(edgeselection)));    
        end
        edgeselection = edgeselection(:,1:maxcols);
        colnumb = colnumb(:,ones(1,maxcols));
        
        indexvec = uint32(double(edgeselection) + (double(colnumb)-1)*nedges);
        indexvec(find(edgeselection==0)) = [];
        
        visedgesfromr(indexvec) = uint8(2*ones(size(indexvec)));
    end
    if ~isempty(colnumb_ind)
        % Select all the edges that are connected to these planes
	
        edgeselection = edgesatplane(rownumb,1:maxncornersperplane);
        if length(rownumb) > 1
            maxcols = sum(sign(sum(edgeselection)));
        else
            maxcols = sum(sign(edgeselection));    
        end
        edgeselection = edgeselection(:,1:maxcols);
        colnumb_ind = colnumb_ind(:,ones(1,maxcols));
        
        indexvec = uint32(edgeselection + (colnumb_ind-1)*nedges);
        indexvec(find(edgeselection==0)) = [];
        
        visedgesfromr(indexvec) = uint8(3*ones(size(indexvec)));
    end
end

%       0   R can never see the edge because it is behind the edge-planes.
%           NB! If the closwedang > pi, then it is enough if R is behind
%           one plane. If the closwedang <= pi, then R must be behind both
%           planes in order not to see the edge.

ivec = find(visplanesfromr==0);

if ~isempty(ivec)
 	if nreceivers < 256
        colnumb = uint8(ceil(ivec/nplanes));             % This is the receiver number
	elseif nreceivers < 65536
        colnumb = uint16(ceil(ivec/nplanes));             % This is the receiver number
	else
        colnumb = uint32(ceil(ivec/nplanes));             % This is the receiver number    
	end
	if nplanes < 256
        rownumb = uint8(ivec - (double(colnumb)-1)*nplanes);     % This is the plane number
	elseif nplanes < 65536
        rownumb = uint16(ivec - (double(colnumb)-1)*nplanes);     % This is the plane number
	else
        rownumb = uint32(ivec - (double(colnumb)-1)*nplanes);     % This is the plane number
	end
    clear ivec
    
    % Select all the edges that are connected to these planes

    if nedges < 256
        edgeselection = uint8(edgesatplane(rownumb,1:maxncornersperplane));        
    elseif nedges < 65536
        edgeselection = uint16(edgesatplane(rownumb,1:maxncornersperplane));
    else
        edgeselection = uint32(edgesatplane(rownumb,1:maxncornersperplane));
    end
    if length(rownumb) > 1
        maxcols = sum(sign(sum(double(edgeselection))));
    else
        maxcols = sum(sign(double(edgeselection)));    
    end
    edgeselection = edgeselection(:,1:maxcols);
    colnumb = colnumb(:,ones(1,maxcols));
    ntot2 = length(rownumb)*maxcols;
    
    indexvec = uint32(double(edgeselection) + (double(colnumb)-1)*nedges);
    indexvec = reshape(indexvec,ntot2,1);
    edgeselection = reshape(edgeselection,ntot2,1);

    iv2 = find(edgeselection==0);    
    indexvec(iv2) = [];
    edgeselection(iv2) = [];
    ntot2 = ntot2 - length(iv2);
    
    [indexvec,sortvec] = sort(indexvec);
    edgeselection = edgeselection(sortvec);
    clear sortvec
    
    markdoubles = [edgeselection(1:ntot2-1)==edgeselection(2:ntot2)];
    markdoubles = sign([markdoubles;0] + [0;markdoubles]);
    iv2 = uint32(find(markdoubles));
    clear markdoubles
    visedgesfromr(indexvec(iv2)) = 0;
    indexvec(iv2) = [];
    edgeselection(iv2) = [];
    iv2 = uint32(find(closwedangvec(edgeselection)<=pi));
    indexvec(iv2) = [];
    visedgesfromr(indexvec) = 0;
    clear indexvec
end

clear edgeselection iv2

%##################################################################
%
% Mask out edges that should be switched off

mask = (ones(nedges,1));
mask(offedges) = mask(offedges) & 0;

visedgesfromr = uint8(double(visedgesfromr).*mask(:,onesvec1R));

clear onesvec1R

%##################################################################
%##################################################################
%##################################################################
%
%		CHECK OBSTRUCTION OF R-TO-EDGE PATHS
%
%   vispartedgesfromr   0-2^nedgesubs-1, [nedges,nrec]
%                       telling which segments of an edge that
%                       are visible
%
% A. If visedgesfromr = 6 then check aligned-plane obstructions
%           If OK,     set visedgesfromr = 1 (check ordinary obst.)
%           If not OK, set visedgesfromr = 0
%
% B. If visedgesfromr = 3 then check within-plane obstructions
%           Set visedgesfromr = 0 and
%           set vispartedgesfromr = 0-full
%
% C. If visedgesfromr = 2 then
%           Set visedgesfromr = 0 and
%           set vispartedgesfromr = Full
%
% D. If visedgesfromr = 1 then check ordinary obstructions
%           Set visedgesfromr = 0 and
%           set vispartedgesfromr = 0-full

edgedividervec = [0:nedgesubs-1].';
weightvec = 2.^edgedividervec;
maxvisibilityval = 2^nedgesubs-1;

if SHOWTEXT >= 3
    if lower(typeofcoords(1)) == 'r'
    	disp(['         Checking for obstructing planes between R and edges'])
    else
    	disp(['         Checking for obstructing planes between S and edges'])
    end
end

if nedgesubs<8
    vispartedgesfromr = uint8(zeros(size(visedgesfromr)));
elseif nedgesubs<16
    vispartedgesfromr = uint16(zeros(size(visedgesfromr)));
else
    vispartedgesfromr = uint32(zeros(size(visedgesfromr)));
end

%##################################################################
%
% The receiver related quantities

iv = uint32(find(visedgesfromr>=5));
if ~isempty(iv)
    disp('      We check aligned-edge obstructions')

    % All edges must be checked vs all other edges?
    % Make subdivision of edges. If a line from R to edge segments
    % pass through a plane that is constructed by an edge and perpendicular to
    % the studied plane, then that edge is obstructing. 
    
    combnumbers = double(iv);
	recnumbers = ceil(combnumbers/nedges);
    if nedges < 256
		edgenumbers = uint8(combnumbers - (recnumbers-1)*nedges);
    elseif nedges < 65536
		edgenumbers = uint16(combnumbers - (recnumbers-1)*nedges);
    else
   		edgenumbers = uint32(combnumbers - (recnumbers-1)*nedges);
    end
    
    nedgesperreceiver = histc(recnumbers,[1:nreceivers]);
    iv1 = find(nedgesperreceiver >= 2);
    if isempty(iv1)
        visedgesfromr(iv)    = ones(size(iv));
    else
        disp(['      Several edges are in-plane with the receiver'])
        disp(['      If possible, move the receivers a little bit'])
        [uniquerecs,exampleindex] = unique(recnumbers);
        for jj = 1:length(uniquerecs)
           problemrec = uniquerecs(jj);
           iv1 = find(recnumbers==problemrec);
           disp(' ')
           disp('      Receiver with coordinates')
           disp(['      ',num2str(pointcoords(problemrec,1)),'  ',num2str(pointcoords(problemrec,2)),'  ',num2str(pointcoords(problemrec,3))])
           disp('      is aligned with edges')
           numvec = [int2str(edgenumbers(iv1(1))),' ',int2str(edgenumbers(iv1(2)))];
           for kk = 3:length(iv1)
              numvec = [numvec,' ', int2str(edgenumbers(iv1(kk)))];
           end
           disp(['      ',numvec])                        
        end
        error('ERROR: Obstruction check of aligned-edges not implemented yet!')        
    end
 end

iv = uint32(find(visedgesfromr==3));
if ~isempty(iv)
    disp('      We check within-plane obstructions')

    visedgesfromr(iv)      = zeros(size(iv));
    error('ERROR: Obstruction check of within-same-plane-edges not implemented yet!')        
end

iv = uint32(find(visedgesfromr==2));
if ~isempty(iv)
    disp('      Edge fully visible')
    visedgesfromr(iv)      = zeros(size(iv));
    vispartedgesfromr(iv)  =  maxvisibilityval*ones(size(iv));
    
end

iv = uint32(find(visedgesfromr==4));
if ~isempty(iv)
    disp('      Edge in totabs plane')
    visedgesfromr(iv)      = zeros(size(iv));
    vispartedgesfromr(iv)  =  zeros(size(iv));
    
end

% Below is the main part of the obstruction test. In the matrix
% visedgesfromr, values 1 indicate that edge-receiver combos should be
% tested. visedgesfromr has the size [nedges,nreceivers].

iv = uint32(find(visedgesfromr==1));
if ~isempty(iv)
            
    visedgesfromr(iv) = zeros(size(iv));
    combnumbers = double(iv);
    clear iv
    
    % combnumbers is a list of the "combined" values for edge and receiver. 
    
    if nreceivers < 256
    	recnumbers = uint8(ceil(combnumbers/nedges));
    elseif nreceivers < 65536
    	recnumbers = uint16(ceil(combnumbers/nedges));
    else
    	recnumbers = uint32(ceil(combnumbers/nedges));
    end
    if nedges < 256
		edgenumbers = uint8(combnumbers - (double(recnumbers)-1)*nedges);
    elseif nedges < 65536
		edgenumbers = uint16(combnumbers - (double(recnumbers)-1)*nedges);            
    else
		edgenumbers = uint32(combnumbers - (double(recnumbers)-1)*nedges);            
    end
    combnumbers = uint32(combnumbers);
    ncombs = length(recnumbers);
    ntot = ncombs*nplanes;
    
    % Mark, in a big matrix, which planes can at all obstruct.
    % Planes that can obstruct must be seen by the receiver or the edge but
    % not both - because then they would be on the same side!
    % NB! Also totabs planes can obstruct.
    % The big matrix will have nedges*nplanes rows. The vector planenumb
    % will have values such as [1 1 1 1... 2 2 2 2 .....]

    iv = [1:ncombs*nplanes].';                
    % Plane number is given by the col no.
    if nplanes < 256
        planenumb = uint8(ceil(iv/ncombs));
    elseif nplanes < 65536
        planenumb = uint16(ceil(iv/ncombs));
    else
        planenumb = uint32(ceil(iv/ncombs));        
    end
    indexvec = planenumb;
    
    % The rownumber will simply have values [1 2 3... N 1 2 3 ... N 1 2 3
    % ... N] where N is the number of combinations.
    if ncombs < 256
        rownumber = uint8(iv - (double(planenumb)-1)*ncombs);
    elseif ncombs < 65536
        rownumber = uint16(iv - (double(planenumb)-1)*ncombs);
    else
        rownumber = uint32(iv - (double(planenumb)-1)*ncombs);
    end
    clear planenumb iv 
    
    % The peak in the memory need might be after the two indexvec lines
    % below

    % indexvec2 will point to the right location in an [nplanes,nedges] matrix
    % indexvec will point to the right location in an [nplanes,nreceivers] matrix
    indexvec2 = uint32(double(indexvec) + (double(edgenumbers(rownumber))-1)*nplanes);
    indexvec = uint32(double(indexvec) + (double(recnumbers(rownumber))-1)*nplanes);
    clear rownumber

    % The big matrix checkplane, size [1,ncombs*nplanes], will contain the
    % value 1 if a plane should be checked for obstruction.
    % The index order of checkplane is:
    % [Plane1&comb1 Plane1&comb2 Plane1&comb3 ...] where comb1 is the first
    % combination of receiver and edge in the recnumbers & edgenumbers
    % lists.
    % Only planes that are seen by the S/R *or* the edge, but not both
    % should be checked for obstruction!!
    % We remove combinations where:
    %   ... S/R is aligned with a plane (visplanesfromr == 1)
    %   ... edge belongs to a plane or is aligned with a plane (edgeseesplane <0)
    %   ... S/R is behind and edge is behind a plane (visplanesfromr == 0 & edgeseesplane == 0)
    %   ... S/R is in front of and edge is in front of a plane ( (visplanesfromr == 2 | visplanesfromr == 3) & edgeseesplane == 1)
    %
    % Comment 050116 (PS):
    %   We would actually need a matrix called planeseesedge. Here we use
    %   edgeseesplane instead, and it is not completely clear if that works
    %   fine. 
    
    checkplane = (visplanesfromr(indexvec)~=1) & (edgeseesplane(indexvec2)>=0) & not(visplanesfromr(indexvec)==0 & edgeseesplane(indexvec2)==0 ) & not( (visplanesfromr(indexvec)==2 | visplanesfromr(indexvec)==3) & edgeseesplane(indexvec2)==1 );
    if size(checkplane,1) > size(checkplane,2)
        checkplane = checkplane.';    
    end

    if size(checkplane,1) == 1 | size(checkplane,2) == 1
        checkplane = reshape(checkplane,ncombs,nplanes);
    end

    clear indexvec indexvec2 edgeseesplane

    % If there are some R-edge combos that have no planes to check
    % obstruction for, we can mark those combos ("combsnottocheck") as fully visible.
    %
    % The remaining ones ("combstocheck") still need to be checked.
    
    [n1,n2] = size(checkplane);
    if n1 < 65536
        combstocheck = uint16(find( sum(checkplane.') ));        
    elseif n1 < 4e9
        combstocheck = uint32(find( sum(checkplane.') ));                    
    end
    combsnottocheck = [1:ncombs];
    combsnottocheck(combstocheck) = [];
    if ~isempty(combsnottocheck)
        vispartedgesfromr(combnumbers(combsnottocheck)) = maxvisibilityval*ones(size(combsnottocheck));
    end

    if ~isempty(combstocheck)
        
        checkplane = checkplane(combstocheck,:);
        recnumbers = recnumbers(combstocheck);
        maxrec = max(recnumbers);
        if maxrec < 256
            recnumbers = uint8(recnumbers);        
        elseif maxrec < 65536
            recnumbers = uint16(recnumbers);                
        else
            recnumbers = uint32(recnumbers);                
        end
        edgenumbers = edgenumbers(combstocheck);
        combnumbers = combnumbers(combstocheck);
        ncombs = length(combstocheck);

        % Now, checkplane is a matrix of size [ncombs,nplanes] where each
        % row corresponds to one path (from one receiver to one edge) that needs
        % an obstruction check. For that row, checkplane has the value 1 for
        % each plane that needs to be checked.
        %    
        % Expand all edges into their edge segment/subdivision
        % We treat all the little edge subdivisions as separate edges
	
        nposs = ncombs*nedgesubs;        
        
        expandedrecnumbers = recnumbers(:,onesvec1ES);
        clear recnumbers
        expandedrecnumbers = reshape(expandedrecnumbers.',nposs,1);
        expandedcombnumbers = combnumbers(:,onesvec1ES);
        clear combnumbers
        expandedcombnumbers = reshape(expandedcombnumbers.',nposs,1);
        if nposs < 65536
            okcombs = uint16(zeros(size(expandedrecnumbers)));
        elseif nposs < 4e9
            okcombs = uint32(zeros(size(expandedrecnumbers)));                
        end
        
         [tocoords,expandededgeweightlist,expandededgenumbers] = EDB2getedgepoints(edgestartcoords(edgenumbers,:),...
             edgeendcoords(edgenumbers,:),edgelengthvec(edgenumbers,:),nedgesubs,1);
         expandededgenumbers = edgenumbers(expandededgenumbers);

		[nonobstructedpaths,nobstructions,edgehits,cornerhits] = EDB2checkobstr_pointtoedge(pointcoords,expandedrecnumbers,tocoords,reshape(repmat(checkplane.',[nedgesubs,1]),nplanes,nposs).',planeseesplane,...
            planeeqs,planenvecs,minvals,maxvals,planecorners,corners,ncornersperplanevec,rearsideplane);        
        clear checkplane planeseesplane

        expandedcombnumbers = expandedcombnumbers(nonobstructedpaths);
        expandededgeweightlist = expandededgeweightlist(nonobstructedpaths);
        expandedrecnumbers = expandedrecnumbers(nonobstructedpaths);
        expandededgenumbers = expandededgenumbers(nonobstructedpaths);

        % Pack all non-obstructed edge segments together and add their weights together
        
        test = [double(expandedrecnumbers) double(expandededgenumbers)];
		ncombs = length(expandedrecnumbers);
		dtest = diff([0;prod(   double(test(1:ncombs-1,:)==test(2:ncombs,:)).'  ).']);
		ivremove = find(dtest==1);
		
		while ~isempty(ivremove)
            expandededgeweightlist(ivremove+1) = double(expandededgeweightlist(ivremove+1)) + double(expandededgeweightlist(ivremove));
            expandededgeweightlist(ivremove) = [];
            expandedrecnumbers(ivremove) = [];
            expandededgenumbers(ivremove,:) = [];
		
            test = [double(expandedrecnumbers) double(expandededgenumbers)];
            ncombs = length(expandedrecnumbers);
            dtest = diff([0;prod(   double(test(1:ncombs-1,:)==test(2:ncombs,:)).'  ).']);
            ivremove = find(dtest==1);
		end
         
        indexvec = uint32(nedges*(double(expandedrecnumbers)-1)+double(expandededgenumbers));
        vispartedgesfromr(indexvec)  =  expandededgeweightlist;
        clear indexvec                
    end
end

%----------------------------------------------------------------------------
%
%		SAVE THE VARIABLES
%
%--------------------------------------------------------------------------

if typeofcoords=='r'
    receivers = pointcoords;
    Varlist = ['  visplanesfromr  receivers '];
    Varlist = [Varlist,'  vispartedgesfromr  routsidemodel '];
    eval(['save ',outputfile,Varlist])
    return
else
    sources = pointcoords;
    visplanesfroms = visplanesfromr;
    vispartedgesfroms = vispartedgesfromr;
    soutsidemodel = routsidemodel;
    Varlist = ['  visplanesfroms  sources '];
    Varlist = [Varlist,'  vispartedgesfroms  soutsidemodel '];
    eval(['save ',outputfile,Varlist])
    return
end
