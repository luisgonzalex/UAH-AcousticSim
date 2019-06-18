function ISEStreefile = EDB2findISEStree(eddatafile,S,isou,specorder,difforder,visplanesfromS,vispartedgesfromS,nedgesubs)
% EDB2findISEStree - Constructs a list ("an ISES tree") of all possible specular-diffraction combinations.
% EDB2findISEStree constructs a list ("an ISES tree") of all possible specular-diffraction combinations that
% are possible for the given source. This list is based on the visibility
% matrices from source to planes/edges and between planes/edges. Visibility
% checks are employed wherever possibble. NB! Visibility checks that depend
% on the receiver position are not used at this stage.
% Planes are numbered 1-nplanes, and edges have numbers
% [nplanes+1:nplanes+nedges].
%
% Input parameters:
% 	eddatafile                  See a description in EDB2edgeo, EDB2srgeo and EDB2mainISES
% 	S
% 	isou
% 	specorder
% 	difforder
% 	visplanesfromS
% 	vispartedgesfromS
% 	nedgesubs
%
% GLOBAL parameters
%   SHOWTEXT JJ JJnumbofchars   See EDB2mainISES
%   SUPPRESSFILES               If this global parameter has the value 1
%                               then the results from this function will be returned in a struct
%                               rather than in a file.
%
% Output parameters:
%   ISEStreefile                The name of the output file which contains
%                               the output data below.
%
% Global output data:
% 	POTENTIALISES               A matrix, [npossiblecombs,specorder], of all the 
% 	                            possible specular-diffraction combos that are possible for the
% 	                            given source. The combos are denoted by
% 	                            plane numbers and edge numbers, but edge
% 	                            numbers have a constant number (=nplanes)
% 	                            added to them, i.e., if there are 20 planes
% 	                            in the model, edge number 1 will be denoted
% 	                            by number 21.
% 	ORIGINSFROM                 A list, [npossiblecombs,1], where the value
% 	                            in row N states which other row in ISCOORDS that
% 	                            the image source in row N originates from.
% 	ISCOORDS                    A matrix, [npossiblecombs,3], containing
% 	                            the image source coordinates, where this is applicable
%                               i.e., where the combo starts with a specular reflection. 
% 	ISESVISIBILITY              A list, [npossiblecombs,1], containing the visibility of
%                               the first edge in a sequence, seen from the source. 
% 	IVNSPECMATRIX               A matrix, [nspeccombs,specorder], where each column contains
%                               the row numbers in POTENTIALISES that contain combos with
%                               1 or 2 or... or "specorder" purely specular reflections.
% 	REFLORDER                   A list, [npossiblecombs,1], containing the
% 	                            order of reflection for each row of POTENTIALISES
% 	IVNDIFFMATRIX               A matrix, [ndiffcombs,specorder], where each column contains
%                               the row numbers in POTENTIALISES that contain combos with
%                               1 or 2 or... or "specorder" diffractions.
% Output data in the ISEStreefile:
% 	lengthNspecmatrix           A list, [1,specorder], with the number of
% 	                            entries in each column of IVNSPECMATRIX.
% 	lengthNdiffmatrix           A list, [1,specorder], with the number of
% 	                            entries in each column of IVNDIFFMATRIX.
% 	singlediffcol               A list, [nfirstorderdiff,1], containing the column number
%                               that the first-order diffracted component can be found in (in
% 	                            POTENTIALISES).
% 	startindicessinglediff      A list, [specorder,1], where the value in position N contains
%                               which row in IVNDIFFMATRIX(:,1) that has the first combination
%                               of order N and with one diffraction component.
% 	endindicessinglediff        A list, [specorder,1], where the value in position N contains
%                               which row in IVNDIFFMATRIX(:,1) that has the last combination
%                               of order N and with one diffraction 
% 	ndecimaldivider             See below.
% 	PointertoIRcombs            A sparse list which for position N contains a row number
%                               (in POTENTIALISES) that has a combination
%                               that ends with a specular reflection in
%                               plane N, after a diffraction. A row with one of the specular
%                               reflections M,N after a diffraction can be
%                               found in position M*ndecimaldivider+N etc
%                               for higher orders.
% 	IRoriginsfrom               A list, [npossiblecombs,1], where the value in position N
%                               states which other row number (in POTENTIALISES) that the
%                               image receiver in row N origins from.
%
% Uses functions EDB2strpend EDB2findis EDB2getedgepoints EDB2chkISvisible EDB2checkobstrpaths
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
% Peter Svensson (svensson@iet.ntnu.no) 20130813
%
% ISEStreefile = EDB2findISEStree(eddatafile,S,isou,specorder,difforder,visplanesfroms,vispartedgesfroms,nedgesubs);

% 20100210 Functioning version
% 20130813 Fixed a bug which gave errors for specular order > 4 in models with only two
%		   walls.

global SHOWTEXT JJ JJnumbofchars
global POTENTIALISES ISCOORDS IVNDIFFMATRIX
global IVNSPECMATRIX ORIGINSFROM ISESVISIBILITY REFLORDER SUPPRESSFILES

GEOMACC = 1e-10;

eval(['load ',eddatafile])
clear cornerinfrontofplane

[filepath,filestem,fileext] = fileparts(eddatafile);
ISEStreefile = [[filepath,filesep],EDB2strpend(filestem,'_eddata'),'_',int2str(isou),'_ISEStree.mat'];

ncorners = size(corners,1);
nplanes = length(planeisthin);
nedges = length(closwedangvec);
nhighval = nplanes + nedges;

onesvec1 = ones(1,nedgesubs);
obstructtestneeded = (sum(canplaneobstruct)~=0);
onesvec3 = ones(1,3);

doconetrace = 0;
if SHOWTEXT >= 3
    disp('WARNING: doconetrace is switched off')    
end

%--------------------------------------------------------------------------
% We modify edgeseesplane because in EDB2edgeo, edgeseesplane was set to 1
% for totabs planes (this was necessary for EDB2srgeo).

% Set all edges belonging to planes with reflfactors = 0 to edgeseesplane = 0.

listofabsorbingplanes = find(reflfactors == 0);

if ~isempty(listofabsorbingplanes)
	edgeseesplane(listofabsorbingplanes,:) = zeros(length(listofabsorbingplanes),nedges);
end

% If we should do cone tracing, then we need plane midpoints and plane
% sizes.

if doconetrace == 1
    planemidpoints = zeros(nplanes,3);
    maxdisttocorner = zeros(nplanes,1);

    iv4planes = find(ncornersperplanevec==4);
    if any(iv4planes)
        xvalues = corners(planecorners(iv4planes,1:4).',1);
        xvalues = reshape(xvalues,4,length(xvalues)/4);
        yvalues = corners(planecorners(iv4planes,1:4).',2);
        yvalues = reshape(yvalues,4,length(yvalues)/4);
        zvalues = corners(planecorners(iv4planes,1:4).',3);
        zvalues = reshape(zvalues,4,length(zvalues)/4);
        planemidpoints(iv4planes,:) = [mean(xvalues).' mean(yvalues).' mean(zvalues).'];
        
        dist1 = sqrt(sum((corners(planecorners(iv4planes,1),:) - planemidpoints(iv4planes,:)).'.^2)).';
        dist2 = sqrt(sum((corners(planecorners(iv4planes,2),:) - planemidpoints(iv4planes,:)).'.^2)).';
        dist3 = sqrt(sum((corners(planecorners(iv4planes,3),:) - planemidpoints(iv4planes,:)).'.^2)).';
        dist4 = sqrt(sum((corners(planecorners(iv4planes,4),:) - planemidpoints(iv4planes,:)).'.^2)).';
        maxdisttocorner(iv4planes) = max([dist1 dist2 dist3 dist4].');
        
    end

    iv3planes = find(ncornersperplanevec==3);
    if any(iv3planes)
        xvalues = corners(planecorners(iv3planes,1:3).',1);
        xvalues = reshape(xvalues,3,length(xvalues)/3);
        yvalues = corners(planecorners(iv3planes,1:3).',2);
        yvalues = reshape(yvalues,3,length(yvalues)/3);
        zvalues = corners(planecorners(iv3planes,1:3).',3);
        zvalues = reshape(zvalues,3,length(zvalues)/3);
        planemidpoints(iv3planes,:) = [mean(xvalues).' mean(yvalues).' mean(zvalues).'];

        dist1 = sqrt(sum((corners(planecorners(iv3planes,1),:) - planemidpoints(iv3planes,:)).'.^2)).';
        dist2 = sqrt(sum((corners(planecorners(iv3planes,2),:) - planemidpoints(iv3planes,:)).'.^2)).';
        dist3 = sqrt(sum((corners(planecorners(iv3planes,3),:) - planemidpoints(iv3planes,:)).'.^2)).';
        maxdisttocorner(iv3planes) = max([dist1 dist2 dist3].');
    
    end

    iv5planes = find(ncornersperplanevec==5);
    if any(iv5planes)
        xvalues = corners(planecorners(iv5planes,1:5).',1);
        xvalues = reshape(xvalues,5,length(xvalues)/5);
        yvalues = corners(planecorners(iv5planes,1:5).',2);
        yvalues = reshape(yvalues,5,length(yvalues)/5);
        zvalues = corners(planecorners(iv5planes,1:5).',3);
        zvalues = reshape(zvalues,5,length(zvalues)/5);
        planemidpoints(iv5planes,:) = [mean(xvalues).' mean(yvalues).' mean(zvalues).'];
    
        dist1 = sqrt(sum((corners(planecorners(iv5planes,1),:) - planemidpoints(iv5planes,:)).'.^2)).';
        dist2 = sqrt(sum((corners(planecorners(iv5planes,2),:) - planemidpoints(iv5planes,:)).'.^2)).';
        dist3 = sqrt(sum((corners(planecorners(iv5planes,3),:) - planemidpoints(iv5planes,:)).'.^2)).';
        dist4 = sqrt(sum((corners(planecorners(iv5planes,4),:) - planemidpoints(iv5planes,:)).'.^2)).';
        dist5 = sqrt(sum((corners(planecorners(iv5planes,5),:) - planemidpoints(iv5planes,:)).'.^2)).';
        maxdisttocorner(iv5planes) = max([dist1 dist2 dist3 dist4 dist5].');

    end

    iv6planes = find(ncornersperplanevec==6);
    if any(iv6planes)
        xvalues = corners(planecorners(iv6planes,1:6).',1);
        xvalues = reshape(xvalues,6,length(xvalues)/6);
        yvalues = corners(planecorners(iv6planes,1:6).',2);
        yvalues = reshape(yvalues,6,length(yvalues)/6);
        zvalues = corners(planecorners(iv6planes,1:6).',3);
        zvalues = reshape(zvalues,6,length(zvalues)/6);
        planemidpoints(iv6planes,:) = [mean(xvalues).' mean(yvalues).' mean(zvalues).'];

        dist1 = sqrt(sum((corners(planecorners(iv6planes,1),:) - planemidpoints(iv6planes,:)).'.^2)).';
        dist2 = sqrt(sum((corners(planecorners(iv6planes,2),:) - planemidpoints(iv6planes,:)).'.^2)).';
        dist3 = sqrt(sum((corners(planecorners(iv6planes,3),:) - planemidpoints(iv6planes,:)).'.^2)).';
        dist4 = sqrt(sum((corners(planecorners(iv6planes,4),:) - planemidpoints(iv6planes,:)).'.^2)).';
        dist5 = sqrt(sum((corners(planecorners(iv6planes,5),:) - planemidpoints(iv6planes,:)).'.^2)).';
        dist6 = sqrt(sum((corners(planecorners(iv6planes,6),:) - planemidpoints(iv6planes,:)).'.^2)).';
        maxdisttocorner(iv6planes) = max([dist1 dist2 dist3 dist4 dist5 dist6].');
    
    end

end

%--------------------------------------------------------------------------
% Combine planeseesplane with edgeseesplane
% and visplanesfroms with vispartedgesfroms
% so that edges appear as planes

if difforder > 0
    visPLANESfroms = [visplanesfromS;2*sign(double(vispartedgesfromS))];
    % source sees planes that have visplanesfroms:
    %   2 => in front of reflective planes
    %   4 => inside plane which is reflective
    %   old version:  sousees = (visPLANESfroms==1 | visPLANESfroms==-2);
    sousees = (visPLANESfroms==2 | visPLANESfroms==4);

    if exist('edgeseespartialedge') ~= 1
        edgeseesedge = int8(zeros(nedges,nedges));    
    else
        if ~isempty(edgeseespartialedge)

            edgeseesedge = int8(full(abs(double(edgeseespartialedge))>0));
            
        else
            edgeseesedge = int8(zeros(nedges,nedges));    
        end
    end

    PLANEseesPLANE = [[planeseesplane edgeseesplane];[edgeseesplane.' edgeseesedge]];
    
    clear edgeseesplane edgeseesedge
else
    visPLANESfroms = visplanesfromS;
    PLANEseesPLANE = planeseesplane;
    sousees = (visPLANESfroms==2 | visPLANESfroms==4);

    clear visplanesfromS
end

startindices = zeros(specorder,1);
endindices = zeros(specorder,1);

%--------------------------------------------------------------------------
% Construct a list of which planes a plane sees so that a search isn't
% needed later.

nPLANES = size(PLANEseesPLANE,1);
listofvisPLANES = zeros(nPLANES,nPLANES-1);    
nvisPLANES = zeros(nPLANES,1);
for ii = 1:nPLANES
    visPLANES = find(PLANEseesPLANE(ii,:)==1);
    nvisPLANES(ii) = length(visPLANES);
    listofvisPLANES(ii,1:nvisPLANES(ii)) = visPLANES;
end

%##################################################################
%##################################################################
%##################################################################
%
%         First order
%
%------------------------------------------------------------------

startindices(1) = 1;

possiblefirsts  =  find( sousees );
npossiblefirsts = length(possiblefirsts);

if nhighval < 255
    POTENTIALISES = uint8( [[possiblefirsts ] zeros(npossiblefirsts,specorder-1)] );
else
    POTENTIALISES = uint16( [[possiblefirsts ] zeros(npossiblefirsts,specorder-1)] );
end

ORIGINSFROM = uint32(zeros(npossiblefirsts,1));
if nedgesubs < 8 
    ISESVISIBILITY = uint8(zeros(npossiblefirsts,1));
elseif nedgesubs < 16
    ISESVISIBILITY = uint16(zeros(npossiblefirsts,1));
else    
    ISESVISIBILITY = uint32(zeros(npossiblefirsts,1));
end
iv = find(possiblefirsts>nplanes);
if ~isempty(iv)
    ISESVISIBILITY(iv) = vispartedgesfromS(possiblefirsts(iv)-nplanes);    
end

endindices(1) = npossiblefirsts;

% Compute the IS coords

ivdiff = [startindices(1):endindices(1)];
ivspec = find(POTENTIALISES(ivdiff,1)<=nplanes);
ivdiff(ivspec) = [];

ISCOORDS = zeros(npossiblefirsts,3);
ISCOORDS(ivspec,:) = EDB2findis(S,POTENTIALISES(ivspec,1),planeeqs,1,onesvec3);

%##################################################################
%##################################################################
%##################################################################
%
% Second order
%
%---------------------------------------------------------------------------------------------

startindices(2) = startindices(1) + length(possiblefirsts);

if specorder > 1
    if SHOWTEXT >= 3
        disp(['     Order number 2'])
    end
    
    if nhighval < 255
        startplanes = uint8(possiblefirsts);
    else
        startplanes = uint16(possiblefirsts);
    end
        
	addtocomb = startplanes;
	nadds = size(addtocomb,1);
    if nadds > 0

       	maxnumberofvispla = max(sum(PLANEseesPLANE(:,startplanes)==1));
	
        addtocomb = (reshape(addtocomb(:,ones(1,maxnumberofvispla)).',length(addtocomb)*maxnumberofvispla,1));
		addtocomb = [addtocomb zeros(size(addtocomb))];
        
        addtoISESVISIBILITY = reshape(ISESVISIBILITY(:,ones(1,maxnumberofvispla)).',nadds*maxnumberofvispla,1);
        
        startpos = [[1:maxnumberofvispla:length(addtocomb)]  length(addtocomb)+1  ].';
		
		for ii = 1:length(startpos)-1
            if nvisPLANES(addtocomb(startpos(ii))) > 0
    			possibleplanes = listofvisPLANES(addtocomb(startpos(ii)),1:nvisPLANES(addtocomb(startpos(ii)))).';
			    nposs = length(possibleplanes);
			    addtocomb( startpos(ii):startpos(ii)+nposs-1,2) = possibleplanes;
            end
            
		end
	
        addtooriginsfrom = uint32([1:startindices(2)-1].');
		addtooriginsfrom = addtooriginsfrom(:,ones(1,maxnumberofvispla));
		addtooriginsfrom = reshape(addtooriginsfrom.',maxnumberofvispla*(startindices(2)-1),1);
		
		cutout = find(addtocomb(:,2)==0);
		addtocomb(cutout,:) = [];
		addtooriginsfrom(cutout) = [];
        addtoISESVISIBILITY(cutout) = [];
        clear cutout
        nadds = size(addtocomb,1);
		    
        if nadds > 0

            %--------------------------------------------------------------------------
			% Try to get rid of some combinations that we know can not be propagated
                        
			% Find those combinations of specular-specular where the IS is behind the second reflecting plane
			% because they can then never ever be propagated, or valid
			
			if difforder > 0
                ivss = uint32(find(addtocomb(:,1)<=nplanes & addtocomb(:,2)<=nplanes));
			else
                ivss = uint32([1:nadds].');
			end

			imsoucheck = dot((ISCOORDS(addtooriginsfrom(ivss),1:3) - corners(planecorners(addtocomb(ivss,2),1),:)).',planenvecs(addtocomb(ivss,2),:).').';
			imsounotOK = find(imsoucheck < GEOMACC);
			addtocomb(ivss(imsounotOK),:) = [];
			addtooriginsfrom(ivss(imsounotOK),:) = [];
            addtoISESVISIBILITY(ivss(imsounotOK)) = [];

            nadds = size(addtocomb,1);

            if nadds > 0 & doconetrace == 1

                if difforder > 0
                    ivss = uint32(find(addtocomb(:,1)<=nplanes & addtocomb(:,2)<=nplanes));
				else
                    ivss = uint32([1:nadds].');
				end
                
                beamdirection1 = planemidpoints(addtocomb(ivss,1),:) - ISCOORDS(addtocomb(ivss,1),:);
                beamlength = sqrt( sum( beamdirection1.'.^2 ) ).';
                beamdirection1  = beamdirection1./beamlength(:,ones(1,3));                
                maxprojradius = maxdisttocorner(addtocomb(ivss,1));
                % The maxdisttocorner below should be multiplied
                % by sine of the angle between the plane normal vector and the
                % beam direction.
                iv = find(beamlength<maxprojradius);
                if ~isempty(iv)
                    beamlength(iv) = beamlength(iv)*0+ eps*1e6;    
                end
                beamangle1 = atan(maxprojradius./beamlength);
                
                % Now we calculate the beam-widths to the second reflection
                % plane, seen from the first-order image source. If this
                % beamwidth is outside the previously calculated beam of
                % reflection plane 1, then plane 2 can never be seen
                % through plane 1.

                beamdirection2 = planemidpoints(addtocomb(ivss,2),:) - ISCOORDS(addtocomb(ivss,1),:);
                beamlength = sqrt( sum( beamdirection2.'.^2 ) ).';
                beamdirection2  = beamdirection2./beamlength(:,ones(1,3));                
                maxprojradius = maxdisttocorner(addtocomb(ivss,2));
                % The maxdisttocorner below should be multiplied
                % by sine of the angle between the plane normal vector and the
                % beam direction.
                iv = find(beamlength<maxprojradius);
                if ~isempty(iv)
                    beamlength(iv) = beamlength(iv)*0+ eps*10;    
                end
                beamangle2 = atan(maxprojradius./beamlength);
                anglebetweenbeams = acos(sum( (beamdirection1.*beamdirection2).' ).');
                ivcutout = find( (beamangle1 + beamangle2 < anglebetweenbeams) & (beamangle1 <1.56) & (beamangle2 < 1.56));
                if ~isempty(ivcutout)
					addtocomb(ivss(ivcutout),:) = [];
					addtooriginsfrom(ivss(ivcutout),:) = [];
                    addtoISESVISIBILITY(ivss(ivcutout)) = [];
                    nadds = size(addtocomb,1);
                end
                
            end
            
			% Combinations of diffraction-specular don't need to be checked because
			% 'edgeseesplane' should have taken care of this.
			%%%%ivds = find(addtocomb(:,1)>nplanes  & addtocomb(:,2)<=nplanes);
			
			%--------------------------------------------------------------------------
			% Combinations of specular-diffraction: check if the IS is behind both of
			% the two planes that the reflecting edge connects to.
            %
            % Bug found 20080711 PS: The visibility test must be split into two, depending on the wedge angle!! 
            % If the open wedge angle > pi, then it is correct to check if the IS behind both of the planes. 
            % If the open wedge angle < pi, then it should be checked if the IS behind one of the planes!! 
            
            if difforder > 0 & nadds > 0

                % First we check the sd combinations where the open wedge angle > pi
                % For those cases, we can exclude combinations where the IS
                % is behind both planes
                
  				ivsd = uint32(find(addtocomb(:,1)<=nplanes & addtocomb(:,2)>nplanes));
                ivsd1 = ivsd( find(   closwedangvec( double(addtocomb(ivsd,2))-nplanes ) < pi ));

                if ~isempty(ivsd1)
                    imsoucheck1 = dot((ISCOORDS(addtooriginsfrom(ivsd1),1:3) - corners(planecorners(    planesatedge(double(addtocomb(ivsd1,2))-nplanes,1)     ,1),:)).',planenvecs(    planesatedge(double(addtocomb(ivsd1,2))-nplanes,1)    ,:).').';
                    imsoucheck2 = dot((ISCOORDS(addtooriginsfrom(ivsd1),1:3) - corners(planecorners(    planesatedge(double(addtocomb(ivsd1,2))-nplanes,2)     ,1),:)).',planenvecs(    planesatedge(double(addtocomb(ivsd1,2))-nplanes,2)    ,:).').';
                    GEOMACC = 1e-10;
                    imsounotOK = find(imsoucheck1 < GEOMACC & imsoucheck2 < GEOMACC);
                    clear imsoucheck1 imsoucheck2

                    addtocomb(ivsd1(imsounotOK),:) = [];
                    addtooriginsfrom(ivsd1(imsounotOK),:) = [];
                    addtoISESVISIBILITY(ivsd1(imsounotOK)) = [];
                  end
                
                % Second we check the sd combinations where the open wedge angle < pi
                % For those cases, we can exclude combinations where the IS
                % is behind one of the two planes
                
  				ivsd = uint32(find(addtocomb(:,1)<=nplanes & addtocomb(:,2)>nplanes));
                if ~isempty(ivsd)
                    ivsd1 = ivsd( find(   closwedangvec( double(addtocomb(ivsd,2))-nplanes ) > pi ));

                    if ~isempty(ivsd1)
                        imsoucheck1 = dot((ISCOORDS(addtooriginsfrom(ivsd1),1:3) - corners(planecorners(    planesatedge(double(addtocomb(ivsd1,2))-nplanes,1)     ,1),:)).',planenvecs(    planesatedge(double(addtocomb(ivsd1,2))-nplanes,1)    ,:).').';
                        imsoucheck2 = dot((ISCOORDS(addtooriginsfrom(ivsd1),1:3) - corners(planecorners(    planesatedge(double(addtocomb(ivsd1,2))-nplanes,2)     ,1),:)).',planenvecs(    planesatedge(double(addtocomb(ivsd1,2))-nplanes,2)    ,:).').';
                        GEOMACC = 1e-10;
                        imsounotOK = find(imsoucheck1 < GEOMACC | imsoucheck2 < GEOMACC);
                        clear imsoucheck1 imsoucheck2

                        addtocomb(ivsd1(imsounotOK),:) = [];
                        addtooriginsfrom(ivsd1(imsounotOK),:) = [];
                        addtoISESVISIBILITY(ivsd1(imsounotOK)) = [];
                    end
                end
                
				nadds = size(addtocomb,1);
                
                ivsd = uint32(find(addtocomb(:,1)<=nplanes & addtocomb(:,2)>nplanes));
                
                if ~isempty(ivsd)
		
                    %----------------------------------------------------------------------
                    % Combinations of specular-diffraction that are not visible at
                    % all should be removed and then they are not propagated either.
				
                    if nedges < 255
                        possibleedges = uint8(double(addtocomb(ivsd,2)) - nplanes);
                    else
                        possibleedges = uint16(double(addtocomb(ivsd,2)) - nplanes);
                    end
                     
                    possiblecombs = addtocomb(ivsd,1);
                    reftoISCOORDS = addtooriginsfrom(ivsd);
				
                    % Expand to take the edge subdivisions into account
                    
                    nposs = length(ivsd);
                    nposs = nposs*nedgesubs;        % We treat all the little edge subdivisions as separate edges
				
                    expandedposscombs = possiblecombs(:,onesvec1);
                    clear possiblecombs
                    expandedposscombs = reshape(expandedposscombs.',nposs,1);			
                    
                    expandedreftoISCOORDS = reftoISCOORDS(:,onesvec1);
                    expandedreftoISCOORDS = reshape(expandedreftoISCOORDS.',nposs,1);
                    expandedpossedges = possibleedges(:,onesvec1);
                    expandedpossedges = reshape(expandedpossedges.',nposs,1);
                    expandedivsd = ivsd(:,onesvec1);
                    expandedivsd = reshape(expandedivsd.',nposs,1);
                    
                    if SHOWTEXT >= 3
						disp(['         ',int2str(nposs),' IS+edge segments found initially '])
					end
                    
                    if nposs > 0
				
                        fromcoords = ISCOORDS(expandedreftoISCOORDS,:);
                        [tocoords,edgeweightlist,edgenumberlist] = EDB2getedgepoints(edgestartcoords(possibleedges,:),edgeendcoords(possibleedges,:),edgelengthvec(possibleedges,:),nedgesubs);
                        clear possibleedges
                        
                         [hitplanes,reflpoints,edgehits,edgehitpoints,cornerhits,cornerhitpoints] = EDB2chkISvisible(fromcoords,tocoords,planeeqs(expandedposscombs(:,1),4),planenvecs(expandedposscombs(:,1),:),minvals(expandedposscombs(:,1),:),...
 						    maxvals(expandedposscombs(:,1),:),planecorners(expandedposscombs(:,1),:),corners,ncornersperplanevec(expandedposscombs(:,1)));
                        if ~isempty(edgehits) | ~isempty(cornerhits)
                            disp('WARNING! An edgehit or cornerhit occurred during the visibility test but this is not')
                            disp('         handled correctly yet.')
                        end

                        eval(['reflpoints',JJ(1,1),' = reflpoints;'])
				
                        expandedivsd          = expandedivsd(hitplanes);
                        expandedposscombs     = expandedposscombs(hitplanes,:);
                        expandedreftoISCOORDS = expandedreftoISCOORDS(hitplanes);
                        expandedpossedges     = expandedpossedges(hitplanes);
                        edgeweightlist        = edgeweightlist(hitplanes);
                        toedgecoords          = tocoords(hitplanes,:);
				
                        nposs = length(expandedivsd);
				
                    end
                    
                    if SHOWTEXT >= 3
						disp(['         ',int2str(nposs),' IS+edge segments visible '])
					end
			
                    if obstructtestneeded & nposs > 0
                        
                        for jj = 1:2
                            if nposs > 0
                                if jj==1
                                    fromcoords = S;    
                                    startplanes = [];    
                                else
                                    startplanes = expandedposscombs(:,jj-1);
                                    eval(['fromcoords = reflpoints',JJ(jj-1,1:JJnumbofchars(jj-1)),';'])
                                end
                                if jj == 2
                                    tocoords = toedgecoords;
                                    endplanes = [planesatedge(expandedpossedges,1) planesatedge(expandedpossedges,2)];
                                else
                                    eval(['tocoords = reflpoints',JJ(jj,1:JJnumbofchars(jj)),';'])    
                                    endplanes = expandedposscombs(:,jj);    
                                end

                                [nonobstructedpaths,nobstructions] = EDB2checkobstrpaths(fromcoords,tocoords,startplanes,endplanes,canplaneobstruct,planeseesplane,...
                                    planeeqs,planenvecs,minvals,maxvals,planecorners,corners,ncornersperplanevec,rearsideplane);
                                
                                if nobstructions > 0
                                    expandedivsd          = expandedivsd(nonobstructedpaths);
                                    expandedposscombs     = expandedposscombs(nonobstructedpaths,:);
                                    expandedreftoISCOORDS = expandedreftoISCOORDS(nonobstructedpaths);
                                    expandedpossedges     = expandedpossedges(nonobstructedpaths);
                                    edgeweightlist        = edgeweightlist(nonobstructedpaths);
                                    toedgecoords          = tocoords(nonobstructedpaths,:);
                                    nposs = length(expandedivsd);
                                    
                                    for kk = 1:1, %norder
                                        eval(['reflpoints',JJ(kk,1:JJnumbofchars(kk)),' = reflpoints',JJ(kk,1:JJnumbofchars(kk)),'(nonobstructedpaths,:);'])    
                                    end
                                    
                                end
                            end
                 
                        end
                    end        
                    
                    if SHOWTEXT >= 3
					    disp(['         ',int2str(nposs),' IS+edge segments survived the obstruction test'])
                    end
                    
                    % There are repetitions of each combination since each edge segment has
                    % been treated separately. Pack them together now
			
                    test = [expandedposscombs expandedpossedges];
                    ncombs = length(expandedpossedges);
                    dtest = diff([0;prod(   double(test(1:ncombs-1,:)==test(2:ncombs,:)).'  ).']);
                    ivremove = find(dtest==1);
                    
                    while ~isempty(ivremove) & ~isempty(edgeweightlist),
                        edgeweightlist(ivremove+1) = double(edgeweightlist(ivremove+1)) + double(edgeweightlist(ivremove));
                        edgeweightlist(ivremove) = [];
                        expandedpossedges(ivremove) = [];
                        expandedposscombs(ivremove,:) = [];
                        expandedivsd(ivremove) = [];
                        nposs = length(expandedivsd);   
                        
                        test = [expandedposscombs expandedpossedges];
                        ncombs = length(expandedpossedges);
                        dtest = diff([0;prod(   double(test(1:ncombs-1,:)==test(2:ncombs,:)).'  ).']);
                        ivremove = find(dtest==1);
			
                    end
				       
                    if SHOWTEXT >= 3
						disp(['         ',int2str(nposs),' IS+edge combinations after edge segments are packed together'])
                    end
                    combstoremove = setdiff(ivsd,expandedivsd);
                    addtocomb(combstoremove,:) = [];
                    addtooriginsfrom(combstoremove) = [];  
                    addtoISESVISIBILITY(combstoremove) = [];
                    nadds = length(addtooriginsfrom);
                end
                
            end
            
			% Combinations of diffraction-diffraction don't need to be checked because
			% 'edgeseesedge' should have taken care of this.
			
            %----------------------------------------------------------------------
            % Add the new combinations to the master list
		
            if nadds > 0
                POTENTIALISES = [POTENTIALISES;[addtocomb zeros(nadds,specorder-2)]];
				ORIGINSFROM = [ORIGINSFROM;addtooriginsfrom];
			                
                if difforder > 0
                    ivtemp = find(addtocomb(:,1)<=nplanes & addtocomb(:,2)>nplanes);
                    if ~isempty(ivtemp),        
                        addtoISESVISIBILITY(ivtemp) = edgeweightlist;
                    end
                    ivtemp = find(addtocomb(:,1)>nplanes);
                end
                ISESVISIBILITY = [ISESVISIBILITY;addtoISESVISIBILITY];                

				endindices(2) = length(ORIGINSFROM);
			
				% Compute the IS coords of the combinations with only specular reflections
				
				ISCOORDStoadd = zeros(nadds,3);
				if difforder > 0
                    ivss = uint32(find(addtocomb(:,1)<=nplanes & addtocomb(:,2)<=nplanes));
				else
                    ivss = uint32([1:nadds].');    
				end
				soutoreflect = ISCOORDS(ORIGINSFROM(double(ivss)+endindices(1)),1:3);
				ISCOORDStoadd(ivss,:) = EDB2findis(soutoreflect,POTENTIALISES(double(ivss)+endindices(1),2),planeeqs,size(soutoreflect,1),onesvec3);
                clear soutoreflect
				ISCOORDS = [ISCOORDS;ISCOORDStoadd];
		
            end
        end
    end
end

%##################################################################
%##################################################################
%##################################################################
%
% Third order and higher
%
%---------------------------------------------------------------------------------------------

for ordernum = 3:specorder
    if SHOWTEXT >= 3
        disp(['     Order number ',int2str(ordernum)])
    end
        
    % Fixed bug 2013_08_13 PS: the first line below gave errors when a
    % geometry had only two planes.
    % startindices(ordernum) = startindices(ordernum-1) + length(addtocomb);
    startindices(ordernum) = startindices(ordernum-1) + size(addtocomb,1);

    % The lines below could use
    % planestoprop = unique(addtocomb(:,ordernum-1));
    
    planelist = sort(addtocomb(:,ordernum-1));
    
    if ~isempty(planelist),            
        planestoprop = planelist([1;find(diff(double(planelist))~=0)+1]);
    else
        planestoprop = [];    
    end
    clear planelist

	maxnumberofvispla = max(sum(PLANEseesPLANE(:,planestoprop)==1));
    oldaddtocomb = addtocomb;

	nadds = size(addtocomb,1);
    if nhighval < 255
        addtocomb = uint8(zeros(nadds*maxnumberofvispla,ordernum));
    else
        addtocomb = uint16(zeros(nadds*maxnumberofvispla,ordernum));
    end

    onesvec2 = ones(1,maxnumberofvispla);
	for ii = 1:ordernum-1
		temp = oldaddtocomb(:,ii);
		addtocomb(:,ii) = reshape(temp(:,onesvec2).',length(temp)*maxnumberofvispla,1);
	end
    nrows = size(addtocomb,1);
	clear oldaddtocomb temp
    
    if nrows > 0
       	startpos = [[1:maxnumberofvispla:nrows]  nrows+1  ].';
    else
        startpos = 1;    
    end

	for ii = 1:length(startpos)-1
            if nvisPLANES(addtocomb(startpos(ii),ordernum-1)) > 0
    			possibleplanes = listofvisPLANES(addtocomb(startpos(ii),ordernum-1),1:nvisPLANES(addtocomb(startpos(ii),ordernum-1))).';
			    nposs = length(possibleplanes);
			    addtocomb( startpos(ii):startpos(ii)+nposs-1,ordernum) = possibleplanes;
            end
	end

    addtooriginsfrom = uint32([1:nadds].' + startindices(ordernum-1)-1);
    addtooriginsfrom = addtooriginsfrom(:,ones(1,maxnumberofvispla));
    addtooriginsfrom = reshape(addtooriginsfrom.',maxnumberofvispla*nadds,1);

    addtoISESVISIBILITY = reshape(addtoISESVISIBILITY(:,ones(1,maxnumberofvispla)).',nadds*maxnumberofvispla,1);
    
	clear startpos
	cutout = find(addtocomb(:,ordernum)==0);
	addtocomb(cutout,:) = [];
    addtooriginsfrom(cutout) = [];
    addtoISESVISIBILITY(cutout) = [];
    nadds = size(addtocomb,1);
    
	if ordernum == specorder
		clear cutout
	end

    %--------------------------------------------------------------------------
    % Try to get rid of some combinations that we know can not be propagated

    % Find those combinations of specular-specular where the IS is behind the second reflecting plane
    % because they can then never ever be propagated.

    if difforder > 0
        ivss = uint32(find(prod( double(addtocomb<=nplanes).' ).'));
    else
        ivss = uint32([1:nadds].');
    end
    
    imsoucheck = dot((ISCOORDS(addtooriginsfrom(ivss),1:3) - corners(planecorners(addtocomb(ivss,ordernum),1),:)).',planenvecs(addtocomb(ivss,ordernum),:).').';
	GEOMACC = 1e-10;
	imsounotOK = uint32(find(imsoucheck < GEOMACC));
    clear imsoucheck
	addtocomb(ivss(imsounotOK),:) = [];
	addtooriginsfrom(ivss(imsounotOK),:) = [];
    addtoISESVISIBILITY(ivss(imsounotOK)) = [];
    clear imsounotOK

    % Addition 20100210 PS: we should check all combinations with dss in the
    % three last orders of addtocomb: if the edge in the first 'd' belongs to the
    % plane of the last 's' *and* the two 's' planes form a 90 degree
    % corner, then these combos should be removed.
    % 
    % As a first step, we check if there are any interior-90-degree
    % corners/edges of the model, because only then can these combos occur.
    
    if ordernum == 3
       deg90edges = find( abs(closwedangvec-3*pi/2) < GEOMACC );
       if any(deg90edges)
           deg90planepairs = planesatedge(deg90edges,:);    
       end
    end    
    
    if difforder > 0 & any(deg90edges)
        ivdss = uint32(find(    double(addtocomb(:,ordernum-2)>nplanes).*double(addtocomb(:,ordernum-1)<=nplanes).*double(addtocomb(:,ordernum)<=nplanes)      ));
    else
        ivdss = [];
    end

    if ~isempty(ivdss)
        
        % First we check if the 'd' (the edge) belongs to the last 's'
        A = addtocomb(ivdss,ordernum-2:ordernum);

        ivsubset = find( (planesatedge(A(:,1)-nplanes,1) == A(:,3)) | (planesatedge(A(:,1)-nplanes,2) == A(:,3)) );
        ivdss = ivdss(ivsubset);

        if ~isempty(ivdss)
            % Second, we check if the two 's' planes form a 90 deg corner
            planesare90degpairs = ismember( A(:,2),deg90planepairs(:,1) ).*ismember( A(:,3),deg90planepairs(:,2) ) + ismember( A(:,2),deg90planepairs(:,2) ).*ismember( A(:,3),deg90planepairs(:,1) );
            ivsubset = find(planesare90degpairs);
            ivdss = ivdss(ivsubset);

            if ~isempty(ivdss)
                addtocomb(ivdss,:) = [];
                addtooriginsfrom(ivdss,:) = [];
                addtoISESVISIBILITY(ivdss) = [];
                clear ivdss ivsubset planesare90degpairs                
            end
        end
    end    
    
    nadds = size(addtocomb,1);
    
    if nadds > 0 & doconetrace == 1

        if difforder > 0
            ivss = uint32(find(prod( double(addtocomb<=nplanes).' ).'));
        else
            ivss = uint32([1:nadds].');
        end
        
        beamdirection1 = planemidpoints(addtocomb(ivss,ordernum-1),:) - ISCOORDS(addtooriginsfrom(ivss),1:3);
        beamlength = sqrt( sum( beamdirection1.'.^2 ) ).';
        beamdirection1  = beamdirection1./beamlength(:,ones(1,3));                
        maxprojradius = maxdisttocorner(addtocomb(ivss,ordernum-1));
        iv = find(beamlength<maxprojradius);
        if ~isempty(iv)
            beamlength(iv) = beamlength(iv)*0+ eps*1e6;    
        end
        beamangle1 = atan(maxprojradius./beamlength);
                
        % Now we calculate the beam-widths to the second reflection
        % plane, seen from the first-order image source. If this
        % beamwidth is outside the previously calculated beam of
        % reflection plane 1, then plane 2 can never be seen
        % through plane 1.

        beamdirection2 = planemidpoints(addtocomb(ivss,ordernum),:) - ISCOORDS(addtooriginsfrom(ivss),1:3);
        beamlength = sqrt( sum( beamdirection2.'.^2 ) ).';
        beamdirection2  = beamdirection2./beamlength(:,ones(1,3));                
        maxprojradius = maxdisttocorner(addtocomb(ivss,ordernum));
        iv = find(beamlength<maxprojradius);
        if ~isempty(iv)
            beamlength(iv) = beamlength(iv)*0+ eps*10;    
        end
        beamangle2 = atan(maxprojradius./beamlength);
       clear maxprojradius beamlength
       
       ivcutout = find( (beamangle1 + beamangle2 < acos(sum( (beamdirection1.*beamdirection2).' ).')) & (beamangle1 <1.56) & (beamangle2 < 1.56));
       clear beamdirection1 beamdirection2 beamangle1 beamangle2
       
        if ~isempty(ivcutout)
			addtocomb(ivss(ivcutout),:) = [];
			addtooriginsfrom(ivss(ivcutout),:) = [];
            addtoISESVISIBILITY(ivss(ivcutout)) = [];
            nadds = size(addtocomb,1);
        end
                
    end
    
	% Combinations with all-specular plus diffraction as the last one:
    % check if the IS is behind both of the two planes that the reflecting edge connects to.
    %
    % Bug found 080711: The visibility test must be split into two, depending on the wedge angle!! 
    % If the open wedge angle > pi, then it is correct to check if the IS behind both of the planes. 
    % If the open wedge angle < pi, then it should be checked if the IS behind one of the planes!! 

    if difforder > 0
        
        ivsd = uint32(find(prod( double(addtocomb(:,1:ordernum-1)<=nplanes).' ).' & (addtocomb(:,ordernum)>nplanes)));

        if ~isempty(ivsd)

             % First we check the sssssd combinations where the open wedge angle > pi
             % For those cases, we can exclude combinations where the IS
             % is behind both planes
                
             ivsd1 = ivsd( find(   closwedangvec( double(addtocomb(ivsd,ordernum))-nplanes ) < pi ));

             if ~isempty(ivsd1)
                imsoucheck1 = dot((ISCOORDS(addtooriginsfrom(ivsd1),1:3) - corners(planecorners(    planesatedge(double(addtocomb(ivsd1,ordernum))-nplanes,1)     ,1),:)).',planenvecs(    planesatedge(double(addtocomb(ivsd1,ordernum))-nplanes,1)    ,:).').';
                imsoucheck2 = dot((ISCOORDS(addtooriginsfrom(ivsd1),1:3) - corners(planecorners(    planesatedge(double(addtocomb(ivsd1,ordernum))-nplanes,2)     ,1),:)).',planenvecs(    planesatedge(double(addtocomb(ivsd1,ordernum))-nplanes,2)    ,:).').';
                GEOMACC = 1e-10;
                imsounotOK = find(imsoucheck1 < GEOMACC & imsoucheck2 < GEOMACC);
                clear imsoucheck1 imsoucheck2
                addtocomb(ivsd1(imsounotOK),:) = [];
                addtooriginsfrom(ivsd1(imsounotOK),:) = [];
                addtoISESVISIBILITY(ivsd1(imsounotOK)) = [];
                clear imsounotOK
             end
             
             % Second we check the sssssd combinations where the open wedge angle < pi
             % For those cases, we can exclude combinations where the IS
             % is behind one of the two planes
                
             ivsd = uint32(find(prod( double(addtocomb(:,1:ordernum-1)<=nplanes).' ).' & (addtocomb(:,ordernum)>nplanes)));
             if ~isempty(ivsd)
                 ivsd1 = ivsd( find(   closwedangvec( double(addtocomb(ivsd,ordernum))-nplanes ) > pi ));

                 if ~isempty(ivsd1)
                    imsoucheck1 = dot((ISCOORDS(addtooriginsfrom(ivsd1),1:3) - corners(planecorners(    planesatedge(double(addtocomb(ivsd1,ordernum))-nplanes,1)     ,1),:)).',planenvecs(    planesatedge(double(addtocomb(ivsd1,ordernum))-nplanes,1)    ,:).').';
                    imsoucheck2 = dot((ISCOORDS(addtooriginsfrom(ivsd1),1:3) - corners(planecorners(    planesatedge(double(addtocomb(ivsd1,ordernum))-nplanes,2)     ,1),:)).',planenvecs(    planesatedge(double(addtocomb(ivsd1,ordernum))-nplanes,2)    ,:).').';
                    GEOMACC = 1e-10;
                    imsounotOK = find(imsoucheck1 < GEOMACC | imsoucheck2 < GEOMACC);
                    clear imsoucheck1 imsoucheck2
                    addtocomb(ivsd1(imsounotOK),:) = [];
                    addtooriginsfrom(ivsd1(imsounotOK),:) = [];
                    addtoISESVISIBILITY(ivsd1(imsounotOK)) = [];
                    clear imsounotOK
                 end
             end
             
            % For combinations including dsd, if the two diff. edges are
            % aligned and the intermediate specular reflections are
            % caused by planes that are perpendicular to the edges
            % remove these combinations
            
            if difforder >= 2
                for kk = 1:ordernum-2
                    ncombs = size(addtocomb,1); 
                    matchvec = ones(ncombs,1);
                    nprespecs = kk-1;
                    nmidspecs = ordernum-1-kk;
                    matchvec = matchvec.*(addtocomb(:,ordernum)>nplanes).*(addtocomb(:,kk)>nplanes);
                    if nmidspecs == 1
                        matchvec = matchvec.*(addtocomb(:,2)<=nplanes);
                    elseif nmidspecs >= 2
                        matchvec = matchvec.*prod( double(addtocomb(:,kk+1:ordernum-1)<=nplanes).' ).';
                    end
                    ivdsd = uint32(find( matchvec ));
                    if ~isempty(ivdsd)
                        edge1 = double(addtocomb(ivdsd,kk)) - nplanes;
                        edge2 = double(addtocomb(ivdsd,ordernum)) - nplanes;
                        midspecs = double(addtocomb(ivdsd,kk+1:ordernum-1));
                        
                        lookupindmat1 = (edge1-1)*nedges + edge2;
                        matrixcolnumbers = (edge1-1)*nplanes;
                        lookupindmat2 = matrixcolnumbers(:,ones(1,nmidspecs)) + midspecs;                    
                        
                        if nmidspecs == 1
                            ivinvalid = find(edgealignedwithedge(lookupindmat1) & edgeperptoplane(lookupindmat2));
                        else
                            ivinvalid = find(edgealignedwithedge(lookupindmat1) & prod(edgeperptoplane(lookupindmat2).').');
                        end
            			addtocomb(ivdsd(ivinvalid),:) = [];
				        addtooriginsfrom(ivdsd(ivinvalid),:) = [];
                        addtoISESVISIBILITY(ivdsd(ivinvalid)) = [];
                    end
                 end
            end            
            
            % We need to find the valid combinations again after we have removed a
            % number of combinations
	
            ivsd = uint32(find(prod( double(addtocomb(:,1:ordernum-1)<=nplanes).' ).' & (addtocomb(:,ordernum)>nplanes)));
	
            %----------------------------------------------------------------------
            % Combinations of specular-diffraction that are not visible at
            % all should be removed and then they are not propagated either.
		
            if nhighval < 255
                possibleedges = uint8(double(addtocomb(ivsd,ordernum)) - nplanes);
                possiblecombs = uint8(double(addtocomb(ivsd,1:ordernum-1)));
            else
                possibleedges = uint16(double(addtocomb(ivsd,ordernum)) - nplanes);
                possiblecombs = uint16(double(addtocomb(ivsd,1:ordernum-1)));
            end
            reftoISCOORDS = addtooriginsfrom(ivsd);
		
            % Expand to take the edge subdivisions into account
            
            nposs = length(ivsd);
            nposs = nposs*nedgesubs;        % We treat all the little edge subdivisions as separate edges

            expandedposscombs = reshape(repmat(possiblecombs.',[nedgesubs,1]),ordernum-1,nposs).';
            clear possiblecombs                        
            expandedreftoISCOORDS = reftoISCOORDS(:,onesvec1);
            expandedreftoISCOORDS = reshape(expandedreftoISCOORDS.',nposs,1);
            expandedpossedges = possibleedges(:,onesvec1);
            expandedpossedges = reshape(expandedpossedges.',nposs,1);
            expandedivsd = ivsd(:,onesvec1);
            expandedivsd = reshape(expandedivsd.',nposs,1);
            
            if SHOWTEXT >= 3
				disp(['         ',int2str(nposs),' IS+edge segments found initially '])
			end
	
            % Go through, iteratively, and check if the path from S to edge passes
            % through all reflection planes along the way
            
            for jj = ordernum-1:-1:1
	
                if nposs > 0
		
                    if jj == ordernum-1
                        fromcoords = ISCOORDS(expandedreftoISCOORDS,:);
                        [tocoords,edgeweightlist,edgenumberlist] = EDB2getedgepoints(edgestartcoords(possibleedges,:),edgeendcoords(possibleedges,:),edgelengthvec(possibleedges,:),nedgesubs);
                        clear possibleedges
                    else
                        eval(['tocoords = reflpoints',JJ(jj+1,1:JJnumbofchars(jj+1)),';'])    
                        ivlist = ORIGINSFROM(expandedreftoISCOORDS);
                        for kk = jj:ordernum-3
                            ivlist = ORIGINSFROM(ivlist);                            
                        end
                        fromcoords = ISCOORDS(ivlist,:);
                                           
                    end

                    [hitplanes,reflpoints,edgehits,edgehitpoints,cornerhits,cornerhitpoints] = EDB2chkISvisible(fromcoords,tocoords,planeeqs(expandedposscombs(:,jj),4),planenvecs(expandedposscombs(:,jj),:),minvals(expandedposscombs(:,jj),:),...
					    maxvals(expandedposscombs(:,jj),:),planecorners(expandedposscombs(:,jj),:),corners,ncornersperplanevec(expandedposscombs(:,jj)));
                    if ~isempty(edgehits) | ~isempty(cornerhits)
                        disp('WARNING! An edgehit or cornerhit occurred during the visibility test but this is not')
                        disp('         handled correctly yet.')
                    end
                    eval(['reflpoints',JJ(jj,1:JJnumbofchars(jj)),' = reflpoints;'])
			
                    expandedivsd          = expandedivsd(hitplanes);
                    expandedposscombs     = expandedposscombs(hitplanes,:);
                    expandedreftoISCOORDS = expandedreftoISCOORDS(hitplanes);
                    expandedpossedges     = expandedpossedges(hitplanes);
                    edgeweightlist        = edgeweightlist(hitplanes);
                    toedgecoords          = tocoords(hitplanes,:);
                    if jj < ordernum-1
                        for kk = jj+1:ordernum-1
                           eval(['reflpoints',JJ(kk,1:JJnumbofchars(kk)),' = reflpoints',JJ(kk,1:JJnumbofchars(kk)),'(hitplanes,:);'])
                        end
                    end
		
                    nposs = length(expandedivsd);
                    
                end
                if SHOWTEXT >= 3
				    disp(['         ',int2str(nposs),' IS+edge segments survived the visibility test in refl plane ',JJ(jj,1:JJnumbofchars(jj))])
                end
                
            end
	
            if obstructtestneeded & nposs > 0
                for jj = 1:ordernum
                    if nposs > 0
                        
                        if jj==1
                            fromcoords = S;    
                            startplanes = [];    
                        else
                            startplanes = expandedposscombs(:,jj-1);
                            eval(['fromcoords = reflpoints',JJ(jj-1,1:JJnumbofchars(jj-1)),';'])
                        end
                        if jj == ordernum
                            tocoords = toedgecoords;
                            endplanes = [planesatedge(expandedpossedges,1) planesatedge(expandedpossedges,2)];
                        else
                            eval(['tocoords = reflpoints',JJ(jj,1:JJnumbofchars(jj)),';'])    
                            endplanes = expandedposscombs(:,jj);    
                        end
             
                        [nonobstructedpaths,nobstructions] = EDB2checkobstrpaths(fromcoords,tocoords,startplanes,endplanes,canplaneobstruct,planeseesplane,...
                            planeeqs,planenvecs,minvals,maxvals,planecorners,corners,ncornersperplanevec,rearsideplane);
		
                        if nobstructions > 0
                            expandedivsd          = expandedivsd(nonobstructedpaths);
                            expandedposscombs     = expandedposscombs(nonobstructedpaths,:);
                            expandedreftoISCOORDS = expandedreftoISCOORDS(nonobstructedpaths);
                            expandedpossedges     = expandedpossedges(nonobstructedpaths);
                            edgeweightlist        = edgeweightlist(nonobstructedpaths);
                            toedgecoords          = tocoords(nonobstructedpaths,:);
                            nposs = length(expandedivsd);
                            for kk = 1:ordernum-1
                                eval(['reflpoints',JJ(kk,1:JJnumbofchars(kk)),' = reflpoints',JJ(kk,1:JJnumbofchars(kk)),'(nonobstructedpaths,:);'])    
                            end
                            
                        end
                    end
                end
            end        
        
            if SHOWTEXT >= 3
		        disp(['         ',int2str(nposs),' IS+edge segments survived the obstruction test'])
            end

            % There are repetitions of each combination since each edge segment has
            % been treated separately. Pack them together now
	
            test = [expandedposscombs expandedpossedges];
            ncombs = length(expandedpossedges);
            dtest = diff([0;prod(   double(test(1:ncombs-1,:)==test(2:ncombs,:)).'  ).']);
            ivremove = find(dtest==1);
            
            while ~isempty(ivremove)
           
                edgeweightlist(ivremove+1) = double(edgeweightlist(ivremove+1)) + double(edgeweightlist(ivremove));
                edgeweightlist(ivremove) = [];
                expandedpossedges(ivremove) = [];
                expandedposscombs(ivremove,:) = [];
                expandedivsd(ivremove) = [];
                nposs = length(expandedivsd);   
                
                test = [expandedposscombs expandedpossedges];
                ncombs = length(expandedpossedges);
                dtest = diff([0;prod(   double(test(1:ncombs-1,:)==test(2:ncombs,:)).'  ).']);
                ivremove = find(dtest==1);
	
            end
            	       
            if SHOWTEXT >= 3
				disp(['         ',int2str(nposs),' IS+edge segments after edge segments are packed together'])
            end
	
            combstoremove = setdiff(ivsd,expandedivsd);
            addtocomb(combstoremove,:) = [];
            addtooriginsfrom(combstoremove) = [];  
            addtoISESVISIBILITY(combstoremove) = [];
            nadds = length(addtooriginsfrom);
        end        
    end   

    POTENTIALISES = [POTENTIALISES;[addtocomb zeros(nadds,specorder-ordernum)]];
    ORIGINSFROM = [ORIGINSFROM;addtooriginsfrom];

    ivtemp = find(prod( double(addtocomb(:,1:ordernum-1)<=nplanes).' ).' & (addtocomb(:,ordernum)>nplanes));
    if difforder > 0
        if ~isempty(ivtemp)
            addtoISESVISIBILITY(ivtemp) = edgeweightlist;
        end
    end
    ISESVISIBILITY = [ISESVISIBILITY;addtoISESVISIBILITY];
    
    endindices(ordernum) = length(ORIGINSFROM);

    % Compute the IS coords

    ISCOORDStoadd = zeros(nadds,3);
    if difforder > 0
        ivss = uint32(find(prod( double(addtocomb<=nplanes).' ).'));
    else
        ivss = uint32([1:nadds].');    
    end
    soutoreflect = ISCOORDS(ORIGINSFROM(double(ivss)+endindices(ordernum-1)),1:3);
    ISCOORDStoadd(ivss,:) = EDB2findis(soutoreflect,POTENTIALISES(double(ivss)+endindices(ordernum-1),ordernum),planeeqs,size(soutoreflect,1),onesvec3);
    clear soutoreflect
    ISCOORDS = [ISCOORDS;ISCOORDStoadd];
    clear ISCOORDStoadd
    
end

ntot = size(POTENTIALISES,1);

%-------------------------------------------------------
% Add the direct sound to the start of the list

POTENTIALISES = [zeros(1,specorder);POTENTIALISES];
ORIGINSFROM = [0;ORIGINSFROM];
startindices = startindices+1;
endindices = endindices+1;
ORIGINSFROM = uint32(double(ORIGINSFROM)+1);

ISCOORDS = [S;ISCOORDS];
ISESVISIBILITY = [0;ISESVISIBILITY];

%-------------------------------------------------------
% Trim the list if there are zeros-only in the last column(s)

[n1,n2] = size(POTENTIALISES);
if n2 > 1
    columnstokeep = find(sum(POTENTIALISES)~=0);
    POTENTIALISES = POTENTIALISES(:,columnstokeep);    
    [n1,n2new] = size(POTENTIALISES);    
    if n2new < n2
        specorderold = specorder;
        specorder = n2new;
    end
end

%-------------------------------------------------------
% Create index lists

% First, the total reflection order

[n1,n2] = size(POTENTIALISES);

lengthNspecmatrix = [];
lengthNdiffmatrix = [];
singlediffcol = [];
startindicessinglediff = [];
endindicessinglediff = [];

if n1 > 0 & n2 > 0
if n2 > 1
    REFLORDER = uint8(sum(POTENTIALISES.'>0).');    
else
    REFLORDER = uint8(POTENTIALISES>0);
end

if difforder > 0
    if specorder > 1
        ivss = uint32(find(prod( double(POTENTIALISES<=nplanes).' ).'));
    else
        ivss = uint32(find( POTENTIALISES<=nplanes ));    
    end

    ivother = uint32([1:length(ORIGINSFROM)].');
    ivother(ivss) = [];
    ivss(1) = [];
    IVNDIFFMATRIX = uint32(zeros(2,specorder));
    lengthNdiffmatrix = zeros(1,specorder);
    nrowsdiff = 2;
    IVNSPECMATRIX = uint32(zeros(2,specorder));
    lengthNspecmatrix = zeros(1,specorder);
    nrowsspec = 2;
    
    for ii = 1:specorder
        if specorder > 1
            ivreftoNdiff = uint32(find(sum( double(POTENTIALISES(ivother,:)> nplanes).' ).'==ii));
        else
            ivreftoNdiff = uint32(find(( double(POTENTIALISES(ivother,:)> nplanes).' ).'==ii));
        end
        ivreftoNspec = uint32(find(REFLORDER(ivss)==ii));

        ivNdiff = ivother(ivreftoNdiff);
        if ~isempty(ivNdiff)
            lengthNdiffmatrix(ii) = length(ivNdiff);
            if lengthNdiffmatrix(ii) > nrowsdiff
               IVNDIFFMATRIX = [IVNDIFFMATRIX;uint32(zeros(lengthNdiffmatrix(ii)-nrowsdiff,specorder))];
               nrowsdiff = lengthNdiffmatrix(ii);
            end
            IVNDIFFMATRIX(1: lengthNdiffmatrix(ii),ii) =  ivNdiff;     
        end
        ivother(ivreftoNdiff) = [];
        
        ivNspec = ivss(ivreftoNspec);
        if ~isempty(ivNspec)
            lengthNspecmatrix(ii) = length(ivNspec);
            if lengthNspecmatrix(ii) > nrowsspec
                IVNSPECMATRIX = [IVNSPECMATRIX;uint32(zeros(lengthNspecmatrix(ii)-nrowsspec,specorder))];
                nrowsspec = lengthNspecmatrix(ii);
            end
            IVNSPECMATRIX(1: lengthNspecmatrix(ii),ii) =  ivNspec;      
        end        
    end
     
    % Determine which column the single diffraction is in
    
    test = POTENTIALISES(IVNDIFFMATRIX(1:lengthNdiffmatrix(1),1),:)>nplanes;
       
    colweight = [1:specorder];
    nsingles = length(IVNDIFFMATRIX(1:lengthNdiffmatrix(1),1));
    if specorder > 1
        singlediffcol = uint8(sum( (test.*colweight(ones(nsingles,1),:)).').');
    else
        singlediffcol = uint8(ones(nsingles,1));        
    end
        
    startindicessinglediff = zeros(specorder,1);
    endindicessinglediff = zeros(specorder,1);
    startindicessinglediff(1) = 1;
    for ii = 1:specorder-1
        iv = find(IVNDIFFMATRIX(1:lengthNdiffmatrix(1),1)<=endindices(ii));
        if ~isempty(iv)
            endindicessinglediff(ii) = iv(length(iv));
        else
            endindicessinglediff(ii) = 0;
        end
        startindicessinglediff(ii+1) = endindicessinglediff(ii)+1;
    end
    endindicessinglediff(specorder) = lengthNdiffmatrix(1);
else
    ivss = uint32([2:endindices(length(endindices))].');
    IVNSPECMATRIX = uint32(zeros(2,specorder));
    lengthNspecmatrix = zeros(1,specorder);
    nrowsspec = 2;

    for ii = 1:specorder
        ivreftoNspec = uint32(find(REFLORDER(ivss)==ii));

        ivNspec = ivss(ivreftoNspec);
        lengthNspecmatrix(ii) = length(ivNspec);
        if lengthNspecmatrix(ii) > nrowsspec
           IVNSPECMATRIX = [IVNSPECMATRIX;uint32(zeros(lengthNspecmatrix(ii)-nrowsspec,specorder))];
           nrowsspec = lengthNspecmatrix(ii);
        end
        IVNSPECMATRIX(1: lengthNspecmatrix(ii),ii) =  ivNspec;      
    
    end
    
    IVNDIFFMATRIX = [];
    lengthNdiffmatrix = [];

    singlediffcol = [];
    startindicessinglediff = [];
    endindicessinglediff = [];
end

end

%-------------------------------------------------------
% Create a pointer list so that it is possible to find a combination
% that ends with edge-plane2, given an edge-plane1-plane2 combination
%
% NB!! The list points to the original POTENTIALISES (and related lists)

if difforder > 0 & specorder > 1

    if SHOWTEXT >= 3
        disp(['       Building a pointer list for image receiver combinations'])    
    end

    % First select only those that have a single diff combo and that
    % start with the diffraction
    
    ndecimaldivider = (nplanes+2);
    PointertoIRcombs = sparse(zeros(ndecimaldivider^(specorder-1),1));
    IRoriginsfrom = uint32(zeros(size(ORIGINSFROM)));
    
	for ii = 1:specorder-1
        if SHOWTEXT >= 3
            disp(['          order ',int2str(ii)])    
        end
                
        if ii == 1
            ivrange = uint32([startindicessinglediff(2):endindicessinglediff(2)]);
            masterivlistselect = IVNDIFFMATRIX(ivrange,1);
            masterivlistselect = masterivlistselect(find(singlediffcol(ivrange)==1));

            % IRoriginsfrom should be given the value zero for these
            % first-order specular combos (after one diff) so we don't
            % bother assigning any value
            
            [B,IA,JA] = unique(POTENTIALISES(masterivlistselect,2));
            PointertoIRcombs(B) = masterivlistselect(IA);
            
            % Now we check if there any active wall numbers that don't
            % occur in an edge-spec combination. For those combinations
            % we can point to a specular combos instead, that ends with the
            % same plane number.
            
            wallist = POTENTIALISES(IVNSPECMATRIX(1:lengthNspecmatrix(1),1),1);
            ivnotincludedyet = find(~ismember(wallist,B));
            if ~isempty(ivnotincludedyet)
                PointertoIRcombs(wallist(ivnotincludedyet)) = IVNSPECMATRIX(ivnotincludedyet,1);
            end
            
        elseif ii >= 2
            ivrange = uint32([startindicessinglediff(ii+1):endindicessinglediff(ii+1)]);
            masterivlistselect = IVNDIFFMATRIX(ivrange,1);
            masterivlistselect = masterivlistselect(find(singlediffcol(ivrange)==1));
 
            ivlist = 0;
            for jj = 3:ii+1
                ivlist = ivlist + double(POTENTIALISES(masterivlistselect,jj))*ndecimaldivider^(ii+1-jj);
            end
            
            IRoriginsfrom(masterivlistselect) = uint32(full(PointertoIRcombs(ivlist)));

            ivlist = ivlist + double(POTENTIALISES(masterivlistselect,2))*ndecimaldivider^(ii-1);
            A = uint32(ivlist);
            [B,IA,JA] = unique(A);
            PointertoIRcombs(B) = masterivlistselect(IA);

            % Now we check if there any active wall numbers that don't
            % occur as spec1 in an edge-spec1-spec2-spec3 combination.

            wallist = 0;
            for jj = 1:ii
                wallist = wallist + double(POTENTIALISES(IVNSPECMATRIX(1:lengthNspecmatrix(ii),ii),jj))*ndecimaldivider^(ii-jj);
            end

            ivnotincludedyet = find(~ismember(wallist,B));

            if ~isempty(ivnotincludedyet)
                PointertoIRcombs(wallist(ivnotincludedyet)) = IVNSPECMATRIX(ivnotincludedyet,ii);
            end
            
        end
	
	end
        
else
    ndecimaldivider = 0;
    PointertoIRcombs = [];
    IRoriginsfrom = [];
end

%-------------------------------------------------------
% Save the output data

maxval = max(IRoriginsfrom);
if maxval < 256
    IRoriginsfrom = uint8(IRoriginsfrom);
elseif maxval < 65536
    IRoriginsfrom = uint16(IRoriginsfrom);    
end
maxval = max(ORIGINSFROM);
if maxval < 256
    ORIGINSFROM = uint8(ORIGINSFROM);
elseif maxval < 65536
    ORIGINSFROM = uint16(ORIGINSFROM);    
end

if SUPPRESSFILES == 0
    % The complete variable set could be saved if one wanted to.
    Varlist = [' POTENTIALISES ORIGINSFROM ISCOORDS  ISESVISIBILITY IVNSPECMATRIX lengthNspecmatrix IVNDIFFMATRIX lengthNdiffmatrix '];
    Varlist = [Varlist,' singlediffcol REFLORDER startindicessinglediff endindicessinglediff ndecimaldivider PointertoIRcombs IRoriginsfrom'];		

    % Varlist = [' lengthNspecmatrix lengthNdiffmatrix '];
    % Varlist = [Varlist,' singlediffcol startindicessinglediff endindicessinglediff ndecimaldivider PointertoIRcombs IRoriginsfrom'];		

    eval(['save ',ISEStreefile,Varlist])
else
   ISEStreefile = struct('POTENTIALISES',POTENTIALISES,'ORIGINSFROM',ORIGINSFROM,'ISCOORDS',ISCOORDS,'ISESVISIBILITY',ISESVISIBILITY,'IVNSPECMATRIX',IVNSPECMATRIX,...
       'lengthNspecmatrix',lengthNspecmatrix,'IVNDIFFMATRIX',IVNDIFFMATRIX,'lengthNdiffmatrix',lengthNdiffmatrix,'singlediffcol',singlediffcol,...
       'REFLORDER',REFLORDER,'startindicessinglediff',startindicessinglediff,'endindicessinglediff',endindicessinglediff,'ndecimaldivider',ndecimaldivider,...
       'PointertoIRcombs',PointertoIRcombs,'IRoriginsfrom',IRoriginsfrom);    
end
    
    
    
    