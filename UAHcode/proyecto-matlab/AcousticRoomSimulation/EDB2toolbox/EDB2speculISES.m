function 	[validISlist,validIScoords,allreflpoints,listguide,listofreflorder] = EDB2speculISES(eddatafile,...
    S,R,lengthNspecmatrix,specorder,visplanesfromoner)
% EDB2speculISES - Finds the valid specular reflections by checking visibility and obstruction.
% Finds the valid specular reflections, given a list of
% potential IS combinations, by checking visibility and obstruction.
%
% Input parameters:
%   eddatafile          See a description in EDB2edgeo, EDB2srgeo and EDB2mainISES
%   S
%   R
%   lengthNspecmatrix
%   specorder
%   visplanesfromoner
%   POTENTIALISES (global)
%   ISCOORDS (global)
%   ORIGINSFROM (global)
%   IVNSPECMATRIX (global)
%
% GLOBAL parameters:
%   SHOWTEXT JJ JJnumbofchars   See EDB2mainISES
%
% Output parameters
%   validISlist         Matrix [nreflections,specorder] of all reflections, with one row for each reflection.
%				        The reflecting plane numbers are given, one for each column.
%   validIScoords       Matrix [nreflections,3] of all image source coordinates.
%   allreflpoints       Matrix [nreflections,(specorder-1)*3] of all hit points in planes.
%   listguide           Matrix [nactiveorders,3] which for each row gives
%                       1. the number of valid reflections for one reflection order
%                       2. the first, and 3. the last
%                       row in validISlist etc which contain reflections of that order.
%                       Which reflection order each row refers to is given
%                       by listofreflorder.
%   listofreflorder     List [nactiveorder,3] which for each row gives the specular
%                       reflection order that the corresponding row in listguide refers to.
%
% Uses functions EDB2chkISvisible EDB2checkobstrpaths
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
% Peter Svensson (svensson@iet.ntnu.no) 20041109
%
% [validISlist,validIScoords,allreflpoints,listguide,listofreflorder] = EDB2speculISES(eddatafile,...
%     S,R,lengthNspecmatrix,specorder,visplanesfromoner);

global SHOWTEXT JJ JJnumbofchars
global POTENTIALISES ISCOORDS ORIGINSFROM IVNSPECMATRIX

eval(['load ',eddatafile])

[ncorners,slask] = size(corners);
[nplanes,ncornersperplane] = size(planecorners);

[n1,n2] = size(POTENTIALISES);
if n2 < specorder
    disp(['   WARNING: The ISEStree-file does not contain high enough order to support'])
    disp(['            spec. refl. order ',int2str(specorder)])
    specorder_maxpossible = n2;
else
    specorder_maxpossible = specorder;
end

%----------------------------------------------
% Find all IS that can be propagated to higher order.
% These must come via planes that S can see, and will be stored in xis1.

thinlist = find(planeisthin==1);

%   ###########################################
%   #                                         #
%   #         First order or higher           #
%   #                                         #
%   ###########################################

validISlist = [];
validIScoords = [];
allreflpoints = [];
listguide        = zeros(specorder,3);
listofreflorder = zeros(specorder,1);
listguide(1,2)   = 1;

obstructtestneeded = (sum(canplaneobstruct)~=0);

for norder = 1:specorder_maxpossible

    if SHOWTEXT >= 3
		disp(['      Reflection order ',int2str(norder)])
	end

    % Start with all the theoretically possible IS
    masterivlist = IVNSPECMATRIX(1:lengthNspecmatrix(norder),norder);
    possiblecombs = POTENTIALISES(masterivlist,1:norder);  
    
    % Pick out only those IS that come via a last refl plane that the
    % receiver is in front of
    okcombs = find(visplanesfromoner(possiblecombs(:,norder))==2);

    masterivlist = masterivlist(okcombs);
    possiblecombs = POTENTIALISES(masterivlist,1:norder); 
    possibleIS = ISCOORDS(masterivlist,:);
    nis = length(masterivlist);
 
    if SHOWTEXT >= 3
		disp(['         ',int2str(nis),' IS found initially for order ',int2str(norder)])
        if SHOWTEXT >= 5
           possiblecombs    
        end
	end

    % Go through iteratively the reflection points, and check if they are
    % inside the reflection plane. Start with the last reflection plane and
    % work backwards.
    for jj = norder:-1:1
        if nis > 0

            if jj == norder
                fromcoords = possibleIS;
                tocoords = R; 
            else    
                
                eval(['tocoords = reflpoints',JJ(jj+1,1:JJnumbofchars(jj+1)),';'])    
                ivreflist = ORIGINSFROM(masterivlist);
                for kk = jj:norder-2
                    ivreflist = ORIGINSFROM(ivreflist);
                end
                fromcoords = ISCOORDS(ivreflist,:);
            end
            
            [hitplanes,reflpoints,edgehits,edgehitpoints,cornerhits,cornerhitpoints] = EDB2chkISvisible(fromcoords,tocoords,planeeqs(possiblecombs(:,jj),4),planenvecs(possiblecombs(:,jj),:),minvals(possiblecombs(:,jj),:),...
				maxvals(possiblecombs(:,jj),:),planecorners(possiblecombs(:,jj),:),corners,ncornersperplanevec(possiblecombs(:,jj)));
            if ~isempty(edgehits) | ~isempty(cornerhits)
                disp('WARNING! An edgehit or cornerhit occurred during the visibility test but this is not')
                disp('         handled correctly yet.')
            end
            eval(['reflpoints',JJ(jj,1:JJnumbofchars(jj)),' = reflpoints;'])
            
            masterivlist = masterivlist(hitplanes);
            possiblecombs = POTENTIALISES(masterivlist,1:norder);  
            possibleIS = ISCOORDS(masterivlist,:);
            if jj < norder
                for kk = jj+1:norder
                   eval(['reflpoints',JJ(kk,1:JJnumbofchars(kk)),' = reflpoints',JJ(kk,1:JJnumbofchars(kk)),'(hitplanes,:);'])
                end
            end
            nis = length(masterivlist);
	
			if SHOWTEXT >= 3
				disp(['         ',int2str(nis),' IS survived the visibility test in refl plane number ',int2str(jj)])
                if SHOWTEXT >= 5
                    possiblecombs    
                end
            end
        end    

    end

    if obstructtestneeded & nis > 0
        % Check obstructions for all the paths: S -> plane1 -> plane2 -> ...
        % -> planeN -> R

        for jj = 1:norder+1
 
            if nis > 0
                if jj==1
                    fromcoords = S;    
                    startplanes = [];    
                else
                    startplanes = possiblecombs(:,jj-1);
                    eval(['fromcoords = reflpoints',JJ(jj-1,1:JJnumbofchars(jj-1)),';'])
                end
                if jj == norder+1
                    tocoords = R;
                    endplanes = [];    
                else
                    eval(['tocoords = reflpoints',JJ(jj,1:JJnumbofchars(jj)),';'])    
                    endplanes = possiblecombs(:,jj);    
                end
                [nonobstructedpaths,nobstructions] = EDB2checkobstrpaths(fromcoords,tocoords,startplanes,endplanes,canplaneobstruct,planeseesplane,...
                    planeeqs,planenvecs,minvals,maxvals,planecorners,corners,ncornersperplanevec,rearsideplane);

                if nobstructions > 0
                    masterivlist = masterivlist(nonobstructedpaths);
                    possiblecombs = POTENTIALISES(masterivlist,1:norder);  
                    possibleIS = ISCOORDS(masterivlist,:);
                    for kk = 1:norder
                        eval(['reflpoints',JJ(kk,1:JJnumbofchars(kk)),' = reflpoints',JJ(kk,1:JJnumbofchars(kk)),'(nonobstructedpaths,:);'])    
                    end
                    nis = length(masterivlist);
                end
                if SHOWTEXT >= 3
			        disp(['         ',int2str(nis),' IS survived the obstruction test for path ',int2str(jj)])
                    if SHOWTEXT >= 5
                       possiblecombs    
                    end
            	end
            end
        
        end
    
    end
    if nis > 0
        if ~isempty(validISlist)
            [n1,n2] = size(validISlist);
            [n3,n4] = size(possiblecombs);
            validISlist = [[validISlist zeros(n1,n4-n2)];(possiblecombs)];
        else
            validISlist = possiblecombs;
        end
        validIScoords = [validIScoords;possibleIS];

        newestreflpoints = [];
        for kk = 1:norder
            eval(['newestreflpoints = [newestreflpoints reflpoints',JJ(kk,1:JJnumbofchars(kk)),'];'])    
        end
        
        if ~isempty(allreflpoints)
            [n1,n2] = size(allreflpoints);
            [n3,n4] = size(newestreflpoints);
            allreflpoints = [[allreflpoints zeros(n1,n4-n2)];newestreflpoints];    
        else
            allreflpoints = newestreflpoints;    
        end
    end

    listguide(norder,1) = nis;
    listguide(norder,3) = listguide(norder,2)+nis-1;
    if norder < specorder
        listguide(norder+1,2) = listguide(norder,3)+1;    
    end
    listofreflorder(norder) = norder;
    
end

iv = find(listguide(:,3)==0);
listguide(iv,2) = zeros(size(iv));

[n1,n2] = size(validISlist);
if n2 < specorder
    validISlist = [validISlist zeros(n1,specorder-n2)];    
end

[n1,n2] = size(allreflpoints);
if n2 < specorder*3
    allreflpoints = [allreflpoints zeros(n1,specorder*3-n2)];    
end

iv = find(listguide(:,1)==0);
listguide(iv,:) = [];
listofreflorder(iv) = [];
