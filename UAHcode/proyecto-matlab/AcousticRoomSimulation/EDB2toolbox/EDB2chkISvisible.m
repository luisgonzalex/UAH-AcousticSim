function [hitplanes,hitpoints,edgehits,edgehitpoints,cornerhits,cornerhitpoints] = EDB2chkISvisible(ISlist,R,planeeqs_lastvalue,planenvecs,minvals,maxvals,planecorners,corners,ncornersperplanevec)
% EDB2chkISvisible - Checks if paths from a set of IS to a set of R pass through their refl. planes. 
% EDB2chkISvisible checks if the paths between a number of IS and a single R, or a number of R, 
% pass through their respective reflecting planes.
%
% Input parameters:
%   ISlist      Matrix, [nIS,3], of the nIS image source coordinates
%   R			Receiver coordinates, [1,3] or [nIS,3]
%   planeeqs_lastvalue  List [nIS,1] of the fourth value of the plane
%               equations.
%   planenvecs,minvals, maxvals, planecorners, corners, ncornersperplanevec
%               Data that should have been taken from the corresponding
%               variables in the eddatafile,see EDB2edgeo for more information.
%               NB!! The matrices planeeqs, minvals, maxvals
%               have been rearranged so that they have nIS rows, and each
%               row contain the data for one specific plane: the reflecting
%               plane for the IS.
%
% Output parameters:
%	 hitplanes		List, [nhits,1], of the planes that have been hit. The values
%	                in this list refer to the index numbers of the input list ISlist
%    hitpoints      List, [nhits,3] of the hitpoint coordinates, for the planes that have been hit
%    edgehits       List, [nedgehits,1] of the planes that were hit right at an edge.
%                   These combinations were marked as hit in the list hitplanes, 
%                   but the extra information in edgehits can possibly be used.
%    edgehitpoints  List, [nedgehits,3] of the hitpoint coordinates for the
%                   planes that were hit at an edge.
%    cornerhits     List, [ncornerhits,1] of the planes that were hit right at a corner.
%                   These combinations were marked as hit in the list hitplanes, 
%                   but the extra information in edgehits can possibly be used.
%    cornerhitpoints  List, [ncornerhits,3] of the hitpoint coordinates for the
%                   planes that were hit at a corner.
%
% Uses the function EDB2poinpla.
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
% Peter Svensson (svensson@iet.ntnu.no) 20050116
% 
% [hitplanes,hitpoints,edgehits,edgehitpoints,cornerhits,cornerhitpoints] = EDB2chkISvisible(ISlist,R,planeeqs_lastvalue,planenvecs,minvals,maxvals,planecorners,corners,ncornersperplanevec)

nsou = size(ISlist,1);
nplanes = size(planeeqs_lastvalue,1);

if size(planeeqs_lastvalue,2) > 1
    error(['Fix the call to EDB2chkISvisible!!!'])    
end

[nrec,slask] = size(R);
if nrec == 1
    onesvec = uint8(ones(nplanes,1));
    R = R(onesvec,:);
    clear onesvec    
end

planenvecs = planenvecs.';

%----------------------------------------------------------------
% First we check if the IS and R are on the same side of the
% respective planes. Then the path can not pass through the plane.

tempmatrix = corners(planecorners(:,1),:);
sseesplanes = sum(  ((ISlist - tempmatrix).').*planenvecs   ).';
rseesplanes = sum(  ((R - tempmatrix).').*planenvecs   ).';
iv1 = uint32(find( (sseesplanes>= (-eps*30) )~=(rseesplanes>= (-eps*30)  ) ));
clear sseesplanes rseesplanes


%--------------------------------------------------------------------------
% Next tests (if there were any planes that had IS and R at different sides)
% 1. Are the vectors IS -> R parallel to the planes?
% 2. Are the hitpoints inside the planes?

if ~isempty(iv1)
	npossible = length(iv1);
    % Below tempmatrix is the direction vector
	tempmatrix = R - ISlist;
    clear R
    dirveclengths = sqrt(sum(tempmatrix.'.^2)).' + eps*100;
    tempmatrix = tempmatrix./dirveclengths(:,ones(1,3));

	% If test == 0, then the vector S -> R runs along the
	% plane.
	
	iv1 = iv1(find(sum(   planenvecs(:,iv1).*(tempmatrix(iv1,:).')   ).'~=0));
    
	if ~isempty(iv1)

		% The last test is if the hitpoint is inside the plane
	
        % Step 1: calculate the hit point using u = the line parameter (with
        % unit meters) from the source towards the receiver.
		udir = ( planeeqs_lastvalue(iv1) - sum( planenvecs(:,iv1).*(ISlist(iv1,:).') ).' );        
        clear planeeqs_lastvalue
        udir = udir./( sum(   planenvecs(:,iv1).*(tempmatrix(iv1,:).')   ).');
        
        % Step 2: check that the hit point is at a shorter distance than
        % the receiver point, and that the hit point is not in the opposite
        % direction than the receiver.

        iv2 = find( udir<(dirveclengths(iv1)*1.001) & udir>=0 );
        clear dirveclengths
        if length(iv2) < length(iv1)
            iv1 = iv1(iv2);
            udir = udir(iv2);
            clear iv2
        end
        
        if ~isempty(iv1)
            % Step 3: calculate the actual xyz coordinates of the hit points
            % Now tempmatrix gets the values of the hit points
            tempmatrix = ISlist(iv1,:) + udir(:,ones(1,3)).*tempmatrix(iv1,:);       
        
            clear ISlist udir

            [hitvec,edgehitvec,cornerhitvec] = EDB2poinpla(tempmatrix,iv1,minvals,maxvals,planecorners,corners,ncornersperplanevec,planenvecs.');
        
    	   	hitplanes = [];
            hitpoints = [];
            edgehits = [];
         	edgehitpoints = [];
            cornerhits = [];
            cornerhitpoints = [];
            if any(any(hitvec)) ~=0
                ivhit = find(hitvec==1);
			    hitplanes = iv1( ivhit );  
             hitpoints = tempmatrix(ivhit,:);
            end
            if ~isempty(edgehitvec)
                ivhit = find(edgehitvec==1);
                edgehits      = iv1( ivhit );
                edgehitpoints = tempmatrix(ivhit,:); 
            end
            if ~isempty(cornerhitvec)
                ivhit = find(cornerhitvec==1);
                cornerhits  = iv1( ivhit );
                cornerhitpoints = tempmatrix(ivhit,:);    
            end
        else
    		hitplanes = [];
            hitpoints = [];
            edgehits = [];
            edgehitpoints = [];
            cornerhits = [];
            cornerhitpoints = [];        
        end
	else
		hitplanes = [];
        hitpoints = [];
        edgehits = [];
        edgehitpoints = [];
        cornerhits = [];
        cornerhitpoints = [];
	end
else    
	hitplanes = [];   
    hitpoints = [];
    edgehits = [];
    edgehitpoints = [];
    cornerhits = [];
    cornerhitpoints = [];
end
