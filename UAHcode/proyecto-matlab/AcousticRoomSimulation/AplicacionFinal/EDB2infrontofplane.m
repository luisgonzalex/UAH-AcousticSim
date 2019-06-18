function pointinfront = EDB2infrontofplane(pointcoords,planenvecs,planecoco,backupplanecoco,cornernumb,planenumb)
% EDB2infrontofplane - Checks if a point is in front of, in plane with, or behind a list of planes. 
% Alternatively, a single plane is checked versus a list of points.
% Alternatively, it checks n points versus n planes.
%
% Input parameters:
%	pointcoords				A list of points, or a single point, as [x1 y1 z1;x2 y2 z2;...]
%	planenvecs				A single planenvec, or a list of planenvecs, as
%								[nx1,ny1,nz1;nx2,ny2,nz2;...]
%	planecoco				A single, or a list of, points belonging to the plane
%	backupplanecoco			A single, or a list of, points belonging to the plane, but different
%								points than in planecoco
%   cornernumb              (optional) A list of the values [1:ncorners] for huge problems.
%                           If this list is used, the other input parameters pointcoords etc must be expanded.
%   planenumb              (optional) A list of the values [1:nplanes] for huge problems.
%                           If this list is used, the other input parameters pointcoords etc must be expanded.
%
% Output parameters:
%	pointinfront			A list of values, +1 (in front) 0 (in plane) or -1 (behind)
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
% Peter Svensson (svensson@iet.ntnu.no) 20090103
%
% pointinfront = EDB2infrontofplane(pointcoords,planenvecs,planecoco,backupplanecoco,cornernumb,planenumb);

geomacc = 1e-6;
batchsize = 2.5e6;

[nplanes,slask] = size(planenvecs);
[npoints,slask] = size(pointcoords);

if npoints == 1

    onesvec = ones(nplanes,1);
	bigpointcoords = pointcoords(onesvec,:);

    vectors = bigpointcoords - planecoco;
	lvectors = sqrt(sum( (vectors.').^2 )).';
	iv = find( lvectors < geomacc );
	if ~isempty(iv)
		vectors(iv,:) = bigpointcoords(iv,:) - backupplanecoco(iv,:);
		if length(iv) > 1
			lvectors(iv) = sqrt(sum( (vectors(iv,:).').^2 )).';
		else
			lvectors(iv) = norm(vectors(iv,:));		
		end	
	end

	pointrelplane = sum((vectors.').*(planenvecs.'));
    pointinfront = double(pointrelplane > geomacc);
	iv = find(pointrelplane < -geomacc);
	if ~isempty(iv)
		pointinfront(iv) = -1*ones(size(iv));
	end
elseif npoints > 1 & nplanes==1
 	onesvec = ones(npoints,1);
	bigplanenvecs = planenvecs(onesvec,:);
	bigplanecoco = planecoco(onesvec,:);
	vectors = pointcoords - bigplanecoco;
	lvectors = sqrt(sum( (vectors.').^2 )).';
	iv = find( lvectors < geomacc );
	if ~isempty(iv)
		bigbackupplanecoco = backupplanecoco(onesvec,:);
		vectors(iv,:) = pointcoords(iv,:) - bigbackupplanecoco(iv,:);
	end
	pointrelplane = sum((vectors.').*(bigplanenvecs.'));
	pointinfront = double(pointrelplane > geomacc);
	iv = find(pointrelplane < -geomacc);
	if ~isempty(iv)
		pointinfront(iv) = -1*ones(size(iv));
	end
else
    if nargin == 4
		vectors = pointcoords - planecoco;
		lvectors = sqrt(sum( (vectors.').^2 )).';
		iv = find( lvectors < geomacc );
		if ~isempty(iv)
			vectors(iv,:) = pointcoords(iv,:) - backupplanecoco(iv,:);
		end
		pointrelplane = sum((vectors.').*(planenvecs.'));
		pointinfront = double(pointrelplane > geomacc);
		iv = find(pointrelplane < -geomacc);
		if ~isempty(iv)
			pointinfront(iv) = -1*ones(size(iv));
    	end
    elseif nargin == 6
        nplanes = length(planenumb);
        ncorners = length(cornernumb);
        corners = pointcoords; clear pointcoords
        planecorners = planecoco; clear planecoco

        planenumb = reshape(planenumb(:,ones(1,ncorners)),nplanes*ncorners,1);

        cornernumb = reshape(cornernumb(ones(nplanes,1),:),nplanes*ncorners,1);
        
        nproblemsize = size(cornernumb,1);
        
        if nproblemsize <= batchsize        
            vectors = corners(cornernumb,:) - corners(planecorners(planenumb,1),:);
            iv = find( sqrt(sum( (vectors.').^2 )).' < geomacc );
        else
           nbatches = ceil(nproblemsize/batchsize);           
           vectors = zeros(nproblemsize,3);
           lengthvec = zeros(nproblemsize,1);
           for ii = 1:nbatches
               ivbatch = [1+(ii-1)*batchsize : min([ii*batchsize nproblemsize])];
               vectors(ivbatch,:) = corners(cornernumb(ivbatch),:) - corners(planecorners(planenumb(ivbatch),1),:);               
               lengthvec(ivbatch) = sqrt(sum( (vectors(ivbatch,:).').^2 ));
           end
            
            iv = find( lengthvec.' < geomacc );
           
        end        

        if ~isempty(iv)
            backupplanecoco = corners(planecorners(planenumb(iv),2),:);
            cornernumb = cornernumb(iv);
            vectors(iv,:) = corners(cornernumb,:) - backupplanecoco;
            clear backupplanecoco   corners
        end
            
        % Here vectors doubles so that vectors gets the value of
        % pointrelplane. For the batches case, only column 1 is reused.

        if nproblemsize <= batchsize        
            vectors = sum((vectors.').*(planenvecs(planenumb,:).'));            
        else
           for ii = 1:nbatches
               ivbatch = [1+(ii-1)*batchsize : min([ii*batchsize nproblemsize])];
               vectors(ivbatch,1) = sum((vectors(ivbatch,:).').*(planenvecs(planenumb(ivbatch),:).'));                
           end
           vectors = vectors(:,1);
        end
        clear planenumb planenvecs
        
		pointinfront = double(vectors > geomacc);
		iv = find(vectors < -geomacc);
		if ~isempty(iv)
			pointinfront(iv) = -1*ones(size(iv));
    	end       
    end
end

pointinfront = pointinfront.';
