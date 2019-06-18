function [edgepointcoords,weightvec,edgenumberlist] = EDB2getedgepoints(edgestartcoords,edgeendcoords,edgelengths,nedgesubs,returnedgenumbers)
% EDB2getedgepoints - Calculates a number of edge coordinates.
%
% Input parameters:
%	edgestartcoords         Matrix, [nedges,3], of the startpoint coordinates of nedges edges.
%	edgeendcoords           Matrix, [nedges,3], of the endpoint coordinates of nedges edges.
%   edgelengths             List, [nedges,1], of edge lengths.
%	nedgesubs               The desired number of subdivisions per edge.
%   returnedgenumbers (optional)    If this optional parameter is given the
%                           value 1, the output parameter edgenumberlist will be non-empty.
%
%   NB! The input parameters edgestartcoords, edgeendcoords and edgelengths
%   should have been passed on directly from the eddata file.
%   The parameter nedgesubs should have been taken from the setup file.
%
% Output parameters:
%	edgepointcoords         Matrix, [nedges*nedgesubs,3], of the midpoints
%	                        of the edge subdivision segments.
%	weightvec               List, [nedges*nedgesubs,3], with the weights
%	                        (1,2,4,8,...) per edge subdivision segment.
%	edgenumberlist          List, [nedges*nedgesubs,1], of which edge
%	                        number each subdivision segment originated from.
%                           See also the input parameter returnedgenumbers.
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
% Peter Svensson (svensson@iet.ntnu.no) 20030503
%
% [edgepointcoords,weightvec,edgenumberlist] = EDB2getedgepoints(edgestartcoords,edgeendcoords,edgelengths,nedgesubs,returnedgenumbers);

if nargin < 5
    returnedgenumbers = 0;
end

[nedges,slask] = size(edgestartcoords);

if nedgesubs == 1
    edgepointcoords = (edgestartcoords+edgeendcoords)/2;
    weightvec = uint8(ones(nedges,1));
	if returnedgenumbers == 1
		if nedges < 256
           edgenumberlist = uint8([1:nedges]);
		elseif nedges < 65536
           edgenumberlist = uint16([1:nedges]);
		else
           edgenumberlist = uint32([1:nedges]);
		end
	else
        edgenumberlist = [];    
	end    
    return
end

ntot = nedges*nedgesubs;

edgestart = zeros(ntot,3);
edgeend   = zeros(ntot,3);

onesvec1 = ones(nedgesubs,1);
onesvec2 = ones(1,3);
onesvec3 = ones(1,nedges);

%-------------------------------------------------
% Start points

A = edgestartcoords(:,1).';
edgestart(:,1) = reshape(A(onesvec1,:),ntot,1);

A = edgestartcoords(:,2).';
edgestart(:,2) = reshape(A(onesvec1,:),ntot,1);

A = edgestartcoords(:,3).';
edgestart(:,3) = reshape(A(onesvec1,:),ntot,1);

%-------------------------------------------------
% End points

A = edgeendcoords(:,1).';
edgeend(:,1) = reshape(A(onesvec1,:),ntot,1);

A = edgeendcoords(:,2).';
edgeend(:,2) = reshape(A(onesvec1,:),ntot,1);

A = edgeendcoords(:,3).';
edgeend(:,3) = reshape(A(onesvec1,:),ntot,1);

%------------------------------------------------------------
% Spread the edge points along the edge, including the end
% points, but nudge them in a bit.

displacement = 1e-1;

edgedividervec = [displacement 1:nedgesubs-2 nedgesubs-1-displacement].'/(nedgesubs-1);
edgedividervec = reshape(edgedividervec(:,onesvec3),ntot,1);

edgepointcoords = edgestart + (edgeend - edgestart).*edgedividervec(:,onesvec2);

%------------------------------------------------------------
% Extra output data

if nedgesubs < 8
   weightvec = uint8(2.^[0:nedgesubs-1].');
elseif nedgesubs < 16
   weightvec = uint16(2.^[0:nedgesubs-1].');
else
   weightvec = uint32(2.^[0:nedgesubs-1].');
end

weightvec = reshape(weightvec(:,onesvec3),ntot,1);

if returnedgenumbers == 1
	if nedges < 256
       edgenumberlist = uint8([1:nedges]);
	elseif nedges < 65536
       edgenumberlist = uint16([1:nedges]);
	else
       edgenumberlist = uint32([1:nedges]);
	end
    edgenumberlist = reshape(edgenumberlist(onesvec1,:),ntot,1);
else
    edgenumberlist = [];    
end
