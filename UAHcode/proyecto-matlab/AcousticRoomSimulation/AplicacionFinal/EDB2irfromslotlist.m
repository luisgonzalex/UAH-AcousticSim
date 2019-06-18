function ir = EDB2irfromslotlist(slotnumbervec,ampvec,method)
% EDB2irfromslotlist - Creates an IR from a list of sample numbers and amplitudes
%
% Input parameters:
%	slotnumbervec	List of sample numbers which don't need to be integers.
%                   For non-integers, each pulse is divided into two
%                   neighboring sample slots. This gives an accurate phase
%                   response, but a high-frequency roll-off.
%	ampvec			List of corresponding amplitudes
%   method          (optional) 'oneslot' or 'twoslot'. Default is 'twoslot'.
%                       'oneslot' means that each pulse is placed in the
%                       nearest sample slot.
%                       'twoslot' means that each pulse is divided between
%                       two neighbour sample slots.
%
% Output parameters:
%   ir              A list containing the impulse response that is created
%                   by a number of pulses in the input lists. 
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
% Peter Svensson (svensson@iet.ntnu.no) 20041201
%
% ir = EDB2irfromslotlist(slotnumbervec,ampvec,method);

if nargin < 3
    method = 'twoslots';    
end
if method(1) == 'o' | method(1) == 'O'
    slotnumbervec = round(slotnumbervec);   
end

[n1,n2] = size(slotnumbervec);
if n2 ~= 1
	slotnumbervec = reshape(slotnumbervec,n1*n2,1);
end
[n1,n2] = size(ampvec);
if n2 ~= 1
	ampvec = reshape(ampvec,n1*n2,1);
end

% B1 is the sampleslot numbers in integer number
B1 = floor(slotnumbervec);

% slotnumbervec is a value between 0 and 1 stating how much of dir that should be added
% to the first sample slot, i.e. the one given in B1. The rest of dir should
% be added to the following sample slot.
slotnumbervec = slotnumbervec - B1;

% Sort the values in B1, because the values occur many times.
% The slotnumbervec and ampvec values must be sorted in the same order.
% The isteps vactor will contain the occurences of a step in slotnumbervec.

[B1,sortvec] = sort(B1);
isteps = find(diff(B1) ~= 0) +1;

slotnumbervec = slotnumbervec(sortvec);
ampvec = ampvec(sortvec);

nelems = length(B1);
nir = B1(nelems) + 2;
ir = zeros( nir,1 );

% Now we can sum all the contributions that should end up in the same
% sample slot. There is a for loop over the steps given in isteps, but
% the first and last values are treated outside the for loop.

slotnumber = B1(1);

if isempty(isteps) | nelems == 1
	ir(slotnumber) = ir(slotnumber) + sum( ampvec.*(1-slotnumbervec) );
	ir(slotnumber+1) = ir(slotnumber+1) + sum( ampvec.*slotnumbervec );	
else
	ir(slotnumber) = ir(slotnumber) + sum( ampvec(1:isteps(1)-1).*(1-slotnumbervec(1:isteps(1)-1)) );
	ir(slotnumber+1) = ir(slotnumber+1) + sum( ampvec(1:isteps(1)-1).*slotnumbervec(1:isteps(1)-1) );

	slotnumbers = B1( isteps );
	for jj = 1:length(isteps)-1
		ir(slotnumbers(jj))   = ir(slotnumbers(jj)) + ...
		sum(ampvec(isteps(jj):isteps(jj+1)-1).*(1-slotnumbervec(isteps(jj):isteps(jj+1)-1)) );
		ir(slotnumbers(jj)+1) = ir(slotnumbers(jj)+1) + ...
		sum( ampvec(isteps(jj):isteps(jj+1)-1).*slotnumbervec(isteps(jj):isteps(jj+1)-1) );
	end
	
	startindex = isteps(length(isteps));
	slotnumber = B1( startindex );
	ir(slotnumber) = ir(slotnumber) + ...
		sum( ampvec(startindex:nelems).*(1-slotnumbervec(startindex:nelems)) );
	ir(slotnumber+1) = ir(slotnumber+1) + ...
		sum( ampvec(startindex:nelems).*slotnumbervec(startindex:nelems) );
end
