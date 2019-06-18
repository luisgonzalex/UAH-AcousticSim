function ivmatrix = EDB2creindexmatrix(ndivvec)
% EDB2creindeixmatrix creates a matrix with index numbers.
%
% Input parameters:
%   ndivvec     A matrix, [1,specorder], with the maximum counter
%               number for each dimension.
%
% Output parameters:
%   ivmatrix    A matrix, [prod(ndivvec),specorder] of all possible
%               combinations of the integers
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
% Peter Svensson (svensson@iet.ntnu.no) 20050505
%
% ivmatrix = EDB2creindexmatrix(ndivvec)

n = length(ndivvec);

maxval = max(ndivvec);
if maxval > 2^32
    error(['ERROR: This version of EDB2creindexmatrix can not create such large matrics'])    
end

if n == 2
	iv1 = uint32([1:ndivvec(1)].');
	iv2 = uint32([1:ndivvec(2)]);
	iv1 = iv1(:,uint8(ones(1,ndivvec(2))));
	iv2 = iv2(uint8(ones(ndivvec(1),1)),:);

    ivmatrix = [reshape(iv1.',prod(ndivvec),1) reshape(iv2.',prod(ndivvec),1)];
elseif n >= 3
    ivmatrix = EDB2creindexmatrix(ndivvec(2:n));
    ivmatrix = repmat(ivmatrix,[ndivvec(1),1]);
    ivfirstcol = uint32([1:ndivvec(1)].');
    ivfirstcol = ivfirstcol(:,uint8(ones(1,prod(ndivvec(2:n)))));
    ivfirstcol = reshape(ivfirstcol.',prod(ndivvec),1);
    ivmatrix = [ivfirstcol ivmatrix];
end
