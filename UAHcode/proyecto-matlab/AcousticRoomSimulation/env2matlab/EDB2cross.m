function c = EDB2cross(a,b)
% EDB2cross - Stripped down version of Matlab's built in function cross.
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
%  c = EDB2cross(a,b)

%CROSS  Vector cross product.
%   C = CROSS(A,B) returns the cross product of the vectors
%   A and B.  That is, C = A x B.  A and B must be 3 element
%   vectors.
%
%   C = CROSS(A,B) returns the cross product of A and B along the
%   first dimension of length 3.
%

%   Clay M. Thompson
%   updated 12-21-94, Denise Chen
%   Copyright 1984-2002 The MathWorks, Inc.
%   $Revision: 1.1 $  $Date: 2016/07/27 12:10:06 $

% Special case: A and B are vectors
%rowvec = 0;
%if ndims(a)==2 & ndims(b)==2
%    if size(a,1)==1, a = a(:); rowvec = 1; end
%    if size(b,1)==1, b = b(:); rowvec = 1; end
%end;

% Calculate cross product
c = [a(2,:).*b(3,:)-a(3,:).*b(2,:)
     a(3,:).*b(1,:)-a(1,:).*b(3,:)
     a(1,:).*b(2,:)-a(2,:).*b(1,:)];
c = reshape(c,size(a));
