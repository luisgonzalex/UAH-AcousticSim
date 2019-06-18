function [rs,thetas,zs] = EDB2coordtrans1(xsou,xwedge,nvec1)
% EDB2coordtrans1 - Transforms one set of cartesian coordinates to edge-related cylindrical coordinates.
% The cyl. coord. system is defined so that:
%     A z-axis is placed along the edge, from the given endpoint 1 to the given
%     endpoint 2.
%     The origo of the cyl. syst. will be edge endpoint 1.
%     The theta-angles of the cyl. coord. syst. will refer to the
%     reference plane of the edge.
% The ref. plane of the edge is described by its normal vector.
% NB! The order of the edge points is important!!
% The vector going from xwedge(1,:) to xwedge(2,:) must be
% oriented so that if the RH thumb is along this vector, the tips
% of the fingers "come out of" the open face of plane1, i.e. where nvec1
% is the normal vector.
%
% Input parameters:
%	xsou		Matrix, [n1,3] of cartesian coordinates of n1 points. 
%	xwedge 		Matrix, [2,3], with the cartesian coordinates of the two 
%				wedge end points: [xw1 yw1 zw1;xw2 yw2 zw2].
%   nvec1       List, [1,3], with the normal vector of the reference plane
%               of the edge.
%
% Output parameters:
%	rs, thetas, zs		cyl. coord. of the points in xsou
%
% Uses the function EDB2cross
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
% Peter Svensson (svensson@iet.ntnu.no) 20061118
%
% [rs,thetas,zs] = EDB2coordtrans1(xsou,xwedge,nvec1)

npoints = size(xsou,1);
if npoints == 1
    xneworigo = xwedge(1,:);

    xknown1 = xwedge(2,:) - xneworigo;
    xknown1 = xknown1 / norm(xknown1);

    A = [nvec1(2)*xknown1(3)-nvec1(3)*xknown1(2) ; nvec1(3)*xknown1(1)-nvec1(1)*xknown1(3) ; nvec1(1)*xknown1(2)-nvec1(2)*xknown1(1)];
    A = inv([xknown1.' A nvec1.']);

    xsou =    (A([2 3 1],:)*( xsou.' - xneworigo.' )).';

    rs = norm(xsou(1:2));
    zs = xsou(:,3);
    thetas = 0;
    if rs > 0
        thetas = real( acos( xsou(1)./rs ).*( xsou(2) ~= 0) );
        thetas = thetas + pi*( (xsou(2)==0) & xsou(1) < 0 );
        thetas = thetas.*( xsou(2) >=0 ) + (2*pi - thetas).*( xsou(2) < 0 );
    end
    
else
    xneworigo = xwedge(1,:);

    xknown1 = xwedge(2,:) - xneworigo;
    xknown1 = xknown1 / sqrt( sum( xknown1.^2 ));


    A = [0 1 0;0 0 1;1 0 0]*inv([xknown1.' EDB2cross(nvec1.',xknown1.') nvec1.']);

    [npoints,slask] = size(xsou);
    xsou =    (A*( xsou.' - xneworigo(ones(npoints,1),:).' )).';

    rs = sqrt( sum(xsou(:,1:2).'.^2) ).';
    zs = xsou(:,3);
    thetas = zeros(npoints,1);
    iv = find(rs>0);
    if ~isempty(iv)
        thetas(iv) = real( acos( xsou(iv,1)./rs(iv) ).*( xsou(iv,2) ~= 0) );
        thetas(iv) = thetas(iv) + pi*( (xsou(iv,2)==0) & xsou(iv,1) < 0 );
        thetas(iv) = thetas(iv).*( xsou(iv,2) >=0 ) + (2*pi - thetas(iv)).*( xsou(iv,2) < 0 );
    end

end