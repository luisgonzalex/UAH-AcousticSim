function outval = EDB2betaoverml(zvec,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec)
% EDB2betaoverml - Integrand function which is called for num. int. of ED IR.
% EDB2betaoverml is copied directly
% into the numerical integration function EDB2quadstep for faster operation
% but EDB2betaoverml is also called directly by EDbwedge1st_int.
%
% Uses no special functions.
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
% Peter Svensson 20040323 (svensson@iet.ntnu.no)
%
% outval = EDB2betaoverml(zvec,rs,rr,zs,zr,ny,sinnyfivec,cosnyfivec);

ml = sqrt( (zvec-zs).^2 + rs.^2 ).*sqrt( (zvec-zr).^2 + rr.^2 );
	
%------------------------------------------------
% We would have liked to use the following code:
%       y = (ml + (zvec-zs).*(zvec-zr))./(rs*rr);
%       expnyeta = real((real(sqrt(y.^2-1)) + y).^ny);
%       coshnyeta = 0.5*( expnyeta + 1./expnyeta);
% but to save three vectors, zvec will be re-used for
% y and expnyeta and coshnyeta

zvec = (ml + (zvec-zs).*(zvec-zr))./(rs*rr);
zvec = real((real(sqrt(zvec.^2-1) + zvec)).^ny);
zvec = 0.5*( zvec + 1./zvec);

outval =          sinnyfivec(1)./(zvec - cosnyfivec(1));
outval = outval + sinnyfivec(2)./(zvec - cosnyfivec(2));
outval = outval + sinnyfivec(3)./(zvec - cosnyfivec(3));
outval = outval + sinnyfivec(4)./(zvec - cosnyfivec(4));         

outval = outval./ml;
