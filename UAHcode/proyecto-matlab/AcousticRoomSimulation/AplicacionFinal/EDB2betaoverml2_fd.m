function integrandvec = EDB2betaoverml2_fd(zvec1,zpos2,k,cylS,cylR,zwedge1,zwedge2,edge1refnvec,edge2refnvec,...
    nyveclist,pathalongplane,Rstart,bc,cair)
% EDB2betaoverml2_fd - Integrand function which is called for num. int. of
% 2nd order ED TF.
%
% Input parameters:
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
%
% Peter Svensson (svensson@iet.ntnu.no) 100718
%
% integrandvec = EDB2betaoverml2_fd(zvec1,zpos2,k,cylS,cylR,cylE2_r1,cylE1_r2,...
% nyveclist,edgelengthlist,dzvec,method,pathalongplane,BigB,Rstart,bc,cair,fs);

multfac = 1/(pathalongplane + 1);

%-----------------------------------------------------------------------
% Extract the individual parameters of the source and the receiver

rS     = cylS(1);
thetaS = cylS(2);
zS     = cylS(3);

rR     = cylR(1);
thetaR = cylR(2);
zR     = cylR(3);

%-----------------------------------------------------------------------
% First we derive the value of the point on edge 2, "zpos2" in edge 1's 
% coordinate system.

edgelength1 = norm( zwedge1(2,:) - zwedge1(1,:) );
edgelength2 = norm( zwedge2(2,:) - zwedge2(1,:) );

zpos2_global = zpos2/edgelength2*(zwedge2(2,:)-zwedge2(1,:)) + zwedge2(1,:);

[rpos2_re1,thetapos2_re1,zpos2_re1,slask,slask,slask] = EDB2coordtrans2(zpos2_global,[],zwedge1,edge1refnvec);

zvec1 = zvec1(:);
zstart_global = zwedge1(1,:);
zend_global = zwedge1(2,:);
wedgevector = zend_global-zstart_global;

zpos1_global = zvec1(:,ones(1,3))/edgelength1.*wedgevector(ones(length(zvec1),1),:) + zstart_global(ones(length(zvec1),1),:);
[rvec1_re2,thetavec1_re2,zvec1_re2,slask,slask,slask] = EDB2coordtrans2(zpos1_global,[],zwedge2,edge2refnvec);

% S2E
S2Edist = sqrt( rS^2 + ( zvec1 - zS ).^2);

% E2R
E2Rdist = sqrt( rR^2 + ( zpos2 - zR ).^2);

thetaedge1_re2 = thetavec1_re2;
thetaedge2_re1 = thetapos2_re1;

	%-----------------------------------------------------------------------
	% In directivity function 2, DF2, we need the quantity
	%
	%   (  ( ze2_re2 - ze1_re2 ) * ( ze2_re2 - zr_re2 ) + n*l )/re1_re2/rr

    % First, E2E is relative to edge 2
    E2Edist = sqrt( (zpos2 - zvec1_re2).^2 + rvec1_re2.^2 );

	% B2 will get the coshnyeta values for DF2
	B2 =  ( ( zpos2 - zvec1_re2 ).*( zpos2 - zR ) + E2Edist.*E2Rdist )./rvec1_re2/rR;
	B2 = ( sqrt( B2.^2-1) + B2 ).^nyveclist(2);
	B2 = real( B2 + 1./B2)/2;
	
	%-----------------------------------------------------------------------
	% In directivity function 1, DF1, we need the quantity
	%
	%   (  ( ze1_re1  -zs_re1 )*( ze1_re1 - ze2_re1 ) + n*m )/rs/re2_re1

    % Now, we need E2E relative to edge 1
    E2Edist = sqrt( (zvec1 - zpos2_re1).^2 + rpos2_re1.^2 );
    
	% B1 will get the coshnyeta values for DF1
	
%  	B1 = (( B4 - zS ).*(  B4 - zedge2_re1 ) + S2Edist.*E2Edist)/rS./redge2_re1;
%  	B1 = (( B4 - zS ).*(  zvec1 - zedge2_re1 ) + S2Edist.*E2Edist)/rS./redge2_re1;
%  	B1 = (( B4 - zS ).*(  zvec1 - zedge2_re1 ) + S2Edist.*E2Edist)/rS./rpos2_re1;
 	B1 = (( zvec1 - zS ).*(  zvec1 - zpos2_re1 ) + S2Edist.*E2Edist)/rS./rpos2_re1;
 	B1 = ( real(sqrt( B1.^2-1)) + B1 ).^nyveclist(1);
 	B1 = ( B1 + 1./B1)/2;

%     if min(B1) < 1,
%         disp(['   Min B1 = ',num2str(min(B1))])  
%         pause
%     end
    
%     clear zedge2_re1 % redge2_re1
	
% Calculate DF1.*DF2

if pathalongplane == 0,

    % Note that here we use E2E re. edge 1!
    %B3 = multfac*(-nyveclist(1)/4/pi*dzvec(1))*(-nyveclist(2)/4/pi*dzvec(2))*...
%     B3 = multfac*nyveclist(1)*nyveclist(2)/16/pi^2*dzvec(1)*dzvec(2)*...
% 	              ( sin(nyveclist(1)*(pi + thetaS + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetaS + thetaedge2_re1   ))) + ...
% 				    sin(nyveclist(1)*(pi + thetaS - thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetaS - thetaedge2_re1   ))) + ...
% 				    sin(nyveclist(1)*(pi - thetaS + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi - thetaS + thetaedge2_re1   ))) + ...
% 					sin(nyveclist(1)*(pi - thetaS - thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi - thetaS - thetaedge2_re1   ))))./S2Edist./E2Edist;
    B3 = multfac*nyveclist(1)*nyveclist(2)/16/pi^2*...
	              ( sin(nyveclist(1)*(pi + thetaS + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetaS + thetaedge2_re1   ))) + ...
				    sin(nyveclist(1)*(pi + thetaS - thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetaS - thetaedge2_re1   ))) + ...
				    sin(nyveclist(1)*(pi - thetaS + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi - thetaS + thetaedge2_re1   ))) + ...
					sin(nyveclist(1)*(pi - thetaS - thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi - thetaS - thetaedge2_re1   ))))./S2Edist./E2Edist;

    B3 = 	  B3.*( sin(nyveclist(2)*(pi + thetaedge1_re2 + thetaR))./(B2-cos(nyveclist(2)*(pi + thetaedge1_re2 + thetaR))) + ...
				    sin(nyveclist(2)*(pi + thetaedge1_re2 - thetaR))./(B2-cos(nyveclist(2)*(pi + thetaedge1_re2 - thetaR))) + ...
				    sin(nyveclist(2)*(pi - thetaedge1_re2 + thetaR))./(B2-cos(nyveclist(2)*(pi - thetaedge1_re2 + thetaR))) + ...
					sin(nyveclist(2)*(pi - thetaedge1_re2 - thetaR))./(B2-cos(nyveclist(2)*(pi - thetaedge1_re2 - thetaR))))./E2Rdist;
                
    else

    %B3 = multfac*(-nyveclist(1)/4/pi*dzvec(1))*(-nyveclist(2)/4/pi*dzvec(2))*...
    if thetaS == 0,
%         B3 = multfac*nyveclist(1)*nyveclist(2)/2/pi^2*dzvec(1)*dzvec(2)*...
% 		              ( sin(nyveclist(1)*(pi + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetaedge2_re1   ))))./S2Edist./E2Edist;
        B3 = multfac*nyveclist(1)*nyveclist(2)/2/pi^2*...
		              ( sin(nyveclist(1)*(pi + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetaedge2_re1   ))))./S2Edist./E2Edist;
    else
        
%         B3 = multfac*nyveclist(1)*nyveclist(2)/4/pi^2*dzvec(1)*dzvec(2)*...
% 		              ( sin(nyveclist(1)*(pi + thetaS + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetaS + thetaedge2_re1   ))) + ...
% 					    sin(nyveclist(1)*(pi - thetaS + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi - thetaS + thetaedge2_re1   ))))./S2Edist./E2Edist;
        B3 = multfac*nyveclist(1)*nyveclist(2)/4/pi^2*...
		              ( sin(nyveclist(1)*(pi + thetaS + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetaS + thetaedge2_re1   ))) + ...
					    sin(nyveclist(1)*(pi - thetaS + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi - thetaS + thetaedge2_re1   ))))./S2Edist./E2Edist;
                    
    end    

    
    B3 = 	      B3.*( sin(nyveclist(2)*(pi + thetaedge1_re2 + thetaR))./(B2-cos(nyveclist(2)*(pi + thetaedge1_re2 + thetaR))) + ...
					    sin(nyveclist(2)*(pi + thetaedge1_re2 - thetaR))./(B2-cos(nyveclist(2)*(pi + thetaedge1_re2 - thetaR))))./E2Rdist;

end

integrandvec = B3.*exp(-1i*k*(S2Edist+E2Edist+E2Rdist-Rstart));
