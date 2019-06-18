function integrandvec = EDB2wedge3rd_intcorefd(zvec1,zpos2,zpos3,k,cylS,cylR,wedgeparams,edgelengths,edge1refnvec,edge2refnvec,edge3refnvec,...
    nyveclist,pathalongplane,R_irstart,bc,cair);
%
% Peter Svensson (svensson@tele.ntnu.no) 050117
%   svensson@iet.ntnu.no 13 March 2011
%
% tfvalue = EDB2wedge3rd_intcorefd(zvec1,zpos2,zpos3,k,cylS,cylR,wedgeparams,edgelengths,edge1refnvec,edge2refnvec,edge3refnvec,...
% nyveclist,pathalongplane,BigB,R_irstart,bc,cair);

% % t00 = clock;

multfac = 1./(pathalongplane + 1);

%-----------------------------------------------------------------------
% Extract the individual parameters of the source and the receiver

rS     = cylS(1);
thetaS = cylS(2);
zS     = cylS(3);

rR     = cylR(1);
thetaR = cylR(2);
zR     = cylR(3);

%-----------------------------------------------------------------------
% We extract the edge coordinates relative to the other edges.
% First, calculate the relative position along each edge.
% Then, use interpolation to find the r, theta, z values from the edge end
% points, which we know:
% wedgeparams = [cylE2_r1;cylE1_r2;cylE3_r2;cylE2_r3]

% Edge 3 (single point) relative to edge 2

E3re2 = wedgeparams(5:6,:);
edge3relpos = zpos3/edgelengths(3);
rpos3_re2 = edge3relpos*(E3re2(2,1)-E3re2(1,1)) + E3re2(1,1);
thetapos3_re2 = edge3relpos*(E3re2(2,2)-E3re2(1,2)) + E3re2(1,2);
zpos3_re2 = edge3relpos*(E3re2(2,3)-E3re2(1,3)) + E3re2(1,3);

% Edge 2 (single point) relative to edge 3

E2re3 = wedgeparams(7:8,:);
edge2relpos = zpos2/edgelengths(2);
rpos2_re3     = edge2relpos*(E2re3(2,1)-E2re3(1,1)) + E2re3(1,1);
thetapos2_re3 = edge2relpos*(E2re3(2,2)-E2re3(1,2)) + E2re3(1,2);
zpos2_re3     = edge2relpos*(E2re3(2,3)-E2re3(1,3)) + E2re3(1,3);

% Edge 2 (single point) relative to edge 1

E2re1 = wedgeparams(1:2,:);
rpos2_re1     = edge2relpos*(E2re1(2,1)-E2re1(1,1)) + E2re1(1,1);
thetapos2_re1 = edge2relpos*(E2re1(2,2)-E2re1(1,2)) + E2re1(1,2);
zpos2_re1     = edge2relpos*(E2re1(2,3)-E2re1(1,3)) + E2re1(1,3);

% Edge 1 (vector) relative to edge 2

E1re2 = wedgeparams(3:4,:);
edge1relpos = zvec1/edgelengths(1);
rvec1_re2     = edge1relpos*(E1re2(2,1)-E1re2(1,1)) + E1re2(1,1);
thetavec1_re2 = edge1relpos*(E1re2(2,2)-E1re2(1,2)) + E1re2(1,2);
zvec1_re2     = edge1relpos*(E1re2(2,3)-E1re2(1,3)) + E1re2(1,3);

% % % edgelength1 = norm( zwedge1(2,:) - zwedge1(1,:) );
% % % edgelength2 = norm( zwedge2(2,:) - zwedge2(1,:) );
% % % edgelength3 = norm( zwedge3(2,:) - zwedge3(1,:) );
% % % 
% % % t02 = etime(clock,t00);
% % % 
% % % zpos2_global = zpos2/edgelength2*(zwedge2(2,:)-zwedge2(1,:)) + zwedge2(1,:);
% % % [rpos2_re1,thetapos2_re1,zpos2_re1,slask,slask,slask] = EDB2coordtrans(zpos2_global,[],zwedge1,edge1refnvec);
% % % 
% % % zpos3_global = zpos3/edgelength3*(zwedge3(2,:)-zwedge3(1,:)) + zwedge3(1,:);
% % % [rpos3_re2,thetapos3_re2,zpos3_re2,slask,slask,slask] = EDB2coordtrans(zpos3_global,[],zwedge2,edge2refnvec);
% % % [rpos2_re3,thetapos2_re3,zpos2_re3,slask,slask,slask] = EDB2coordtrans(zpos2_global,[],zwedge3,edge3refnvec);
% % % 
% % % zvec1 = zvec1(:);
% % % zstart_global = zwedge1(1,:);
% % % zend_global = zwedge1(2,:);
% % % wedgevector = zend_global-zstart_global;
% % % 
% % % zpos1_global = zvec1(:,ones(1,3))/edgelength1.*wedgevector(ones(length(zvec1),1),:) + zstart_global(ones(length(zvec1),1),:);
% % % [rvec1_re2,thetavec1_re2,zvec1_re2,slask,slask,slask] = EDB2coordtrans(zpos1_global,[],zwedge2,edge2refnvec);

%----------------------------------------------------

% S2E
S2Edist = sqrt( rS^2 + ( zvec1 - zS ).^2);

% E2R
E2Rdist = sqrt( rR^2 + ( zpos3 - zR ).^2);

% thetaedge1_re2 = thetavec1_re2;
% thetaedge2_re1 = thetapos2_re1;
% thetaedge2_re3 = thetapos2_re3;
% thetaedge3_re2 = thetapos3_re2;

% t02 = etime(clock,t00);

	%-----------------------------------------------------------------------
	% In directivity function 3, DF3, we need the quantity
	%
	%   (  ( ze3_re3 - ze2_re3 ) * ( ze3_re3 - zr_re3 ) + n*l )/re2_re3/rr

    % E2Edist2 is relative to edge 3 (paths onto edge 3, from edge 2)
    E2Edist2 = sqrt( (zpos3 - zpos2_re3).^2 + rpos2_re3.^2 );

	% B2 will get the coshnyeta values for DF2
	B3 =  ( ( zpos3 - zpos2_re3 ).*( zpos3 - zR ) + E2Edist2.*E2Rdist )./rpos2_re3/rR;
	B3 = ( sqrt( B3.^2-1) + B3 ).^nyveclist(3);
	B3 = real( B3 + 1./B3)/2;

	%-----------------------------------------------------------------------
	% In directivity function 2, DF2, we need the quantity
	%
	%   (  ( ze2_re2 - ze1_re2 ) * ( ze2_re2 - zr_re2 ) + n*l )/re1_re2/rr

    % E2Edist1 is relative to edge 2 (paths onto edge 2, from edge 1)
    % E2Edist2 is relative to edge 2 (paths from edge 2, to edge 3)
    E2Edist1 = sqrt( (zpos2 - zvec1_re2).^2 + rvec1_re2.^2 );
    E2Edist2mod = sqrt( (zpos2 - zpos3_re2).^2 + rpos3_re2.^2 );

	% B2 will get the coshnyeta values for DF2
	B2 =  ( ( zpos2 - zvec1_re2 ).*( zpos2 - zpos3_re2 ) + E2Edist1.*E2Edist2 )./rvec1_re2/rpos3_re2;
	B2 = ( sqrt( B2.^2-1) + B2 ).^nyveclist(2);
	B2 = real( B2 + 1./B2)/2;
	
	%-----------------------------------------------------------------------
	% In directivity function 1, DF1, we need the quantity
	%
	%   (  ( ze1_re1 - zs_re1 )*( ze1_re1 - ze2_re1 ) + n*m )/rs/re2_re1

    % E2E is relative to edge 1 (paths from edge 1, to edge 2)
    E2Edist1mod = sqrt( (zvec1 - zpos2_re1).^2 + rpos2_re1.^2 );
 
	% B1 will get the coshnyeta values for DF1
	B1 = (( zvec1 - zS ).*(  zvec1 - zpos2_re1 ) + S2Edist.*E2Edist1)/rS./rpos2_re1;
 	B1 = ( real(sqrt( B1.^2-1)) + B1 ).^nyveclist(1);
 	B1 = ( B1 + 1./B1)/2;

%     if min(B1) < 1,
%         disp(['   Min B1 = ',num2str(min(B1))])  
%         pause
%     end
    
%     clear zedge2_re1 % redge2_re1
	
% Calculate DF1.*DF2.*DF3

% % t02 = etime(clock,t00);

if pathalongplane == 0,

    % Note that here we use E2E re. edge 1!
    %B3 = multfac*(-nyveclist(1)/4/pi*dzvec(1))*(-nyveclist(2)/4/pi*dzvec(2))*...
%     B3 = multfac*nyveclist(1)*nyveclist(2)/16/pi^2*dzvec(1)*dzvec(2)*...
% 	              ( sin(nyveclist(1)*(pi + thetaS + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetaS + thetaedge2_re1   ))) + ...
% 				    sin(nyveclist(1)*(pi + thetaS - thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetaS - thetaedge2_re1   ))) + ...
% 				    sin(nyveclist(1)*(pi - thetaS + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi - thetaS + thetaedge2_re1   ))) + ...
% 					sin(nyveclist(1)*(pi - thetaS - thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi - thetaS - thetaedge2_re1   ))))./S2Edist./E2Edist;
error('Not done yet!')
B4 = multfac*nyveclist(1)*nyveclist(2)/16/pi^2*...
	              ( sin(nyveclist(1)*(pi + thetaS + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetaS + thetaedge2_re1   ))) + ...
				    sin(nyveclist(1)*(pi + thetaS - thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetaS - thetaedge2_re1   ))) + ...
				    sin(nyveclist(1)*(pi - thetaS + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi - thetaS + thetaedge2_re1   ))) + ...
					sin(nyveclist(1)*(pi - thetaS - thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi - thetaS - thetaedge2_re1   ))))./S2Edist./E2Edist1;

    B4 = 	  B4.*( sin(nyveclist(2)*(pi + thetaedge1_re2 + thetaR))./(B2-cos(nyveclist(2)*(pi + thetaedge1_re2 + thetaR))) + ...
				    sin(nyveclist(2)*(pi + thetaedge1_re2 - thetaR))./(B2-cos(nyveclist(2)*(pi + thetaedge1_re2 - thetaR))) + ...
				    sin(nyveclist(2)*(pi - thetaedge1_re2 + thetaR))./(B2-cos(nyveclist(2)*(pi - thetaedge1_re2 + thetaR))) + ...
					sin(nyveclist(2)*(pi - thetaedge1_re2 - thetaR))./(B2-cos(nyveclist(2)*(pi - thetaedge1_re2 - thetaR))))./E2Rdist;
                
    else,

    if thetaS == 0,
        B4 = -prod(multfac)*nyveclist(1)*nyveclist(2)*nyveclist(3)/4/pi^3*...
		              ( sin(nyveclist(1)*(pi + thetapos2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetapos2_re1   ))))./S2Edist./E2Edist1;
% 		              ( sin(nyveclist(1)*(pi + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetaedge2_re1   ))))./S2Edist./E2Edist1;
    else,        
        B4 = -prod(multfac)*nyveclist(1)*nyveclist(2)*nyveclist(3)/8/pi^3*...
		              ( sin(nyveclist(1)*(pi + thetaS + thetapos2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetaS + thetapos2_re1   ))) + ...
					    sin(nyveclist(1)*(pi - thetaS + thetapos2_re1  ))./(B1- cos(nyveclist(1)*(pi - thetaS + thetapos2_re1   ))))./S2Edist./E2Edist1;
% 		              ( sin(nyveclist(1)*(pi + thetaS + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi + thetaS + thetaedge2_re1   ))) + ...
% 					    sin(nyveclist(1)*(pi - thetaS + thetaedge2_re1  ))./(B1- cos(nyveclist(1)*(pi - thetaS + thetaedge2_re1   ))))./S2Edist./E2Edist1;
                    
    end    

    
    B4 = 	      B4.*( sin(nyveclist(2)*(pi + thetavec1_re2 + thetapos3_re2))./(B2-cos(nyveclist(2)*(pi + thetavec1_re2 + thetapos3_re2))) + ...
					    sin(nyveclist(2)*(pi + thetavec1_re2 - thetapos3_re2))./(B2-cos(nyveclist(2)*(pi + thetavec1_re2 - thetapos3_re2))))./E2Edist2;
    B4 = 	      B4.*( sin(nyveclist(3)*(pi + thetapos2_re3 + thetaR))./(B3-cos(nyveclist(3)*(pi + thetapos2_re3 + thetaR))) + ...
					    sin(nyveclist(3)*(pi + thetapos2_re3 - thetaR))./(B3-cos(nyveclist(3)*(pi + thetapos2_re3 - thetaR))))./E2Rdist;
%     B4 = 	      B4.*( sin(nyveclist(2)*(pi + thetaedge1_re2 + thetaedge3_re2))./(B2-cos(nyveclist(2)*(pi + thetaedge1_re2 + thetaedge3_re2))) + ...
% 					    sin(nyveclist(2)*(pi + thetaedge1_re2 - thetaedge3_re2))./(B2-cos(nyveclist(2)*(pi + thetaedge1_re2 - thetaedge3_re2))))./E2Edist2;
%     B4 = 	      B4.*( sin(nyveclist(3)*(pi + thetaedge2_re3 + thetaR))./(B3-cos(nyveclist(3)*(pi + thetaedge2_re3 + thetaR))) + ...
% 					    sin(nyveclist(3)*(pi + thetaedge2_re3 - thetaR))./(B3-cos(nyveclist(3)*(pi + thetaedge2_re3 - thetaR))))./E2Rdist;

end

integrandvec = B4.*exp(-j*k*(S2Edist + E2Edist1 + E2Edist2 + E2Rdist - R_irstart));

% % % t01 = etime(clock,t00);
% % % disp(['   Times = ',num2str(t01),' s and .',num2str(t02),' Length of vec is ',int2str(length(zvec1))])
% % % % disp(['   Times = ',num2str(t01),' s. Length of vec is ',int2str(length(zvec1))])
% % % pause
% % % 

% disp([num2str(B4)])
