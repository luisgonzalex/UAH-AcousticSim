function [hitvec,edgehit,cornerhit] = EDB2poinpla(xpoints,planelist,minvals,maxvals,planecorners,corners,ncornersperplanevec,planenvecs)
% EDB2poinpla - Detects if one or more points are inside a number of finite planes. 
% If one point is given as input, it will be checked against all
% planes in a list of planes. If N points are given as inputs, a list of
% planes should have the N planes, and each point will be checked against
% its corresponding plane.
%
% Input parameters:
%   xpoints         Matrix, [N,3], of coordinates for the N points to check
%   planelist       List, [nplanes,1], of planes to check the N points
%                   against. If N~= 1, then nplanes must be equal to N.
%                   NB! 
%   minvals, maxvals, planecorners, corners, ncornersperplanevec, planenvecs
%                   Data that should have been taken from the corresponding
%                   variables in the eddatafile,see EDB2edgeo for more information.
%                   NB!! All of these matrices except corners
%                   have been rearranged so that they have N rows, and each
%                   row contain the data for one specific plane, the one
%                   that the inside-check should be done for.
%
% Output parameters:
%   hitvec			List, [N,1], with 1 or 0 indicating whether a point is
%                   inside or outside the plane given by planelist.
%   edgehit         List, [N,1], with 1 or 0 indicating if a hit was right at the
%                   edge of a plane. These hits were not marked in hitvec.
%   cornerhit       List, [N,1], with 1 or 0 indicating if a hit was right at the
%                   corner of a plane. These hits were not marked in hitvec.
%
% Uses no special functions
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
% Peter Svensson (svensson@iet.ntnu.no) 20130813
% 
% [hitvec,edgehit,cornerhit] = EDB2poinpla(xpoints,planelist,minvals,maxvals,planecorners,corners,ncornersperplanevec,planenvecs);

% 20100204 Functioning version
% 20111201 Fixed a bug: when the ray (from the hitpoint in the positive x-direction) passed exactly
%		 through a vertex, then there was an error. For those rare cases an extra ray is shot in a random direction.
% 20130813 The random-ray version didn't seem to work properly so the
%          function was rewritten quite much. At the same time, the edgehit and
%          cornerhit was implemented correctly.

global SHOWTEXT

geomacc = 1e-12;

npoints = size(xpoints,1);
nplanestotest = length(planelist);
if npoints == 1
	xpoints = xpoints(ones(nplanestotest,1),:);
end

planecorners = [planecorners planecorners(:,1)];

%------------------------------------------------------------
% First test: are the points inside the cubic boxes?
%
% NB! The values in possibleones tell which entries in the lists xpoint
%     and planelist that passed the first test.

if nplanestotest <= 65535
    possibleones = uint16(1:nplanestotest);
else
    possibleones = uint32(1:nplanestotest);    
end
possibleones = possibleones(:);

iv = uint32(find(xpoints(:,1) > maxvals(planelist,1)));
possibleones(iv) = [];
iv = uint32(find(xpoints(possibleones,1) < minvals(planelist(possibleones),1)));
possibleones(iv) = [];
iv = uint32(find(xpoints(possibleones,2) > maxvals(planelist(possibleones),2)));
possibleones(iv) = [];
iv = uint32(find(xpoints(possibleones,2) < minvals(planelist(possibleones),2)));
possibleones(iv) = [];
iv = uint32(find(xpoints(possibleones,3) > maxvals(planelist(possibleones),3)));
clear maxvals
possibleones(iv) = [];
iv = uint32(find(xpoints(possibleones,3) < minvals(planelist(possibleones),3)));
possibleones(iv) = [];
clear minvals iv

hitvec = uint8(zeros(nplanestotest,1));
edgehit   = uint8(zeros(nplanestotest,1));
cornerhit = uint8(zeros(nplanestotest,1));

nposs = length(possibleones);

if SHOWTEXT >= 4
    disp(['         Of the ',int2str(npoints),' points,'])
    disp(['         ',int2str(length(possibleones)),' survived the cube test:'])      
end



%------------------------------------------------------------
% Second test: project onto two dimensions.
% Start by finding which dimension of the planes that has the strongest
% normal vector component. That dimension should be tossed.
%
% Easiest way to handle: make three subsets:
%   Combos that should be projected onto xy
%   Combos that should be projected onto xz
%   Combos that should be projected onto yz

if nposs>0,    

    A = abs((planenvecs(planelist(possibleones),:)));
    maxA = max(A.').';
    markwhereismax = (A == maxA(:,[1 1 1]));
    iv = find(markwhereismax(:,1).*markwhereismax(:,2)~=0);
    if ~isempty(iv)
        markwhereismax(iv,2) = 0;
    end
    iv = find(markwhereismax(:,1).*markwhereismax(:,3)~=0);
    if ~isempty(iv)
        markwhereismax(iv,1) = 0;
    end
    iv = find(markwhereismax(:,2).*markwhereismax(:,3)~=0);
    if ~isempty(iv)
        markwhereismax(iv,2) = 0;
    end
    colno = [1 2 3];
    finalcolno = markwhereismax.*colno(ones(nposs,1),:);
    finalcolno = sum(finalcolno.').';

    yzsubset = find(finalcolno==1);
    xzsubset = find(finalcolno==2);
    xysubset = find(finalcolno==3);
    
    %---------------------------------------------
    % First the xysubset
    %
    % Create a ray that starts in xpoint and extends parallel to the x-axis
    % in the positive x-direction, that is:
    %       y = xpoints(:,2):
    %       xstart = xpoints(:,1);
    
    if ~isempty(xysubset)
        if SHOWTEXT >= 4
            disp(['            ',int2str(length(xysubset)),' xy-projected points:'])
        end

        numberofedgestocheck = ncornersperplanevec(planelist(possibleones(xysubset)));    
        edgenumbers = unique(numberofedgestocheck);

        yray = xpoints(possibleones(xysubset),2);
        xstart = xpoints(possibleones(xysubset),1);
        edgecrossings = zeros(size(xysubset));
        loweredgehits = zeros(size(xysubset));
        upperedgehits = zeros(size(xysubset));        
        
        for ii = 1:max(edgenumbers)
            y1shift = corners(planecorners(planelist(possibleones(xysubset)),ii),2) - yray;
            y2shift = corners(planecorners(planelist(possibleones(xysubset)),ii+1),2) - yray;
            x1shift = corners(planecorners(planelist(possibleones(xysubset)),ii),1) - xstart;
            x2shift = corners(planecorners(planelist(possibleones(xysubset)),ii+1),1) - xstart;

            signy1 = sign(y1shift);
            signy2 = sign(y2shift);
            signx1 = sign(x1shift);
            signx2 = sign(x2shift);
            
            catmp_or_pm = find(signy1.*signy2 == -1);
            catzz = find(signy1 ==  0 & signy2 == 0);
            catzp = find(signy1 ==  0 & signy2 == +1);
            catzm = find(signy1 ==  0 & signy2 == -1);
            catpz = find(signy1 == +1 & signy2 ==  0);
            catmz = find(signy1 == -1 & signy2 ==  0);
            
            if ~isempty(catmp_or_pm)
                tedgecrossing = y1shift(catmp_or_pm)./(y1shift(catmp_or_pm)-y2shift(catmp_or_pm));
                xsign = sign( x1shift(catmp_or_pm) + tedgecrossing.*(x2shift(catmp_or_pm)-x1shift(catmp_or_pm)) );
                iv = find(xsign == 0);
                if ~isempty(iv)
                   edgehit( possibleones(xysubset(catmp_or_pm(iv))) ) = 1; 
                end
                iv = find(xsign == +1);
                if ~isempty(iv)
                   edgecrossings( catmp_or_pm(iv) ) = edgecrossings( catmp_or_pm(iv) ) + (ii <= numberofedgestocheck(catmp_or_pm(iv)));
                end                
            end
            
            if ~isempty(catpz)
               iv = find(signx2(catpz) == 0);
               if ~isempty(iv)
                  cornerhit( possibleones(xysubset(catpz(iv))) ) = 1;
               end
               iv = find(signx2(catpz) == 1);
               if ~isempty(iv)
                   upperedgehits( catpz(iv) ) = upperedgehits( catpz(iv) ) + 1*(ii <= numberofedgestocheck(catpz(iv)));
               end                
            end
            if ~isempty(catmz)
               iv = find(signx2(catmz) == 0);
               if ~isempty(iv)
                  cornerhit( possibleones(xysubset(catmz(iv))) ) = 1; 
               end
               iv = find(signx2(catmz) == 1);
               if ~isempty(iv)
                   loweredgehits( catmz(iv) ) = loweredgehits( catmz(iv) ) + 1*(ii <= numberofedgestocheck(catmz(iv)));
               end                
            end
            
            if ~isempty(catzp)
               iv = find(signx1(catzp) == 0);
               if ~isempty(iv)
                  cornerhit( possibleones(xysubset(catzp(iv))) ) = 1; 
               end
               iv = find(signx1(catzp) == 1);
               if ~isempty(iv)
                   upperedgehits( catzp(iv) ) = upperedgehits( catzp(iv) ) + 1*(ii <= numberofedgestocheck(catzp(iv)));
               end                
            end
            if ~isempty(catzm)
               iv = find(signx1(catzm) == 0);
               if ~isempty(iv)
                  cornerhit( possibleones(xysubset(catzm(iv))) ) = 1; 
               end
               iv = find(signx1(catzm) == 1);
               if ~isempty(iv)
                   loweredgehits( catzm(iv) ) = loweredgehits( catzm(iv) ) + 1*(ii <= numberofedgestocheck(catzm(iv))); 
               end                
            end
            
            if ~isempty(catzz)
               iv = find(signx1(catzz).*signx2(catzz) == -1);
               if ~isempty(iv)
                  edgehit( possibleones(xysubset(catzz(iv))) ) = 1;                
               end               
               iv = find(signx1(catzz).*signx2(catzz) == 0);
               if ~isempty(iv)
                  cornerhit( possibleones(xysubset(catzz(iv))) ) = 1;                
               end                                              
            end            
        end
         
        hitvec(possibleones(xysubset)) = rem(edgecrossings,2);        
        hitvec(possibleones(xysubset)) = hitvec(possibleones(xysubset)) + uint8(((loweredgehits==1).*(upperedgehits==1)));
        hitvec = hitvec.*rem(hitvec,2);
        
        if SHOWTEXT >= 4
            disp(['               ',int2str(sum((edgecrossings==1))),' survived the xyplane projections test:'])  
        end
        
    end
        
    %---------------------------------------------
    % Then the xzsubset
    %
    % Create a ray that starts in xpoint and extends parallel to the x-axis
    % in the positive x-direction, that is:
    %       z = xpoints(:,3):
    %       xstart = xpoints(:,1);

    if ~isempty(xzsubset)
        if SHOWTEXT >= 4
            disp(['            ',int2str(length(xzsubset)),' xz-projected points:'])
        end

        numberofedgestocheck = ncornersperplanevec(planelist(possibleones(xzsubset)));    
        edgenumbers = unique(numberofedgestocheck);

        zray = xpoints(possibleones(xzsubset),3);
        xstart = xpoints(possibleones(xzsubset),1);
        edgecrossings = zeros(size(xzsubset));
        loweredgehits = zeros(size(xzsubset));
        upperedgehits = zeros(size(xzsubset));        

        for ii = 1:max(edgenumbers)
            z1shift = corners(planecorners(planelist(possibleones(xzsubset)),ii),3) - zray;
            z2shift = corners(planecorners(planelist(possibleones(xzsubset)),ii+1),3) - zray;
            x1shift = corners(planecorners(planelist(possibleones(xzsubset)),ii),1) - xstart;
            x2shift = corners(planecorners(planelist(possibleones(xzsubset)),ii+1),1) - xstart;

            signz1 = sign(z1shift);
            signz2 = sign(z2shift);
            signx1 = sign(x1shift);
            signx2 = sign(x2shift);
            
            catmp_or_pm = find(signz1.*signz2 == -1);
            catzz       = find(signz1 ==  0 & signz2 == 0);
            catzp       = find(signz1 ==  0 & signz2 == +1);
            catzm       = find(signz1 ==  0 & signz2 == -1);
            catpz       = find(signz1 == +1 & signz2 ==  0);
            catmz       = find(signz1 == -1 & signz2 ==  0);
        
            if ~isempty(catmp_or_pm)
                tedgecrossing = z1shift(catmp_or_pm)./(z1shift(catmp_or_pm)-z2shift(catmp_or_pm));
                xsign = sign( x1shift(catmp_or_pm) + tedgecrossing.*(x2shift(catmp_or_pm)-x1shift(catmp_or_pm)) );
                iv = find(xsign == 0);
                if ~isempty(iv)
                   edgehit( possibleones(xzsubset(catmp_or_pm(iv))) ) = 1; 
                end
                iv = find(xsign == +1);
                if ~isempty(iv)
                   edgecrossings( catmp_or_pm(iv) ) = edgecrossings( catmp_or_pm(iv) ) + (ii <= numberofedgestocheck(catmp_or_pm(iv)));
                end                
            end
        
            if ~isempty(catpz)
               iv = find(signx2(catpz) == 0);
               if ~isempty(iv)
                  cornerhit( possibleones(xzsubset(catpz(iv))) ) = 1;
               end
               iv = find(signx2(catpz) == 1);
               if ~isempty(iv)
                   upperedgehits( catpz(iv) ) = upperedgehits( catpz(iv) ) + 1*(ii <= numberofedgestocheck(catpz(iv)));
               end                
            end
            if ~isempty(catmz)
               iv = find(signx2(catmz) == 0);
               if ~isempty(iv)
                  cornerhit( possibleones(xzsubset(catmz(iv))) ) = 1; 
               end
               iv = find(signx2(catmz) == 1);
               if ~isempty(iv)
                   loweredgehits( catmz(iv) ) = loweredgehits( catmz(iv) ) + 1*(ii <= numberofedgestocheck(catmz(iv)));
               end                
            end
            
            if ~isempty(catzp)
               iv = find(signx1(catzp) == 0);
               if ~isempty(iv)
                  cornerhit( possibleones(xzsubset(catzp(iv))) ) = 1; 
               end
               iv = find(signx1(catzp) == 1);
               if ~isempty(iv)
                   upperedgehits( catzp(iv) ) = upperedgehits( catzp(iv) ) + 1*(ii <= numberofedgestocheck(catzp(iv)));
               end                
            end
            if ~isempty(catzm)
               iv = find(signx1(catzm) == 0);
               if ~isempty(iv)
                  cornerhit( possibleones(xzsubset(catzm(iv))) ) = 1; 
               end
               iv = find(signx1(catzm) == 1);
               if ~isempty(iv)
                   loweredgehits( catzm(iv) ) = loweredgehits( catzm(iv) ) + 1*(ii <= numberofedgestocheck(catzm(iv))); 
               end                
            end

            if ~isempty(catzz)
               iv = find(signx1(catzz).*signx2(catzz) == -1);
               if ~isempty(iv)
                  edgehit( possibleones(xzsubset(catzz(iv))) ) = 1;                
               end               
               iv = find(signx1(catzz).*signx2(catzz) == 0);
               if ~isempty(iv)
                  cornerhit( possibleones(xzsubset(catzz(iv))) ) = 1;                
               end                                              
            end            
        
        end
        
        hitvec(possibleones(xzsubset)) = rem(edgecrossings,2);
        hitvec(possibleones(xzsubset)) = hitvec(possibleones(xzsubset)) + uint8(((loweredgehits==1).*(upperedgehits==1)));
        hitvec = hitvec.*rem(hitvec,2);

        if SHOWTEXT >= 4
            disp(['               ',int2str(sum((edgecrossings==1))),' survived the xzplane projections test:'])  
        end
    end
        
        
    %---------------------------------------------
    % Third the yzsubset
    %
    % Create a ray that starts in xpoint and extends parallel to the y-axis
    % in the positive y-direction, that is:
    %       z = xpoints(:,3):
    %       ystart = xpoints(:,2);

    if ~isempty(yzsubset)
        if SHOWTEXT >= 4
            disp(['            ',int2str(length(yzsubset)),' yz-projected points:'])
        end

        numberofedgestocheck = ncornersperplanevec(planelist(possibleones(yzsubset)));    
        edgenumbers = unique(numberofedgestocheck);

        zray = xpoints(possibleones(yzsubset),3);
        ystart = xpoints(possibleones(yzsubset),2);
        edgecrossings = zeros(size(yzsubset));
        loweredgehits = zeros(size(yzsubset));
        upperedgehits = zeros(size(yzsubset));        

        % Use a parametric representation for each edge:
        % y_edge = y_1 + t*(y_2 - y_1)
        % z_edge = z_1 + t*(z_2 - z_1)
        % Find t by setting z_ray = z_edge

        for ii = 1:max(edgenumbers)
            z1shift = corners(planecorners(planelist(possibleones(yzsubset)),ii),3) - zray;
            z2shift = corners(planecorners(planelist(possibleones(yzsubset)),ii+1),3) - zray;
            y1shift = corners(planecorners(planelist(possibleones(yzsubset)),ii),2) - ystart;
            y2shift = corners(planecorners(planelist(possibleones(yzsubset)),ii+1),2) - ystart;

            signz1 = sign(z1shift);
            signz2 = sign(z2shift);
            signy1 = sign(y1shift);
            signy2 = sign(y2shift);
            
            catmp_or_pm = find(signz1.*signz2 == -1);
            catzz       = find(signz1 ==  0 & signz2 == 0);
            catzp       = find(signz1 ==  0 & signz2 == +1);
            catzm       = find(signz1 ==  0 & signz2 == -1);
            catpz       = find(signz1 == +1 & signz2 ==  0);
            catmz       = find(signz1 == -1 & signz2 ==  0);
            
            if ~isempty(catmp_or_pm)
                tedgecrossing = z1shift(catmp_or_pm)./(z1shift(catmp_or_pm)-z2shift(catmp_or_pm));
                ysign = sign( y1shift(catmp_or_pm) + tedgecrossing.*(y2shift(catmp_or_pm)-y1shift(catmp_or_pm)) );
                iv = find(ysign == 0);
                if ~isempty(iv)
                   edgehit( possibleones(yzsubset(catmp_or_pm(iv))) ) = 1; 
                end
                iv = find(ysign == +1);
                if ~isempty(iv)
                   edgecrossings( catmp_or_pm(iv) ) = edgecrossings( catmp_or_pm(iv) ) + (ii <= numberofedgestocheck(catmp_or_pm(iv)));
                end                
            end
        
            if ~isempty(catpz)
               iv = find(signy2(catpz) == 0);
               if ~isempty(iv)
                  cornerhit( possibleones(yzsubset(catpz(iv))) ) = 1;
               end
               iv = find(signy2(catpz) == 1);
               if ~isempty(iv)
                   upperedgehits( catpz(iv) ) = upperedgehits( catpz(iv) ) + 1*(ii <= numberofedgestocheck(catpz(iv)));
               end                
            end
            if ~isempty(catmz)
               iv = find(signy2(catmz) == 0);
               if ~isempty(iv)
                  cornerhit( possibleones(yzsubset(catmz(iv))) ) = 1; 
               end
               iv = find(signy2(catmz) == 1);
               if ~isempty(iv)
                   loweredgehits( catmz(iv) ) = loweredgehits( catmz(iv) ) + 1*(ii <= numberofedgestocheck(catmz(iv)));
               end                
            end
 
            if ~isempty(catzp)
               iv = find(signy1(catzp) == 0);
               if ~isempty(iv)
                  cornerhit( possibleones(yzsubset(catzp(iv))) ) = 1; 
               end
               iv = find(signy1(catzp) == 1);
               if ~isempty(iv)
                   upperedgehits( catzp(iv) ) = upperedgehits( catzp(iv) ) + 1*(ii <= numberofedgestocheck(catzp(iv)));
               end                
            end
            if ~isempty(catzm)
               iv = find(signy1(catzm) == 0);
               if ~isempty(iv)
                  cornerhit( possibleones(yzsubset(catzm(iv))) ) = 1; 
               end
               iv = find(signy1(catzm) == 1);
               if ~isempty(iv)
                   loweredgehits( catzm(iv) ) = loweredgehits( catzm(iv) ) + 1*(ii <= numberofedgestocheck(catzm(iv))); 
               end                
            end

             if ~isempty(catzz)
               iv = find(signy1(catzz).*signy2(catzz) == -1);
               if ~isempty(iv)
                  edgehit( possibleones(yzsubset(catzz(iv))) ) = 1;                
               end               
               iv = find(signy1(catzz).*signy2(catzz) == 0);
               if ~isempty(iv)
                  cornerhit( possibleones(yzsubset(catzz(iv))) ) = 1;                
               end                                              
            end                       
            
        end
        
        hitvec(possibleones(yzsubset)) = rem(edgecrossings,2);
        hitvec(possibleones(yzsubset)) = hitvec(possibleones(yzsubset)) + uint8(((loweredgehits==1).*(upperedgehits==1)));
        hitvec = hitvec.*rem(hitvec,2);

        if SHOWTEXT >= 4
            disp(['               ',int2str(sum((edgecrossings==1))),' survived the yzplane projections test:'])      
        end
    end
        
end

hitvec = hitvec.*(1-edgehit);
hitvec = hitvec.*(1-cornerhit);
