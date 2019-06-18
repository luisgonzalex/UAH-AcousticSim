function edirfile = EDB2makeirs(edpathsfile,specorder,Rstart,...
                        EDcalcmethod,edgestartcoords,edgeendcoords,edgenvecs,...
                        edgelengthvec,planeeqs,approxplanemidpoints,reflfactors,closwedangvec,planesatedge,elemsize,reftoshortlistE,re1sho,re2sho,...
                        thetae1sho,thetae2sho,ze1sho,ze2sho,edgeseespartialedge,edgeplaneperptoplane1,desiredirfile,guiderowstouse,includedirectsound,...
                        saveindividualdiffirs)
% EDB2makeirs - Constructs impulse responses from a list of paths in a file.
%
% Input parameters:
% 	edpathsfile         The name of the file containing the plane and edge data
% 	specorder           The maximum order of reflections that is calculated.
% 	Rstart              The reference distance that the time zero of the impulse
% 	                    response refers to.
% 	EDcalcmethod        'n' or 'v' or 'k', meaning:
%                       'n': the new method by Svensson et al
%                       'v': Vanderkooys method
%                       'k': the Kirchhoff diffraction approximation
% 	edgestartcoords     A matrix [nedges,3] of the startpoint coordinates
% 	of
% 	                    the edges.
% 	edgeendcoords       A matrix [nedges,3] of the endpoint coordinates of
% 	                    the edges.
% 	edgenvecs           A matrix [nedges,3] of the normal vectors of
% 	                    the reference plane of each edge.
% 	edgelengthvec       A list [nedge,1] of the lenghts of each edge.
%   planeeqs            A matrix, [nplanes,4], of all the plane equations.
%   approxplanemidpoints   A matrix, [nplanes,3], of all plane midpoint
%                       coordinates. It may be an approximate location of
%                       the midpoint.
% 	reflfactors         A list [nplanes,1] of the reflection factors of
% 	                    each plane.
% 	closwedangvec       A list [nedge,1] of the close-wedge angle of each
% 	edge.
% 	planesatedge        A matrix [nedge,2] of the two planes of each edge
% 	                    with the reference plane first.
% 	elemsize            A list [1,difforder] with the relative edge element
% 	                    sizes that will be used for each order. The first
% 	                    value, for diffraction order 1, is redundant and
% 	                    not used. For higher orders, a value of, e.g., 0.5
% 	                    means that edges will be subdivided into elements
% 	                    such that the wavelength at half the sampling
% 	                    frequency is 0.5 times the element size. A lower
% 	                    value gives faster but less accurate results.
% 	reftoshortlistE     A matrix, [nedges,nedges] with a pointer to the
% 	                    short lists that contain edge-to-edge data.
% 	re1sho              A short list, [nshort,1] of the cylindrical radii
%                       of each edge startpoint relative to all other
%                       edges.
% 	re2sho              A short list, [nshort,1] of the cylindrical radii
%                       of each edge endpoint relative to all other
%                       edges.
% 	thetae1sho          A short list, [nshort,1] of the theta angle
%                       of each edge startpoint relative to all other
%                       edges.
% 	thetae2sho          A short list, [nshort,1] of the theta angle
%                       of each edge endpoint relative to all other
%                       edges.
% 	ze1sho              A short list, [nshort,1] of the z-value
%                       of each edge startpoint relative to all other
%                       edges.
% 	ze2sho              A short list, [nshort,1] of the z-value
%                       of each edge endpoint relative to all other
%                       edges.
% 	edgeseespartialedge A matrix, [nedges,nedges], with edge-to-edge
% 	                    visibility values. The values are 0 (invisible)
%                       up to (2^(nedgesubs)-1)^2 (full visibility).
%                       A pos. value indicates that the edge-to-edge path
%                       runs along a plane; a neg. avlue indicates that it
%                       does not run along a plane.
%   edgeplaneperptoplane1   A matrix, [nplanes,nedges], with 0 or 1. A
%                       value 1 indicates that one of the edge's two
%                       defining planes is perpendicular to another plane.
% 	desiredirfile       The file name that will be given to the
% 	                    output file. 
%   guiderowstouse      (optional) A list of values, indicating which rows in the
%                       mainlistguide and mainlistguidepattern that should
%                       be used. This way it is possible to select only
%                       diffraction, for instance. If this list is not
%                       specified, all components will be used.
%   includedirectsound  (optional) 0 or 1, indicating whether the direct
%                       sound should be included or not. Default: 1 
%   saveindividualdiffirs (optional) Vector of length 2 with values 0 or 1
%						Pos 1 indicating whether
%                       	0 - all orders of diffraction IR:s will be summed
%                       	in a single vector 'irdiff', or
%                       	1 - each order of diffraction IR:s will be placed
%                       	in individual columns in a matrix 'irdiff'.
%						Pos 2 indicating whether
%							0 - only the total diffraction ir will be saved in the file.
%							1 - all individual diffraction irs will be saved in a large 
%							matrix alldiffirs.
%                       Default: [0 0].
%   FSAMP, CAIR, SHOWTEXT   Global parameters.
%
% Output parameters:
%   edirfile            The filename used for the output file that contains
%                       the impulse responses.
%
% Data in the output file:
%   irdirect            The IR containing the direct sound, size [nirlength,1].
%   irgeom              The IR containing the specular reflections, size [nirlength,1].
%   irdiff              The IR containing the diffracted components, size [nirlength,1].
%   irtot               The IR containing the sum of all components, size [nirlength,1].
%	alldiffirs			(optional) Matrix with all individual first-order diffraction IRs
%						on individual lines.
%	alldiffirsextradata	(optional) Vector with combination number (form reflpaths) that
%						matches alldiffirs.
%   FSAMP Rstart CAIR   Same as the input parameters.
%
% Uses functions  EDB2irfromslotlist EDB2calcdist EDB2coordtrans1 EDB2coordtrans2
% 				  EDB2wedge1st_int EDB2wedge2nd
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
% Peter Svensson (svensson@iet.ntnu.no) 20061208
%
% edirfile = EDB2makeirs(edpathsfile,specorder,...
%            Rstart,EDcalcmethod,edgestartcoords,edgeendcoords,edgenvecs,...
%            edgelengthvec,planeeqs,approxplanemidpoints,reflfactors,closwedangvec,planesatedge,elemsize,reftoshortlistE,re1sho,re2sho,...
%            thetae1sho,thetae2sho,ze1sho,ze2sho,edgeseespartialedge,edgeplaneperptoplane1,desiredirfile,guiderowstouse,includedirectsound,...
%            saveindividualdiffirs);

global SHOWTEXT FSAMP CAIR BIGEDGESTEPMATRIX

global IRDIFFVEC

eval(['load ',edpathsfile])

%--------------------------------------------------------------------------
% We look for the optional vector multfactors in the edpathsfile.
% If it was there, then those values will be used to switch off, or boost
% any reflection path.
% Those values will only be used for diffraction paths.

if exist('multfactors') == 0
    multfactors = ones(size(reflpaths,1),1);    
else
    if size(multfactors,1) ~= size(reflpaths,1)
        error(['The edpaths file contained a vector multfactors which does not have the same length as reflpaths.'])    
    end
end

%--------------------------------------------------------------------------

edirfile = desiredirfile;

Sdirection = [1 0 0];
maxnedgesorplanes = max(max(reflpaths));
if ~isempty(maxnedgesorplanes)
    multfac = 10^ceil(log10(double(maxnedgesorplanes)+1));
else
    multfac = 0;
end

if nargin < 27
    saveindividualdiffirs = [0 0];
end
if nargin < 26
    includedirectsound = 1;
end
if nargin < 25
    guiderowstouse = [];
end
if isempty(guiderowstouse)
    usesubset = 0;    
    if SHOWTEXT > 2
        disp('Using all components')
    end
else
    usesubset = 1;        
    if SHOWTEXT > 2
        disp('Using only some components')
    end
end

nyvec = pi./(2*pi - closwedangvec);    
onesvec1 = ones(1,3);

souspecboost = 1;
if ~isempty(Sinsideplanenumber)
    if reflfactors(Sinsideplanenumber(1)) ~= 0
        souspecboost = 2;
    end
end

recspecboost = 1;
if ~isempty(Rinsideplanenumber)
    if reflfactors(Rinsideplanenumber(1)) ~= 0
        recspecboost = 2;
    end
end

if isempty(re1sho)
    userwantsdiff2 = 0;
else
    userwantsdiff2 = 1;    
end

%-------------------------------------------------------

if exist('mainlistguide') ~= 1
    mainlistguide = [];    
else
    mainlistguide = double(mainlistguide);    
end

[ncomponents,ncols] = size(reflpaths);
[nrowsinlist,slask] = size(mainlistguide);
if ncols > 1
    difforderinlist = sum(mainlistguidepattern.'=='d').';
else
    difforderinlist = (mainlistguidepattern=='d');    
end

lastNdiffrow = zeros(1,specorder);
for ii = 1:specorder
    iv = find(difforderinlist <= ii );
    if ~isempty(iv)
        lastNdiffrow(ii) = iv(end);
    end
end

nrefltogetherwdiff = sum(mainlistguidepattern.'=='d').'.*(  sum(mainlistguidepattern.'=='s').');

if SHOWTEXT >= 3
	disp(['   Constructing IR. ',int2str(ncomponents),' components.'])
end

irdirect = [0 0].';
irgeom = [0 0].';
irdiff = [0 0].';
irtot = [0 0].';

Varlist = [' irdirect irgeom irdiff irtot FSAMP Rstart CAIR EDcalcmethod'];

if ncomponents == 0
    eval(['save ',edirfile,Varlist])
    return
end

%##############################################################################
%##############################################################################
%##############################################################################
%##############################################################################
%
%  Diffraction once, possibly with pre- and post-specular reflections
%
%##############################################################################

directsoundonboundary = 0;
specreflonboundary = zeros(size(planeeqs,1),1);

if firstdiffrow ~= 0
    
    % First we remove the 'd' for all the rows in the mainlistguidepattern that
    % we do not want to use.
    
    if usesubset == 1
        singdiffrows = find(difforderinlist==1);
        rowstoclear = setdiff(singdiffrows,guiderowstouse);
        if ~isempty(rowstoclear)
            mainlistguidepattern(rowstoclear,:) = mainlistguidepattern(rowstoclear,:)*0;    
        end
    end

    % Here we can use the single-diffraction calculation for all
    % combinations. The source should either be the original source or
    % an image source. The receiver should either be the original
    % receiver or an image receiver.
    
    % First we create lists that specify whether there are specular
    % reflections before the edge diffraction, and after the edge
    % diffraction. Also, a list which gives the column where the edge
    % diffraction is present can be used to efficiently extract the
    % edge number. 
    % NB! The lists prespeclist, postspeclist and singdiffcol all have
    % the (short) length of the mainlistguide. 
    % The data in these lists will then be used to create the long
    % lists prespecpresent and postspecpresent which simply are 0 or 1.
    % A matrix called edgemask will have the same number of columns as
    % the reflpaths matrix. For each row, there will be a 1 in the
    % column where the edge diffraction is so we can efficiently find
    % the edge number.

    multmatrix = [1:specorder];
    if size(mainlistguidepattern,2) > specorder
        multmatrix = [1:size(mainlistguidepattern,2)];
    end
    multmatrix = multmatrix( ones(nrowsinlist,1),:  );
    
    if ncols > 1
        singdiffcol = sum( (multmatrix.*(mainlistguidepattern=='d')).' ).';
        singdiffcol = singdiffcol.*(sum(mainlistguidepattern.'=='d').'<=1);
        nrefl = sum( (mainlistguidepattern=='s' | mainlistguidepattern=='d').' ).';
    else
        singdiffcol = mainlistguidepattern=='d';
        nrefl = singdiffcol | mainlistguidepattern=='s';
    end
    
    prespeclist = singdiffcol-1;
    prespeclist = prespeclist.*(prespeclist>0);
    postspeclist = (nrefl-singdiffcol).*(singdiffcol~=0);
    
    prespecpresent = zeros(ncomponents,1);
    postspecpresent = zeros(ncomponents,1);
    edgemask = zeros(ncomponents,specorder);

    diffrowsthatareOK = [firstdiffrow:lastNdiffrow(1)];
    if usesubset == 1
        diffrowsthatareOK = intersect(diffrowsthatareOK,guiderowstouse);
    end
    for ii = diffrowsthatareOK,    %firstdiffrow:lastNdiffrow(1)
        iv = [(mainlistguide(ii,2)):(mainlistguide(ii,3))];
        onesvec2 = ones(length(iv),1);
        edgemask(iv,singdiffcol(ii)) = onesvec2;            
    	prespecpresent(iv) = (prespeclist(ii)>0)*onesvec2;
        postspecpresent(iv) = (postspeclist(ii)>0)*onesvec2;
    end
    prespecpresent = prespecpresent(:,onesvec1);
    postspecpresent = postspecpresent(:,onesvec1);

    if ncols > 1
        edgenumb = sum( (edgemask.*double(reflpaths(:,1:specorder))).').';
    else
        edgenumb = edgemask.*double(reflpaths);    
    end

    nedgeirs = size(edgenumb,1);
    nnonzeroedgeirs = sum(edgenumb(:,1)>0);

    ircounter = 0;
    edgelist = unique(edgenumb);
    for ii = 1: length(edgelist)
        edge = edgelist(ii);
        if edge~= 0
           iv = find(edgenumb==edge);
           ncombs = length(iv);
           
           if ncombs > 0 & any(multfactors(iv))
               if SHOWTEXT >= 4
                    disp(['   Edge ',int2str(edge),': ',int2str(ncombs),' combinations'])    
               end               
                 
                onesvec3 = ones(ncombs,1);
	
                IS = full(specextradata(iv,1:3));
                IR = full(specextradata(iv,4:6));                
            
                edgestart = full(edgeextradata(iv,1));
                edgeend   = full(edgeextradata(iv,2));
                
                % Calculate the cyl coordinates of all IS
                % 
                % Calculate the cyl coord of all IR
	
                edgecoords = [edgestartcoords(edge,:);edgeendcoords(edge,:)];
            	[rs,thetas,zs,rr,thetar,zr] = EDB2coordtrans2(IS,IR,edgecoords,edgenvecs(edge,:));
	
				bc = reflfactors(planesatedge(edge,:)).';
				if prod(double(bc==[1 1])) ~= 1
					disp('   Non-rigid wedge')
				end
	
                for jj = 1:ncombs
                    if multfactors(iv(jj)) > 0
                        if EDcalcmethod(1) == 'n'

                           [irnew,slask,singularterm] = EDB2wedge1st_int(FSAMP,closwedangvec(edge),rs(jj),thetas(jj),zs(jj),rr(jj),thetar(jj),zr(jj),...
                             edgelengthvec(edge)*[edgestart(jj) edgeend(jj)],EDcalcmethod,Rstart,bc);                  
                         
                              if any(singularterm)
                                  if singularterm(2) | singularterm(3)
                                      directsoundonboundary = 1;
                                elseif singularterm(1)
                                    if specorder == 1
                                       specreflonboundary(planesatedge(edge,2)) = 1;
                                   else
                                       disp('WARNING! A specular refl. of order > 1 is half obscured but this is not handled yet!');    
                                   end
                                elseif singularterm(4)
                                      if specorder == 1
                                           specreflonboundary(planesatedge(edge,1)) = 1;
                                       else
                                           disp('WARNING! A specular refl. of order > 1 is half obscured but this is not handled yet!');    
                                       end
                                end
                              end
		
                        else    % ...   if EDcalcmethod(1) == 'n'
                           
                              [irnew,slask,singularterm] = EDB2wedge1stcombo(FSAMP,closwedangvec(edge),rs(jj),thetas(jj),zs(jj),rr(jj),thetar(jj),zr(jj),...
                                  edgelengthvec(edge)*[edgestart(jj) edgeend(jj)],EDcalcmethod,Rstart,bc);                      
                              
                        end      % ...   if EDcalcmethod(1) == 'n'
                    
                        % Decide whether the IR should be boosted or not
                        %
                        % The factors souspecboost and recspecboost have
                        % values other than 1 if the source or receiver is directly
                        % at a plane, for the boosting of the direct sound, and
                        % specular reflections.
                        %
                        % This boost factor should be used also for edge
                        % diffraction if:
                        %   1.  There are some specular reflections between the 
                        %       source and the edge
                        %   or
                        %   2.  There are no specular reflections between the
                        %       source and the edge, but the source/receiver is at
                        %       a plane which does not belong to the edge.
                        
                        if souspecboost ~= 1 | recspecboost ~= 1
                        
                            boostfactor = 1;
                            if prespecpresent(iv(jj)) == 1
                                boostfactor = souspecboost;    
                            else
                                if ~isempty(Sinsideplanenumber)
                                    if Sinsideplanenumber(1)~=planesatedge(edge,1) & Sinsideplanenumber(1)~=planesatedge(edge,2), 
                                        boostfactor = souspecboost;
                                    end
                                end
                            end
                            if postspecpresent(iv(jj)) == 1
                                boostfactor = boostfactor*recspecboost;    
                            else
                                if ~isempty(Rinsideplanenumber)
                                    if Rinsideplanenumber(1)~=planesatedge(edge,1) & Rinsideplanenumber(1)~=planesatedge(edge,2), 
                                        boostfactor = boostfactor*recspecboost;
                                    end
                                end
                            end
                            if boostfactor ~= 1
                                irnew = irnew*boostfactor;    
                            end
                        
                        end    % ....    if souspecboost ~= 1 | recspecboost ~= 1

                        ndiff = length(irdiff);
                        nnew = length(irnew);
                        
                        ircounter = ircounter + 1;                        
                        
                        if nnew > ndiff
                           irdiff = [irdiff;zeros(nnew-ndiff,1)];
                        end
                        irdiff(1:nnew) = irdiff(1:nnew) + irnew*double(multfactors(iv(jj)));

						if saveindividualdiffirs(2) == 1
							if exist('alldiffirs','var') == 0
								alldiffirs = sparse(zeros(nnonzeroedgeirs,nnew));
								alldiffirs(1,:) = irnew*double(multfactors(iv(jj))).';
								alldiffirsextradata = zeros(nnonzeroedgeirs,1);
								alldiffirsextradata(1) = iv(jj);
							else
								if nnew > ndiff
									alldiffirs = [alldiffirs zeros(nnonzeroedgeirs,nnew-ndiff)];
								end
								alldiffirs(ircounter,1:nnew) = irnew*double(multfactors(iv(jj))).';							
								alldiffirsextradata(ircounter) = iv(jj); 
							end						
						
						end


                        if SHOWTEXT >= 7
                            figure(1)
                            plot(irnew)
                            figure(2)
                            plot(irdiff)
                            save ~/Documents/Temp/irdiff2.mat
                            pause
                        end  

                    end   %  ...  if multfactors(iv) > 0
                    
                end  % .....                for jj = 1:ncombs

            end   %....    if ncombs > 0 & any(multfactors(iv))
            
        end    % ....         if edge~= 0
        
    end  % ....     for ii = 1: length(edgelist)
    
end   % ...if firstdiffrow ~= 0

nspecreflonboundary = sum(specreflonboundary);

%##############################################################################
%##############################################################################
%##############################################################################
%##############################################################################
%
%      The direct sound
%
%##############################################################################

if usesubset == 1 & directsoundrow == 1
    if any(guiderowstouse == 1) == 0
        directsoundrow = 0;    
    end
end

if directsoundrow == 1 & includedirectsound == 1
    dist = norm(R-S);
    slotnumberfrac = (dist-Rstart)/CAIR*FSAMP+1;
    amp = souspecboost*recspecboost./dist;
    if directsoundonboundary
        amp = amp/2;    
        if SHOWTEXT >= 4
            disp('HALVING DIRECT SOUND')
        end
    end
    slotnumber = floor(slotnumberfrac);
    if any(slotnumber<1)
        error('ERROR: Rstart is set to too low a value!')
    end
    slotnumberfrac = slotnumberfrac - slotnumber;
    irdirect = zeros( slotnumber+2,1 );

	irdirect(slotnumber)   = amp.*(1-slotnumberfrac) ;
	irdirect(slotnumber+1) = amp.*slotnumberfrac;	

	if ncomponents == 1
        nnew = length(irdirect);
        if nnew > 2
          irdiff = [irdiff;zeros(nnew-2,1)];
          irtot =  [irtot;zeros(nnew-2,1)];
          irgeom = [irgeom;zeros(nnew-2,1)];
        end
        irtot(1:nnew) = irtot(1:nnew) + irdirect;

        irtot = sparse(irtot);
        irgeom = sparse(irgeom);
        irdirect = sparse(irdirect);
        irdiff = sparse(irdiff);
        
        eval(['save ',edirfile,Varlist])
        return
	end    
else
    if directsoundonboundary == 1
        
        dist = norm(R-S);
        slotnumberfrac = (dist-Rstart)/CAIR*FSAMP+1;
        amp = souspecboost*recspecboost./dist;
        if directsoundonboundary
            amp = amp/2;    
            if SHOWTEXT >= 4
                disp('HALVING DIRECT SOUND')
            end
        end
        slotnumber = floor(slotnumberfrac);
        if any(slotnumber<1)
            error('ERROR: Rstart is set to too low a value!')
        end
        slotnumberfrac = slotnumberfrac - slotnumber;
        irdirect = zeros( slotnumber+2,1 );
	
		irdirect(slotnumber)   = amp.*(1-slotnumberfrac) ;
		irdirect(slotnumber+1) = amp.*slotnumberfrac;	
        
    end
    if SHOWTEXT >= 4
        disp(['      No direct sound'])    
    end
end

%##############################################################################
%##############################################################################
%##############################################################################
%##############################################################################
%
%  All-specular			s, ss, sss, etc
%
%##############################################################################

if allspecrows(1) ~= 0
    if usesubset == 1
        specrowsthatareOK = intersect(allspecrows,guiderowstouse);
    else
        specrowsthatareOK = [allspecrows(1):allspecrows(2)];
    end
    if ~isempty(specrowsthatareOK)
        if usesubset == 1
            temp = mainlistguide(specrowsthatareOK,2:3);
            ivspec = [double(mainlistguide(specrowsthatareOK(1),2)):double(mainlistguide(specrowsthatareOK(1),3))];
            for ii = 2:size(temp,1);
                ivspec = [ivspec  [double(mainlistguide(specrowsthatareOK(ii),2)):double(mainlistguide(specrowsthatareOK(ii),3))]];
            end
        else
            ivspec = [mainlistguide(allspecrows(1),2):mainlistguide(allspecrows(2),3)];        
        end
        
        if nspecreflonboundary > 0
            listofspecreflonboundary = find(specreflonboundary);
            specrefltoadd = setdiff(listofspecreflonboundary,ivspec)
            if ~isempty(specrefltoadd)
                error(['ERROR: We need to add some code here, to reinsert pruned spec. refl.'])    
            end
        end
        
        if SHOWTEXT >= 4
            disp(['      ',int2str(length(ivspec)),' specular reflections'])    
        end
	   
        dist = EDB2calcdist(full(specextradata(ivspec,1:3)),R);
    
		if nspecreflonboundary > 0
			specamp = ones(size(souspecboost));
			amplitudeshouldbehalf = ismember(reflpaths(ivspec,1),listofspecreflonboundary);
			specamp = specamp - amplitudeshouldbehalf*0.5;
			irnew = specamp*souspecboost*recspecboost*EDB2irfromslotlist((dist-Rstart)/CAIR*FSAMP+1,1./dist);		                
		else
			irnew = souspecboost*recspecboost*EDB2irfromslotlist((dist-Rstart)/CAIR*FSAMP+1,1./dist);		
		end
		ngeom = length(irgeom);
		nnew = length(irnew);
		if nnew > ngeom
		   irgeom = [irgeom;zeros(nnew-ngeom,1)];
		end
		irgeom(1:nnew) = irgeom(1:nnew) + irnew;
		
    else
        if nspecreflonboundary > 0
            listofspecreflonboundary = find(specreflonboundary);
            error(['ERROR: We need to add some code here, to reinsert pruned spec. refl., pos. 2'])    
        end        
    end    
else
    if nspecreflonboundary > 0
        listofspecreflonboundary = find(specreflonboundary);
        
        [xis] = EDB2findis(S,listofspecreflonboundary,planeeqs,1,[1 1 1]);
        
        disp(['WARNING: We have some new code here, to reinsert pruned spec. refl., pos. 3'])
        dist = EDB2calcdist(xis,R);
        irnew = 0.5*souspecboost*recspecboost*EDB2irfromslotlist((dist-Rstart)/CAIR*FSAMP+1,1./dist);		                
        ngeom = length(irgeom);
        nnew = length(irnew);
        if nnew > ngeom
           irgeom = [irgeom;zeros(nnew-ngeom,1)];
        end
        irgeom(1:nnew) = irgeom(1:nnew) + irnew;
        
    end        
end

%##############################################################################
%##############################################################################
%##############################################################################
%##############################################################################
%
%  Multiple diffraction, possibly with pre- and post-specular reflections
%
%##############################################################################

if userwantsdiff2 == 1
    
    JJ = setstr(32*ones(specorder,1));
    for jj=1:specorder
        jjstr = int2str(jj);
        JJ(jj,1:length(jjstr)) = jjstr;
    end
    [n1,n2] = size(JJ);

    
	for Ndifforder = 2:specorder
        if SHOWTEXT >= 3
            disp(['   Diffraction order ',int2str(Ndifforder)])    
        end
	
        if any(difforderinlist==Ndifforder) & elemsize(Ndifforder) > 0
	
            % Calculate some general parameters that are shared for all
            % N-diffraction calculations
            
            divmin = CAIR/(FSAMP*elemsize(Ndifforder));
            ndivvec = ceil(abs( edgelengthvec.' )/divmin);
            dzvec = (edgelengthvec.')./ndivvec;
	
            ncylrows = 4*(Ndifforder-1);
            
            % Here we can use the double-diffraction calculation for all
            % combinations. The source should either be the original source or
            % an image source. The receiver should either be the original
            % receiver or an image receiver.
            
            noffset = lastNdiffrow(Ndifforder-1);
            ivNdiff = [noffset+1:lastNdiffrow(Ndifforder)];
            ndiffcombs = length(ivNdiff);
            zerosvec1 = zeros(ndiffcombs,1);
            
            ndoubrows = length(ivNdiff);        
            previousrow = lastNdiffrow(Ndifforder-1);
            if previousrow > 0
                noffsetlonglist = mainlistguide(previousrow,3);
            else
                noffsetlonglist = 0;    
            end
            nremaining = ncomponents - (mainlistguide(lastNdiffrow(Ndifforder),3));
            nlonglist = ncomponents-noffsetlonglist-nremaining;
            ivlonglist = [mainlistguide(ivNdiff(1),2):mainlistguide(ivNdiff(end),3)].';
	
            diffcols = zerosvec1(:,ones(1,Ndifforder));
            for ii = ivNdiff(1):ivNdiff(end)
                diffcols(ii-ivNdiff(1)+1,:) = find(mainlistguidepattern(ii,:)=='d');
            end
            
            nreflorder = sum( (mainlistguidepattern(ivNdiff(1):ivNdiff(end),:)=='d'|mainlistguidepattern(ivNdiff(1):ivNdiff(end),:)=='s').').';
            nprespecs = diffcols(:,1)-1;
            nmidspecs = diffcols(:,2)-diffcols(:,1)-1;
            npostspecs = nreflorder-diffcols(:,Ndifforder);
	
            if ndiffcombs > 1
                ivreftolonglist = ivlonglist(:,ones(1,ndiffcombs));
                comppattern = mainlistguide(ivNdiff,2).';        
                comppattern = comppattern(ones(nlonglist,1),:);
            
                rowinpatternlist = sum( (ivreftolonglist>=comppattern).' ).';
            else
                rowinpatternlist = ones(size(ivlonglist));    
            end
            
            % Construct a long matrix with the edge numbers extracted from the
            % reflpaths matrix.
	
            longnmidspec = zeros(nlonglist,1);
            iv = find(nmidspecs>0);
            for ii = 1:length(iv)
                ivreftolonglist = [mainlistguide(ivNdiff(iv(ii)),2):mainlistguide(ivNdiff(iv(ii)),3)] - noffsetlonglist;
                longnmidspec(ivreftolonglist) = nmidspecs(iv(ii));
            end
            
            %------------------------------------------------------------------
            % edgepattern will be a matrix which, for each reflection combination, contains
            % the 2,3,4,... edge numbers that are involved in each path
            % regardless of there are specular reflections before, in the
            % middle, or after.
            
            edgepattern = zeros(nlonglist,Ndifforder);
            countvec = [1:nlonglist].';
            reflpathscut = reflpaths(ivlonglist,:);            
            for ii = 1:Ndifforder
                ivrefcol = countvec + (diffcols(rowinpatternlist,ii)-1)*nlonglist;
                edgepattern(:,ii) = reflpathscut(ivrefcol);
            end
	
            %------------------------------------------------------------------
            % midspecpattern will be a matrix which, for each reflection combination, contains
            % the 1,2,3,4,... specular reflection numbers (i.e., plane numbers) that are involved
            % in each path between diffraction 1 and diffraction 2.
            % First, the reflection component immediately after diffraction 1
            % is selected for every reflection. Then there is a masking which
            % zeroes the midspecpattern value for all the combinations that
            % don't have any midspec.
            % NB! For combinations like ssssdd, the first step will point to a
            % column which is outside the matrix. This is fixed by maximizing
            % the column number.
            
            midspecpattern = zeros(nlonglist,max([specorder-Ndifforder 1]));
            maxpointvalue = nlonglist*max([specorder-Ndifforder 1]);
            for ii = 1:specorder-Ndifforder
                ivrefcol = countvec + (diffcols(rowinpatternlist,1)+(ii-1))*nlonglist;
                ivrefcol = mod(ivrefcol-1,maxpointvalue)+1;
                midspecpattern(:,ii) = double(reflpathscut(ivrefcol)).*(longnmidspec>=ii);
            end
            diff1col = diffcols(rowinpatternlist,1);
	
            % Many combinations include the same edge-pair/N-let, so we extract the
            % unique edge pairs/N-lets, and go through them in the ii-loop
            
            [edgeshortlist,ivec,jvec] = unique(edgepattern,'rows');
            lastndivcomb = zeros(1,Ndifforder);
            
            for ii = 1: length(ivec)
                
                if SHOWTEXT >= 3
                    if round(ii/ceil(length(ivec)/10))*ceil(length(ivec)/10) == ii
                        disp(['      Combination no. ',int2str(ii),' of ',int2str(length(ivec))]) 
                    end
                end
                
                iv = find(jvec==ii);
                ncombs = length(iv);
                onesvec3 = ones(ncombs,1);
	
                if SHOWTEXT >= 4
                    numvec = int2str(edgeshortlist(ii,1));
                    for ll = 2:Ndifforder
                        numvec = [numvec,' ',  int2str(edgeshortlist(ii,ll))];
                    end
                    disp(['   Edges ',numvec])
                        
                end
                
                if Ndifforder >= 2 & any(multfactors(ivlonglist(iv)))
	
                    newndivvec = ndivvec(edgeshortlist(ii,:));
                    
                    if any(lastndivcomb~=newndivvec)
                        if Ndifforder == 2
                            ivmatrix = EDB2creindexmatrix(newndivvec);
                        else
                            ivmatrix = EDB2creindexmatrix(newndivvec(2:end));                        
                        end
                        [nedgeelcombs,slask] = size(ivmatrix);
                        if Ndifforder == 2
                            BIGEDGESTEPMATRIX = (double(ivmatrix)-0.5)./newndivvec(uint8(ones(nedgeelcombs,1)),:);
                        else
                            BIGEDGESTEPMATRIX = (double(ivmatrix)-0.5)./newndivvec(uint8(ones(nedgeelcombs,1)),2:end);
                        end
                        clear ivmatrix
                        lastndivcomb = newndivvec;
	
                    end
                end            

                pathalongplane = (edgeseespartialedge(edgeshortlist(ii,2),edgeshortlist(ii,1))>0);
                for jj = 3:Ndifforder
                    pathalongplane = [pathalongplane,(edgeseespartialedge(edgeshortlist(ii,jj),edgeshortlist(ii,jj-1))>0)];                    
                end
                
                if any(multfactors(ivlonglist(iv)))
	
                    IS = full(specextradata(ivlonglist(iv),1:3));
                    IR = full(specextradata(ivlonglist(iv),4:6));
                    
                    firstedgestart = full(edgeextradata(ivlonglist(iv),1));
                    firstedgeend   = full(edgeextradata(ivlonglist(iv),2));
		
                    lastedgestart = full(edgeextradata(ivlonglist(iv),3));
                    lastedgeend   = full(edgeextradata(ivlonglist(iv),4));
		
                    % Calculate the cyl coordinates of all IS/S and IR/R
		
                    cylS = zeros(ncombs,3);
                    edgecoords = [edgestartcoords(edgeshortlist(ii,1),:);edgeendcoords(edgeshortlist(ii,1),:)];
                	[cylS(:,1),cylS(:,2),cylS(:,3)] = EDB2coordtrans1(IS,edgecoords,edgenvecs(edgeshortlist(ii,1),:));               
		
                    cylR = zeros(ncombs,3);
                    edgecoords = [edgestartcoords(edgeshortlist(ii,Ndifforder),:);edgeendcoords(edgeshortlist(ii,Ndifforder),:)];
                	[cylR(:,1),cylR(:,2),cylR(:,3)] = EDB2coordtrans1(IR,edgecoords,edgenvecs(edgeshortlist(ii,Ndifforder),:));     
                        
                    bc = ones(Ndifforder,2);    % Check real reflfactors!!!    
                    bc = reshape(bc.',2*Ndifforder,1);
		
                    % Pick out the edge-to-edge coordinates for this specific
                    % edge pair/N-let.                        
                                                
                      if ~isempty(reftoshortlistE),                
                            if Ndifforder >= 2
                                index1 = reftoshortlistE(edgeshortlist(ii,2),edgeshortlist(ii,1));
                                cylE2_r1 = [re1sho(index1,:) thetae1sho(index1,:) ze1sho(index1,:);re2sho(index1,:) thetae2sho(index1,:) ze2sho(index1,:)];
                                index2 = reftoshortlistE(edgeshortlist(ii,1),edgeshortlist(ii,2));
                                cylE1_r2 = [re1sho(index2,:) thetae1sho(index2,:) ze1sho(index2,:);re2sho(index2,:) thetae2sho(index2,:) ze2sho(index2,:)];
                            end
                            if Ndifforder >= 3
                                index1 = reftoshortlistE(edgeshortlist(ii,3),edgeshortlist(ii,2));
                                cylE3_r2 = [re1sho(index1,:) thetae1sho(index1,:) ze1sho(index1,:);re2sho(index1,:) thetae2sho(index1,:) ze2sho(index1,:)];
                                index2 = reftoshortlistE(edgeshortlist(ii,2),edgeshortlist(ii,3));
                                cylE2_r3 = [re1sho(index2,:) thetae1sho(index2,:) ze1sho(index2,:);re2sho(index2,:) thetae2sho(index2,:) ze2sho(index2,:)];
                            end
                            if Ndifforder >= 4
                                index1 = reftoshortlistE(edgeshortlist(ii,4),edgeshortlist(ii,3));
                                cylE4_r3 = [re1sho(index1,:) thetae1sho(index1,:) ze1sho(index1,:);re2sho(index1,:) thetae2sho(index1,:) ze2sho(index1,:)];
                                index2 = reftoshortlistE(edgeshortlist(ii,3),edgeshortlist(ii,4));
                                cylE3_r4 = [re1sho(index2,:) thetae1sho(index2,:) ze1sho(index2,:);re2sho(index2,:) thetae2sho(index2,:) ze2sho(index2,:)];
                            end
                            if Ndifforder >= 5
                                index1 = reftoshortlistE(edgeshortlist(ii,5),edgeshortlist(ii,4));
                                cylE5_r4 = [re1sho(index1,:) thetae1sho(index1,:) ze1sho(index1,:);re2sho(index1,:) thetae2sho(index1,:) ze2sho(index1,:)];
                                index2 = reftoshortlistE(edgeshortlist(ii,4),edgeshortlist(ii,5));
                                cylE4_r5 = [re1sho(index2,:) thetae1sho(index2,:) ze1sho(index2,:);re2sho(index2,:) thetae2sho(index2,:) ze2sho(index2,:)];
                            end
                            if Ndifforder >= 6
                                index1 = reftoshortlistE(edgeshortlist(ii,6),edgeshortlist(ii,5));
                                cylE6_r5 = [re1sho(index1,:) thetae1sho(index1,:) ze1sho(index1,:);re2sho(index1,:) thetae2sho(index1,:) ze2sho(index1,:)];
                                index2 = reftoshortlistE(edgeshortlist(ii,5),edgeshortlist(ii,6));
                                cylE5_r6 = [re1sho(index2,:) thetae1sho(index2,:) ze1sho(index2,:);re2sho(index2,:) thetae2sho(index2,:) ze2sho(index2,:)];
                            end
                            if Ndifforder >= 7
                                index1 = reftoshortlistE(edgeshortlist(ii,7),edgeshortlist(ii,6));
                                cylE7_r6 = [re1sho(index1,:) thetae1sho(index1,:) ze1sho(index1,:);re2sho(index1,:) thetae2sho(index1,:) ze2sho(index1,:)];
                                index2 = reftoshortlistE(edgeshortlist(ii,6),edgeshortlist(ii,7));
                                cylE6_r7 = [re1sho(index2,:) thetae1sho(index2,:) ze1sho(index2,:);re2sho(index2,:) thetae2sho(index2,:) ze2sho(index2,:)];
                            end
                            if Ndifforder >= 8
                                index1 = reftoshortlistE(edgeshortlist(ii,8),edgeshortlist(ii,7));
                                cylE8_r7 = [re1sho(index1,:) thetae1sho(index1,:) ze1sho(index1,:);re2sho(index1,:) thetae2sho(index1,:) ze2sho(index1,:)];
                                index2 = reftoshortlistE(edgeshortlist(ii,7),edgeshortlist(ii,8));
                                cylE7_r8 = [re1sho(index2,:) thetae1sho(index2,:) ze1sho(index2,:);re2sho(index2,:) thetae2sho(index2,:) ze2sho(index2,:)];
                            end
                            
                       else
                            error('STRANGE TO END UP HERE??')
                            if Ndifforder == 2
                                cylE1_r2 = [0 0 edgelengthvec(edgeshortlist(ii,jj));0 0 edgelengthvec(edgeshortlist(ii,jj))]; 
                                cylE2_r1 = [0 0 edgelengthvec(edgeshortlist(ii,jj));0 0 edgelengthvec(edgeshortlist(ii,jj))];
                            else
                                error(['ERROR: Geometries with a single edge can not handle difforder >= 3'])
                            end
                       end
                        
                    for jj = 1:ncombs               

                        if multfactors(ivlonglist(iv(jj)))
                            
                            if Ndifforder == 2
                                
                                cylE1_r2frac = cylE1_r2;
                                e1length = edgelengthvec(edgeshortlist(ii,1));
                                cylE2_r1frac = cylE2_r1;
                                e2length = edgelengthvec(edgeshortlist(ii,2));
                                
                                if firstedgestart(jj) ~= 0
                                    cylE1_r2frac(1,3) = cylE1_r2frac(1,3) + e1length*firstedgestart(jj);
                                end
                                if firstedgeend(jj) ~= 1
                                    cylE1_r2frac(2,3) = cylE1_r2frac(1,3) + e1length*(1-firstedgeend(jj));
                                end
                                if lastedgestart(jj) ~= 0
                                    cylE2_r1frac(1,3) = cylE2_r1frac(1,3) + e2length*lastedgestart(jj);
                                end
                                if lastedgeend(jj) ~= 1
                                    cylE2_r1frac(2,3) = cylE2_r1frac(1,3)+ e2length*(1-lastedgeend(jj));
                                end
             
                                if midspecpattern(iv(jj))~=0
                                    
                                    % For the cases with specular reflections
                                    % in-between a double diffraction we must
                                    % mirror the two edges to calculate the
                                    % relative-to-edge cylindrical coordinates.
                                    
                                    nspec = sum(midspecpattern(iv(jj),:)>0);
			
                                    % Mirror edge 2 in all the specular reflection
                                    % planes, in reversed order
                                    edgestartmirr = edgestartcoords(edgeshortlist(ii,2),:);
                                    edgeendmirr =   edgeendcoords(edgeshortlist(ii,2),:);
                                    edgevector = edgeendmirr - edgestartmirr;
                                    if lastedgestart(jj) ~= 0
                                        edgestartmirr = edgestartmirr + edgevector*lastedgestart(jj);
                                    else
                                        edgestartmirr = edgestartmirr + edgevector*1e-6;
                                    end
                                    if lastedgeend(jj) ~= 1
                                        edgeendmirr =   edgestartmirr + edgevector*(1-lastedgeend(jj));
                                    else
                                        edgeendmirr =   edgestartmirr + edgevector*(1-1e-6);
                                    end
                                    edgerefcoords = [edgestartcoords(edgeshortlist(ii,1),:);edgeendcoords(edgeshortlist(ii,1),:)];
                                    edgerefnvec =   edgenvecs(edgeshortlist(ii,1),:);
			
                                    % If we have a specular reflection in a plane
                                    % which is perpendicular to the edge plane, we
                                    % should nudge the mirrored edge out a bit so
                                    % that there is no 0/(2*pi) mistake
                                    if nspec == 1 & ( edgeplaneperptoplane1(midspecpattern(iv(jj),1),edgeshortlist(ii,1)) | edgeplaneperptoplane1(midspecpattern(iv(jj),1),edgeshortlist(ii,2)) )
                                        vectowardsmidpoint = approxplanemidpoints(midspecpattern(iv(jj),1),:) - edgestartmirr;
                                        edgestartmirr = edgestartmirr + vectowardsmidpoint*1e-10;
                                        vectowardsmidpoint = approxplanemidpoints(midspecpattern(iv(jj),1),:) - edgeendmirr;
                                        edgeendmirr = edgeendmirr + vectowardsmidpoint*1e-10;
                                    end
                                    xis = [edgestartmirr;edgeendmirr];
                                    for kk = nspec:-1:1
                                        xis = EDB2findis(xis,[midspecpattern(iv(jj),kk);midspecpattern(iv(jj),kk)],planeeqs,2,onesvec1);
                                    end
                                    [rstart,thetastart,zstart,rend,thetaend,zend] = EDB2coordtrans2(xis(1,:),xis(2,:),edgerefcoords,edgerefnvec);
                                    cylE2mirr_r1 = [rstart thetastart zstart;rend thetaend zend];
			
                                    % Mirror edge 1 in all the specular reflection
                                    % planes, in forward order
                                    edgestartmirr = edgestartcoords(edgeshortlist(ii,1),:);
                                    edgeendmirr =   edgeendcoords(edgeshortlist(ii,1),:);
                                    edgevector = edgeendmirr - edgestartmirr;
                                    if firstedgestart(jj) ~= 0
                                        edgestartmirr = edgestartmirr + edgevector*firstedgestart(jj);
                                    else
                                        edgestartmirr = edgestartmirr + edgevector*1e-6;
                                    end
                                    if firstedgeend(jj) ~= 1
                                        edgeendmirr =   edgestartmirr + edgevector*(1-firstedgeend(jj));
                                    else
                                        edgeendmirr =   edgestartmirr + edgevector*(1-1e-6);                               
                                    end
                                    edgerefcoords = [edgestartcoords(edgeshortlist(ii,2),:);edgeendcoords(edgeshortlist(ii,2),:)];
                                    edgerefnvec =   edgenvecs(edgeshortlist(ii,2),:);
			
                                    % Normally, when there is a specular reflection
                                    % in-between, the diffraction path will not be
                                    % along a plane, unless: 
                                    % 1. The same edge is involved twice, with a
                                    % reflection in a perpendicular plane
                                    % in-between.
                                    % 2. See below
                                    pathalongplanewmidspec = 0;
                                    if edgeshortlist(ii,1) == edgeshortlist(ii,2)
                                        if thetastart == 0 | thetastart == (2*pi-closwedangvec(edgeshortlist(ii,2)))
                                            pathalongplanewmidspec = 1;
                                        end
                                    end
                                    
                                    % If we have a specular reflection in a plane
                                    % which is perpendicular to the edge plane, we
                                    % should nudge the mirrored edge out a bit so
                                    % that there is no 0/(2*pi) mistake
                                    %
                                    % We also have case 2 here (see above) for when
                                    % we could have a pathalongplanewmidspec
                                    % 2. When two different edges have a perpendicular reflection
                                    % in-between 
			
                                    if nspec == 1 & ( edgeplaneperptoplane1(midspecpattern(iv(jj),1),edgeshortlist(ii,1)) | edgeplaneperptoplane1(midspecpattern(iv(jj),1),edgeshortlist(ii,2)) )
                                        vectowardsmidpoint = approxplanemidpoints(midspecpattern(iv(jj),1),:) - edgestartmirr;
                                        edgestartmirr = edgestartmirr + vectowardsmidpoint*1e-10;
                                        vectowardsmidpoint = approxplanemidpoints(midspecpattern(iv(jj),1),:) - edgeendmirr;
                                        edgeendmirr = edgeendmirr + vectowardsmidpoint*1e-10;
                                        pathalongplanewmidspec = 1;
                                    end
                                    xis = [edgestartmirr;edgeendmirr];
                                    for kk = 1:nspec
                                        xis = EDB2findis(xis,[midspecpattern(iv(jj),kk);midspecpattern(iv(jj),kk)],planeeqs,2,onesvec1);
                                    end
                                    [rstart,thetastart,zstart,rend,thetaend,zend] = EDB2coordtrans2(xis(1,:),xis(2,:),edgerefcoords,edgerefnvec);
                                    cylE1mirr_r2 = [rstart thetastart zstart;rend thetaend zend];
                                    
                                    [irnew,slask] = EDB2wedge2nd(cylS(jj,:),cylR(jj,:),cylE2mirr_r1,cylE1mirr_r2,...
                                        nyvec(edgeshortlist(ii,:)),[edgelengthvec(edgeshortlist(ii,1))*[firstedgestart(jj) firstedgeend(jj)];edgelengthvec(edgeshortlist(ii,2))*[lastedgestart(jj) lastedgeend(jj)]],dzvec(edgeshortlist(ii,:)),...
                                        EDcalcmethod,pathalongplanewmidspec,Rstart,bc,CAIR,FSAMP);
                                    
                                else   %  ....   if midspecpattern(iv(jj))~=0		
                                    
									[irnew,slask] = EDB2wedge2nd(cylS(jj,:),cylR(jj,:),cylE2_r1frac,cylE1_r2frac,...
										nyvec(edgeshortlist(ii,:)),[edgelengthvec(edgeshortlist(ii,1))*[firstedgestart(jj) firstedgeend(jj)];edgelengthvec(edgeshortlist(ii,2))*[lastedgestart(jj) lastedgeend(jj)]],dzvec(edgeshortlist(ii,:)),...
										EDcalcmethod,pathalongplane,Rstart,bc,CAIR,FSAMP);
                                end   %  ....   if midspecpattern(iv(jj))~=0
			
                                if SHOWTEXT >= 6
                                    figure(1)
                                    plot(irnew)
                                   figure(2)
                                    plot(irdiff)
                                    pause
                                end  
                                IRDIFFVEC = [IRDIFFVEC;sum(irnew)];
                                			                
                            elseif Ndifforder == 3,    %   if Ndifforder == 2			
			
								for kk = 1:newndivvec(1)
									if SHOWTEXT >= 5
										disp(['   ',int2str(kk),' of ',int2str(newndivvec(1))]) 
									end
									BIGEDGE1stvalue = (kk-0.5)./newndivvec(1);
									wedgeparams = [cylE2_r1;cylE1_r2;cylE3_r2;cylE2_r3];
									
									[irnewpartition,slask] = EDB2wedgeN(cylS(jj,:),cylR(jj,:),wedgeparams,ncylrows,...
										nyvec(edgeshortlist(ii,:)),edgelengthvec(edgeshortlist(ii,:)).',...
										dzvec(edgeshortlist(ii,:)),EDcalcmethod,pathalongplane,nedgeelcombs,Rstart,bc,CAIR,FSAMP,BIGEDGE1stvalue);
									irnewpartition = real(irnewpartition);
																			
									if kk == 1
										irnew = irnewpartition;    
									else
										lengthaddition = length(irnewpartition);
										lengthaccum = length(irnew);
										if lengthaddition > lengthaccum
											irnew = [irnew;zeros(lengthaddition-lengthaccum,1)];    
										end
										irnew(1:lengthaddition) = irnew(1:lengthaddition) + irnewpartition; 
									end
								end
			
                                if SHOWTEXT >= 6
                                    sum(irnew)
                                    plot(irnew)
                                    pause
                                end
                            end
                            if Ndifforder == 4
                                    wedgeparams = [cylE2_r1;cylE1_r2;cylE3_r2;cylE2_r3;cylE4_r3;cylE3_r4];                                
                            elseif Ndifforder == 5
                                    wedgeparams = [cylE2_r1;cylE1_r2;cylE3_r2;cylE2_r3;cylE4_r3;cylE3_r4;cylE5_r4;cylE4_r5];
                            elseif Ndifforder == 6
                                    wedgeparams = [cylE2_r1;cylE1_r2;cylE3_r2;cylE2_r3;cylE4_r3;cylE3_r4;cylE5_r4;cylE4_r5;cylE6_r5;cylE5_r6];
                            elseif Ndifforder == 7
                                    wedgeparams = [cylE2_r1;cylE1_r2;cylE3_r2;cylE2_r3;cylE4_r3;cylE3_r4;cylE5_r4;cylE4_r5;cylE6_r5;cylE5_r6;cylE7_r6;cylE6_r7];
                            elseif Ndifforder == 8,                                
                                    wedgeparams = [cylE2_r1;cylE1_r2;cylE3_r2;cylE2_r3;cylE4_r3;cylE3_r4;cylE5_r4;cylE4_r5;cylE6_r5;cylE5_r6;cylE7_r6;cylE6_r7;cylE8_r7;cylE7_r8];
                            end
                            
                            if Ndifforder >= 4,    
                                
                                for kk = 1:newndivvec(1)
                                    if SHOWTEXT >= 5
                                        disp(['   ',int2str(kk),' of ',int2str(newndivvec(1))]) 
                                    end
                                    BIGEDGE1stvalue = (kk-0.5)./newndivvec(1);

                                    [irnewpartition,slask] = EDB2wedgeN(cylS(jj,:),cylR(jj,:),wedgeparams,ncylrows,...
                                        nyvec(edgeshortlist(ii,:)),edgelengthvec(edgeshortlist(ii,:)).',...
                                        dzvec(edgeshortlist(ii,:)),EDcalcmethod,pathalongplane,nedgeelcombs,Rstart,bc,CAIR,FSAMP,BIGEDGE1stvalue);
                                    irnewpartition = real(irnewpartition);
                                        
                                    if kk == 1
                                        irnew = irnewpartition;
                                    else
                                        lengthaddition = length(irnewpartition);
                                        lengthaccum = length(irnew);
                                        if lengthaddition > lengthaccum
                                            irnew = [irnew;zeros(lengthaddition-lengthaccum,1)];    
                                        end
                                        irnew(1:lengthaddition) = irnew(1:lengthaddition) + irnewpartition; 
                                    end
                                end
                            end
                                                            
                            if SHOWTEXT >= 5
                                varname = ['allirs',int2str(Ndifforder)];
                                if exist(varname) == 0
                                        eval([varname,' = irnew;'])
                                else
                                        lengthaddition = length(irnew);
                                        eval(['lengthaccum = size(',varname,',1);'])
                                        if lengthaddition > lengthaccum
		                                    eval([varname,' = [',varname,';zeros(lengthaddition-lengthaccum,size(',varname,',2))];'])
                                        end
                                        eval([varname,' = [',varname,' zeros(size(',varname,',1),1)];'])
                                        eval([varname,'(1:length(irnew),end) = irnew;'])
                                end
                             end
                            
                            % Decide whether the IR should be boosted or not
			
                            boostfactor = 1;
                            if souspecboost ~= 1 | recspecboost ~= 1
                            
                                if prespecpresent(iv(jj)) == 1
                                    boostfactor = souspecboost;    
                                else
                                    if ~isempty(Sinsideplanenumber)
                                        if Sinsideplanenumber(1)~=planesatedge(edgeshortlist(ii,1),1) & Sinsideplanenumber(1)~=planesatedge(edgeshortlist(ii,1),2), 
                                            boostfactor = souspecboost;
                                        end
                                    end    
			
                                end
                                if postspecpresent(iv(jj)) == 1
                                    boostfactor = boostfactor*recspecboost;    
                                else
                                    if ~isempty(Rinsideplanenumber)
                                        if Rinsideplanenumber(1)~=planesatedge(edgeshortlist(ii,Ndifforder),1) & Rinsideplanenumber(1)~=planesatedge(edgeshortlist(ii,Ndifforder),2), 
                                            boostfactor = boostfactor*recspecboost;
                                        end
                                    end
                                end
                            
                            end   % ...   if souspecboost ~= 1 | recspecboost ~= 1
                            
                            % For thin plates, we must have a boost factor!
                            % This is because there will be multiple equivalent
                            % combinations passing on the rear side of the thin plate
                            
                            if all( nyvec(edgeshortlist(ii,:)) == 0.5 )
                                boostfactor = boostfactor*2^(Ndifforder-1);                        
                            end
			
                            if boostfactor ~= 1
                                irnew = irnew*boostfactor;    
                            end
                            
                            irnew = irnew*double(multfactors(ivlonglist(iv(jj))));
                            
                            [ndiff,ncolsdiff] = size(irdiff);
                            nnew = size(irnew,1);
                            if saveindividualdiffirs(1) == 0
                                if nnew > ndiff
                                   irdiff = [irdiff;zeros(nnew-ndiff,1)];
                                end
                                irdiff(1:nnew) = irdiff(1:nnew) + irnew;
                            else
                                if Ndifforder > ncolsdiff
                                    irdiff = [irdiff zeros(ndiff,Ndifforder-ncolsdiff)];    
                                    ncolsdiff = ncolsdiff+1;
                                end
                                if nnew > ndiff
                                   irdiff = [irdiff;zeros(nnew-ndiff,Ndifforder)];
                                end

                                irdiff(1:nnew,Ndifforder) = irdiff(1:nnew,Ndifforder) + irnew;
                            end
                            
                        end    % ...  if multfactors(ivlonglist(iv(jj)))
                        
                    end    % ....for jj = 1:ncombs

                else   %  ...    if any(multfactors(ivlonglist(iv)))
                    if SHOWTEXT >= 4
                        disp(['      Combination not computed because of repetition'])                         
                    end
                    
                end    % ...    if any(multfactors(ivlonglist(iv)))
                
            end   % .... for ii = 1: length(ivec)
        
        end   % ...     if any(difforderinlist==Ndifforder)
    
        if size(irdiff,2) < Ndifforder & saveindividualdiffirs(1) == 1
            irdiff = [irdiff zeros(size(irdiff,1),Ndifforder-size(irdiff,2))];    
        end
	end    % ... for Ndifforder = 2:specorder
end     % .... if userwantsdiff2 == 1

ntot = length(irtot);
ndirect = length(irdirect);
ngeom = size(irgeom,1);
ndiff = size(irdiff,1);

if ndirect > ntot
   irtot = [irtot;zeros(ndirect-ntot,1)];
    ntot = length(irtot);
end
irtot(1:ndirect) = irtot(1:ndirect) + irdirect;

if ngeom > ntot
   irtot = [irtot;zeros(ngeom-ntot,1)];
    ntot = length(irtot);
end
irtot(1:ngeom) = irtot(1:ngeom) + irgeom;

if ndiff > ntot
   irtot = [irtot;zeros(ndiff-ntot,1)];
    ntot = length(irtot);
end
if saveindividualdiffirs(1) == 0 | userwantsdiff2 == 0
    irtot(1:ndiff) = irtot(1:ndiff) + irdiff;
else
     irtot(1:ndiff) = irtot(1:ndiff) + sum(irdiff.').';    
end

%-------------------------------------------------------
% Make the IRs the same length

nmax = max([ length(irtot) size(irgeom,1) size(irdiff,1) length(irdirect)]);

if length(irtot) < nmax
	irtot = [irtot;zeros(nmax-length(irtot),1)];
end
if length(irdirect) < nmax
	irdirect = [irdirect;zeros(nmax-length(irdirect),1)];
end
if length(irgeom) < nmax
	irgeom = [irgeom;zeros(nmax-size(irgeom,1),size(irgeom,2))];
end
if length(irdiff) < nmax
	irdiff = [irdiff;zeros(nmax-size(irdiff,1),size(irdiff,2))];
end

%-------------------------------------------------------
% Save the IRs

irtot = sparse(irtot);
irgeom = sparse(irgeom);
irdirect = sparse(irdirect);
irdiff = sparse(irdiff);
Varlist = [' irdirect irgeom irdiff irtot FSAMP Rstart CAIR EDcalcmethod'];

if length(saveindividualdiffirs) > 1
	if saveindividualdiffirs(2) == 1
		Varlist = [Varlist,' alldiffirs alldiffirsextradata'];
	end
end

if SHOWTEXT >= 5
    if specorder >= 2
        if exist('allirs2') == 0
            allirs2 = [];    
        end
        Varlist = [Varlist,' allirs2'];
    end       
    if specorder >= 3
        if exist('allirs3') == 0
            allirs3 = [];    
        end
        Varlist = [Varlist,' allirs3'];
    end       
    if specorder >= 4
        if exist('allirs4') == 0
            allirs4 = [];    
        end
        Varlist = [Varlist,' allirs4'];
    end       
    if specorder >= 5
        if exist('allirs5') == 0
            allirs5 = [];    
        end
        Varlist = [Varlist,' allirs5'];
    end       
    if specorder >= 6
        if exist('allirs6') == 0
            allirs6 = [];    
        end
        Varlist = [Varlist,' allirs6'];
    end       
    if specorder >= 7
        if exist('allirs7') == 0
            allirs7 = [];    
        end
        Varlist = [Varlist,' allirs7'];
    end       
    if specorder >= 8
        if exist('allirs8') == 0
            allirs8 = [];    
        end
        Varlist = [Varlist,' allirs8'];
    end       
end

eval(['save ',edirfile,Varlist])
