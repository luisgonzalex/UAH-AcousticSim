function edpathsfile = EDB2findpathsISESx(eddatafile,lengthNspecmatrix,lengthNdiffmatrix,...
    singlediffcol,startindicessinglediff,endindicessinglediff,...
    ndecimaldivider,PointertoIRcombs,IRoriginsfrom,S,R,isou,irec,...
    directsound,specorder,difforder,nedgesubs,visplanesfromS,visplanesfromR,...
    vispartedgesfromS,vispartedgesfromR,desiredfilename)
% EDB2findpathsISESx - Finds all possible paths that include direct sound, specular, diffraction.
% Finds all possible paths of types direct sound, specular, diffraction
% combinations of specular and diffraction.
%
% Input parameters:
% 	eddatafile          Taken directly from the setup-file or created automatically
%   lengthNspecmatrix, lengthNdiffmatrix
% 	singlediffcol, startindicessinglediff, endindicessinglediff
% 	ndecimaldivider, PointertoIRcombs, IRoriginsfrom
%                       Taken directly from the ISEStree-file
% 	S, R, isou, irec    Coordinates and counter number of source and receiver, taken
%                       from the setup file.
% 	directsound, specorder, difforder
%                       Taken directly from the setup-file.
% 	nedgesubs           Taken directly from the setup-file.
% 	visplanesfromS, visplanesfromR, vispartedgesfromS, vispartedgesfromR
%                       Taken directly from the srdatafile.
% 	desiredfilename     (optional) The desired name of the output file.
%   ISCOORDS (global)   The matrix ISCOORDS from the ISEStreefile
%   POTENTIALISES (global)  The matrix POTENTIALISES from the ISEStreefile
%   ORIGINSFROM (global)
%   ISESVISIBILITY (global)
%   IVNDIFFMATRIX (global)  The list IVNDIFFMATRIX from the ISEStreefile
%   IVNSPECMATRIX (global)
%   REFLORDER (global)
%
% See descriptions in EDB2findISEStree, EDB2srgeo and EDB2mainISES
%
% Output parameters:
%   edpathsfile             The output file name, where the output data is
%                           stored. If the optional input parameter desiredfilename, 
%                           is given, then edpathsfile will simple repeat
%                           this name. If desiredfilename was not specified
%                           a file name is created by extracting the file
%                           stem from the input parameter eddatafile, and
%                           with '_edpaths' added.
%
% Output data in the edpathsfile:
% 	pathtypevec             A matrix, [ncombs,specorder], describing the
% 	                        type of reflection, using 'f' for the direct
% 	                        sound, 's' for specular reflection and 'd' for
% 	                        diffraction.
% 	reflpaths               A matrix, [ncombs,specorder], giving the plane 
%                           and edge numbers involved in each refl. comb.
% 	specextradata           A sparse matrix, [ncombs,(specorder+1)*3]
% 	                        containing the coordinates of the image source
% 	                        (col 1-3), the image receiver (col 4-6)
%                           and the coordinates of all specular hit points
%                           (only for purely specular combinations).
% 	edgeextradata           A sparse matrix, [ncombs,(specorder*2)], with
% 	                        the visibility of the involved edge, denoted
%                           by two values between 0 and 1.
% 	S                       The coordinates of the source.
% 	R                       The coordinates of the receiver.
% 	mainlistguide           A matrix, [nposscombs,3], which for each row gives
%                           1. the number of rows in "reflpaths" and
%                           "pathtypevec" that have one type of combination.
%                           2. the first row and 3. the last row
% 	mainlistguidepattern    A matrix, [nposscombs,specorder] which for each
% 	                        row gives the description of each possible combination
%                           using 'f', 's' and 'd'.
%   directsoundrow          The value 0 or 1 indicating whether the first row in
%                           mainlistguide & mainlistguidepattern & reflpaths & pathtypevec 
%                           contains the direct sound or not.
%   allspecrows             A list, [1,2], containing the first and last row numbers
%                           of mainlistguide & mainlistguidepattern that contain specular
%                           reflections. If there are no specular reflections, allspecrows
%                           contains [0 0].
%   firstdiffrow            A number indicating which is the first row of mainlistguide &
%                           mainlistguidepattern that contain a diffraction component. If
%                           there are no diffraction components, firstdiffrow = 0.
%   Sinsideplanenumber      The plane number that the source is directly at. If the
%                           source is not placed directly at a plane, then
%                           souinsideplanenumber = [];
%   Rinsideplanenumber      The plane number that the receiver is directly at. If the
%                           receiver is not placed directly at a plane, then
%                           recinsideplanenumber = [];
%
% Uses the functions EDB2strpend EDB2directsound
% EDB2speculISES EDB2diffISESx EDB2diff2ISES EDB2diffNISES
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
% Peter Svensson (svensson@iet.ntnu.no) 20050404
%
% edpathsfile = EDB2findpathsISES(eddatafile,...
%     lengthNspecmatrix,lengthNdiffmatrix,singlediffcol,startindicessinglediff,...
%     endindicessinglediff,ndecimaldivider,PointertoIRcombs,IRoriginsfrom,S,R,isou,irec,...
%     directsound,specorder,difforder,nedgesubs,visplanesfromS,visplanesfromR,...
%     vispartedgesfromS,vispartedgesfromR,desiredfilename);

global SHOWTEXT 
global POTENTIALISES ISCOORDS IVNDIFFMATRIX
global IVNSPECMATRIX ORIGINSFROM ISESVISIBILITY REFLORDER

eval(['load ',eddatafile])
clear cornerinfrontofplane edgeseesplane planeseesplane

if nargin < 22
    [filepath,filestem,fileext] = fileparts(eddatafile);
    edpathsfile = [[filepath,filesep],EDB2strpend(filestem,'_eddata'),'_',int2str(isou),'_',int2str(irec),'_edpaths.mat'];    
else
    edpathsfile = desiredfilename;
end

nplanes = length(planeisthin);

pathtypevec = [];
reflpaths = [];
specextradata = [];
edgeextradata = [];
mainlistguide = [];
mainlistguidepattern = [];

%-------------------------------------------------------
%		#############################
%		##     DIRECT SOUND        ##
%		#############################

if directsound == 1
	if SHOWTEXT >= 2
		disp('   Direct sound')
	end
    dirsoundok = EDB2directsound(eddatafile,S,R,visplanesfromS,visplanesfromR);
    if dirsoundok == 1
        pathtypevec = [pathtypevec;['f',zeros(1,specorder-1)]];
        reflpaths = [reflpaths;zeros(1,specorder)];
        if specorder > 0
            specextradata = [specextradata;[S R zeros(1,(specorder-1)*3)]];
        else
            specextradata = [specextradata;[S R]];   
        end
        edgeextradata = [edgeextradata;zeros(1,specorder*2)];    
        mainlistguide = [1 1 1];   
        mainlistguidepattern = ['f' ' '*ones(1,specorder-1)];
        directsoundrow = 1;
    else
        directsoundrow = 0;    
    end
else
    directsoundrow = 0;
end

%-------------------------------------------------------
%		#############################
%		##        SPECULAR         ##
%		#############################

if specorder >= 1
	if SHOWTEXT >= 2
		disp('   SPECULAR reflections')
	end
	[validISlist,validIScoords,allreflpoints,listguide,listofreflorder] = EDB2speculISES(eddatafile,...
        S,R,lengthNspecmatrix,specorder,visplanesfromR);
    [nvalidreflorders,slask] = size(listguide);
    for ii = 1:nvalidreflorders
        ncombs = listguide(ii,1);
        norder = listofreflorder(ii);
        
        if ncombs > 0
            pathtypevec = [pathtypevec;['s'*ones(ncombs,norder) zeros(ncombs,specorder-norder)]];
            reflpaths = [reflpaths;validISlist(listguide(ii,2):listguide(ii,3),:)];
            if size(pathtypevec,2) > size(reflpaths,2)
                reflpaths = [reflpaths zeros(size(reflpaths,1),size(pathtypevec,2)-size(reflpaths,2))];    
            end
            specextradata = [specextradata;[validIScoords(listguide(ii,2):listguide(ii,3),1:3) allreflpoints(listguide(ii,2):listguide(ii,3),:)]];
            if size(specextradata,2) < (size(pathtypevec,2)+1)*3
                specextradata = [specextradata zeros(size(specextradata,1),(size(pathtypevec,2)+1)*3-size(specextradata,2))];
            end
            edgeextradata = [edgeextradata;zeros(ncombs,specorder*2)];    
        end            

        mainlistguidepattern = [mainlistguidepattern;['s'*ones(1,norder) ' '*ones(1,specorder-norder)]];
    
    end

    [nprevious,slask] = size(mainlistguide);
    if nprevious > 0
        listguide(:,2:3) = listguide(:,2:3)+mainlistguide(nprevious,3);
    end
    mainlistguide = [mainlistguide;listguide];
    if nvalidreflorders == 0
        allspecrows = [0 0];
    else        
        allspecrows = [min([1 nvalidreflorders]) nvalidreflorders]+directsoundrow;
    end
    firstdiffrow = nvalidreflorders + directsoundrow + 1;
   
else
    allspecrows = [0 0];        
    firstdiffrow = directsoundrow + 1;
end

%-------------------------------------------------------
%		####################################
%		##      SINGLE DIFFRACTION        ##
%		####################################

if difforder >=1 & ~isempty(lengthNdiffmatrix)
	if SHOWTEXT >= 2
		disp('   SINGLE DIFFRACTION combined with any order specular')
	end
    [edgedifflist,startandendpoints,prespeclist,postspeclist,validISEDcoords,validEDIRcoords,listguide,listoforders,...
        bigedgeweightlist] = EDB2diffISESx(eddatafile,S,R,...
        IVNDIFFMATRIX(1:lengthNdiffmatrix(1),1),singlediffcol,startindicessinglediff,endindicessinglediff,...
        specorder,visplanesfromR,vispartedgesfromS,vispartedgesfromR,nedgesubs,ndecimaldivider,PointertoIRcombs,IRoriginsfrom);    
    [nvalidcombs,slask] = size(listguide);
    nprespecs = listoforders(:,2) - 1;
    npostspecs = listoforders(:,1) - nprespecs - 1;
    
    [nrows,nsubsegmentcols] = size(startandendpoints);
    nsubsegments = uint8( (startandendpoints(:,2:2:nsubsegmentcols))~=0 );

    if nrows > 1 & nsubsegmentcols > 2
        nsubsegments = sum(nsubsegments.').';
    end
    nmaxsubsegments = max(nsubsegments);

    expandedlistguide = listguide;
    for ii = 1:nvalidcombs
        ncombs = listguide(ii,1);
        iv = [listguide(ii,2):listguide(ii,3)];

        pathtypevec   = [pathtypevec;  ['s'*ones(ncombs,nprespecs(ii))       'd'*ones(ncombs,1)         's'*ones(ncombs,npostspecs(ii))       zeros(ncombs,specorder-listoforders(ii)) ]];
        reflpaths     = [reflpaths;    [prespeclist(iv,1:nprespecs(ii))   edgedifflist(iv)         postspeclist(iv,1:npostspecs(ii))  zeros(ncombs,specorder-listoforders(ii)) ]];  
        specextradata = [specextradata;[  validISEDcoords(iv,1:3)         validEDIRcoords(iv,1:3)    zeros(ncombs,(specorder-1)*3)            ]];
        edgeextradata = [edgeextradata;[startandendpoints(iv,1:2)           zeros(ncombs,(specorder-1)*2)                         ]];    
        mainlistguidepattern = [mainlistguidepattern;['s'*ones(1,nprespecs(ii))       'd'         's'*ones(1,npostspecs(ii))       zeros(1,specorder-listoforders(ii))]];
        if nmaxsubsegments > 2
        for jj = 2:nmaxsubsegments
            ivsub = find(nsubsegments(iv)>=jj);
            ncombssubset = length(ivsub);
            if ncombssubset > 0
                pathtypevec   = [pathtypevec;  ['s'*ones(ncombssubset,nprespecs(ii))       'd'*ones(ncombssubset,1)         's'*ones(ncombssubset,npostspecs(ii))       zeros(ncombssubset,specorder-listoforders(ii)) ]];
                reflpaths     = [reflpaths;    [prespeclist(iv(ivsub),1:nprespecs(ii))   edgedifflist(iv(ivsub))         postspeclist(iv(ivsub),1:npostspecs(ii))  zeros(ncombssubset,specorder-listoforders(ii)) ]];  
                specextradata = [specextradata;[  validISEDcoords(iv(ivsub),1:3)         validEDIRcoords(iv(ivsub),1:3)    zeros(ncombssubset,(specorder-1)*3)            ]];
                edgeextradata = [edgeextradata;[startandendpoints(iv(ivsub),(jj-1)*2+1:(jj-1)*2+2)           zeros(ncombssubset,(specorder-1)*2)                         ]];
                expandedlistguide(ii,1) = expandedlistguide(ii,1) + ncombssubset;
                expandedlistguide(ii,3) = expandedlistguide(ii,3) + ncombssubset;
                for kk = ii+1:nvalidcombs
                    expandedlistguide(kk,2:3) = expandedlistguide(kk,2:3) + ncombssubset; 
                end
            end    
        end
        end
    end
    [n1,n2] = size(mainlistguide);
    if n1 > 0
        listoffset = mainlistguide(n1,3);
    else
        listoffset = 0;    
    end
    if size(expandedlistguide,1) > 0
        mainlistguide = [mainlistguide;[expandedlistguide(:,1) expandedlistguide(:,2:3)+listoffset]];
    end
end

%-------------------------------------------------------
%		################################################
%		##      DOUBLE AND HIGHER DIFFRACTION         ##
%		################################################

for kk = 2:min([specorder difforder])

    if ~isempty(lengthNdiffmatrix)
        
    if SHOWTEXT >= 2
		disp(['   DIFFRACTION, ORDER ',int2str(kk),' combined with any order specular'])
	end

    nmaxsubsegments = 0;
    
    % We will use the spec-edge-spec combs that were found valid for
    % first-order diffraction, because if higher-order diffraction combs
    % end with the same edge-spec-spec-... comb we know the visibility
    % already.
    % First, pick out only the combs that have postspecs. Second, keep
    % only a unique set because there will be many repeated combinations.
    
   
    if kk == 2
        iv = find(sum(postspeclist.')>0 );
        edgedifflist = edgedifflist(iv,:);
        postspeclist = postspeclist(iv,:);
        bigedgeweightlist = bigedgeweightlist(iv,:);
        validEDIRcoords = validEDIRcoords(iv,:);

        patternthatisOK = [edgedifflist postspeclist];
        [patternthatisOK,iv,slask] = unique(patternthatisOK,'rows');
        edgedifflistin = edgedifflist(iv,:);
        postspeclistin = postspeclist(iv,:);
        bigedgeweightlistin = bigedgeweightlist(iv,:);
        validEDIRcoordsin = validEDIRcoords(iv,:);

    end

    if kk == 2
        maxrownumber = max(lengthNdiffmatrix(1:2));
        [edgedifflist,startandendpoints,prespeclist,midspeclist,postspeclist,validISEDcoords,validEDIRcoords,listguide,listofallspecs] = EDB2diff2ISES(eddatafile,S,R,...
            IVNDIFFMATRIX(1:maxrownumber,1:2),lengthNdiffmatrix(1:2),...
            specorder,visplanesfromR,vispartedgesfromS,vispartedgesfromR,nedgesubs,ndecimaldivider,edgedifflistin,postspeclistin,...
            bigedgeweightlistin,validEDIRcoordsin,edgeplaneperptoplane1,edgeplaneperptoplane2);

        [nrows,nsubsegmentcols] = size(startandendpoints);
        nsubsegments = uint8((startandendpoints(:,2:4:nsubsegmentcols))~=0);
        if nrows > 1 & nsubsegmentcols > 4
            nsubsegments = sum(nsubsegments.').';
        end
        nmaxsubsegments = max(nsubsegments);

    elseif kk >= 3
        maxrownumber = max(lengthNdiffmatrix(1:kk));

        [edgedifflist,startandendpoints,prespeclist,postspeclist,validISEDcoords,validEDIRcoords,listguide,listoforders] = EDB2diffNISES(eddatafile,S,R,...
            IVNDIFFMATRIX(1:maxrownumber,1:kk),lengthNdiffmatrix(1:kk),kk,...
            specorder,visplanesfromR,vispartedgesfromS,vispartedgesfromR,nedgesubs,ndecimaldivider,edgedifflistin,postspeclistin,...
            bigedgeweightlistin,validEDIRcoordsin);
    end

    [nvalidcombs,slask] = size(listguide); 
    if kk >= 3
        nprespecs = listoforders(:,2) - 1;
        npostspecs = listoforders(:,1) - nprespecs - kk;
    else
        nprespecs = listofallspecs(:,1);    
        nmidspecs = listofallspecs(:,2);    
        npostspecs = listofallspecs(:,3);    
        listoforders = nprespecs+nmidspecs+npostspecs+2;
    end
    
    for ii = 1:nvalidcombs
        ncombs = listguide(ii,1);
        i1 = listguide(ii,2);     i2 = listguide(ii,3);
        if kk == 2
            pathtypevec   = [pathtypevec;  ['s'*ones(ncombs,nprespecs(ii))  'd'*ones(ncombs,1)   's'*ones(ncombs,nmidspecs(ii))  'd'*ones(ncombs,1)    's'*ones(ncombs,npostspecs(ii))       zeros(ncombs,specorder-listoforders(ii)) ]];
        else
            pathtypevec   = [pathtypevec;  ['s'*ones(ncombs,nprespecs(ii))  'd'*ones(ncombs,kk)   's'*ones(ncombs,npostspecs(ii))       zeros(ncombs,specorder-listoforders(ii)) ]];            
        end
        
        if kk == 2        
            reflpaths     = [reflpaths;    [prespeclist(i1:i2,1:nprespecs(ii))   edgedifflist(i1:i2,1)  midspeclist(i1:i2,1:nmidspecs(ii)) edgedifflist(i1:i2,2)         postspeclist(i1:i2,1:npostspecs(ii))  zeros(ncombs,specorder-listoforders(ii)) ]];  
        else
            reflpaths     = [reflpaths;    [prespeclist(i1:i2,1:nprespecs(ii))   edgedifflist(i1:i2,1:kk)         postspeclist(i1:i2,1:npostspecs(ii))  zeros(ncombs,specorder-listoforders(ii)) ]];              
        end
        
        specextradata = [specextradata;[  validISEDcoords(i1:i2,1:3)         validEDIRcoords(i1:i2,1:3)    zeros(ncombs,(specorder-1)*3)            ]];
        edgeextradata = [edgeextradata;[startandendpoints(i1:i2,1:kk*2)           zeros(ncombs,(specorder-kk)*2)                         ]];    
        if kk == 2
            mainlistguidepattern = [mainlistguidepattern;['s'*ones(1,nprespecs(ii))   'd' 's'*ones(1,nmidspecs(ii)) 'd'         's'*ones(1,npostspecs(ii))       zeros(1,specorder-listoforders(ii))]];
        else
            mainlistguidepattern = [mainlistguidepattern;['s'*ones(1,nprespecs(ii))   'd'*ones(1,kk)  's'*ones(1,npostspecs(ii))       zeros(1,specorder-listoforders(ii))]];            
        end
        
        for jj = 2:nmaxsubsegments
            ivsub = find(nsubsegments(iv)>=jj);
            ncombssubset = length(ivsub);
            if ncombssubset > 0
                pathtypevec   = [pathtypevec;  ['s'*ones(ncombssubset,nprespecs(ii))       'd'*ones(ncombssubset,1) 's'*ones(ncombssubset,nmidspecs(ii)) 'd'*ones(ncombssubset,1)        's'*ones(ncombssubset,npostspecs(ii))       zeros(ncombssubset,specorder-listoforders(ii)) ]];
                reflpaths     = [reflpaths;    [prespeclist(iv(ivsub),1:nprespecs(ii))   edgedifflist(iv(ivsub),1)  midspeclist(iv(ivsub),1:nmidspecs(ii))  edgedifflist(iv(ivsub),2)         postspeclist(iv(ivsub),1:npostspecs(ii))  zeros(ncombssubset,specorder-listoforders(ii)) ]];  
                specextradata = [specextradata;[  validISEDcoords(iv(ivsub),1:3)         validEDIRcoords(iv(ivsub),1:3)    zeros(ncombssubset,(specorder-1)*3)            ]];
                edgeextradata = [edgeextradata;[startandendpoints(iv(ivsub),(jj-1)*4+1:(jj-1)*4+4)           zeros(ncombssubset,(specorder-1)*4)                         ]];
                expandedlistguide(ii,1) = expandedlistguide(ii,1) + ncombssubset;
                expandedlistguide(ii,3) = expandedlistguide(ii,3) + ncombssubset;
                for kk = ii+1:nvalidcombs
                    expandedlistguide(kk,2:3) = expandedlistguide(kk,2:3) + ncombssubset; 
                end
            end    
         end
         
     end

    [n1,n2] = size(mainlistguide);
    if n1 > 0
        listoffset = mainlistguide(n1,3);
    else
        listoffset = 0;    
    end
    mainlistguide = [mainlistguide;[listguide(:,1) listguide(:,2:3)+listoffset]];

    end

end    

[n1,n2] = size(mainlistguide);
if firstdiffrow > n1
    firstdiffrow = 0;    
end

%-------------------------------------------------------
% Is the source or receiver directly at a plane?

Sinsideplanenumber = find(visplanesfromS==4);

Rinsideplanenumber = find(visplanesfromR==4);

%-------------------------------------------------------
% Get rid of empty columns

if ~isempty(pathtypevec)
    checksum = sum(pathtypevec);
    ncols = length(checksum);
    while checksum(ncols) == 0
        pathtypevec(:,ncols) = [];
        checksum = sum(pathtypevec);
        ncols = length(checksum);
    end
end
    
%-------------------------------------------------------
% Save the output data

specextradata = sparse(specextradata);
edgeextradata = sparse(edgeextradata);
pathtypevec = uint8(pathtypevec);
[ncombs,slask] = size(reflpaths);
if ncombs+1 < 256
    mainlistguide = uint8(mainlistguide);
elseif ncombs+1 < 65536
    mainlistguide = uint16(mainlistguide);
else
    mainlistguide = uint32(mainlistguide);
end

Varlist = [' pathtypevec reflpaths specextradata edgeextradata S R mainlistguide mainlistguidepattern directsoundrow allspecrows firstdiffrow Sinsideplanenumber Rinsideplanenumber'];		

eval(['save ',edpathsfile,Varlist])
