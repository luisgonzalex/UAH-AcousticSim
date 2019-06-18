function [edgedifflist,startandendpoints,prespeclist,midspeclist,postspeclist,validIScoords,validIRcoords,listguide,...
    listofallspecs] = EDB2diff2ISES(eddatafile,S,R,ivNdiffmatrix,...
    lengthNdiffmatrix,specorder,visplanesfromR,vispartedgesfromS,vispartedgesfromR,nedgesubs,...
    ndecimaldivider,edgedifflistin,postspeclistin,bigedgeweightlistin,validEDIRcoords,edgeplaneperptoplane1,edgeplaneperptoplane2)
% EDB2diff2ISES - Gives list of paths that includes a 2nd-order diff. combination.
% EDB2diff2ISES gives the list of possible diffraction paths that includes one
% second-order diffraction path, and possibly specular reflections before
% and after.
%
% Input parameters:
%       eddatafile,S,R,ivNdiffmatrix,...
%       lengthNdiffmatrix,specorder,visplanesfromR,vispartedgesfromS,vispartedgesfromR,nedgesubs,...
%       ndecimaldivider,edgeplaneperptoplane1,edgeplaneperptoplane2
%                       All these are taken from the ISEStreefile or the
%                       setup file.
%   edgedifflistin      List [nd1combs,1] of edges involved in all
%                       first-order diffraction combs that have already
%                       been found
%   postspeclistin      Matrix [nd1combs,specorder-1] of all specular reflections 
%                       following the diffraction for all first-order diff combs
%                       that have already been found
%   bigedgeweightlistin List [nd1combs,1] of visibility for edges involved in all
%                       first-order diffraction combs that have already been found
%   validEDIRcoords     Matrix [nd1combs,3] of image receiver coordinates for all 
%                       first-order diff combs that have already been found
% GLOBAL parameters:
%   POTENTIALISES,ISCOORDS,ORIGINSFROM,ISESVISIBILITY,REFLORDER    See EDB2findISEStree
%   SHOWTEXT JJ JJnumbofchars   See EDB2mainISES
%
% Output parameters:
%   edgedifflist        List [ncombs,2] of the edge numbers involved in each
%                       spec-diff-diff-spec combination.
%   startandendpoints   Matrix [ncombs,4] of the relative start and end
%                       points of each edge. The values, [0,1], indicate
%                       which part of the two edges that are visible.
%   prespeclist         Matrix [ncombs,specorder-2] of the specular
%                       reflections that precede every diffraction.
%   midspeclist         Matrix [ncombs,specorder-2] of the specular
%                       reflections inbetween the two diffractions.
%   postspeclist        Matrix [ncombs,specorder-2] of the specular
%                       reflections that follow every diffraction.
%   validIScoords       Matrix [ncombs,3] of the image source for each
%                       multiple-spec that precedes the diffraction. If
%                       there is no spec refl before the diffraction, the
%                       value [0 0 0] is given.
%   validIRcoords       Matrix [ncombs,3] of the image receiver for each
%                       multiple-spec that follows the diffraction. If
%                       there is no spec refl after the diffraction, the
%                       value [0 0 0] is given.
%   listguide           Matrix [nuniquecombs,3] which for each row gives
%                       1. The number of examples in edgefdifflist etc that
%                          are the same type of spec-diff-diff-spec comb.
%                       2. The first row number and 3. The last row number.
%   listofallspecs      Matrix [nuniquecombs,3] which for each row gives
%                       1. The number of pre-specular reflections for the spec-diff-spec-diff-spec comb
%                          in the same row in listguide.
%                       2. The number of mid-specular reflections for the spec-diff-spec-diff-spec comb
%                          in the same row in listguide.
%                       3. The number of post-specular reflections for the spec-diff-spec-diff-spec comb
%                          in the same row in listguide.
%
% Uses functions  EDB2findis EDB2getedgepoints EDB2chkISvisible EDB2checkobstrpaths
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
% Peter Svensson (svensson@iet.ntnu.no) 20050202
%
% [edgedifflist,startandendpoints,prespeclist,midspeclist,postspeclist,validIScoords,validIRcoords,listguide,...
%   listofallspecs] = EDB2diff2ISES(eddatafile,S,R,...
%   ivNdiffmatrix,lengthNdiffmatrix,specorder,visplanesfromR,vispartedgesfromS,...
%   vispartedgesfromR,nedgesubs,ndecimaldivider,edgedifflistin,postspeclistin,bigedgeweightlist,validEDIRcoords,edgeplaneperptoplane1,edgeplaneperptoplane2)

global SHOWTEXT JJ JJnumbofchars
global POTENTIALISES ISCOORDS ORIGINSFROM ISESVISIBILITY REFLORDER

eval(['load ',eddatafile])

[nedges,slask] = size(planesatedge);
[nplanes,slask] = size(planecorners);
multfac = 10^(ceil(log10(max([nedges nplanes]))));

edgedifflist = [];
prespeclist = [];
midspeclist = [];
postspeclist = [];
startandendpoints = [];
bigedgeweightlist = [];
validIScoords = [];
validIRcoords = [];

[n1,n2] = size(POTENTIALISES);
if n2 < specorder
    specorder = n2;    
end

maxvisibilityvalue = 2^nedgesubs-1;
zerosvec1 = zeros(1,specorder-2);
zerosvec2 = zeros(1,3);
listguide = zeros(specorder*2-1,3);
bitmultvec = 2.^[0:nedgesubs-1];

obstructtestneeded = (sum(canplaneobstruct)~=0);
onesvec = ones(1,nedgesubs);
onesvec3 = ones(1,3);

npostspecs = sum(double(postspeclistin.'>0)).';
planesatedge = double(planesatedge);

if specorder >= 3
    coplanarsviaflatedge = sparse(zeros(nplanes,nplanes));
    listofflatedges = find(closwedangvec==pi);
    if ~isempty(listofflatedges)
        ivreftomatrix = planesatedge(listofflatedges,1) + (planesatedge(listofflatedges,2)-1)*nplanes;
        coplanarsviaflatedge(ivreftomatrix) = ones(size(ivreftomatrix));
        coplanarsviaflatedge = sign(coplanarsviaflatedge + coplanarsviaflatedge.');
    end
end

if specorder >= 4
    cateyeplanecombs = sparse(zeros(nplanes,nplanes));
    listof90edges = find(closwedangvec==3*pi/2);
    if ~isempty(listof90edges)
        ivreftomatrix = planesatedge(listof90edges,1) + (planesatedge(listof90edges,2)-1)*nplanes;
        cateyeplanecombs(ivreftomatrix) = ones(size(ivreftomatrix));
        cateyeplanecombs = sign(cateyeplanecombs + cateyeplanecombs.');
    end
    
end

%   ###########################################
%   #                                         #
%   #         S - edge - edge - R cases       #
%   #                                         #
%   ###########################################
%
% Possible edges for S-E-E-R are seen (at least partly) by the source and by the receiver.
%
% The visibility by the source is taken care of by the findISEStree
% so we must check which ones are visible by the receiver.
% Also, the active edge segment(s) must be selected but this is done
% further down together with the S-spec-spec-edge-R cases

ivNdiff = ivNdiffmatrix(1:lengthNdiffmatrix(2),2);
ivsinglediff = ivNdiffmatrix(1:lengthNdiffmatrix(1),1);

ivSEER = ivNdiff(find(REFLORDER(ivNdiff)==2));

possibleedgepairs = double(POTENTIALISES(ivSEER,1:2))-nplanes;
ivnotvisiblefromr = find(vispartedgesfromR(possibleedgepairs(:,2))==0);
if ~isempty(ivnotvisiblefromr)
    possibleedgepairs(ivnotvisiblefromr,:) = [];
end

edgedifflist      = [edgedifflist;possibleedgepairs];
bigedgeweightlist = [bigedgeweightlist;[vispartedgesfromS(edgedifflist(:,1)) vispartedgesfromR(edgedifflist(:,2))]];

[nedgesadded,slask] = size(edgedifflist);
zerosvec3 = zeros(nedgesadded,1);
ndiffonly = nedgesadded;

prespeclist =  [prespeclist;   zerosvec3(:,ones(1,max(specorder-2,1)))];
midspeclist =  [midspeclist;  zerosvec3(:,ones(1,max(specorder-2,1)))];
postspeclist = [postspeclist;  zerosvec3(:,ones(1,max(specorder-2,1)))];
validIScoords = [validIScoords;S(ones(nedgesadded,1),:)];
validIRcoords = [validIRcoords;R(ones(nedgesadded,1),:)];

if SHOWTEXT >= 3
	disp(['         ',int2str(nedgesadded),' double diff valid'])
end

%   ###########################################
%   #                                         #
%   #      S - spec - edge - edge - R cases   #
%   #                                         #
%   #         Prespec cases                   #
%   #                                         #
%   ###########################################
%
% Possible edges for S-spec-E-E-R are seen (at least partly) by the receiver.
%
% The visibility doesn't need to be checked since the source-to-edge paths
% were checked in the ISEStree, and the visibility from the receiver also
% has been checked.

% The vector ivmultidiff will always refer to the original data vector
% i.e. POTENTIALISES, ORIGINSFROM, REFLORDER etc
%
% We should remove the combinations that involve an edge which the
% receiver can not see, but since the edge number is in different columns
% we do that in the for loop.

% The ii-loop will go through: spec-diff, spec-spec-diff
% spec-spec-spec-diff etc

for ii = 1:specorder-2
    
    if SHOWTEXT >= 3
        disp(['      Checking for ',JJ(ii,1:JJnumbofchars(ii)),' spec refl before the double edge diff'])    
    end

    % Select the combinations where the reflection order == ii+2
    % which means ii specular reflections before the diffraction
    
    iv = find(REFLORDER(ivNdiff) == ii+2 & POTENTIALISES(ivNdiff,ii+1)>nplanes & POTENTIALISES(ivNdiff,ii+2)>nplanes);
    masterivlist = ivNdiff(iv);
    possibleedgepairs = double(POTENTIALISES(masterivlist,ii+1:ii+2)) - nplanes;
    
    % Keep only combinations for which the receiver can see the edge
    
    ivnotvisiblefromr = find(vispartedgesfromR(possibleedgepairs(:,2))==0);
    if ~isempty(ivnotvisiblefromr)
        masterivlist(ivnotvisiblefromr) = [];
        possibleedgepairs(ivnotvisiblefromr,:) = [];   
    end        

    possiblecombs = POTENTIALISES(masterivlist,1:ii);
    
    edgeweightlist = [ISESVISIBILITY(masterivlist) vispartedgesfromR(possibleedgepairs(:,2))];
        
    nposs = length(masterivlist);

    if SHOWTEXT >= 3
 		disp(['         ',int2str(nposs),' IS+edge pairs valid'])
    end

    edgedifflist = [edgedifflist;possibleedgepairs];
    prespeclist = [prespeclist;[possiblecombs zeros(nposs,specorder-2-ii)]];        
    midspeclist = [midspeclist;zerosvec1(ones(nposs,1),:)];
    postspeclist = [postspeclist;zerosvec1(ones(nposs,1),:)];
    bigedgeweightlist = [bigedgeweightlist;edgeweightlist];    
   
    % NB! It is correct below that the indices for the ISCOORDS should be
    % ORIGINSFROM(ORIGINSFROM(masterivlist)), rather than masterivlist.
    % The combinations in POTENTIALISES(masterivlist,:) all have
    % spec-spec-...-diff-diff combinations and then
    % ISCOORDS(masterivlist,:) are zeros since a comb. that
    % ends with a diff has no image source. 
    % Also, two recursive references are needed since we need to get back
    % through the two last diffractions to reach the last specular
    % reflection.

    validIScoords = [validIScoords;ISCOORDS(ORIGINSFROM(ORIGINSFROM(masterivlist)),:)];
    validIRcoords = [validIRcoords;R(ones(nposs,1),:)];
    
end

%   #######################################################################
%   #######################################################################
%   #######################################################################
%   #######################################################################
%   #######################################################################

%   #############################################
%   #                                           #
%   #         S - edge - edge -spec - R cases   #
%   #                                           #
%   #         Postspec cases                    #
%   #                                           #
%   #############################################
%
% Possible edges for S-E-E-spec-spec-R are seen (at least partly) by the receiver.
%
% For some of the combos, we have the IR coords and the edge visibility
% from the single-diffraction run. For potential combos that were not
% included there, they are either invisible/obstructed or were not tested.
% We can figure out which ones were tested because they can be found in the
% POTENTIALISES under edge-spec-spec but not in the final valid list.

% The vector masterivlist will always refer to the original data vector
% i.e. PotentialIR, OriginsfromIR, reflorderIR etc
%
% First we pick out those indices where there was a single diffraction, but
% skip those with only diffraction (because we dealt with them already).
% Also, select only those where the diffraction is the first in the sequence
% of reflections.

for ii = 1:specorder-2
    
    if SHOWTEXT >= 3
        disp(['      Checking for ',JJ(ii,1:JJnumbofchars(ii)),' spec refl after the double edge diff'])    
    end

    % Select the combinations where the reflection order == ii+2
    % (which means ii specular reflections in addition to the diffraction)
    % and where the first two columns in POTENTIALISES contain edges.

    iv = find(REFLORDER(ivNdiff) == ii+2 & POTENTIALISES(ivNdiff,1)>nplanes & POTENTIALISES(ivNdiff,2)>nplanes);
    masterivlist = ivNdiff(iv);
    possibleedgepairs = double(POTENTIALISES(masterivlist,1:2)) - nplanes;
    possiblecombs = POTENTIALISES(masterivlist,3:2+ii);
    possibleweights = ISESVISIBILITY(masterivlist);
    
    % Compare with those combinations that were found OK
    % in the first-order diffraction step (EDB2diffISES), 
    % and that were input matrices to EDB2diff2ISES.
    % The index vector ivOK refers to these input matrices
    % (edgedifflistin,posspeclistin,validEDIRcoords,bigedgeweightlistin) 
    % npostspecs is a list that was calculated inside EDB2diff2ISES but
    % refers to the input matrices.
    
    ivOK = find(npostspecs==ii);
    if ~isempty(ivOK)
        patternOK = [edgedifflistin(ivOK) postspeclistin(ivOK,1:ii)];
    
        % Find out which ones, of all the possible first-order diffraction combos
        % in POTENTIALISES, that were indeed tested and found
        % invisible/obstructed in EDB2diffISES.
    
        ivallcombs = ivsinglediff(find( POTENTIALISES(ivsinglediff,1)>nplanes & REFLORDER(ivsinglediff) == ii+1));
        patternALL = [double(POTENTIALISES(ivallcombs,1))-nplanes double(POTENTIALISES(ivallcombs,2:1+ii))];
        if ~isempty(patternOK) & ~isempty(patternALL)
            patternNOTOK = setdiff(patternALL,patternOK,'rows');
        else
            if isempty(patternOK)
                patternNOTOK = patternALL;   
            else  % Then patternALL must be empty
                patternNOTOK = [];
            end
        end
    
        % Now, the patterns in patternNOTOK can be removed from
        % masterivlist.

        patterntocompare = [possibleedgepairs(:,2) possiblecombs(:,1:ii)];
    
        ivtocancel = find(ismember(patterntocompare,patternNOTOK,'rows'));
        masterivlist(ivtocancel) = [];
        possibleedgepairs(ivtocancel,:) = [];
        possiblecombs(ivtocancel,:) = [];
        possibleweights(ivtocancel,:) = [];
        patterntocompare(ivtocancel,:) = [];
    
        [ivcompletelyOK,ivreftoindata] = ismember(patterntocompare,patternOK,'rows');
        ivmustbechecked = find(ivcompletelyOK==0);
        ivcompletelyOK = find(ivcompletelyOK);    
    
    if ~isempty(ivmustbechecked)
        masterlisttocheckmore = masterivlist(ivmustbechecked);
        ntocheckmore = length(masterlisttocheckmore);
        
        %----------------------------------------------
        % Must carry out a visibility and obstruction check for the special
        % combinations here.
        % These combinations have a postspec-combination that hasn't been
        % encountered in the single diffraction cases, so no visibility
        % test has been made for these.
        
        lastedgenumbers = double(POTENTIALISES(masterlisttocheckmore,2))-nplanes;
        newIRcoords = R;
        reflplanesexpand = zeros(ntocheckmore*nedgesubs,ii);
        for jj = 1:ii
            reflplanes = POTENTIALISES(masterlisttocheckmore,3+ii-jj);
            reflplanesexpand(:,jj) = reshape(reflplanes(:,onesvec).',ntocheckmore*nedgesubs,1);
            newIRcoords = EDB2findis(newIRcoords,reflplanes,planeeqs,1,onesvec3);
            newIRcoordsexpand = reshape(repmat(newIRcoords.',nedgesubs,1),3,ntocheckmore*nedgesubs).';
            eval(['newIRcoords',JJ(jj,1:JJnumbofchars(jj)),' = newIRcoordsexpand;'])                    
        end
        [toedgecoords,edgeweightlist,edgenumberlist] = EDB2getedgepoints(edgestartcoords(lastedgenumbers,:),edgeendcoords(lastedgenumbers,:),edgelengthvec(lastedgenumbers,:),nedgesubs);
        tocoords = toedgecoords;
        lastedgenumbers = lastedgenumbers(:,onesvec);
        lastedgenumbers = reshape(lastedgenumbers.',ntocheckmore*nedgesubs,1);
        masterlisttocheckmore = masterlisttocheckmore(:,onesvec);
        masterlisttocheckmore = reshape(masterlisttocheckmore.',ntocheckmore*nedgesubs,1);
        
        ntocheckmore = length(masterlisttocheckmore);
        if SHOWTEXT >= 3
            disp(['         ',int2str(ntocheckmore),' special edge+edge+IR combinations to check'])    
        end

        for jj = ii:-1:1
            if length(masterlisttocheckmore) > 0
                eval(['fromcoords = newIRcoords',JJ(jj,1:JJnumbofchars(jj)),';']);
                if jj < ii
                    eval(['tocoords = reflpoints',JJ(jj+1,1:JJnumbofchars(jj+1)),';'])    
                end
                
                [hitplanes,reflpoints,edgehits,edgehitpoints,cornerhits,cornerhitpoints] = EDB2chkISvisible(fromcoords,tocoords,planeeqs(reflplanesexpand(:,jj),4),planenvecs(reflplanesexpand(:,jj),:),minvals(reflplanesexpand(:,jj),:),...
				    maxvals(reflplanesexpand(:,jj),:),planecorners(reflplanesexpand(:,jj),:),corners,ncornersperplanevec(reflplanesexpand(:,jj)));
                if ~isempty(edgehits) | ~isempty(cornerhits)
                    disp('WARNING! An edgehit or cornerhit occurred during the visibility test but this is not')
                    disp('         handled correctly yet.')
                end
                eval(['reflpoints',JJ(jj,1:JJnumbofchars(jj)),' = reflpoints;'])

                masterlisttocheckmore = masterlisttocheckmore(hitplanes);
                edgeweightlist = edgeweightlist(hitplanes);
                lastedgenumbers = lastedgenumbers(hitplanes);
                reflplanesexpand = reflplanesexpand(hitplanes,:);
                toedgecoords = toedgecoords(hitplanes,:);

                for kk = 1:ii
                    eval(['newIRcoords',JJ(kk,1:JJnumbofchars(kk)),' = newIRcoords',JJ(kk,1:JJnumbofchars(kk)),'(hitplanes,:);']);
                    if kk > jj
                        eval(['reflpoints',JJ(kk,1:JJnumbofchars(kk)),' = reflpoints',JJ(kk,1:JJnumbofchars(kk)),'(hitplanes,:);']);
                    end
                end
                ntocheckmore = length(masterlisttocheckmore);
            
            end
            
            if SHOWTEXT >= 3
                disp(['         ',int2str(ntocheckmore),' of the special edge+edge+IR combinations survived the visibility test in refl plane ',int2str(jj)])
            end
        end
        
        % Obstruction test of all the involved paths: R ->
        % reflplane1 -> reflplane2 -> ... -> last edge
        
        if obstructtestneeded & ntocheckmore > 0

            for jj = 1:ii+1
                if ntocheckmore > 0
                    if jj == 1
                        fromcoords = R;
                        startplanes = [];
                    else
                        eval(['fromcoords = reflpoints',JJ(jj-1,1:JJnumbofchars(jj-1)),';']) 
                        startplanes = reflplanesexpand(:,jj-1);
                    end
                    if jj == ii+1,                    
                        tocoords = toedgecoords;
                        endplanes = [planesatedge(lastedgenumbers,1) planesatedge(lastedgenumbers,2)];
                    else
                        eval(['tocoords = reflpoints',JJ(jj,1:JJnumbofchars(jj)),';'])    
                        endplanes = reflplanesexpand(:,jj);    
                    end
                    
                    [nonobstructedpaths,nobstructions] = EDB2checkobstrpaths(fromcoords,tocoords,startplanes,endplanes,canplaneobstruct,planeseesplane,...
                        planeeqs,planenvecs,minvals,maxvals,planecorners,corners,ncornersperplanevec,rearsideplane);

                    if nobstructions > 0
                        masterlisttocheckmore = masterlisttocheckmore(nonobstructedpaths);
                        edgeweightlist = edgeweightlist(nonobstructedpaths);
                        lastedgenumbers = lastedgenumbers(nonobstructedpaths);
                        reflplanesexpand = reflplanesexpand(nonobstructedpaths,:);
                        toedgecoords = toedgecoords(nonobstructedpaths,:);
                        for kk = 1:ii
                            eval(['reflpoints',JJ(kk,1:JJnumbofchars(kk)),' = reflpoints',JJ(kk,1:JJnumbofchars(kk)),'(nonobstructedpaths,:);']);
                            eval(['newIRcoords',JJ(kk,1:JJnumbofchars(kk)),' = newIRcoords',JJ(kk,1:JJnumbofchars(kk)),'(nonobstructedpaths,:);']);
                        end
                        
                    end
                    ntocheckmore = length(masterlisttocheckmore);
                    
                end
                
            end
            if SHOWTEXT >= 3
                disp(['         ',int2str(ntocheckmore),' of the special edge+edge+IR combinations survived the obstruction test'])
            end
            
        end
        
        % Add the found special combinations to the outdata list
    
        edgedifflist = [edgedifflist;double(POTENTIALISES(masterlisttocheckmore,1:2))-nplanes];
        prespeclist = [prespeclist;zerosvec1(ones(ntocheckmore,1),:)];
        midspeclist = [midspeclist;zerosvec1(ones(ntocheckmore,1),:)];
        postspeclist = [postspeclist;[reflplanesexpand zeros(ntocheckmore,specorder-2-ii)]];
        bigedgeweightlist = [bigedgeweightlist;[ ISESVISIBILITY(masterlisttocheckmore) edgeweightlist]];    
                
        eval(['validIRcoords = [validIRcoords;newIRcoords',JJ(ii,1:JJnumbofchars(ii)),'];']);
        validIScoords = [validIScoords;S(ones(ntocheckmore,1),:)];
        
    end
    
    masterivlist = masterivlist(ivcompletelyOK);    
    possibleedgepairs = possibleedgepairs(ivcompletelyOK,:);
    possiblecombs = possiblecombs(ivcompletelyOK,:);
    possibleweights = possibleweights(ivcompletelyOK,:);
        
    nposs = length(ivcompletelyOK);

    if SHOWTEXT >= 3
		disp(['         ',int2str(nposs),' Edge+edge+IR segments (non-special combinations) survived the obstruction test'])
    end
    
    % Add the found "standard" combinations to the outdata list
    
    edgedifflist = [edgedifflist;possibleedgepairs];
    prespeclist = [prespeclist;zerosvec1(ones(nposs,1),:)];
    midspeclist = [midspeclist;zerosvec1(ones(nposs,1),:)];
    postspeclist = [postspeclist;[possiblecombs zeros(nposs,specorder-2-ii)]];
    bigedgeweightlist = [bigedgeweightlist;[possibleweights bigedgeweightlistin(ivOK(ivreftoindata(ivcompletelyOK)))]];    
    validIRcoords = [validIRcoords;validEDIRcoords(ivOK(ivreftoindata(ivcompletelyOK)),:)];
    validIScoords = [validIScoords;S(ones(nposs,1),:)];
    
end
end


%   #######################################################################
%   #######################################################################
%   #######################################################################
%   #######################################################################
%   #######################################################################

%   ####################################################
%   #                                                  #
%   #         S - spec - edge - edge - spec - R cases  #
%   #                                                  #
%   #         Pre and postspec cases                   #
%   #                                                  #
%   ####################################################


% The ii- and jj-loops will go through all combinations of specular
% reflection before and after the double diffraction.
%
% Unlike before, we will look in the list of already found combinations
% (edgedifflist etc)

for ii = 1:specorder-3

    for jj = 1:specorder-ii-2
    
        if SHOWTEXT >= 3
            disp(['      Checking for ',JJ(ii,1:JJnumbofchars(ii)),' spec refl before and ',JJ(jj,1:JJnumbofchars(jj)),' spec refl after the double edge diff'])    
        end
	
        % Select the combinations where the reflection order == ii+jj+2
        % (which means ii+jj specular reflections in addition to the
        % diffraction), and where columns ii+1 and ii+2 of POTENTIALISES
        % contain edges.
	
        iv = find(REFLORDER(ivNdiff) == ii+jj+2 & POTENTIALISES(ivNdiff,ii+1)>nplanes & POTENTIALISES(ivNdiff,ii+2)>nplanes);
        masterivlist = ivNdiff(iv);
        possibleedgepairs = double(POTENTIALISES(masterivlist,ii+1:ii+2)) - nplanes;
        possibleprespecs  = POTENTIALISES(masterivlist,1:ii);
        possiblepostspecs = POTENTIALISES(masterivlist,ii+3:ii+2+jj);
        possibleweights = ISESVISIBILITY(masterivlist);

        % Compare with those that have already been found OK
        ivOK = find(npostspecs==jj);
        if ~isempty(ivOK)
            patternOK = [edgedifflistin(ivOK) postspeclistin(ivOK,1:jj)];
        else
            patternOK = [];    
        end
        
        % Find out which ones have been checked and found invisible/obstructed
        ivallcombs = ivsinglediff(find( POTENTIALISES(ivsinglediff,1)>nplanes & REFLORDER(ivsinglediff) == jj+1));
        patternALL = [double(POTENTIALISES(ivallcombs,1))-nplanes double(POTENTIALISES(ivallcombs,2:1+jj))];
        if ~isempty(patternOK) & ~isempty(patternALL)
            patternNOTOK = setdiff(patternALL,patternOK,'rows');
        else
            if isempty(patternOK)
                patternNOTOK = patternALL;   
            else  % Then patternALL must be empty
                patternNOTOK = [];
            end
        end
	        
        patterntocompare = [possibleedgepairs(:,2) possiblepostspecs(:,1:jj)];
        
       ivtocancel = find(ismember(patterntocompare,patternNOTOK,'rows'));
        masterivlist(ivtocancel) = [];
        possibleedgepairs(ivtocancel,:) = [];
        possibleprespecs(ivtocancel,:) = [];
        possiblepostspecs(ivtocancel,:) = [];
        possibleweights(ivtocancel,:) = [];
        patterntocompare(ivtocancel,:) = [];
	
       [ivcompletelyOK,ivreftoindata] = ismember(patterntocompare,patternOK,'rows');
        ivmustbechecked = find(ivcompletelyOK==0);
        ivcompletelyOK = find(ivcompletelyOK);
        if ~isempty(ivmustbechecked)
            masterlisttocheckmore = masterivlist(ivmustbechecked);
            ntocheckmore = length(masterlisttocheckmore);
        
            %----------------------------------------------
            % Must carry out a visibility and obstruction check for the special
            % combinations here.
            %
            % NB! toedgecoords are the coordinates of the last edge in the sequence.
            %     This name is because for post-specular reflections, the propagation
            %	  is viewed from the receiver towards the last edge!

            lastedgenumbers = double(POTENTIALISES(masterlisttocheckmore,ii+2))-nplanes;
            newIRcoords = R;
            reflplanesexpand = zeros(ntocheckmore*nedgesubs,ii);
            for kk = 1:jj
                reflplanes = POTENTIALISES(masterlisttocheckmore,3+ii+jj-kk);
                reflplanesexpand(:,kk) = reshape(reflplanes(:,onesvec).',ntocheckmore*nedgesubs,1);
                newIRcoords = EDB2findis(newIRcoords,reflplanes,planeeqs,1,onesvec3);
                newIRcoordsexpand = reshape(repmat(newIRcoords.',nedgesubs,1),3,ntocheckmore*nedgesubs).';
                eval(['newIRcoords',JJ(kk,1:JJnumbofchars(kk)),' = newIRcoordsexpand;'])                    
            end
            [toedgecoords,edgeweightlist,edgenumberlist] = EDB2getedgepoints(edgestartcoords(lastedgenumbers,:),edgeendcoords(lastedgenumbers,:),edgelengthvec(lastedgenumbers,:),nedgesubs);
            tocoords = toedgecoords;
            lastedgenumbers = lastedgenumbers(:,onesvec);
            lastedgenumbers = reshape(lastedgenumbers.',ntocheckmore*nedgesubs,1);
            masterlisttocheckmore = masterlisttocheckmore(:,onesvec);
            masterlisttocheckmore = reshape(masterlisttocheckmore.',ntocheckmore*nedgesubs,1);
            
            ntocheckmore = length(masterlisttocheckmore);
            if SHOWTEXT >= 3
                disp(['         ',int2str(ntocheckmore),' special IS+edge+edge+IR combinations to check'])    
            end

            for kk = jj:-1:1
                if length(masterlisttocheckmore) > 0
                    eval(['fromcoords = newIRcoords',JJ(kk,1:JJnumbofchars(kk)),';']);
                    if kk < jj
                        eval(['tocoords = reflpoints',JJ(kk+1,1:JJnumbofchars(kk+1)),';'])    
                    end
                    
                    [hitplanes,reflpoints,edgehits,edgehitpoints,cornerhits,cornerhitpoints]  = EDB2chkISvisible(fromcoords,tocoords,planeeqs(reflplanesexpand(:,kk),4),planenvecs(reflplanesexpand(:,kk),:),minvals(reflplanesexpand(:,kk),:),...
					    maxvals(reflplanesexpand(:,kk),:),planecorners(reflplanesexpand(:,kk),:),corners,ncornersperplanevec(reflplanesexpand(:,kk)));
                    if ~isempty(edgehits) | ~isempty(cornerhits)
                        disp('WARNING! An edgehit or cornerhit occurred during the visibility test but this is not')
                        disp('         handled correctly yet.')
                    end
                    eval(['reflpoints',JJ(kk,1:JJnumbofchars(kk)),' = reflpoints;'])

                    masterlisttocheckmore = masterlisttocheckmore(hitplanes);
                    edgeweightlist = edgeweightlist(hitplanes);
                    lastedgenumbers = lastedgenumbers(hitplanes);
                    reflplanesexpand = reflplanesexpand(hitplanes,:);
                    toedgecoords = toedgecoords(hitplanes,:);
	
                    for ll = 1:jj
                        eval(['newIRcoords',JJ(ll,1:JJnumbofchars(ll)),' = newIRcoords',JJ(ll,1:JJnumbofchars(ll)),'(hitplanes,:);']);
                        if ll > kk
                            eval(['reflpoints',JJ(ll,1:JJnumbofchars(ll)),' = reflpoints',JJ(ll,1:JJnumbofchars(ll)),'(hitplanes,:);']);
                        end
                    end
                    ntocheckmore = length(masterlisttocheckmore);
            
                end
                
                if SHOWTEXT >= 3
                    disp(['         ',int2str(ntocheckmore),' of the special IS+edge+edge+IR combinations survived the visibility test in refl plane ',int2str(jj)])
                end
                
            end
            
            % Obstruction test of all the involved paths: R ->
            % reflplane1 -> reflplane2 -> ... -> last edge
            
            if obstructtestneeded & ntocheckmore > 0
	
                for kk = 1:jj+1
                    if ntocheckmore > 0
                        if kk == 1
                            fromcoords = R;
                            startplanes = [];
                        else
                            eval(['fromcoords = reflpoints',JJ(kk-1,1:JJnumbofchars(kk-1)),';']) 
                            startplanes = reflplanesexpand(:,kk-1);
                        end
                        if kk == jj+1,                    
                            tocoords = toedgecoords;
                            endplanes = [planesatedge(lastedgenumbers,1) planesatedge(lastedgenumbers,2)];
                        else
                            eval(['tocoords = reflpoints',JJ(kk,1:JJnumbofchars(kk)),';'])    
                            endplanes = reflplanesexpand(:,kk);    
                        end
                        
                        [nonobstructedpaths,nobstructions] = EDB2checkobstrpaths(fromcoords,tocoords,startplanes,endplanes,canplaneobstruct,planeseesplane,...
                            planeeqs,planenvecs,minvals,maxvals,planecorners,corners,ncornersperplanevec,rearsideplane);
	
                        if nobstructions > 0
                            masterlisttocheckmore = masterlisttocheckmore(nonobstructedpaths);
                            edgeweightlist = edgeweightlist(nonobstructedpaths);
                            lastedgenumbers = lastedgenumbers(nonobstructedpaths);
                            reflplanesexpand = reflplanesexpand(nonobstructedpaths,:);
                            toedgecoords = toedgecoords(nonobstructedpaths,:);
                            for ll = 1:jj
                                eval(['reflpoints',JJ(ll,1:JJnumbofchars(ll)),' = reflpoints',JJ(ll,1:JJnumbofchars(ll)),'(nonobstructedpaths,:);']);
                                eval(['newIRcoords',JJ(ll,1:JJnumbofchars(ll)),' = newIRcoords',JJ(ll,1:JJnumbofchars(ll)),'(nonobstructedpaths,:);']);
                            end
                            
                        end
                        ntocheckmore = length(masterlisttocheckmore);
                        
                    end
                    
                end
                if SHOWTEXT >= 3
                    disp(['         ',int2str(ntocheckmore),' of the special IS+edge+edge+IR combinations survived the obstruction test'])
                end
                
            end
            
            % Add the found special combinations to the outdata list
        
            edgedifflist = [edgedifflist;double(POTENTIALISES(masterlisttocheckmore,ii+1:ii+2))-nplanes];
            prespeclist = [prespeclist;[double(POTENTIALISES(masterlisttocheckmore,1:ii)) zeros(ntocheckmore,specorder-2-ii)]];
            midspeclist = [midspeclist;zerosvec1(ones(ntocheckmore,1),:)];
            postspeclist = [postspeclist;[reflplanesexpand zeros(ntocheckmore,specorder-2-ii)]];
            bigedgeweightlist = [bigedgeweightlist;[ ISESVISIBILITY(masterlisttocheckmore) edgeweightlist]];    
            eval(['validIRcoords = [validIRcoords;newIRcoords',JJ(jj,1:JJnumbofchars(jj)),'];']);
            % NB!! In the same way as earlier, we must a recursive reference
            % method to find the image source of the last specular reflection.
            ivref = ORIGINSFROM(ORIGINSFROM(ORIGINSFROM(masterlisttocheckmore)));
            for kk = 2:jj
                ivref = ORIGINSFROM(ivref);
            end
            validIScoords = [validIScoords;ISCOORDS(ivref,:)];
        end
        
        masterivlist = masterivlist(ivcompletelyOK);    
        possibleedgepairs = possibleedgepairs(ivcompletelyOK,:);
        possibleprespecs = possibleprespecs(ivcompletelyOK,:);
        possiblepostspecs = possiblepostspecs(ivcompletelyOK,:);
        possibleweights = possibleweights(ivcompletelyOK,:);
            
        nposs = length(ivcompletelyOK);
	
        if SHOWTEXT >= 3
			disp(['         ',int2str(nposs),' IS+Edge+edge+IR segments (non-special) survived the obstruction test'])
        end

        % Add the found "standard" combinations to the outdata list
    
        edgedifflist = [edgedifflist;possibleedgepairs];
        prespeclist = [prespeclist;  [possibleprespecs  zeros(nposs,specorder-2-ii)]];
        midspeclist = [midspeclist;zerosvec1(ones(nposs,1),:)];
        postspeclist = [postspeclist;[possiblepostspecs zeros(nposs,specorder-2-jj)]];
        bigedgeweightlist = [bigedgeweightlist;[possibleweights bigedgeweightlistin(ivOK(ivreftoindata(ivcompletelyOK)))]];    
        validIRcoords = [validIRcoords;validEDIRcoords(ivOK(ivreftoindata(ivcompletelyOK)),:)];
        % NB!! In the same way as earlier, we must a recursive reference
        % method to find the image source of the last specular reflection.
        ivref = ORIGINSFROM(ORIGINSFROM(ORIGINSFROM(masterivlist)));
        for kk = 2:jj
            ivref = ORIGINSFROM(ivref);
        end
        validIScoords = [validIScoords;ISCOORDS(ivref,:)];

    end

end

%   #######################################################################
%   #######################################################################
%   #######################################################################
%   #######################################################################
%   #######################################################################

%   #######################################################################
%   #
%   #   Midspec cases, with prespecs and postspecs
%   #
%   #  S - spec - edge - spec - edge - R
%   #
%   #######################################################################

for iipre = 0:specorder-3
    for jjmid = 1:specorder-2-iipre
        for kkpost = 0:specorder-iipre-jjmid-2

            if SHOWTEXT >= 3
                if iipre > 0
                    JJpre = JJ(iipre,1:JJnumbofchars(iipre));
                else
                    JJpre = '0';
                end
                if kkpost > 0
                    JJpost = JJ(kkpost,1:JJnumbofchars(kkpost));
                else
                    JJpost = '0';                    
                end
                disp(['      Checking for ',JJpre,' spec refl before, ',JJ(jjmid,1:JJnumbofchars(jjmid)),' spec refl between the double edge diff, and'])
                disp(['      ',JJpost,' spec refl after the double edge diff.'])
            end
            iv = find(REFLORDER(ivNdiff) == iipre+jjmid+kkpost+2 & POTENTIALISES(ivNdiff,iipre+1)>nplanes & POTENTIALISES(ivNdiff,iipre+jjmid+2)>nplanes);
            masterivlist = ivNdiff(iv);
            possibleedgepairs = double(POTENTIALISES(masterivlist,[iipre+1 iipre+jjmid+2])) - nplanes;
            nfast = 0;

            if kkpost == 0
                possiblepostspecs = [];
                
                % Keep only combinations for which the receiver can see the edge
                ivnotvisiblefromr = find(vispartedgesfromR(possibleedgepairs(:,2))==0);
                if ~isempty(ivnotvisiblefromr)
                    masterivlist(ivnotvisiblefromr) = [];
                    possibleedgepairs(ivnotvisiblefromr,:) = [];
                end        
            else
                possiblepostspecs = POTENTIALISES(masterivlist,iipre+jjmid+3:iipre+jjmid+2+kkpost);

                % Compare with those that have already been found OK
                ivOK = find(npostspecs==kkpost);
                if ~isempty(ivOK)
                    patternOK = [edgedifflistin(ivOK) postspeclistin(ivOK,1:kkpost)];
                else
                    patternOK = [];    
                end
                % Find out which ones have been checked and found invisible/obstructed
                ivallcombs = ivsinglediff(find( POTENTIALISES(ivsinglediff,1)>nplanes & REFLORDER(ivsinglediff) == kkpost+1));
                patternALL = [double(POTENTIALISES(ivallcombs,1))-nplanes double(POTENTIALISES(ivallcombs,2:1+kkpost))];

                if ~isempty(patternOK) & ~isempty(patternALL)
                    patternNOTOK = setdiff(patternALL,patternOK,'rows');
                else
                    if isempty(patternOK)
                        patternNOTOK = patternALL;   
                    else  % Then patternALL must be empty
                        patternNOTOK = [];
                    end
                end

                patterntocompare = [possibleedgepairs(:,2) possiblepostspecs(:,1:kkpost)];
        
                if ~isempty(patternNOTOK)
                   ivtocancel = find(ismember(patterntocompare,patternNOTOK,'rows'));
                    masterivlist(ivtocancel) = [];
                    possibleedgepairs(ivtocancel,:) = [];
                    possiblepostspecs(ivtocancel,:) = [];
                    patterntocompare(ivtocancel,:) = [];
                end
              
               [ivcompletelyOK,ivreftoindata] = ismember(patterntocompare,patternOK,'rows');
                ivmustbechecked = find(ivcompletelyOK==0);
                ivcompletelyOK = find(ivcompletelyOK);
                
                if ~isempty(ivmustbechecked) & SHOWTEXT > 0
                    disp('WARNING!! For midspec and postspec case, all checks have not been implemented yet!!')
                end

                masterivlist = masterivlist(ivcompletelyOK);    
                possibleedgepairs = possibleedgepairs(ivcompletelyOK,:);
                possiblepostspecs = possiblepostspecs(ivcompletelyOK,:);
                    
                nposs = length(ivcompletelyOK);
			                
            end

            if iipre > 0
                possibleprespecs = POTENTIALISES(masterivlist,1:iipre);
            else
                possibleprespecs = [];
            end
            
            % NB!! possiblemidspecs is numbered in reverse order
            % since we view the propagation by starting to mirror the last edge
            % and move towards the first edge. 
            
            possiblemidspecs = POTENTIALISES(masterivlist,iipre+1+jjmid:-1:iipre+2);
                         
            if kkpost > 0
                edgeweightlist = [ISESVISIBILITY(masterivlist) bigedgeweightlistin(ivOK(ivreftoindata(ivcompletelyOK)))];
            else
                edgeweightlist = [ISESVISIBILITY(masterivlist) vispartedgesfromR(possibleedgepairs(:,2))];
            end
            
            % Expand the various lists and matrices to represent the
            % sub-divided edges.

            nposs = length(masterivlist);
            if nposs > 0
                if iipre == 1
                    possibleprespecs = reshape(possibleprespecs(:,onesvec).',nposs*nedgesubs,1);
                elseif iipre > 1
                    possibleprespecs = reshape(repmat(possibleprespecs.',nedgesubs,1),iipre,nposs*nedgesubs).';
                end
                if jjmid == 1
                    possiblemidspecs = reshape(possiblemidspecs(:,onesvec).',nposs*nedgesubs,1);
                elseif jjmid > 1
                    possiblemidspecs = reshape(repmat(possiblemidspecs.',nedgesubs,1),jjmid,nposs*nedgesubs).';
                end
                if kkpost == 1
                    possiblepostspecs = reshape(possiblepostspecs(:,onesvec).',nposs*nedgesubs,1);
                elseif kkpost > 1
                    possiblepostspecs = reshape(repmat(possiblepostspecs.',nedgesubs,1),kkpost,nposs*nedgesubs).';
                end
	
                expandedmasterivlist = reshape(masterivlist(:,onesvec).',nposs*nedgesubs,1);       
                if kkpost > 0
                    expandedivcompletelyOK = reshape(ivcompletelyOK(:,onesvec).',nposs*nedgesubs,1);
                end
                
                edgeweightlist = reshape(repmat(edgeweightlist.',[nedgesubs 1]),2,nposs*nedgesubs).';
                
                for ll = 1:nedgesubs
                    edgeweightlist(ll:nedgesubs:end,1) = double(bitget(edgeweightlist(ll:nedgesubs:end,1),ll))*bitmultvec(ll);
                    edgeweightlist(ll:nedgesubs:end,2) = double(bitget(edgeweightlist(ll:nedgesubs:end,2),ll))*bitmultvec(ll);
                end
                
                %----------------------------------------------
                % Must carry out a visibility and obstruction check for the 
                % edge-spec-edge paths
				%
                % NB! toedgecoords are the coordinates of the first edge in the sequence
                %     and fromedgecoords refers to the last edge, after the mid-specular
                %	  reflections.
                %     This name is because for mid-specular reflections, the propagation
                %	  is viewed from the last edge towards the first edge!
                
                [toedgecoords,firstedgeweightlist,slask]     = EDB2getedgepoints(edgestartcoords(possibleedgepairs(:,2),:),edgeendcoords(possibleedgepairs(:,2),:),edgelengthvec(possibleedgepairs(:,1),:),nedgesubs);
                tocoords = toedgecoords;
                [fromedgecoords,lastedgeweightlist,slask] = EDB2getedgepoints(edgestartcoords(possibleedgepairs(:,1),:),edgeendcoords(possibleedgepairs(:,1),:),edgelengthvec(possibleedgepairs(:,2),:),nedgesubs);
                fromcoords = fromedgecoords;
	
                possibleedgepairs = reshape(repmat(possibleedgepairs.',nedgesubs,1),2,nposs*nedgesubs).';
	
                edgeimagecoords = fromedgecoords;
                for ll = 1:jjmid
                    edgeimagecoords = EDB2findis(edgeimagecoords,possiblemidspecs(:,ll),planeeqs,size(fromedgecoords,1),onesvec3);
                    eval(['bigedgeimagecoords',JJ(ll,1:JJnumbofchars(ll)),' = edgeimagecoords;'])                    
                end
	
                % Some cases do not need to be checked, when jjmid = 2: the
                % cateye cases. For these, we will have doubles (both, e.g.
                % 3-7 and 7-3) and one should be tossed then. The non-doubles
                % can be kept in a "fast lane" so that visibility isn't
                % checked, but obstruction is.
	
                if jjmid == 2
                    specpattern = double(possiblemidspecs);
                    ivreftomatrix = specpattern(:,1) + ( specpattern(:,2)-1)*nplanes;
                    ivcateyes = find( cateyeplanecombs(ivreftomatrix) );
                    
                    if ~isempty(ivcateyes),                    
                        specpattern = specpattern(ivcateyes,:);
                        fliporder = specpattern(:,2)<specpattern(:,1);
                        ivfliporder = find(fliporder);
                        specpattern(ivfliporder,:) = specpattern(ivfliporder,[2 1]);                    
                        [uniquepatterns,ivec,jvec] = unique(specpattern,'rows');
	
                        countcases = histc(jvec,[1:max(jvec)]);
                        ivtossone = find(fliporder & countcases(jvec)==nedgesubs*2);
	
                        ivfastlane = ivcateyes;
                        ivfastlane(ivtossone) = [];
                        if ~isempty(ivfastlane)
                            expandedmasterivlistfast = expandedmasterivlist(ivfastlane);
                            possibleedgepairsfast = possibleedgepairs(ivfastlane,:); 
                            fromedgecoordsfast = fromedgecoords(ivfastlane,:);
                            toedgecoordsfast = toedgecoords(ivfastlane,:);
                            if ~isempty(possibleprespecs)
                                possibleprespecsfast = possibleprespecs(ivfastlane,:);
                            end
                            possiblemidspecsfast = possiblemidspecs(ivfastlane,:);
                            if ~isempty(possiblepostspecs)
                                possiblepostspecsfast = possiblepostspecs(ivfastlane,:);
                            end
                            edgeweightlistfast = edgeweightlist(ivfastlane,:);
                            for mm = 1:jjmid
                                eval(['bigedgeimagecoords',JJ(mm,1:JJnumbofchars(mm)),'fast = bigedgeimagecoords',JJ(mm,1:JJnumbofchars(mm)),'(ivfastlane,:);']);
                            end
                            reflpoints2fast = 0.5*(bigedgeimagecoords2fast-fromedgecoordsfast) + fromedgecoordsfast;
                            reflpoints1fast = reflpoints2fast;
                            nfast = length(edgeweightlistfast);               
                        end
                        
                        if ~isempty(ivtossone)
                            ivtossone = ivcateyes;
                            expandedmasterivlist((ivtossone)) = [];
                            edgeweightlist((ivtossone),:) = [];
                            possibleedgepairs((ivtossone),:) = [];
                            if ~isempty(possibleprespecs)
                                possibleprespecs((ivtossone),:) = [];
                            end
                            possiblemidspecs((ivtossone),:) = [];
                            if ~isempty(possiblepostspecs)
                                possiblepostspecs((ivtossone),:) = [];
                            end
                            fromcoords((ivtossone),:) = [];
                            tocoords((ivtossone),:) = [];
                            for mm = 1:jjmid
                                eval(['bigedgeimagecoords',JJ(mm,1:JJnumbofchars(mm)),'((ivtossone),:) = [];']);
                            end
                            nposs = length(expandedmasterivlist);
                        end
                        
                    end
                    
                end

                nposs = length(expandedmasterivlist);
            end
        
            if SHOWTEXT >= 3
                if jjmid ~= 2
     		        disp(['         ',int2str(nposs),' edge+spec+edge combos found initially:'])
                else
                    disp(['         ',int2str(nposs),' edge+spec+edge combos found initially, + ',int2str(nfast),' cateye combos'])
                end
            end
            
            %--------------------------------------------------------------
            % Check the visibility through all the reflection planes
            %
            
            for ll = 1:jjmid
                eval(['reflpoints',JJ(ll,1:JJnumbofchars(ll)),' = [];'])
            end
            for ll = jjmid:-1:1
                if nposs > 0
                    eval(['fromcoords = bigedgeimagecoords',JJ(ll,1:JJnumbofchars(ll)),';'])
                    if ll < jjmid
                       eval(['tocoords = reflpoints',JJ(ll+1,1:JJnumbofchars(ll+1)),';'])
                    end

                    [hitplanes,reflpoints,edgehits,edgehitpoints,cornerhits,cornerhitpoints] = EDB2chkISvisible(fromcoords,tocoords,planeeqs(possiblemidspecs(:,ll),4),planenvecs(possiblemidspecs(:,ll),:),minvals(possiblemidspecs(:,ll),:),...
     			    maxvals(possiblemidspecs(:,ll),:),planecorners(possiblemidspecs(:,ll),:),corners,ncornersperplanevec(possiblemidspecs(:,ll)));

                    % Make a special treatment for the cases with the
                    % specular reflection point right on an edge since some
                    % of these are special cases:
                    % "edgeplaneperptoplane1" indicates that edge-plane-edge
                    %    travels along the edge's plane and is reflected at a
                    %    90 degree corner (which is an inactive edge).
                    %    These are treated as good hits.
                    % "edgeplaneperptoplane2" indicates that edge-plane-edge
                    %    has a specular reflection right at a flat edge
                    %    between two coplanar planes.
                    %    These come in pairs; one half-hit in the two
                    %    coplanar planes. Only one in each pair should be
                    %    kept.
                    
                    if jjmid == 1 & ~isempty(edgehits)
                        edge1 = possibleedgepairs(edgehits,1);
                        edge2 = possibleedgepairs(edgehits,2);
                        midspec = possiblemidspecs(edgehits,1);
                        ivreftomatrix1 = double(midspec) + double(edge1-1)*nplanes;
                        ivreftomatrix2 = double(midspec) + double(edge2-1)*nplanes;

                        specialhit = edgeplaneperptoplane1(ivreftomatrix1).*edgeplaneperptoplane1(ivreftomatrix1);
                        ivspecial = find(specialhit);
                        if ~isempty(ivspecial)
                            hitplanes = [hitplanes;edgehits(ivspecial)];
                            reflpoints = [reflpoints;edgehitpoints(ivspecial,:)];
                            edgehits(ivspecial) = [];
                            edgehitpoints(ivspecial,:) = [];
                            ivreftomatrix1(ivspecial) = [];
                            ivreftomatrix2(ivspecial) = [];
                        end
                        
                        specialhit = edgeplaneperptoplane2(ivreftomatrix1).*edgeplaneperptoplane2(ivreftomatrix1);
                        ivspecial = find(specialhit);
                        if ~isempty(ivspecial)
                            patternlist = double([possibleedgepairs(edgehits(ivspecial),1) possiblemidspecs(edgehits(ivspecial),1) possibleedgepairs(edgehits(ivspecial),1)]);
                            [uniquepatterns,ivec,jvec] = unique(patternlist,'rows');
                            keeppattern = zeros(size(ivec));
                            for mm = 1:length(ivec)
                                ivreftomatrix = uniquepatterns(mm,2) + (uniquepatterns(mm+1:end,2)-1)*nplanes;
                                coplanarindicator = coplanarsviaflatedge(ivreftomatrix);
                                ivcoplanars = find(uniquepatterns(mm+1:end,1)==uniquepatterns(mm,1) & uniquepatterns(mm+1:end,3)==uniquepatterns(mm,3) & coplanarindicator);
                                if ~isempty(ivcoplanars)
                                    keeppattern(mm) = 1;    
                                end
                            end
                            ivkeepers = find(keeppattern(jvec));
                            hitplanes = [hitplanes;edgehits(ivspecial(ivkeepers))];
                            reflpoints = [reflpoints;edgehitpoints(ivspecial(ivkeepers),:)];

                        end                        
                    end

                    eval(['reflpoints',JJ(ll,1:JJnumbofchars(ll)),' = reflpoints;'])
    
                    expandedmasterivlist = expandedmasterivlist(hitplanes);
                    edgeweightlist = edgeweightlist(hitplanes,:);
                    possibleedgepairs = possibleedgepairs(hitplanes,:);
                    if ~isempty(possibleprespecs)
                        possibleprespecs = possibleprespecs(hitplanes,:);
                    end
                    possiblemidspecs = possiblemidspecs(hitplanes,:);
                    if ~isempty(possiblepostspecs)
                        possiblepostspecs = possiblepostspecs(hitplanes,:);
                    end
                    fromedgecoords = fromedgecoords(hitplanes,:);
                    toedgecoords = toedgecoords(hitplanes,:);
                    if kkpost > 0
                        expandedivcompletelyOK = expandedivcompletelyOK(hitplanes);
                    end
                    
                    for mm = 1:jjmid
                        eval(['bigedgeimagecoords',JJ(mm,1:JJnumbofchars(mm)),' = bigedgeimagecoords',JJ(mm,1:JJnumbofchars(mm)),'(hitplanes,:);']);
                        if mm > ll
                            eval(['reflpoints',JJ(mm,1:JJnumbofchars(mm)),' = reflpoints',JJ(mm,1:JJnumbofchars(mm)),'(hitplanes,:);']);
                        end
                    end
                    nposs = length(expandedmasterivlist);
	
                    if SHOWTEXT >= 3
                        if jjmid ~= 2
     		               disp(['         ',int2str(nposs),' edge+spec+edge combos survived the visibility test in reflection plane ',int2str(ll)])
                        else
                           disp(['         ',int2str(nposs),' edge+spec+edge combos survived the visibility test in reflection plane ',int2str(ll),' + ',int2str(nfast),' cateye combos'])
                        end
                    end
                end                
            end     
            
            %--------------------------------------------------------------
            % Check for obstructions for all the paths, starting from edge number 2
            % towards edge number 1.
            %
            % Reinsert the "fast lane" cases, i.e., the cateye reflections.
            
            if jjmid == 2
                if nfast > 0,                    
                    expandedmasterivlist = [expandedmasterivlist;expandedmasterivlistfast];
                    edgeweightlist = [edgeweightlist;edgeweightlistfast];
                    possibleedgepairs = [possibleedgepairs;possibleedgepairsfast];
                    if ~isempty(possibleprespecs)
                        possibleprespecs = [possibleprespecs;possibleprespecsfast];
                    end
                    possiblemidspecs = [possiblemidspecs;possiblemidspecsfast];
                    if ~isempty(possiblepostspecs)
                        possiblepostspecs = [possiblepostspecs;possiblepostspecsfast];
                    end
                    fromedgecoords = [fromedgecoords;fromedgecoordsfast];
                    toedgecoords   = [toedgecoords;toedgecoordsfast];
                    for mm = 1:jjmid
                        eval(['bigedgeimagecoords',JJ(mm,1:JJnumbofchars(mm)),' = [bigedgeimagecoords',JJ(mm,1:JJnumbofchars(mm)),';bigedgeimagecoords',JJ(mm,1:JJnumbofchars(mm)),'fast];']);
                        eval(['reflpoints',JJ(mm,1:JJnumbofchars(mm)),        ' = [reflpoints',JJ(mm,1:JJnumbofchars(mm)),        ';reflpoints',JJ(mm,1:JJnumbofchars(mm)),'fast];']);
                    end
                    nposs = length(expandedmasterivlist);
                end
            end

            if nposs > 0 & obstructtestneeded

                for ll = 1:jjmid+1
                    if nposs > 0
                        if ll == 1
                            fromcoords = fromedgecoords;
                            startplanes = [planesatedge(possibleedgepairs(:,2),1) planesatedge(possibleedgepairs(:,2),2)];
                        else
                            eval(['fromcoords = reflpoints',JJ(ll-1,1:JJnumbofchars(ll-1)),';']) 
                            startplanes = possiblemidspecs(:,ll-1);
                        end
                        if ll == jjmid+1,                    
                            tocoords = toedgecoords;
                            endplanes = [planesatedge(possibleedgepairs(:,1),1) planesatedge(possibleedgepairs(:,1),2)];
                        else
                            eval(['tocoords = reflpoints',JJ(ll,1:JJnumbofchars(ll)),';'])    
                            endplanes = possiblemidspecs(:,ll);    
                        end
                        [nonobstructedpaths,nobstructions,edgehits,cornerhits] = EDB2checkobstrpaths(fromcoords,tocoords,startplanes,endplanes,canplaneobstruct,planeseesplane,...
                            planeeqs,planenvecs,minvals,maxvals,planecorners,corners,ncornersperplanevec,rearsideplane);

                        if ~isempty(edgehits)
                            nonobstructedpaths = setdiff(nonobstructedpaths,edgehits);
                            nobstructions = nposs - length(nonobstructedpaths);
                        end
    
                        if nobstructions > 0
                            expandedmasterivlist = expandedmasterivlist(nonobstructedpaths);

							edgeweightlist = edgeweightlist(nonobstructedpaths,:);
							possibleedgepairs = possibleedgepairs(nonobstructedpaths,:);
							if ~isempty(possibleprespecs)
								possibleprespecs = possibleprespecs(nonobstructedpaths,:);
							end
							possiblemidspecs = possiblemidspecs(nonobstructedpaths,:);
							if ~isempty(possiblepostspecs)
								possiblepostspecs = possiblepostspecs(nonobstructedpaths,:);
							end
							fromedgecoords = fromedgecoords(nonobstructedpaths,:);
							toedgecoords = toedgecoords(nonobstructedpaths,:);
                            if kkpost > 0
                                expandedivcompletelyOK = expandedivcompletelyOK(nonobstructedpaths,:);
                            end
                            
							for mm = 1:jjmid
								eval(['bigedgeimagecoords',JJ(mm,1:JJnumbofchars(mm)),' = bigedgeimagecoords',JJ(mm,1:JJnumbofchars(mm)),'(nonobstructedpaths,:);']);
								if mm > ll
									eval(['reflpoints',JJ(mm,1:JJnumbofchars(mm)),' = reflpoints',JJ(mm,1:JJnumbofchars(mm)),'(nonobstructedpaths,:);']);
								end
							end
		                            
                        end
                        nposs = length(expandedmasterivlist);
                        
                    end
                end
                
                if SHOWTEXT >= 3
                    if jjmid == 2
                        disp(['         ',int2str(nposs),' edge+spec+edge combos survived the obstruction test (including cateye cases)'])
                    else
                        disp(['         ',int2str(nposs),' edge+spec+edge combos survived the obstruction test'])
                    end
                end
                
            end

            if nposs > 0
                edgedifflist = [edgedifflist;possibleedgepairs];
                if iipre == 0
                    possibleprespecs = zeros(nposs,1);    
                end
                if specorder <= 4
                    if specorder == 3
                        prespeclist = [prespeclist;[possibleprespecs]];
                    else
                       prespeclist = [prespeclist;[possibleprespecs zeros(nposs,1)]];
                   end
                else
                   [n1,n2] = size(prespeclist);
                   [n3,n4] = size(possibleprespecs);
                   if n1 > 0
% Error found 20050202 PS                       
% Old wrong version    prespeclist = [prespeclist;[possibleprespecs zeros(nposs,n4-n2)]];
                       prespeclist = [prespeclist;[possibleprespecs zeros(nposs,n2-n4)]];
                   else
% Error found 20050202 PS                       
% Old wrong version 	prespeclist =    [possibleprespecs zeros(nposs,n4-n2)];
                        prespeclist =    [possibleprespecs zeros(nposs,n2-n4)];
                   end
                end
                if jjmid == specorder-2
                    midspeclist = [midspeclist;possiblemidspecs];
                else
                    midspeclist = [midspeclist;[possiblemidspecs zeros(nposs,specorder-2-jjmid) ]];    
                end
                if kkpost == 0
                    possiblepostspecs = zeros(nposs,1);    
                end
                if specorder <= 4
                    if specorder == 3
                        postspeclist = [postspeclist;[possiblepostspecs]];
                    else
                        postspeclist = [postspeclist;[possiblepostspecs zeros(nposs,1)]];
                    end
                else
                    if kkpost == 0
                        postspeclist = [postspeclist;[possiblepostspecs zeros(nposs,specorder-3)]];
                    else
                        postspeclist = [postspeclist;[possiblepostspecs zeros(nposs,specorder-2-kkpost)]];                    
                    end
                end
                bigedgeweightlist = [bigedgeweightlist;edgeweightlist];  
	
                % NB! It is correct below that the indices for the ISCOORDS should be
                % ORIGINSFROM(ORIGINSFROM(masterivlist)), rather than masterivlist.
                % The combinations in POTENTIALISES(masterivlist,:) all have
                % spec-spec-...-diff-diff combinations and then
                % ISCOORDS(masterivlist,:) are zeros since a comb. that
                % ends with a diff has no image source. 
                % Also, two recursive references are needed since we need to get back
                % through the two last diffractions to reach the last specular
                % reflection.
			
                ivref = ORIGINSFROM(ORIGINSFROM(ORIGINSFROM(expandedmasterivlist)));
                for kk = 2:jjmid+kkpost
                    ivref = ORIGINSFROM(ivref);
                end
                validIScoords = [validIScoords;ISCOORDS(ivref,:)];
                if kkpost == 0
                    validIRcoords = [validIRcoords;R(ones(nposs,1),:)];
                else
                    if jjmid ~= 2
                        validIRcoords = [validIRcoords;validEDIRcoords(ivOK(ivreftoindata(expandedivcompletelyOK)),:)];
                    else
                        error(['ERROR: Calculation of IR coords for midspec = 2 and postspec not implemented yet!']);    
                    end
                end
            end            
        end
    end
end

%   #######################################################################
%   #
%   #   Pack the edge segments together because every little edge segment
%   #   is present as a separate edge
%   #   This can be done for all combinations at once.
%   #
%   #######################################################################

test = [prespeclist edgedifflist(:,1) midspeclist edgedifflist(:,2) postspeclist];

if ~isempty(test)

	[ncombs,slask] = size(edgedifflist);
	dtest = diff([0;prod(   double(test(1:ncombs-1,:)==test(2:ncombs,:)).'  ).']);
	ivremove = find(dtest==1);
	
	while ~isempty(ivremove)
        bigedgeweightlist(ivremove+1,:) = double(bigedgeweightlist(ivremove+1,:)) + double(bigedgeweightlist(ivremove,:));
        bigedgeweightlist(ivremove,:) = [];
        edgedifflist(ivremove,:) = [];
        prespeclist(ivremove,:) = [];
        midspeclist(ivremove,:) = [];
        postspeclist(ivremove,:) = [];
        validIScoords(ivremove,:) = [];
        validIRcoords(ivremove,:) = [];
	
        test = [prespeclist edgedifflist(:,1) midspeclist edgedifflist(:,2) postspeclist];
    	[ncombs,slask] = size(edgedifflist);
        dtest = diff([0;prod(   double(test(1:ncombs-1,:)==test(2:ncombs,:)).'  ).']);
        ivremove = find(dtest==1);

    end
    
end

%   #######################################################################
%   #
%   #   The weights of the visible edge segments should be
%   #   translated into start and end points, together with the visibility
%   #   weights from the receiver.
%   #   This can be done for all combinations at once!
%   #
%   #######################################################################

% As a start we set the start and end values to 0 and 1, i.e. assuming full
% visibility.

startandendpoints = [startandendpoints;...
        [zeros(length(edgedifflist),1),ones(length(edgedifflist),1),...
         zeros(length(edgedifflist),1),ones(length(edgedifflist),1)]];

% Treat edge 1

ivtoaccountfor = [1:size(edgedifflist,1)].';

ivwholeedge1 = find( bigedgeweightlist(:,1) == maxvisibilityvalue);
if ~isempty(ivwholeedge1)
    ivtoaccountfor(ivwholeedge1) = [];
end

if ~isempty(ivtoaccountfor)
    ncombs = length(ivtoaccountfor);
    bitpattern = zeros(ncombs,nedgesubs);
    for ii=1:nedgesubs
        bitpattern(:,ii) = bitget(bigedgeweightlist(ivtoaccountfor,1),ii); 
    end
    dbit1 = diff([zeros(ncombs,1) bitpattern].').';
    dbit2 = [dbit1 -bitpattern(:,nedgesubs)]; 
    
    nsegments = ceil((sum(abs(dbit1.')).')/2);
    ivonesegments = find(nsegments==1);

    if ~isempty(ivonesegments)

        nonesegments = length(ivonesegments);
        multvec = 2.^[0:nedgesubs];
        segstartpos = round(log(sum( ((dbit2(ivonesegments,:)== 1).*multvec(ones(nonesegments,1),:)).').')/log(2))+1;
        segendpos   = round(log(sum( ((dbit2(ivonesegments,:)==-1).*multvec(ones(nonesegments,1),:)).').')/log(2))+1;

        ivmodify = find(segstartpos==1);
        segstartpos(ivmodify) = ones(size(ivmodify))*1.5;
        ivmodify = find(segendpos>nedgesubs);
        segendpos(ivmodify) = ones(size(ivmodify))*(nedgesubs+0.5);

        startandendpoints(ivtoaccountfor(ivonesegments),1) = (segstartpos-1.5)/(nedgesubs-1);
        startandendpoints(ivtoaccountfor(ivonesegments),2) = (segendpos-1.5)/(nedgesubs-1);
                
    end    

    % If we have some two-or-more-subsegments cases, they will be
    % discovered by the if-condition below
    
    if length(ivonesegments) < ncombs
        for nsegmentstocheck = 2:ceil(nedgesubs/2)
                disp(['Checking for ',int2str(nsegmentstocheck),' sub-segments']) 
            ivNsegments = find(nsegments==nsegmentstocheck);
            if ~isempty(ivNsegments)
                [n1,n2] = size(startandendpoints);
                if n2 < 4*nsegmentstocheck
                    startandendpoints = [startandendpoints zeros(n1,4*nsegmentstocheck-n2)];    
                end
                for jj = 1:length(ivNsegments)
                    ivstartbits = find(dbit2(ivNsegments(jj),:) == 1);
                    ivstartbits = (ivstartbits==1)*1.5 + (ivstartbits~=1).*ivstartbits;
                    ivendbits = find(dbit2(ivNsegments(jj),:) == -1);
                    ivendbits = (ivendbits>nedgesubs)*(nedgesubs+0.5) + (ivendbits<=nedgesubs).*ivendbits;
                    
                    for kk = 1:nsegmentstocheck,                                                        
                        startandendpoints(ivtoaccountfor(ivNsegments(jj)),(kk-1)*4+1) = (ivstartbits(kk)-1.5)/(nedgesubs-1);
                        startandendpoints(ivtoaccountfor(ivNsegments(jj)),(kk-1)*4+2) = (ivendbits(kk)-1.5)/(nedgesubs-1);
                    end
                end                
            end
            
        end
    end        
end

% Treat edge 2

ivtoaccountfor = [1:size(edgedifflist,1)].';

ivwholeedge2 = find( bigedgeweightlist(:,2) == maxvisibilityvalue);
if ~isempty(ivwholeedge2)
    ivtoaccountfor(ivwholeedge2) = [];
end

if ~isempty(ivtoaccountfor)
    ncombs = length(ivtoaccountfor);
    bitpattern = zeros(ncombs,nedgesubs);
    for ii=1:nedgesubs
        bitpattern(:,ii) = bitget(bigedgeweightlist(ivtoaccountfor,2),ii); 
    end
    dbit1 = diff([zeros(ncombs,1) bitpattern].').';
    dbit2 = [dbit1 -bitpattern(:,nedgesubs)]; 
    
    nsegments = ceil((sum(abs(dbit1.')).')/2);
    ivonesegments = find(nsegments==1);

    if ~isempty(ivonesegments)

        nonesegments = length(ivonesegments);
        multvec = 2.^[0:nedgesubs];
        segstartpos = round(log(sum( ((dbit2(ivonesegments,:)== 1).*multvec(ones(nonesegments,1),:)).').')/log(2))+1;
        segendpos   = round(log(sum( ((dbit2(ivonesegments,:)==-1).*multvec(ones(nonesegments,1),:)).').')/log(2))+1;

        ivmodify = find(segstartpos==1);
        segstartpos(ivmodify) = ones(size(ivmodify))*1.5;
        ivmodify = find(segendpos>nedgesubs);
        segendpos(ivmodify) = ones(size(ivmodify))*(nedgesubs+0.5);

        startandendpoints(ivtoaccountfor(ivonesegments),3) = (segstartpos-1.5)/(nedgesubs-1);
        startandendpoints(ivtoaccountfor(ivonesegments),4) = (segendpos-1.5)/(nedgesubs-1);
                
    end    

    % If we have some two-or-more-subsegments cases, they will be
    % discovered by the if-condition below
    
    if length(ivonesegments) < ncombs
        for nsegmentstocheck = 2:ceil(nedgesubs/2)
 
            ivNsegments = find(nsegments==nsegmentstocheck);
            if ~isempty(ivNsegments)
                [n1,n2] = size(startandendpoints);
                if n2 < 4*nsegmentstocheck
                    startandendpoints = [startandendpoints zeros(n1,4*nsegmentstocheck-n2)];    
                end
                for jj = 1:length(ivNsegments)
                    ivstartbits = find(dbit2(ivNsegments(jj),:) == 1);
                    ivstartbits = (ivstartbits==1)*1.5 + (ivstartbits~=1).*ivstartbits;
                    ivendbits = find(dbit2(ivNsegments(jj),:) == -1);
                    ivendbits = (ivendbits>nedgesubs)*(nedgesubs+0.5) + (ivendbits<=nedgesubs).*ivendbits;
                    
                    for kk = 1:nsegmentstocheck,                                                        
                        startandendpoints(ivtoaccountfor(ivNsegments(jj)),(kk-1)*4+3) = (ivstartbits(kk)-1.5)/(nedgesubs-1);
                        startandendpoints(ivtoaccountfor(ivNsegments(jj)),(kk-1)*4+4) = (ivendbits(kk)-1.5)/(nedgesubs-1);
                    end
                end                
            end
            
        end
    end        
end

%   #######################################################################
%   #
%   #   Construct a list guide, which will tell which rows have only
%   #   dd, which rows have sdd etc
%   #   Syntax:        dd,sdd,ssdd,sssdd,...,dds,ddss,ddsss,...
%   #       (should also continue with sdds, sddss,...
%   #   
%   #######################################################################

[n1,n2] = size(prespeclist);
if n2 > 1
    nprespecs  = sum(prespeclist.' > 0).';
else
    nprespecs = (prespeclist>0);    
end
[n1,n2] = size(midspeclist);
if n2 > 1
    nmidspecs  = sum(midspeclist.' > 0).';
else
    nmidspecs = (midspeclist>0);    
end
[n1,n2] = size(postspeclist);
if n2 > 1
    npostspecs  = sum(postspeclist.' > 0).';
else
    npostspecs = (postspeclist>0);    
end

[B,ivec,jvec] = unique([nprespecs nmidspecs npostspecs],'rows');
nuniquecombs = length(ivec);
ntotcombs = length(jvec);

listguide = zeros(nuniquecombs,3);
listofallspecs = zeros(nuniquecombs,3);
sortvec = zeros(ntotcombs,1);
for ii = 1:length(ivec)
    ivfindcombs = find(jvec==ii);
    listguide(ii,1) = length(ivfindcombs);
    if ii > 1
        listguide(ii,2) = listguide(ii-1,3)+1;
    else
        listguide(ii,2) = 1;    
    end
    listguide(ii,3) = listguide(ii,2)+listguide(ii,1)-1;
    listofallspecs(ii,:) = [B(ii,:)]; 

    sortvec(listguide(ii,2):listguide(ii,3)) = ivfindcombs;
    
end

prespeclist = prespeclist(sortvec,:);
midspeclist = midspeclist(sortvec,:);
postspeclist = postspeclist(sortvec,:);
validIScoords = validIScoords(sortvec,:);
validIRcoords = validIRcoords(sortvec,:);
edgedifflist = edgedifflist(sortvec,:);
startandendpoints = startandendpoints(sortvec,:);
