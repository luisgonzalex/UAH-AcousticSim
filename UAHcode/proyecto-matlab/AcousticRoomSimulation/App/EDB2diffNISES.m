function [edgedifflist,startandendpoints,prespeclist,postspeclist,validIScoords,validIRcoords,listguide,...
    listoforders] = EDB2diffNISES(eddatafile,S,R,...
    ivNdiffmatrix,lengthNdiffmatrix,Ndifforder,specorder,visplanesfromR,vispartedgesfromS,...
    vispartedgesfromR,nedgesubs,ndecimaldivider,edgedifflistin,postspeclistin,bigedgeweightlistin,validEDIRcoords)
% EDB2diffNISES - Gives a list of paths that includes an N-order diffraction path.
% Gives the list of visible N-diffraction paths, possibly with specular reflections before and after.
%
% Input parameters:
%    eddatafile,S,R,ivsinglediff,...
%    singlediffcol,startindicessinglediff,endindicessinglediff,specorder,visplanesfromR,...
%    vispartedgesfromS,vispartedgesfromR,nedgesubs,ndecimaldivider,PointertoIRcombs,IRoriginsfrom;
%
% GLOBAL parameters:
%   SHOWTEXT JJ JJnumbofchars   See EDB2mainISES
%   POTENTIALISES,ISCOORDS,ORIGINSFROM,ISESVISIBILITY,REFLORDER   See EDB2findISEStree
%
% Output parameters:
%   edgedifflist        List [ncombs,N] of the edge numbers involved in each
%                       spec-diff-...-diff-spec combination.
%   startandendpoints   Matrix [ncombs,4] of the relative start and end
%                       points of each edge. The values, [0,1], indicate
%                       which part of the two edges that are visible.
%   prespeclist         Matrix [ncombs,specorder-N] of the specular
%                       reflections that precede every diffraction.
%   postspeclist        Matrix [ncombs,specorder-N] of the specular
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
%                          are the same type of spec-diff-...-diff-spec comb.
%                       2. The first row number and 3. The last row number.
%   listoforders        Matrix [nuniquecombs,2] which for each row gives
%                       1. The reflection order for the spec-diff-spec comb
%                          in the same row in listguide.
%                       2. The order of the diffraction in the
%                          spec-diff-...-diff-spec comb.
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
% Peter Svensson (svensson@iet.ntnu.no) 20031018
%
% [edgedifflist,startandendpoints,prespeclist,postspeclist,validIScoords,validIRcoords,...
%   listguide] = EDB2diff3ISES(eddatafile,S,R,...
%   ivNdiffmatrix,lengthNdiffmatrix,Ndifforder,specorder,visplanesfromR,vispartedgesfromS,...
%   vispartedgesfromR,nedgesubs,ndecimaldivider,edgedifflistin,postspeclistin,bigedgeweightlist,validEDIRcoords);

global SHOWTEXT JJ JJnumbofchars
global POTENTIALISES ISCOORDS ORIGINSFROM ISESVISIBILITY REFLORDER

eval(['load ',eddatafile])

[nedges,slask] = size(planesatedge);
[nplanes,slask] = size(planecorners);
multfac = 10^(ceil(log10(max([nedges nplanes]))));

edgedifflist = [];
prespeclist = [];
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
zerosvec1 = zeros(1,specorder-Ndifforder);
zerosvec2 = zeros(1,3);
listguide = zeros(specorder*2-1,3);

obstructtestneeded = (sum(canplaneobstruct)~=0);
onesvec = ones(1,nedgesubs);
onesvec3 = ones(1,3);

npostspecs = sum(double(postspeclistin.'>0)).';

%   ###########################################
%   #                                         #
%   #   S - edge - edge - edge -R cases       #
%   #                                         #
%   ###########################################
%
% Possible edges for S-E-E-E-R are seen (at least partly) by the source and by the receiver.
%
% The visibility by the source is taken care of by the findISEStree
% so we must check which ones are visible by the receiver.
% Also, the active edge segment(s) must be selected but this is done
% further down together with the S-spec-spec-edge-R cases

ivNdiff = ivNdiffmatrix(1:lengthNdiffmatrix(Ndifforder),Ndifforder);
ivsinglediff = ivNdiffmatrix(1:lengthNdiffmatrix(1),1);

ivSEEER = ivNdiff(find(REFLORDER(ivNdiff)==Ndifforder));

possibleedgecombs = double(POTENTIALISES(ivSEEER,1:Ndifforder))-nplanes;
ivnotvisiblefromr = find(vispartedgesfromR(possibleedgecombs(:,Ndifforder))==0);
if ~isempty(ivnotvisiblefromr)
    possibleedgecombs(ivnotvisiblefromr,:) = [];
end

edgedifflist      = [edgedifflist;possibleedgecombs];
bigedgeweightlist = [bigedgeweightlist;[vispartedgesfromS(edgedifflist(:,1)) vispartedgesfromR(edgedifflist(:,Ndifforder))]];

[nedgesadded,slask] = size(edgedifflist);
zerosvec3 = zeros(nedgesadded,1);
ndiffonly = nedgesadded;

prespeclist =  [prespeclist;zerosvec3(:,ones(1,max(specorder-Ndifforder,1)))];
postspeclist = [postspeclist;zerosvec3(:,ones(1,max(specorder-Ndifforder,1)))];
validIScoords = [validIScoords;S(ones(nedgesadded,1),:)];
validIRcoords = [validIRcoords;R(ones(nedgesadded,1),:)];

if SHOWTEXT >= 3
	disp(['         ',int2str(nedgesadded),' Nth diff valid'])
end

%   ###########################################
%   ###########################################
%   ###########################################
%   ###########################################
%   #                                         #
%   #         Prespec cases                   #
%   #                                         #
%   ###########################################
%
% Possible edges for S-spec-spec-E-R are seen (at least partly) by the receiver.
%
% The visibility doesn't need to be checked since the source-to-edge paths
% were checked in the ISEStree, and the visibility from the receiver also
% has been checked.

% The vector ivmultidiff will always refer to the original data vector
% i.e. POTENTIALISES, ORIGINSFROM, REFLORDER etc
%
% We should remove the combinations which involve an edge that the
% receiver can not see, but since the edge number is in different columns
% we do that in the for loop.

% The ii-loop will go through: spec-diff, spec-spec-diff
% spec-spec-spec-diff etc

for ii = 1:specorder-Ndifforder
    
    if SHOWTEXT >= 3
        disp(['      Checking for ',JJ(ii,1:JJnumbofchars(ii)),' spec refl before the  Nth edge diff'])    
    end

    % Select the combinations where the reflection order == ii+Ndifforder
    % (which means ii specular reflections in addition to the diffraction)
    % and where the last Ndifforder columns contain edges.
    
    iv = find((REFLORDER(ivNdiff)==ii+Ndifforder) & prod(   double( POTENTIALISES(ivNdiff,ii+1:ii+Ndifforder)>nplanes).'  ).'   );
    masterivlist = ivNdiff(iv);
    possibleedgecombs = double(POTENTIALISES(masterivlist,ii+1:ii+Ndifforder)) - nplanes;
    
    % Keep only combinations for which the receiver can see the edge
    
    ivnotvisiblefromr = find(vispartedgesfromR(possibleedgecombs(:,Ndifforder))==0);
    if ~isempty(ivnotvisiblefromr)
        masterivlist(ivnotvisiblefromr) = [];
        possibleedgecombs = double(POTENTIALISES(masterivlist,ii+1:ii+Ndifforder)) - nplanes;
    end        

    possibleprespecs = POTENTIALISES(masterivlist,1:ii);
    
    reftoIScoords = ORIGINSFROM(masterivlist);
    edgeweightlist = [ISESVISIBILITY(masterivlist) vispartedgesfromR(possibleedgecombs(:,Ndifforder))];
        
    nposs = length(masterivlist);

    if SHOWTEXT >= 3
 		disp(['         ',int2str(nposs),' IS+edge multiples valid'])
    end

    edgedifflist = [edgedifflist;possibleedgecombs];
    prespeclist = [prespeclist;[possibleprespecs zeros(nposs,specorder-Ndifforder-ii)]];        
    postspeclist = [postspeclist;zeros(nposs,specorder-Ndifforder)];
    bigedgeweightlist = [bigedgeweightlist;edgeweightlist];    
    % The ref. to ISCOORDS must be nested the correct number of time to
    % come back to the image source coords of the last spec. refl.
    ivref = ORIGINSFROM(ORIGINSFROM(ORIGINSFROM(masterivlist)));
    for jj = 4:Ndifforder
        ivref = ORIGINSFROM(ivref);    
    end
    validIScoords = [validIScoords;ISCOORDS(ivref,:)];
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
% i.e. PotentialIR, ORIGINSFROMIR, reflorderIR etc
%
% First we pick out those indices where there was a single diffraction, but
% skip those with only diffraction (because we dealt with them already).
% Also, select only those where the diffraction is the first in the sequence
% of reflections.

for ii = 1:specorder-Ndifforder

    if SHOWTEXT >= 3
        disp(['      Checking for ',JJ(ii,1:JJnumbofchars(ii)),' spec refl after the multiple edge diff'])    
    end

    % Select the combinations where the reflection order == ii+Ndifforder
    % (which means ii specular reflections in addition to the diffraction)
    % and where the first Ndifforder columns contain edges.
   
    iv = find((REFLORDER(ivNdiff)==ii+Ndifforder) & prod(   double( POTENTIALISES(ivNdiff,1:Ndifforder)>nplanes).'  ).'   );
    masterivlist = ivNdiff(iv);
    possibleedgecombs = double(POTENTIALISES(masterivlist,1:Ndifforder)) - nplanes;
    possiblepostspecs = POTENTIALISES(masterivlist,Ndifforder+1:specorder);
    possibleweights = ISESVISIBILITY(masterivlist);
    
    % Compare with those that have already been found OK
    ivOK = find(npostspecs==ii);
    if ~isempty(ivOK)
        patternOK = [edgedifflistin(ivOK) postspeclistin(ivOK,1:ii)];
    else
        patternOK = [];    
    end
    % Find out which ones have been checked and found invisible/obstructed
    ivallcombs = ivsinglediff(find( POTENTIALISES(ivsinglediff,1)>nplanes & REFLORDER(ivsinglediff) == ii+1));
    patternALL = [double(POTENTIALISES(ivallcombs,1))-nplanes double(POTENTIALISES(ivallcombs,2:1+ii))];
    if ~isempty(patternOK)
        patternNOTOK = setdiff(patternALL,patternOK,'rows');
    else
        patternNOTOK = patternALL; 
    end
    patterntocompare = [possibleedgecombs(:,Ndifforder) possiblepostspecs(:,1:ii)];
    
    ivtocancel = find(ismember(patterntocompare,patternNOTOK,'rows'));
    masterivlist(ivtocancel) = [];
    possibleedgecombs(ivtocancel,:) = [];
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
        % These combinations have a postspec-combination that hasn't been
        % encountered in the single diffraction cases, so no visibility
        % test has been made for these.
        
        lastedgenumbers = double(POTENTIALISES(masterlisttocheckmore,Ndifforder))-nplanes;
        newIRcoords = R;
        reflplanesexpand = zeros(ntocheckmore*nedgesubs,ii);
        for jj = 1:ii
            reflplanes = POTENTIALISES(masterlisttocheckmore,Ndifforder+1+ii-jj);
            reflplanesexpand(:,jj) = reshape(reflplanes(:,onesvec).',ntocheckmore*nedgesubs,1);
            newIRcoords = EDB2findis(newIRcoords,reflplanes,planeeqs,size(newIRcoords,1),onesvec3);
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
            disp(['         ',int2str(ntocheckmore),' special N-edge - IR combinations to check'])    
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
                    disp('         implemented.')
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
                disp(['         ',int2str(ntocheckmore),' of the special N-edge - IR combinations survived the visibility test in refl plane ',int2str(jj)])
            end
        end

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
                    
                    [nonobstructedpaths,nobstructions,edgehits,cornerhits] = EDB2checkobstrpaths(fromcoords,tocoords,startplanes,endplanes,canplaneobstruct,planeseesplane,...
                        planeeqs,planenvecs,minvals,maxvals,planecorners,corners,ncornersperplanevec,rearsideplane);
                    if ~isempty(edgehits) | ~isempty(cornerhits)
                        disp('WARNING! An edgehit or cornerhit occurred during the obstruction test but this is not')
                        disp('         implemented.')
                    end

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
                disp(['         ',int2str(ntocheckmore),' of the special N-edge - IR combinations survived the obstruction test'])
            end
            
        end

        % Add the found special combinations to the outdata list
    
        edgedifflist = [edgedifflist;double(POTENTIALISES(masterlisttocheckmore,1:Ndifforder))-nplanes];
        prespeclist = [prespeclist;zerosvec1(ones(ntocheckmore,1),:)];
        postspeclist = [postspeclist;[reflplanesexpand zeros(ntocheckmore,specorder-Ndifforder-ii)]];
        bigedgeweightlist = [bigedgeweightlist;[ ISESVISIBILITY(masterlisttocheckmore) edgeweightlist]];    
        eval(['validIRcoords = [validIRcoords;newIRcoords',JJ(ii,1:JJnumbofchars(ii)),'];']);
        validIScoords = [validIScoords;S(ones(ntocheckmore,1),:)];
        
    end
    
    masterivlist = masterivlist(ivcompletelyOK);    
    possibleedgecombs = possibleedgecombs(ivcompletelyOK,:);
    possiblepostspecs = possiblepostspecs(ivcompletelyOK,:);
    possibleweights = possibleweights(ivcompletelyOK,:);
        
    nposs = length(ivcompletelyOK);

    if SHOWTEXT >= 3
		disp(['         ',int2str(nposs),' Edge,Nth-order + IR segments (non-special) survived the obstruction test'])
    end
    
    % Add the found "standard" combinations to the outdata list

    edgedifflist = [edgedifflist;possibleedgecombs];
    prespeclist = [prespeclist;zeros(nposs,specorder-Ndifforder)];
    postspeclist = [postspeclist;possiblepostspecs ];
    bigedgeweightlist = [bigedgeweightlist;[possibleweights bigedgeweightlistin(ivOK(ivreftoindata(ivcompletelyOK)))]];    
    validIRcoords = [validIRcoords;validEDIRcoords(ivOK(ivreftoindata(ivcompletelyOK)),:)];
    validIScoords = [validIScoords;S(ones(nposs,1),:)];
    
end


%   #######################################################################
%   #######################################################################
%   #######################################################################
%   #######################################################################
%   #######################################################################

%   #############################################
%   #                                           #
%   #         S - spec - edge - spec - R cases  #
%   #                                           #s
%   #         Pre and postspec cases            #
%   #                                           #
%   #############################################

% The ii- and jj-loops will go through all combinations of spec before and after
% the Ndifforder diffractions.
%
% Unlike before, we will look in the list of already found combinations
% (edgedifflist etc)


for ii = 1:specorder-Ndifforder-1

    for jj = 1:specorder-ii-Ndifforder
    
        if SHOWTEXT >= 3
            disp(['      Checking for ',JJ(ii,1:JJnumbofchars(ii)),' spec refl before and ',JJ(jj,1:JJnumbofchars(jj)),' spec refl after the N-order edge diff'])    
        end
	
        % Select the combinations where the reflection order ==
        % ii+jj+Ndifforder
        % (which means ii+jj specular reflections in addition to the
        % diffraction), and where columns ii+1:ii+Ndifforder of POTENTIALISES
        % contain edges.
	
        iv = find((REFLORDER(ivNdiff) == ii+jj+Ndifforder) & prod(   double( POTENTIALISES(ivNdiff,ii+1:ii+Ndifforder)>nplanes).'  ).'      );
        masterivlist = ivNdiff(iv);
        possibleedgecombs = double(POTENTIALISES(masterivlist,ii+1:ii+Ndifforder)) - nplanes;
        possibleprespecs  = POTENTIALISES(masterivlist,1:ii);
        possiblepostspecs = POTENTIALISES(masterivlist,ii+Ndifforder+1:ii+Ndifforder+jj);
        possibleweights   = ISESVISIBILITY(masterivlist);

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
        if (   (~isempty(patternOK)) &  (~isempty(patternALL)) )
            patternNOTOK = setdiff(patternALL,patternOK,'rows');
        else
            if (   (~isempty(patternOK)) )
                patternNOTOK = patternOK;
            else
                patternNOTOK = patternALL;
            end
        end
        patterntocompare = [possibleedgecombs(:,Ndifforder) possiblepostspecs(:,1:jj)];

        ivtocancel = find(ismember(patterntocompare,patternNOTOK,'rows'));
        masterivlist(ivtocancel) = [];
        possibleedgecombs(ivtocancel,:) = [];
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

            lastedgenumbers = double(POTENTIALISES(masterlisttocheckmore,ii+Ndifforder))-nplanes;
            newIRcoords = R;
            reflplanesexpand = zeros(ntocheckmore*nedgesubs,ii);
            for kk = 1:jj
                reflplanes = POTENTIALISES(masterlisttocheckmore,Ndifforder+1+ii+jj-kk);
                reflplanesexpand(:,kk) = reshape(reflplanes(:,onesvec).',ntocheckmore*nedgesubs,1);
                newIRcoords = EDB2findis(newIRcoords,reflplanes,planeeqs,size(newIRcoords,1),onesvec3);
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
                disp(['         ',int2str(ntocheckmore),' special IS - N-edge - IR combinations to check'])    
            end
            
            for kk = jj:-1:1
                if length(masterlisttocheckmore) > 0
                    eval(['fromcoords = newIRcoords',JJ(kk,1:JJnumbofchars(kk)),';']);
                    if kk < jj
                        eval(['tocoords = reflpoints',JJ(kk+1,1:JJnumbofchars(kk+1)),';'])    
                    end
                    
                    [hitplanes,reflpoints,edgehits,edgehitpoints,cornerhits,cornerhitpoints] = EDB2chkISvisible(fromcoords,tocoords,planeeqs(reflplanesexpand(:,kk),4),planenvecs(reflplanesexpand(:,kk),:),minvals(reflplanesexpand(:,kk),:),...
					    maxvals(reflplanesexpand(:,kk),:),planecorners(reflplanesexpand(:,kk),:),corners,ncornersperplanevec(reflplanesexpand(:,kk)));
                    if ~isempty(edgehits) | ~isempty(cornerhits)
                        disp('WARNING! An edgehit or cornerhit occurred during the visibility test but this is not')
                        disp('         implemented.')
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
                    disp(['         ',int2str(ntocheckmore),' of the special IS - N-edge - IR combinations survived the visibility test in refl plane ',int2str(kk)])
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
                    disp(['         ',int2str(ntocheckmore),' of the special IS - N-edge - IR combinations survived the obstruction test'])
                end
                
            end
            
            % Add the found special combinations to the outdata list
        
            edgedifflist = [edgedifflist;double(POTENTIALISES(masterlisttocheckmore,ii+1:ii+Ndifforder))-nplanes];
            prespeclist = [prespeclist;[double(POTENTIALISES(masterlisttocheckmore,1:ii)) zeros(ntocheckmore,specorder-Ndifforder-ii)]];
            postspeclist = [postspeclist;[reflplanesexpand zeros(ntocheckmore,specorder-Ndifforder-ii)]];
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
        possibleedgecombs = possibleedgecombs(ivcompletelyOK,:);
        possibleprespecs = possibleprespecs(ivcompletelyOK,:);
        possiblepostspecs = possiblepostspecs(ivcompletelyOK,:);
        possibleweights = possibleweights(ivcompletelyOK,:);
            
        nposs = length(ivcompletelyOK);
	
        if SHOWTEXT >= 3
			disp(['         ',int2str(nposs),' IS + Edge (Nth-order) + IR segments (non-special combinations) survived the obstruction test'])
        end

        % Add the found "standard" combinations to the outdata list
    
        edgedifflist = [edgedifflist;possibleedgecombs];
        prespeclist = [prespeclist;  [possibleprespecs  zeros(nposs,specorder-Ndifforder-ii)]];
        postspeclist = [postspeclist;[possiblepostspecs zeros(nposs,specorder-Ndifforder-jj)]];
        bigedgeweightlist = [bigedgeweightlist;[possibleweights bigedgeweightlistin(ivOK(ivreftoindata(ivcompletelyOK)))]];    
        validIRcoords = [validIRcoords;validEDIRcoords(ivOK(ivreftoindata(ivcompletelyOK)),:)];
        % NB!! In the same way as earlier, we must a recursive reference
        % method to find the image source of the last specular reflection.
        ivref = ORIGINSFROM(ORIGINSFROM(ORIGINSFROM(ORIGINSFROM(ORIGINSFROM(masterivlist)))));
        for kk = 6:jj+Ndifforder
            ivref = ORIGINSFROM(ivref);
        end
        validIScoords = [validIScoords;ISCOORDS(ivref,:)];
        
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

addvec = [];
for ii = 1:Ndifforder
    addvec = [addvec,zeros(length(edgedifflist),1),ones(length(edgedifflist),1)];
end
startandendpoints = [startandendpoints addvec];

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
[n1,n2] = size(postspeclist);
if n2 > 1
    npostspecs  = sum(postspeclist.' > 0).';
else
    npostspecs = (postspeclist>0);    
end

listofdiffcol   = 1 + nprespecs;
listofreflorder = listofdiffcol + npostspecs + Ndifforder - 1;

[B,ivec,jvec] = unique([listofreflorder listofdiffcol],'rows');
nuniquecombs = length(ivec);
ntotcombs = length(jvec);

listguide = zeros(nuniquecombs,3);
listoforders = zeros(nuniquecombs,2);
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
    listoforders(ii,:) = [B(ii,1) B(ii,2)]; 

    sortvec(listguide(ii,2):listguide(ii,3)) = ivfindcombs;
    
end

prespeclist = prespeclist(sortvec,:);
postspeclist = postspeclist(sortvec,:);
validIScoords = validIScoords(sortvec,:);
validIRcoords = validIRcoords(sortvec,:);
edgedifflist = edgedifflist(sortvec,:);
startandendpoints = startandendpoints(sortvec,:);
