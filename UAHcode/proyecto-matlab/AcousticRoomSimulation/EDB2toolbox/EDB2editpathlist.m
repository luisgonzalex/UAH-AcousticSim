function EDB2editpathlist(edpathsfile,useedges,symmetricedges,pathstokeep)
% EDB2editpathlist - Does some semi-automatic editing of an edpathsfile.
% A vector, 'multfactors' is introduced, which will boost or
% keep or switch off paths in the list of paths.
%
% Input parameters:
%   edpathsfile     An edpaths file that should be modified.
%   useedges        (optional) A list of edges that should be used. All
%                   paths that contain an edge which is not this edge will
%                   be switched off.
%   symmetricedges  (optional) A matrix of edge pairs that are symmetrical
%                   around a symmetry plane. If, e.g., edges 2 and 3 are symmetric around the
%                   source/receiver plane, then symmetricedges = [2 3]. Several pairs are
%                   possible; symmetricedges = [2 3;7 8;9 12];
%   pathstokeep     (optional) A matrix of specific diffraction paths to
%                   keep.
% SHOWTEXT (global)
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
% Peter Svensson (svensson@iet.ntnu.no) 20050510
%
% EDB2editpathlist(edpathsfile,useedges,symmetricedges,pathstokeep);

global SHOWTEXT

if isempty(symmetricedges)
    dosymmetry = 0;
else
    dosymmetry = 1;
end

if isempty(pathstokeep)
    dopathpruning = 0;
else
    dopathpruning = 1;
end    

if ~isempty(useedges)
    error(['ERROR: Sorry, the useedges option has not been implemented yet!'])    
end

%-------------------------------------------------------------------------

nsymmetricedgepairs = size(symmetricedges,1);

Varlist = [' S R allspecrows firstdiffrow mainlistguidepattern  mainlistguide reflpaths   pathtypevec'];          
Varlist = [Varlist,' Sinsideplanenumber Rinsideplanenumber    directsoundrow        specextradata  edgeextradata'];         
if dosymmetry | dopathpruning
    Varlist = [Varlist,'  multfactors'];           
end

%-------------------------------------------------------------------------
    
eval(['load ',edpathsfile])
disp(' ')
disp(['Modifying ',edpathsfile])
[ncombs,ncols] = size(reflpaths);
multfactors = ones(ncombs,1);

%-------------------------------------------------------------------------

if dosymmetry

    onesvec = ones(ncombs,1);
    for kkk = 1:ncols
        if SHOWTEXT >= 1
            disp(['   Diffraction order ',int2str(kkk)])    
        end
        irow = sum((mainlistguidepattern=='d').').' == kkk;
        if sum(irow) > 0
            ivselection = [double(mainlistguide(irow,2)):double(mainlistguide(irow,3))].';
            for mmm = 1:nsymmetricedgepairs
                if SHOWTEXT >= 1
                    disp(['   Symmetry pair ',int2str(mmm)])    
                end
                iv = any(reflpaths(ivselection,:).'==symmetricedges(mmm,1));
                if ~isempty(iv)
                    iv = ivselection(iv);
                    for nnn=1:length(iv)
                        reflpathsmatch = reflpaths(iv(nnn),:);
                        if multfactors(iv(nnn)) ~= 0
                            colstoswitch1 = find(reflpathsmatch==symmetricedges(mmm,1));
                            colstoswitch2 = find(reflpathsmatch==symmetricedges(mmm,2));
                            reflpathsmatch(colstoswitch1) = symmetricedges(mmm,2);
                            reflpathsmatch(colstoswitch2) = symmetricedges(mmm,1);
                            
                            reflpathsmatch = reflpathsmatch(onesvec,:);
                            
                            ivmatch = find(sum( (reflpaths==reflpathsmatch).' ).'==ncols);
                            if ~isempty(ivmatch) & multfactors(ivmatch) ~= 0
                                multfactors(iv(nnn)) = 2*multfactors(iv(nnn));    
                                multfactors(ivmatch) = 0;    
                                if SHOWTEXT >= 2
                                    disp(['   Found symmetric edge pair:'])
                                    disp(['   ',int2str(iv(nnn)),': ',int2str(double(reflpaths(iv(nnn),:)))])
                                    disp(['   ',int2str(ivmatch),': ',int2str(double(reflpaths(ivmatch,:)))])
                                end                
                            end
                        end    % ... if multfactors(iv(nnn)) ~= 0,                        
                        
                    end   % ... for nnn=1:length(iv)
                end    % ... if ~isempty(iv)
            end   % for mmm = 1:nsymmetricedgepairs
        end   % ...if sum(irow) > 0
    end     % ... for kkk = 1:ncols
end     % ...if dosymmetry

%-------------------------------------------------------------------------

if dopathpruning
    multfactors = multfactors*0;
    for ii = 1:size(pathstokeep,1)
        pathtomatch = pathstokeep(ii,:);
        pathtomatch = pathtomatch(ones(ncombs,1),:);
        ivmatch = find(sum( (reflpaths==pathtomatch).' ).'==ncols);        
        multfactors(ivmatch) = 1;
    end
end

%-------------------------------------------------------------------------

eval(['save ',edpathsfile,Varlist])
