function EDB2plotpath(eddatafile,reflpathsfile,plotoptions)
% EDB2plotpath - Plots a model which is given in an eddatafile
% and one or more of the paths in an edreflpathsfile.
%
% Input parameters:
%   eddatafile      The input file.
%   reflpathsfile   The file with reflpaths.
%   plotoptions     The row in the reflpathsfile that should be plotted.
%
% Sources and receivers are taken from an sdatafile and an rdatafile, the file name of
% which is assumed to be similar to the eddatafile.
%
% Uses functions EDB2strpend, EDB2myvalueinput
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
% Peter Svensson (svensson@iet.ntnu.no) 20090708
%
% EDB2plotpath(eddatafile,reflpathsfile,plotoptions);

%--------------------------------------------------------------
% Read the input files

[eddatafilepath,eddatafile,fileext] = fileparts(eddatafile);
eddatafile = [eddatafile,fileext];
Filestem = EDB2strpend(eddatafile,'_eddata');

[reflpathsfilepath,reflpathsfile,fileext] = fileparts(reflpathsfile);
reflpathsfile = [reflpathsfile,fileext];

if ~isempty(eddatafilepath)
    eval(['load ',eddatafilepath,filesep,eddatafile])
else
    eval(['load ',eddatafile])    
end

if ~isempty(reflpathsfilepath)
    eval(['load ',reflpathsfilepath,filesep,reflpathsfile])
else
    eval(['load ',reflpathsfile])    
end

plotsources = 1;
plotreceivers = 1;

%--------------------------------------------------------------

ncornersperplanevec = double(ncornersperplanevec);
if plotsources
    sdatafile = [eddatafilepath,Filestem,'_sdata.mat'];
    if exist(sdatafile) == 2
        eval(['load ',sdatafile])
    else
        error(['ERROR: The sdata file named ',sdatafile,' could not be opened'])    
    end
end
if plotreceivers
    rdatafile = [eddatafilepath,Filestem,'_rdata.mat'];
    if exist(rdatafile) == 2
        eval(['load ',rdatafile])
    else
        error(['ERROR: The rdata file named ',rdatafile,' could not be opened'])    
    end
end

%--------------------------------------------------------------

[ncorners,slask] = size(corners);
[nedges,slask] = size(edgecorners);
[nplanes,slask] = size(planenvecs);

planelist = 1:nplanes;

% viewpos = EDB2myvalueinput('From which point do you want to watch the model? (x y z, or az el, with spaces inbetween) ',-1);
viewpos = [1 0 0];
axis equal

hold off
for ii = 1:nedges
	co1 = edgecorners(ii,1);
	co2 = edgecorners(ii,2);
	iv = [co1;co2];
	plot3(corners(iv,1),corners(iv,2),corners(iv,3))
    if ii ==1
        view(viewpos)
        hold
    end
end

if plotsources == 1
    plot3(sources(:,1),sources(:,2),sources(:,3),'*');
end

if plotreceivers == 1
    plot3(receivers(:,1),receivers(:,2),receivers(:,3),'ro');
end

%---------------------------------------------------------------
% Plot the path

pathtype = pathtypevec(plotoptions,:)
reflpath = reflpaths(plotoptions,:)
norder =sum(pathtype ~= 0);
nspecorder = sum(pathtype == 115);
speccols = find(pathtype == 115);
diffcols = find(pathtype == 100);


x1 = sources(1,1);
y1 = sources(1,2);
z1 = sources(1,3);

for ii = 1:nspecorder+1
    x2 = specextradata(plotoptions,1+3*(ii-1));
    y2 = specextradata(plotoptions,2+3*(ii-1));
    z2 = specextradata(plotoptions,3+3*(ii-1));    
    plot3(x2,y2,z2,'o')
    line([x1;x2],[y1;y2],[z1;z2]);
    x1 = x2;
    y1 = y2;
    z1 = z2;
end

x2 = receivers(1,1);
y2 = receivers(1,2);
z2 = receivers(1,3);
line([x1;x2],[y1;y2],[z1;z2]);

% Print the involved plane numbers

reflplanes = reflpath(speccols);

for ii = reflplanes
    midpoint = mean(corners(planecorners(ii,1:ncornersperplanevec(ii)),:));
    endpoint = midpoint + planenvecs(ii,:);
    text(endpoint(1),endpoint(2),endpoint(3),int2str(ii));
end

% Print the involved edge numbers

if ~isempty(diffcols)
    diffedges = reflpath(diffcols);
    for ii = diffedges
         midpoint = mean(corners(edgecorners(ii,1:2),:));
         text(midpoint(1),midpoint(2),midpoint(3),int2str(ii));
    	co1 = edgecorners(ii,1);
        co2 = edgecorners(ii,2);
        iv = [co1;co2];
        h = plot3(corners(iv,1),corners(iv,2),corners(iv,3));
        set(h,'LineWidth',3)
         
         
    end
    
end
