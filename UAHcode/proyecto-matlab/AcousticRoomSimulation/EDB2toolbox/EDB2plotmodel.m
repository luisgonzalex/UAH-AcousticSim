function EDB2plotmodel(eddatafile,plotoptions,plotoptions2,plotoptions3)
% EDB2plotmodel - Plots a model which is given in an eddatafile.
%
% Input parameters:
%   eddatafile (optional)   The input file. If not specified, a file
%                           opening window will be presented.
%   plotoptions (optional)  An integer that can give extra options:
%                           if bit0 = 1 (= 1) => plot sources
%                           if bit1 = 1 (= 2) => plot receivers
%                           if bit2 = 1 (= 4) => plot plane normal vectors
%                           if bit3 = 1 (= 8) => print plane numbers
%                           if bit4 = 1 (=16) => print edge numbers
%                           if bit5 = 1 (=32) => print corner numbers
%                           if bit6 = 1 (=64) => print plane numbers using the
%                                                CAD file numbering
%                           if bit7 = 1 (=128)=> print corner numbers using the
%                                                CAD file numbering
%                           if bit8 = 1 (=256) => offedges are not included
%                           if bit9 = 1 (=512) => only offedges are included
%                           if plotoptions == -1 => plot edge by edge
%                           if plotoptions == -2 => plot plane by plane
%                           Example: the integer 11 = 8+2+1
%                           so sources, receivers, and plane numbers will be plotted.
%   plotoptions2 (optional) A matrix with two columns of 1 or many integers, the first
%                           giving the source numbers to be plotted and the
%                           second giving the receiver numbers to be
%                           plotted. The first in the sequence will be
%                           plotted with a black color.
%                           Alternatively, if plotoptions = -1, then
%                           plotoptions2 can be used to give a list of
%                           edges to plot.
%                           Alternatively, if plotoptions = -2, then
%                           plotoptions2 can be used to give a list of
%                           planes to plot.
%   plotoptions3 (optional) A vector with two or three values, giving the
%                           view position for the 3D plot view.
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
% Peter Svensson (svensson@iet.ntnu.no) 20110704
%
% EDB2plotmodel(eddatafile,plotoptions,plotoptions2,plotoptions3);

%--------------------------------------------------------------
% Read the eddatafile

if nargin == 0 | (nargin >= 1 & isstr(eddatafile) ~= 1)
	[eddatafile,eddatafilepath] = uigetfile('*eddata.mat','Please select the eddatafile');
    [eddatafilepath,temp1,temp2] = fileparts(eddatafilepath);
	if ~isstr(eddatafile) | isempty(eddatafile)
		return
	end
else
	[eddatafilepath,eddatafile,fileext] = fileparts(eddatafile);
    eddatafile = [eddatafile,fileext];
end
Filestem = EDB2strpend(eddatafile,'_eddata');

disp(['eddatafile is ',eddatafilepath,filesep,eddatafile])

if nargin == 1
    if isstr(eddatafile) ~= 1
        plotoptions = eddatafile;
    else
        plotoptions = 0;    
    end
    plotoptions2 = [];
    plotoptions3 = [];
elseif nargin == 0
    plotoptions = 0;    
    plotoptions2 = [];
    plotoptions3 = [];
elseif nargin == 2
    plotoptions2 = [];
    plotoptions3 = [];
elseif nargin == 3
    nplotoptions2 = size(plotoptions2);
    plotoptions3 = [];
else
    nplotoptions2 = size(plotoptions2);    
end

plotedgesequence = 0; 
plotplanesequence = 0;
if plotoptions < 0    
    if plotoptions == -1
        plotedgesequence = 1;
    elseif plotoptions == -2 
        plotplanesequence = 1;
    end
    plotoptions = 0;
end
plotsources = bitget(plotoptions,1);
plotreceivers = bitget(plotoptions,2);
plotnvecs = bitget(plotoptions,3);
plotplnumbers = bitget(plotoptions,4);
plotednumbers = bitget(plotoptions,5);
plotconumbers = bitget(plotoptions,6);
plotplCADnumbers = bitget(plotoptions,7);
plotcoCADnumbers = bitget(plotoptions,8);
plotnooffedges = bitget(plotoptions,9);
plotonlyoffedges = bitget(plotoptions,10);

if plotplCADnumbers, plotplnumbers = 1; end
if plotcoCADnumbers, plotconumbers = 1; end

%--------------------------------------------------------------------------
% Load the needed input files

if ~isempty(eddatafilepath)
    eval(['load ',eddatafilepath,filesep,eddatafile])
else
    eval(['load ',eddatafile])    
end
    
ncornersperplanevec = double(ncornersperplanevec);
if plotsources
    sdatafile = [eddatafilepath,filesep,Filestem,'_sdata_1.mat'];
    if exist(sdatafile) == 2
        allsources = [];
        sfilecounter = 1;
        while exist(sdatafile) == 2
            eval(['load ',sdatafile])
            allsources = [allsources;sources];
            sfilecounter = sfilecounter +1;
            sdatafile = [eddatafilepath,filesep,Filestem,'_sdata_',int2str(sfilecounter),'.mat']            
        end 
        sources = allsources;
    else
        sdatafile = [eddatafilepath,filesep,Filestem,'_sdata.mat'];
        if exist(sdatafile) == 2
            eval(['load ',sdatafile])
        else
            error(['ERROR: The sdata file named ',sdatafile,' could not be opened'])    
        end
    end
end
if plotreceivers
    rdatafile = [eddatafilepath,filesep,Filestem,'_rdata_1.mat'];
    if exist(rdatafile) == 2
        allreceivers = [];
        rfilecounter = 1;
        while exist(rdatafile) == 2
            eval(['load ',rdatafile])
            allreceivers = [allreceivers;receivers];
            rfilecounter = rfilecounter +1;
            rdatafile = [eddatafilepath,filesep,Filestem,'_rdata_',int2str(rfilecounter),'.mat']            
        end 
        receivers = allreceivers;
    else
        rdatafile = [eddatafilepath,filesep,Filestem,'_rdata.mat'];
        if exist(rdatafile) == 2
            eval(['load ',rdatafile])
        else
            error(['ERROR: The rdata file named ',rdatafile,' could not be opened'])    
        end
    end    
end

if plotplCADnumbers | plotcoCADnumbers
    cadgeofile = [eddatafilepath,Filestem,'_cadgeo.mat'];
    if exist(cadgeofile) == 2
        eval(['load ',cadgeofile])
    else
        error(['ERROR: The cadgeo file named ',cadgeofile,' could not be opened'])            
    end
end

%--------------------------------------------------------------

[ncorners,slask] = size(corners);
[nedges,slask] = size(edgecorners);
[nplanes,slask] = size(planenvecs);

planelist = 1:nplanes;

if isempty(plotoptions3)
    viewpos = EDB2myvalueinput('From which point do you want to watch the model? (x y z, or az el, with spaces inbetween) ',-1);
else
    viewpos = plotoptions3;
end

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

if plotedgesequence
    for ii = 1:nedges
    	co1 = edgecorners(ii,1);
        co2 = edgecorners(ii,2);
        iv = [co1;co2];
        h = plot3(corners(iv,1),corners(iv,2),corners(iv,3));
        set(h,'LineWidth',2)
        disp(['Edge no. ',int2str(ii),' Closed wedge angle: ',num2str(closwedangvec(ii)*180/pi)])
        pause
        set(h,'LineWidth',1)    
    
    end
end
if plotplanesequence
    if ~isempty(plotoptions2)
       listofplanes = plotoptions2; 
    else
        listofplanes = [1:nplanes];
    end
    for jj = listofplanes
        edgelist = edgesatplane(jj,:);
        edgelist = edgelist(find(edgelist));
        for ii = edgelist
            co1 = edgecorners(ii,1);
            co2 = edgecorners(ii,2);
            iv = [co1;co2];
            h = plot3(corners(iv,1),corners(iv,2),corners(iv,3));
            set(h,'LineWidth',2)
        end
        disp(['Plane no. ',int2str(jj)])
        pause
        
        % For some reason it doesnt work to replot the same list of edges
        % with a thinner linewidth; even trying out with replotting the
        % edges with another color just replots one of the edges??
        % Therefore the inelegant method of replotting the whole model.
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
        
        
        
    end
end

if plotnvecs
	for ii = 1:nplanes
        midpoint = mean(corners(planecorners(ii,1:(ncornersperplanevec(ii))),:));
        
        endpoint = midpoint + planenvecs(ii,:);
        bothpoints = [midpoint;endpoint];
        plot3(bothpoints(:,1),bothpoints(:,2),bothpoints(:,3));
        plot3(midpoint(1),midpoint(2),midpoint(3),'ro');
	end
end

if plotplnumbers
	for ii = 1:nplanes
        midpoint = mean(corners(planecorners(ii,1:ncornersperplanevec(ii)),:));
        endpoint = midpoint + planenvecs(ii,:)*0.1;
        if plotplCADnumbers == 1
            text(endpoint(1),endpoint(2),endpoint(3),int2str(planenumbers(ii)));
        else
            text(endpoint(1),endpoint(2),endpoint(3),int2str(ii));
        end
    
	end
end

if plotednumbers
    iv = [1:nedges];
    if plotnooffedges, iv(offedges) = []; end
    if plotonlyoffedges, iv = offedges; end
    iv = iv(:).';
    
    for ii = iv
        midpoint = mean(corners(edgecorners(ii,1:2),:));
        text(midpoint(1),midpoint(2),midpoint(3),int2str(ii));
    end
end

if plotconumbers
    for ii = 1:ncorners
        if plotcoCADnumbers == 1
            text(corners(ii,1),corners(ii,2),corners(ii,3),int2str(cornernumbers(ii)))    
        else
            text(corners(ii,1),corners(ii,2),corners(ii,3),int2str(ii))    
        end
    end
    
end

if plotsources == 1
    if isempty(plotoptions2)
        plot3(sources(:,1),sources(:,2),sources(:,3),'*');
    else
        if nplotoptions2 == 1
            plot3(sources(plotoptions2(1),1),sources(plotoptions2(1),2),sources(plotoptions2(1),3),'*')
        else
            plot3(sources(plotoptions2(1,1),1),sources(plotoptions2(1,1),2),sources(plotoptions2(1,1),3),'k*')
            plot3(sources(plotoptions2(2:nplotoptions2,1),1),sources(plotoptions2(2:nplotoptions2,1),2),sources(plotoptions2(2:nplotoptions2,1),3),'*')            
        end
    end
end

if plotreceivers == 1
    if isempty(plotoptions2)
        plot3(receivers(:,1),receivers(:,2),receivers(:,3),'ro');
    else
        if nplotoptions2 == 1
            plot3(receivers(plotoptions2(2),1),receivers(plotoptions2(2),2),receivers(plotoptions2(2),3),'ro');
        else
            plot3(receivers(plotoptions2(1,2),1),receivers(plotoptions2(1,2),2),receivers(plotoptions2(1,2),3),'ko');
            plot3(receivers(plotoptions2(2:nplotoptions2,2),1),receivers(plotoptions2(2:nplotoptions2,2),2),receivers(plotoptions2(2:nplotoptions2,2),3),'ro');
            
        end
    end
end
