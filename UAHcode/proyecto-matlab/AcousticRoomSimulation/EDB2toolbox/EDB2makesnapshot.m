function snapshot = EDB2makesnapshot(firstirfile,xvec,yvec,timestep,filtir,maxvalue,viewpoint,...
    outlinecorners,noutlinepoints,sourcepoint,typeofir,windowsize,ampexp)
% EDB2makesnapshot - Makes a snapshot of a sound field based on IRs.
%
% Input parameters:
%   firstirfile     The first of the IR files
%   xvec, yvec      The ranges of x- and y-values of the receiver positions
%   timestep        The time step to plot
%   filtir          A window to filter with
%   maxvalue        [1,2] The minimum and maximum amplitude in the plot
%   viewpoint       The point to view from [x y z]
%   outlinecorners  A matrix of corners which will be plotted as an outline
%                   in the form [x1 y1;x2 y2;x3 y3;...]
%                   If the matrix has more columns than two, each pair of
%                   column will be plotted as a separate outline, and they
%                   will not be tied together.
%   noutlinepoints  A list with the number of points in each column pair.
%   sourcepoint     The [x y z] coordinates of the source.
%   typeofir        't' (default), 'g' (geom), 'f' (direct sound), 'd'
%                   (diffracted), 'c' (direct sound + geom)
%   windowsize      's', 'm' or 'l' or 'x'
%   ampexp          (optional) If a value is given here, the signal will be
%                   plotted as abs(pressure)^ampexp, so that if ampexp =
%                   0.5, sqrt(abs(pressure)) will be plotted. If no value
%                   is given, or the value zero is given, the linear pressure will be plotted.
%
% Output parameters:
%   snapshot           A snapshot, in the format which can be played by the
%                   command surf(xvec,yvec,snapshot).
%
% Uses functions EDB2strpend
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
% Peter Svensson (svensson@iet.ntnu.no) 20051013
%
% M = EDB2snapshot(firstirfile,xvec,yvec,starttimestep,ntimesteps,filtir,maxvalue,viewpoint
%    outlinecorners,noutlinepoints,sourcepoint,typeofir,windowsize,ampexp);

disp('*************************************')
disp('*')
disp('*   Creating a snapshot....')
disp(' ')


%--------------------------------------------------------------
% Extract the filename stem

if nargin == 0
	[firstirfile,irfilepath] = uigetfile('*ir.mat','Please select the first irfile');
    [irfilepath,temp1,temp2] = fileparts(irfilepath);
	if ~isstr(firstirfile) | isempty(firstirfile)
		return
	end
else
    [irfilepath,firstirfile,fileext] = fileparts(firstirfile);
    irfilepath = [irfilepath,filesep];
    firstirfile = [firstirfile,fileext];
end

if nargin < 11
    typeofir = 't';    
    windowsize = 's';
    ampexp = 0;
else
    typeofir = lower(typeofir(1));    
    if nargin < 12
        windowsize = 's';        
        ampexp = 0;
    else
        windowsize = lower(windowsize(1));                  
        if nargin < 13
            ampexp = 0;
        else
            if ampexp < 0
                error('ERROR: ampexp must be larger than, or equal to, zero')    
            end
        end
    end
end

filestem = EDB2strpend(firstirfile,'_ir');
iv = find(filestem=='_');
iv = iv(length(iv));
firstnumber = str2num(filestem(iv+1:length(filestem)));
filestem = filestem(1:iv);

nx = length(xvec);
ny = length(yvec);
nfiles = nx*ny;

loadcommandbase = ['load ',irfilepath,filestem];
    
%--------------------------------------------------------------
% Read the files, filter and store

Bigir = zeros(nfiles,1);

disp('   Loading the files, starting with:')
for ii = 1:nfiles
    eval([loadcommandbase,sprintf('%.0f',ii),'_ir.mat'])
    if typeofir == 'f'
        irtot = irdirect;    
    elseif typeofir == 'g'
        irtot = irgeom;
    elseif typeofir == 'c'
        irtot = irgeom + irdirect;
    elseif typeofir == 'd',     
        irtot = irdiff;
    end
        
    if any(irtot)
        
        if length(irtot) > timestep
            irtot = irtot(1:timestep);    
        elseif length(irtot) < timestep
            irtot = [irtot;zeros(timestep-length(irtot),1)];            
        end

        irsamplevalue = sum(irtot([timestep:-1:timestep-length(filtir)+1]).*filtir);
    else
        irsamplevalue = 0;
    end
    
    if ampexp == 0
        Bigir(ii) = irsamplevalue;
    else
        Bigir(ii) = (abs(irsamplevalue)).^ampexp;
    end
    
end

disp('   ... files are loaded.')


if length(maxvalue) < 2
    maxvalue = [maxvalue(1) max(Bigir)];    
end

%--------------------------------------------------------------
% Prepare the outline plotting

opengl neverselect

snapshot = reshape(Bigir,ny,nx);

[noutlinerows,ncolumns] = size(outlinecorners);
noutlines = floor(ncolumns/2);
if noutlines == 1
    outlinecorners = [outlinecorners;outlinecorners(1,:)];
    noutlinepoints = noutlinepoints+1;
end

axisvalues = [min(xvec) max(xvec) min(yvec) max(yvec) maxvalue(1) maxvalue(2)];
figure(1)
clf

xyaspectratio = abs(max(xvec)-min(xvec))/abs(max(yvec)-min(yvec));

windowheight = 675;
windowwidth = windowheight*xyaspectratio;

if windowwidth > 800
    scaledownfactor = windowwidth/800;
    windowwidth = windowwidth/scaledownfactor;
    windowheight = windowheight/scaledownfactor;
end

windowpos = [380 80];
if windowsize == 's'
    windowwidth = windowwidth/2;
    windowheight = windowheight/2;
elseif windowsize == 'm'
    windowwidth = windowwidth/1.41;
    windowheight = windowheight/1.41;
elseif windowsize == 'x'
    windowwidth = windowwidth*1.41;
    windowheight = windowheight*1.41;
    windowpos(1) = 100;
end

set(1,'Position',[windowpos(1:2)   windowwidth   windowheight])

snapshot = reshape(Bigir,ny,nx);
surf(xvec,yvec,snapshot);

shading interp
colormap('jet')    
caxis([maxvalue(1) maxvalue(2)])
axis(axisvalues);
axis off
view(viewpoint)

if noutlines > 0
	for kk = 1:noutlines
		for ll = 1:noutlinepoints(kk)-1
			H = line([outlinecorners(ll,(kk-1)*2+1) outlinecorners(ll+1,(kk-1)*2+1)],[outlinecorners(ll,(kk-1)*2+2) outlinecorners(ll+1,(kk-1)*2+2)],[0 0]);
			set(H,'LineWidth',4);
			set(H,'Color',[1 1 0]);
			H = line([outlinecorners(ll,(kk-1)*2+1) outlinecorners(ll+1,(kk-1)*2+1)],[outlinecorners(ll,(kk-1)*2+2) outlinecorners(ll+1,(kk-1)*2+2)],2*[maxvalue(2) maxvalue(2)]);
			set(H,'LineWidth',4);
			set(H,'Color',[1 1 0]);
			H = line([outlinecorners(ll,(kk-1)*2+1) outlinecorners(ll,(kk-1)*2+1)],[outlinecorners(ll,(kk-1)*2+2) outlinecorners(ll,(kk-1)*2+2)],[0 2*maxvalue(2)]);        
			set(H,'LineWidth',4);
			set(H,'Color',[1 1 0]);
		end
		H = line([outlinecorners(ll+1,(kk-1)*2+1) outlinecorners(ll+1,(kk-1)*2+1)],[outlinecorners(ll+1,(kk-1)*2+2) outlinecorners(ll+1,(kk-1)*2+2)],[0 2*maxvalue(2)]);        
			set(H,'LineWidth',4);
			set(H,'Color',[1 1 0]);
	end
end
if ~isempty(sourcepoint)
	H = line([sourcepoint(1)*0.95 sourcepoint(1)*1.05],[sourcepoint(2)*0.95 sourcepoint(2)*1.05],[maxvalue(2) maxvalue(2)]);
			set(H,'LineWidth',4);
			set(H,'Color',[1 0 0]);        
	H = line([sourcepoint(1)*0.95 sourcepoint(1)*1.05],[sourcepoint(2)*1.05 sourcepoint(2)*0.95],[maxvalue(2) maxvalue(2)]);
			set(H,'LineWidth',4);
			set(H,'Color',[1 0 0]);        
end
