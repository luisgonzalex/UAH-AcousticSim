function Bigtf = EDB2plotstackedTFs(firstirfile,nfiles,offset,nfft,typeofir,startsample,endsample)
% EDB2plotstackedTFs - Plots a number of TFs in a stacked way.
%
% Input parameters:
%   firstirfile     The data file containing the first IR file
%   nfiles          The number of files
%   offset          The numerical offset for each file [in dB]
%   nfft            The fft size
%   typeofir        't' (default), 'g' (geom), 'f' (direct sound), 'd'
%                   (diffracted)
%   startsample     (optional) The first sample number that should be
%                   plotted for each IR (after filtering). Default: 1.
%   endsample       (optional) The last sample number that should be
%                   plotted for each IR (after filtering). Default: the
%                   last sample for the longest of all IRs.
%
% Output parameters:
%   Bigtf           Matrix, [nfiles,nlength], of the TFs in the files, one
%                   per row.
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
% Peter Svensson (svensson@iet.ntnu.no) 20050116
%
% EDB2plotstackedTFs(firstirfile,nfiles,offset,nfft,typeofirs,startsample,endsample);

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

if nargin < 5
    typeofir = 't';    
else
    typeofir = lower(typeofir(1));    
end

if nargin < 6
    startsample = 1;    
end
if nargin < 7
    endsample = -1;    
end

filestem = EDB2strpend(firstirfile,'_ir');
iv = find(filestem=='_');
iv = iv(length(iv));
firstnumber = str2num(filestem(iv+1:length(filestem)));
filestem = filestem(1:iv)

%--------------------------------------------------------------
% Read the files

offsetvec = ([1:nfiles].'-1)*offset;

disp('   Loading files')
for ii = nfiles:-1:1
    irfile = [filestem,int2str(ii+firstnumber-1),'_ir.mat'];
    eval(['load ',irfilepath,irfile])
    irtot = real(irtot);
    if typeofir == 'f'
        irtot = irdirect;
    elseif typeofir == 'g'
        irtot = irgeom + irdirect;
    elseif typeofir == 's'
        irtot = irgeom;
    elseif typeofir == 'd'
        irtot = real(irdiff);
    end
        
    if isempty(irtot)
        irtot = [0 0].';    
    end
    
    if endsample == -1
        irtot = irtot(startsample:end);
    else
        if endsample <= length(irtot)
            irtot = irtot(startsample:endsample);        
        else
            irtot = irtot(startsample:end);        
        end
    end
    F = fft(full(irtot),nfft);
    L = 20*log10(abs(F(1:nfft/2)));
    
    if ii == nfiles
        Bigtf = zeros(nfiles,nfft/2);
        BigL = zeros(nfiles,nfft/2);
        [slask,nbigtf] = size(Bigtf);
        BigL = BigL + offsetvec(:,ones(1,nbigtf));
        Bigtf(ii,:) = (F(1:nfft/2)).';
        BigL(ii,:) = BigL(ii,:) + (L).';
    else
        Bigtf(ii,:) = (F(1:nfft/2)).';
        BigL(ii,:) = BigL(ii,:) + (L).';
    end
    
end

fvec = FSAMP/nfft*[0:nfft/2-1];
semilogx(fvec,BigL.','r')
