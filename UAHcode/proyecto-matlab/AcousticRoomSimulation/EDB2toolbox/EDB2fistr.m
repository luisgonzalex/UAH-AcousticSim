function iline = EDB2fistr(longtextstring,iline,linestartvec,CRpositionvec,textstring)
% EDB2fistr - Reads a text file, line by line, until a given textstring is found.
% Reads line by line in a long textstring 'longtextstring'
% until the text in 'textstring' has been found.
%
% Input parameters:
%	longtextstring		The long textstring ot search in
%	iline				The line number to start searching at
%	linestartvec		A vector with the starting positions of all lines
%	CRpositionvec		A vector with the ending positions of all lines (inlcudes CR)
%	textstring			The textstring to search for.
% Output parameters:
%   iline				The line number which contains the text in textstring.
%						If the text doesn't occur, the value -1 is returned.
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
% Peter Svensson (svensson@iet.ntnu.no) 20001010
%
% iline = EDB2fistr(longtextstring,iline,linestartvec,CRpositionvec,textstring);

nlines = length(linestartvec);

foundtext = 0;
foundend = 0;
nchar = length(textstring);

iline = iline - 1;
while foundtext == 0 & foundend == 0
	iline = iline + 1;
	if iline > nlines
		foundend = 1;
		iline = -1;
	else
		Str = longtextstring( linestartvec(iline):CRpositionvec(iline)-1 );	
    	if length(Str) >= nchar, 
    		if Str(1:nchar) == textstring
    			foundtext = 1;
    		end   
   		end
	end
end
