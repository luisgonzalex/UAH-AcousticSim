function outval = EDB2myvalueinput(textstr,nvalues,minval,maxval)
% EDB2myvalueinput - Prints a text on the screen and waits for a numerical to be input.
% Prints the text in textstr, and waits for a numerical to be input.
% It is checked if these are within the range [minval,maxval].
% If one or both is not specified, or given as [], an open range
% is used. It is also checked if exactly nvalues are input. The
% values should be input with spaces inbetween. If a negative
% value of nvalues is input, any number will be accepted.
%
% Input parameters:
%   textstr             A text str which will be printed on the screen.
%   nvalues             The number of numerical values that should be read in from
%                       the keyboard.
%   minval (optional)   The lowest value that should be allowed.
%   maxval (optional)   The highest value that should be allowed.
%
% Output parameters:
%   outval              The nvalues numerical values that should be read in.
%
% Uses the function EDB2extrnums
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
% Peter Svensson (svensson@iet.ntnu.no) 20030602
%
% outval = EDB2myvalueinput(textstr,nvalues,minval,maxval);

if nargin < 4
	maxval = 1e9;
	if nargin < 3
		minval = -1e9;
	end
end
if isempty(maxval)
	maxval = 1e9;	
end
if isempty(minval)
	minval = -1e9;	
end
if maxval < minval
	error(['ERROR: minval > maxval'])
end

foundOKvalue = 0;
while foundOKvalue == 0
	outval = EDB2extrnums( input(textstr,'s') );
	if length(outval) ~= nvalues & nvalues >= 0
		disp('Wrong number of values!');
	else
		if prod( double(outval >= minval) ) & prod( double(outval <= maxval) )
			foundOKvalue = 1;
		else
			disp('One value is outside the allowed range!');
		end
	end
end
