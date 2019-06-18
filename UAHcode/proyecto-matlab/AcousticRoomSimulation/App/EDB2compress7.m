function [reftoshort,sho1,sho2,sho3,sho4,sho5,sho6,sho7,examplecombination] = ...
EDB2compress7(Big1,Big2,Big3,Big4,Big5,Big6,Big7)
% EDB2compress7 - 
% This function looks for combinations of N edges-to-N edges that have identical
% values for 7 [N,N] matrices. The values are rounded to a relative
% accuracy of 1e-12.
%
% Input parameters:
%   Big1, Big2,...,Big 7    Seven matrices, [N,N], of data.
%
% Output parameters:
%   reftoshort              Matrix, [N,N], of pointer values to seven short
%                           lists.
%   sho1, sho2,...,sho7     Seven shortlists, [nshort,1], with values that
%                           are taken from the seven input matrices.
%   examplecombination      Vector of pointers from the short lists, back
%                           to the input matrices.
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
% Peter Svensson (svensson@iet.ntnu.no) 20030503
%
% [reftoshort,sho1,sho2,sho3,sho4,sho5,sho6,sho7,examplecombination] = ...
% EDB2compress7(Big1,Big2,Big3,Big4,Big5,Big6,Big7)

accuracy = 1e-12;

[n1,n2] = size(Big1);

Big1 = reshape(Big1,n1*n2,1);
Big2 = reshape(Big2,n1*n2,1);
Big3 = reshape(Big3,n1*n2,1);
Big4 = reshape(Big4,n1*n2,1);
Big5 = reshape(Big5,n1*n2,1);
Big6 = reshape(Big6,n1*n2,1);
Big7 = reshape(Big7,n1*n2,1);

maxamp = max(abs(Big1));
if maxamp ~= 0
	Big1 = round( Big1/maxamp/accuracy)*accuracy*maxamp;
end
maxamp = max(abs(Big2));
if maxamp ~= 0
	Big2 = round( Big2/maxamp/accuracy)*accuracy*maxamp;
end
maxamp = max(abs(Big3));
if maxamp ~= 0
	Big3 = round( Big3/maxamp/accuracy)*accuracy*maxamp;
end
maxamp = max(abs(Big4));
if maxamp ~= 0
	Big4 = round( Big4/maxamp/accuracy)*accuracy*maxamp;
end
maxamp = max(abs(Big5));
if maxamp ~= 0
	Big5 = round( Big5/maxamp/accuracy)*accuracy*maxamp;
end
maxamp = max(abs(Big6));
if maxamp ~= 0
	Big6 = round( Big6/maxamp/accuracy)*accuracy*maxamp;
end
maxamp = max(abs(Big7));
if maxamp ~= 0
	Big7 = round( Big7/maxamp/accuracy)*accuracy*maxamp;
end

%---------------------------------------
% sort

[Big1,sortvec] = sort(Big1);
Big2 = Big2(sortvec);
Big3 = Big3(sortvec);
Big4 = Big4(sortvec);
Big5 = Big5(sortvec);
Big6 = Big6(sortvec);
Big7 = Big7(sortvec);

[Big2,sortvec2] = sort(Big2);
Big1 = Big1(sortvec2);
Big3 = Big3(sortvec2);
Big4 = Big4(sortvec2);
Big5 = Big5(sortvec2);
Big6 = Big6(sortvec2);
Big7 = Big7(sortvec2);
sortvec = sortvec(sortvec2);

[Big3,sortvec2] = sort(Big3);
Big1 = Big1(sortvec2);
Big2 = Big2(sortvec2);
Big4 = Big4(sortvec2);
Big5 = Big5(sortvec2);
Big6 = Big6(sortvec2);
Big7 = Big7(sortvec2);
sortvec = sortvec(sortvec2);

[Big4,sortvec2] = sort(Big4);
Big1 = Big1(sortvec2);
Big2 = Big2(sortvec2);
Big3 = Big3(sortvec2);
Big5 = Big5(sortvec2);
Big6 = Big6(sortvec2);
Big7 = Big7(sortvec2);
sortvec = sortvec(sortvec2);

[Big5,sortvec2] = sort(Big5);
Big1 = Big1(sortvec2);
Big2 = Big2(sortvec2);
Big3 = Big3(sortvec2);
Big4 = Big4(sortvec2);
Big6 = Big6(sortvec2);
Big7 = Big7(sortvec2);
sortvec = sortvec(sortvec2);

[Big6,sortvec2] = sort(Big6);
Big1 = Big1(sortvec2);
Big2 = Big2(sortvec2);
Big3 = Big3(sortvec2);
Big4 = Big4(sortvec2);
Big5 = Big5(sortvec2);
Big7 = Big7(sortvec2);
sortvec = sortvec(sortvec2);

[Big7,sortvec2] = sort(Big7);
Big1 = Big1(sortvec2);
Big2 = Big2(sortvec2);
Big3 = Big3(sortvec2);
Big4 = Big4(sortvec2);
Big5 = Big5(sortvec2);
Big6 = Big6(sortvec2);
sortvec = sortvec(sortvec2);

ivec1 = find( abs(diff(Big1)) > accuracy );
ivec2 = find( abs(diff(Big2)) > accuracy );
ivec3 = find( abs(diff(Big3)) > accuracy );
ivec4 = find( abs(diff(Big4)) > accuracy );
ivec5= find( abs(diff(Big5)) > accuracy );
ivec6= find( abs(diff(Big6)) > accuracy );
ivec7= find( abs(diff(Big7)) > accuracy );

ivec = [0;sort([ivec1;ivec2;ivec3;ivec4;ivec5;ivec6;ivec7])];

if length(ivec) > 1
	ivec = ivec( find(diff(ivec) ~= 0) + 1 );
	ivec = [1;ivec + 1;n1*n2+1];

	sho1     = Big1( ivec(1:length(ivec)-1));
	sho2     = Big2( ivec(1:length(ivec)-1));
	sho3     = Big3( ivec(1:length(ivec)-1));
	sho4     = Big4( ivec(1:length(ivec)-1));
	sho5     = Big5( ivec(1:length(ivec)-1));
	sho6     = Big6( ivec(1:length(ivec)-1));
	sho7	 = Big7( ivec(1:length(ivec)-1));

    if length(ivec) < 256
    	reftoshort = uint8(zeros( size(Big1) ));
    elseif length(ivec) < 65536
    	reftoshort = uint16(zeros( size(Big1) ));
    else
    	reftoshort = uint32(zeros( size(Big1) ));        
    end
	numberofoccurences = zeros( size(sho1) );
	examplecombination = zeros( size(sho1) );

	for ii = 1:length(ivec)-1
		sameblock = sortvec( [ivec(ii):ivec(ii+1)-1] );
		reftoshort(sameblock) = ii*ones(size(sameblock));
		numberofoccurences(ii) = length( sameblock );
		examplecombination(ii) = sameblock(1);
		patchcomb = sameblock(1);
	end

	reftoshort = reshape(reftoshort,n1,n2);

else % All IRs are the same
	% The extra zero is added to give the right dimensions of the output
	% vector when the shortlists are used.

	reftoshort = ones(n1,n2);
	sho1 = [Big1(1,1);0];
	sho2 = [Big2(1,1);0];
	sho3 = [Big3(1,1);0];
	sho4 = [Big4(1,1);0];
	sho5 = [Big5(1,1);0];
	sho6 = [Big6(1,1);0];
	sho7 = [Big7(1,1);0];
	examplecombination = 1;	
	
end
