% BENCH_Lsp_Kessel.m created on Feb. 18, 2003 by Peter Svensson
 
global FSAMP CAIR RHOAIR SHOWTEXT
FSAMP = 96000;
CAIR = 344;
RHOAIR = 1.21;
SHOWTEXT = 3;	% This parameter determines how much text will be printed on the
 				% screen.

% Filepath: you must specify where the output files should be stored.
%           Note that this example uses UNIX folders but for Windows you 
%           will need something like:
%           Filepath = 'C:\folder\';
Filepath = 'D:\Dropbox\TFM\REPOSITORIO\proyecto-matlab\AcousticRoomSimulation\EDB2toolbox\Salidas';

% Filestem: All your output/result files will get names that start with
%           the name you specify here.
Filestem = 'mioorigin';
 
% CADfile:  This should be the CAD file with the input data.
%			You must include the complete path name.
CADfile = './sinmesa_CAD.CAD';
open_or_closed_model = 'closed';    % Specify if the model is open or closed
int_or_ext_model = 'int';         % Specify if you are interested in the interior or exterior
                                  % of the model.

% Calculation parameters
EDcalcmethod = 'n';     % 'n' = the new method by Svensson et al
directsound = 1;        % 1 if you want to include the direct sound, 0 if you don't
specorder = 2;          % The highest number of specular reflections.
difforder = 2;          % The highest number of edge diffractions. Note that specorder must be >= difforder.
elemsize = [1 0.5];     % This is an accuracy parameter for each order of edge diffraction.
                        % The vector must start with 1. The value 0.5 decides how
                        % small edge elements will be used for second order
                        % diffraction. A higher number gives more accurate
                        % results but takes much longer time.
                        % For third-order: elemsize = [1 0.5 0.25] for
                        % instance.
nedgesubs = 2;          % This is a parameter that decides how many parts each edge is divided into for the 
                        % part-visibility check of edges. A higher number
                        % gets more accurate but takes much longer time.
                        % Minimum = 2.
calcpaths = 1;          % If you want to run the first calculation step (find the paths), set the value 1.
                        % If you have run this part earlier and just want
                        % to change some setting for the second calculation
                        % step, then you can re-use the first step and set
                        % this value to 0.
calcirs = 0;            % If you want to run the second calculation step (construct IRs), set the value 1,
                        % otherwise 0.
 
calctfs = 0; %añadido por mi
% Sources and receivers
sources = [0 -2 0]; % Source coordinates are specified as (example with two sources)
                        %   sources = [0 0 0.0001;0.1 0 0.0001];  
receivers = [0 0 0];   % Receiver coordinates are specified as (example with four receivers)
                        %   receivers = [0 0 10;10*sin(pi/6) 0 10*cos(pi/6);10*sin(2*pi/6) 0 10*cos(2*pi/6);10*sin(3*pi/6) 0 10*cos(3*pi/6)]; 
                        % Also possible:
                        %   fivec = [0:5:180].'*pi/180;
                        %   radius = 10;
                        %   receivers = [radius*sin(fivec) zeros(size(fivec)) radius*cos(fivec)]; 
                        
% Extra parameters
skipcorners = 1000000;  % If you want to not include parts of your model, then you can exclude all corners
                        % with corner numbers higher than this.
Rstart = 9.9;           % All impulse responses will have a lot of zeros at the start, if the distance from source
                        % to reciever is long, and the sampling frequency
                        % is high. By setting Rstart to some non-zero
                        % value, the impulse responses will all start at
                        % the time that corresponds to this distance in
                        % meters. It is important to set this longer than
                        % the minimum that can ever happen since there
                        % might be some cryptic error message otherwise.
