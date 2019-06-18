%API

%infos=InitSimulation('D:\Dropbox\TFM\REPOSITORIO\proyecto-matlab\AcousticRoomSimulation\IdiapRoom\environment.env','apiprueba','C:\Users\Ruth\Documents\borrar\', ...
%'int','closed','n',2,1, ...
%4,0,1.3,96000,344,1.21,1);
infos=InitSimulation('D:\Dropbox\TFM\REPOSITORIO\proyecto-matlab\AcousticRoomSimulation\App\Ispace.cad','apiispace','C:\Users\Ruth\Documents\borrar\', ...
'int','closed','n',2,2, ...
4,0,1.3,96000,344,1.21,1);
% sources=[0 -2 0.1];
% receivers=[0 0 -0.4];
sources=[2 2 2];
receivers=[6.5 5 2];
[raytype,paths]=getpaths(sources,receivers,infos);



