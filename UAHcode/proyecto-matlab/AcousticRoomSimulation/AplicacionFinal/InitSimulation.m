function infos=InitSimulation(ENVfile,Filestem,Filepath, ...
    int_or_ext_model,open_or_closed_model,EDcalcmethod,specorder,difforder, ...
    SHOWTEXT,calcirs,Rstart,FSAMP,CAIR,RHOAIR,separacion_point_edges)

%%ENVfile: path de la ubicación y archivo .env que contiene la definición del modelo (terminar en separador \ o /)
%%Filestem: palabra por la que se desea que comiencen los archivos generados en la simulación
%%Filepath: path donde se desea que se guarden los archivos generados en la simulación (terminar en separador \ o /)
%%open_or_closed_model: 'open' para modelo abierto y 'closed' para modelo cerrado
%%int_or_ext_model: 'int' para modelo interior, 'ext' para modelo exterior
%%EDcalcmethod: 'n' para método Svensson, 'v' para método Vanderkooy, 'k' para aproximación Kirchoff
%%specorder: orden de reflexiones especulares
%%difforder: orden de difracciones, tiene que ser menor o igual que specorder, a menos que specorder sea igual o mayor que 3, en ese caso difforder valdrá 0
%%SHOWTEXT: cantidad de texto que se muestra en la simulación
%%%calcirs:
%%%Rstart:
%%FSAMP,CAIR,RHOAIR solo sirve si se calcula IR


% NO MODIFICABLES
infos.directsound = 1;
infos.elemsize = [1 0.5];
infos.nedgesubs = 2;
infos.calcpaths = 1;
infos.calctfs = 0;
infos.skipcorners = 1000000;


% PARÁMETROS:

[rutapath,name,extension]=fileparts(ENVfile);
if strcmp(extension,'.env')
    ruta=[rutapath filesep];
    env2matlab(ruta,ENVfile);
    matlab2cad();
    path_cad=pwd;   %CUIDADO!!!
    CADfile = [path_cad filesep 'environment.cad'];
    infos.CADfile=CADfile;
else
    infos.CADfile=ENVfile;
end

infos.Filestem=Filestem;
infos.Filepath=Filepath;
infos.int_or_ext_model=int_or_ext_model;
infos.open_or_closed_model=open_or_closed_model;
infos.EDcalcmethod=EDcalcmethod;
infos.specorder=specorder;
infos.difforder=difforder;

infos.calcirs=calcirs;
infos.Rstart=Rstart;
infos.FSAMP=FSAMP;
infos.CAIR=CAIR;
infos.RHOAIR=RHOAIR;
infos.SHOWTEXT=SHOWTEXT;

infos.separacion_point_edges=separacion_point_edges;

