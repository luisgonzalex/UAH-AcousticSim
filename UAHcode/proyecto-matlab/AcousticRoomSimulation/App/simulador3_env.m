function varargout = simulador3_env(varargin)
% SIMULADOR3_ENV MATLAB code for simulador3_env.fig

% Select .env/.cad to simulate, write word of the files of simulation, select directory to save files of simulation
% Select number of sources and receivers and their coordinates
% Select simulation options:
        % MODEL: exterior/interior    open/closed
        % DIFFRACTION METHOD: Svensson/Vanderkooy/Kirchoff
        % TERMINATE: order maximum of specular reflect and edge diffraction
        % how many text display
        
% Generate:
            % info.mat with variables: (...) open_or_closed_model SHOWTEXT specorder difforder int_or_ext_model EDcalcmethod FSAMP CAIR RHOAIR
            % sources_receivers.mat with the number of sources and receivers and their coordinates in variables:
                    % nsources: number of sources
                    % nreceivers: number of receivers
                    % sources: matrix [nsources,3] with the coordinate x, y, z
                    % receivers: matrix [nreceivers,3] with the coordinate x, y, z
                    
% Use functions:
        % EDB2cross.m
        % env2matlab(ruta,envfile): convierte .env to Matlab struct.
        % matlab2cad.m : convierte Matlab struct to .CAD
        % EDB2main : llamada al programa principal de simulación EDT
        % representation_room(): representa el recinto a simular a partir de la estructura de matlab generada en env2matlab()
        % representation_src_rcv(): representa la fuentes y/o receptores indicados
        % ayuda_simulador5_env.m: nueva ventana con información de los parámetros de simulación
        % simulador7_env.m: ventana de visualización de resultados de simulación


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @simulador3_env_OpeningFcn, ...
                   'gui_OutputFcn',  @simulador3_env_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before simulador3_env is made visible.
function simulador3_env_OpeningFcn(hObject, eventdata, handles, varargin)
% Choose default command line output for simulador3_env
handles.output = hObject;

% PARAMETROS DE SIMULACIÓN:
global SHOWTEXT FSAMP CAIR RHOAIR  % FSAMP, CAIR, RHOAIR only for calculate IR
FSAMP=96000;
CAIR = 344;
RHOAIR = 1.21;
SHOWTEXT = 4;

% INITIALIZE VARIABLES:

calcirs = 0;
Rstart = 1.3;
int_or_ext_model = 'int';
open_or_closed_model= 'closed';
specorder = 2;
difforder = 1;
EDcalcmethod = 'n';
save info open_or_closed_model SHOWTEXT specorder difforder int_or_ext_model EDcalcmethod FSAMP CAIR RHOAIR calcirs Rstart

sources=0;
receivers=0;
handles.sources=sources;
handles.receivers=receivers;
nsources=1;
nreceivers=1;
handles.nsources=nsources;
handles.nreceivers=nreceivers;

input_source = 0;
input_receiver = 0;
input_model=0;
handles.input_model=input_model;
handles.input_source=input_source;
handles.input_receiver=input_receiver;

Filestem=handles.dt; % dt calculate at create function Filestem
handles.Filestem=Filestem;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = simulador3_env_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;

% --- Executes on button press in examinar.
function examinar_Callback(hObject, eventdata, handles)
[envfile,ruta] = uigetfile({'*.env;*.cad'},'Please select the .env or .cad file');
envfile = [ruta,envfile];
[rutapath,name,extension]=fileparts(envfile);
set(handles.path_env,'String',envfile);
input_model=1;
if strcmp(extension,'.env')
    extension=1;  % .env
else
    extension=2;  % .cad
end

handles.input_model=input_model;
handles.extension=extension;
handles.ruta=ruta;
handles.envfile=envfile;
guidata(hObject,handles)


% --- Executes on button press in examinar_carpeta.
function examinar_carpeta_Callback(hObject, eventdata, handles)
[direc] = uigetdir('','Choose the directory to save the outfiles');
set(handles.path_store,'string',direc);
Filepath = [direc filesep ];  
save info Filepath '-append' 

% --- Executes on button press in simular.
function simular_Callback(hObject, eventdata, handles)
extension=handles.extension;
nsources=handles.nsources;
nreceivers=handles.nreceivers;
sources=handles.sources;
receivers=handles.receivers;
save sources_receivers nsources nreceivers sources receivers

% VARIABLES to SIMULATE (no se pueden modificar en GUI):
directsound = 1;        % 1 if you want to include the direct sound, 0 if you don't
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
calctfs = 0;            % If you want construct TFs, set the value 1, otherwise 0. (COMENTADO EN EDB2main)
skipcorners = 1000000;  % If you want to not include parts of your model, then you can exclude all corners
                        % with corner numbers higher than this.

save info directsound elemsize nedgesubs calcpaths calctfs skipcorners '-append'

if extension == 1   % .ENV
    ruta=handles.ruta;
    envfile=handles.envfile;
   env2matlab(ruta,envfile); 
    matlab2cad();
    path_cad=pwd;
    name_cad=[filesep 'environment.cad'];
    CADfile = [path_cad name_cad];
    save info CADfile '-append'
    Filestem=handles.Filestem;
    save info Filestem '-append'  
            
elseif extension == 2  % .CAD
    envfile=handles.envfile;
    CADfile = envfile;
    save info CADfile '-append'
    Filestem=handles.Filestem;
    save info Filestem '-append'
    
end

diary('resultados.txt');
api=0;
EDB2main(api);
diary off;

simulador7_env;


% --- Executes during object creation, after setting all properties.
function path_store_CreateFcn(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function path_env_CreateFcn(hObject, eventdata, handles)

function Filestem_Callback(hObject, eventdata, handles)
Filestem = get(hObject,'string');
handles.Filestem=Filestem;
guidata(hObject,handles)


% --- Executes during object creation, after setting all properties.
function Filestem_CreateFcn(hObject, eventdata, handles)
dt=datestr(now,'yyyymmdd_HHMMSS');
set(hObject,'string',dt);
handles.dt=dt;
guidata(hObject,handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in visualizar.
function visualizar_Callback(hObject, eventdata, handles)
input_source=handles.input_source;
input_receiver=handles.input_receiver;
input_model=handles.input_model;
extension=handles.extension;
nsources=handles.nsources;
nreceivers=handles.nreceivers;
sources=handles.sources; 
receivers=handles.receivers; 

if extension == 1 % .env
    ruta=handles.ruta;
    envfile=handles.envfile;
    env2matlab(ruta,envfile);
elseif extension == 2 % .cad
    envfile=handles.envfile;
    CADfile = envfile;
    save info CADfile '-append'
    Filestem=handles.Filestem;
    save info Filestem '-append'
end

if input_model == 1
    figure;
    representation_room(extension); 
    axis equal;
end
if input_source == 1
    if input_receiver == 0
        representation_src_rcv(nsources,0,sources,0);
        axis equal;
    else
        representation_src_rcv(nsources,nreceivers,sources,receivers);
    end
elseif input_receiver == 1
        representation_src_rcv(0,nreceivers,0,receivers);
        axis equal;
end

guidata(hObject,handles)



% --- Executes on selection change in showtext.
function showtext_Callback(hObject, eventdata, handles)
SHOWTEXT = get(hObject,'Value');
save info SHOWTEXT '-append'

% --- Executes during object creation, after setting all properties.
function showtext_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in ayuda.
function ayuda_Callback(hObject, eventdata, handles)
ayuda_simulador5_env;

function nsources_Callback(hObject, eventdata, handles)
nsources=str2num(get(hObject,'String'));
handles.nsources=nsources;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function nsources_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function nreceivers_Callback(hObject, eventdata, handles)
nreceivers=str2num(get(hObject,'String'));
handles.nreceivers=nreceivers;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function nreceivers_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function coor_source_Callback(hObject, eventdata, handles)
sources=str2num(get(hObject,'String'));
nsources=handles.nsources;
if length(sources(:,1)) > nsources
    for i=1:nsources
        coor_src1(i,:)=sources(i,:);
        handles.sources=coor_src1;
    end
end
if length(sources(:,1)) < nsources
    nsources=length(sources);
end
if length(sources(:,1)) == nsources
   handles.sources=sources; 
end
input_source = 1;
handles.input_source=input_source;
handles.nsources=nsources;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function coor_source_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function coor_receiver_Callback(hObject, eventdata, handles)
receivers=str2num(get(hObject,'String'));
nreceivers=handles.nreceivers;
if length(receivers(:,1)) > nreceivers
    for i=1:nreceivers
        coor_receiver1(i,:)=receivers(i,:);
        handles.receivers=coor_receiver1;
    end
end
if length(receivers(:,1)) < nreceivers
    nreceivers=length(receivers);
end
if length(receivers(:,1)) == nreceivers
   handles.receivers=receivers; 
end
input_receiver = 1;
handles.input_receiver=input_receiver;
handles.nreceivers=nreceivers;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function coor_receiver_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function specorder_Callback(hObject, eventdata, handles)
specorder=str2num(get(hObject,'String'));
save info specorder '-append'


% --- Executes during object creation, after setting all properties.
function specorder_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function difforder_Callback(hObject, eventdata, handles)
difforder=str2num(get(hObject,'String'));
save info difforder '-append'


% --- Executes during object creation, after setting all properties.
function difforder_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes when selected object is changed in open_closed.
function open_closed_SelectionChangeFcn(hObject, eventdata, handles)

 if hObject == handles.model_open
    open_or_closed_model='open';
 elseif hObject == handles.model_closed
    open_or_closed_model='closed';
 end
 save info open_or_closed_model '-append'


% --- Executes when selected object is changed in int_ext.
function int_ext_SelectionChangeFcn(hObject, eventdata, handles)
if hObject == handles.model_int
    int_or_ext_model='int';
 elseif hObject == handles.model_ext
    int_or_ext_model='ext';
 end
save info int_or_ext_model '-append'



% --- Executes on button press in ir.
function ir_Callback(hObject, eventdata, handles)
calcirs=get(hObject,'Value');
save info calcirs '-append'




function val_rstart_Callback(hObject, eventdata, handles)
Rstart=str2num(get(hObject,'String'));
save info Rstart '-append'




% --- Executes during object creation, after setting all properties.
function val_rstart_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function freqmuestreo_Callback(hObject, eventdata, handles)
FSAMP = get(hObject,'Value');
save info FSAMP '-append'


% --- Executes during object creation, after setting all properties.
function freqmuestreo_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function velsonido_Callback(hObject, eventdata, handles)
CAIR = get(hObject,'Value');
save info CAIR '-append'


% --- Executes during object creation, after setting all properties.
function velsonido_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function densidadaire_Callback(hObject, eventdata, handles)
RHOAIR = get(hObject,'Value');
save info RHOAIR '-append'


% --- Executes during object creation, after setting all properties.
function densidadaire_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
