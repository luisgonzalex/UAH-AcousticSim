function varargout = simulador6_env(varargin)
% SIMULADOR6_ENV MATLAB code for simulador6_env.fig


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @simulador6_env_OpeningFcn, ...
                   'gui_OutputFcn',  @simulador6_env_OutputFcn, ...
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


% --- Executes just before simulador6_env is made visible.
function simulador6_env_OpeningFcn(hObject, eventdata, handles, varargin)
set(handles.representacion,'Enable', 'off'); %deshabilita boton
set(handles.resultados,'Enable', 'off'); %deshabilita boton

% Choose default command line output for simulador6_env
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = simulador6_env_OutputFcn(hObject, eventdata, handles) 


% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in simular.
function simular_Callback(hObject, eventdata, handles)
%SI RESULTADOS.TXT YA ESTA CREADO NOSE SI VUELVE A ESCRIBIR ENCIMA...
diary('resultados.txt');
EDB2main;
diary off;
set(handles.representacion,'Enable', 'on'); %habilita boton
set(handles.resultados,'Enable', 'on'); %habilita boton

% load('info.mat');
% representation_paths(Filepath,Filestem);

% --- Executes on button press in atras.
function atras_Callback(hObject, eventdata, handles)

simulador5_env;
delete(handles.output);  %cierra la ventana anterior

% --- Executes on button press in representacion.
function representacion_Callback(hObject, eventdata, handles)

simulador7_env;
delete(handles.output);  %cierra la ventana anterior

% representation_room;
% representation_src_rcv;
% load('info.mat');
% representation_paths(Filepath,Filestem);


% --- Executes on button press in resultados.
function resultados_Callback(hObject, eventdata, handles)
!notepad resultados.txt
