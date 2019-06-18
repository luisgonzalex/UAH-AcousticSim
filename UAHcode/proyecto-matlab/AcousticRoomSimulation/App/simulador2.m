function varargout = simulador2(varargin)
% SIMULADOR2 MATLAB code for simulador2.fig
% Select the option .env/.cad

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @simulador2_OpeningFcn, ...
                   'gui_OutputFcn',  @simulador2_OutputFcn, ...
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


% --- Executes just before simulador2 is made visible.
function simulador2_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for simulador2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);




% --- Outputs from this function are returned to the command line.
function varargout = simulador2_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in env.
function env_Callback(hObject, eventdata, handles)

simulador3_env;
delete(handles.output);  %cierra la ventana anterior

% --- Executes on button press in cad.
function cad_Callback(hObject, eventdata, handles)

simulador3_cad;
delete(handles.output);  %cierra la ventana anterior

% --- Executes on button press in atras1.
function atras1_Callback(hObject, eventdata, handles)

simulador1;
delete(handles.output);  %cierra la ventana anterior
