function varargout = ayuda_simulador5_env(varargin)
% AYUDA_SIMULADOR5_ENV MATLAB code for ayuda_simulador5_env.fig   

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ayuda_simulador5_env_OpeningFcn, ...
                   'gui_OutputFcn',  @ayuda_simulador5_env_OutputFcn, ...
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


% --- Executes just before ayuda_simulador5_env is made visible.
function ayuda_simulador5_env_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for ayuda_simulador5_env
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = ayuda_simulador5_env_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in cerrar.
function cerrar_Callback(hObject, eventdata, handles)
delete(handles.output);  %cierra la ventana 
