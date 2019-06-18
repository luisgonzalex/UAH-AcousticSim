function varargout = ayuda_simulador2(varargin)
% AYUDA_SIMULADOR2 MATLAB code for ayuda_simulador2.fig   

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @ayuda_simulador2_OpeningFcn, ...
                   'gui_OutputFcn',  @ayuda_simulador2_OutputFcn, ...
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


% --- Executes just before ayuda_simulador2 is made visible.
function ayuda_simulador2_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for ayuda_simulador2
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% --- Outputs from this function are returned to the command line.
function varargout = ayuda_simulador2_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in cerrar.
function cerrar_Callback(hObject, eventdata, handles)
delete(handles.output);  %cierra la ventana 
