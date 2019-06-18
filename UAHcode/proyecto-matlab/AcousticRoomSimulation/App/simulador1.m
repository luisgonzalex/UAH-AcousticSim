function varargout = simulador1(varargin)
% SIMULADOR1 MATLAB code for simulador1.fig
%    Welcome page to the toolbox 

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @simulador1_OpeningFcn, ...
                   'gui_OutputFcn',  @simulador1_OutputFcn, ...
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


% --- Executes just before simulador1 is made visible.
function simulador1_OpeningFcn(hObject, eventdata, handles, varargin)
imshow(imread('logoUAH.jpg'),'Parent',handles.axes7);
imshow(imread('logoEPS-UAH.jpg'),'Parent',handles.axes8);
imshow(imread('Logo_depeca_azul.jpg'),'Parent',handles.axes9);

% Choose default command line output for simulador1
handles.output = hObject;
% Update handles structure
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = simulador1_OutputFcn(hObject, eventdata, handles) 
% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in entrar.
function entrar_Callback(hObject, eventdata, handles)
simulador3_env;
delete(handles.output); % close this window
