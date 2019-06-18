function varargout = simulador5_env(varargin)
% SIMULADOR5_ENV MATLAB code for simulador5_env.fig
%  Read parameters, store in info.mat

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @simulador5_env_OpeningFcn, ...
                   'gui_OutputFcn',  @simulador5_env_OutputFcn, ...
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


% --- Executes just before simulador5_env is made visible.
function simulador5_env_OpeningFcn(hObject, eventdata, handles, varargin)

% Choose default command line output for simulador5_env
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

global SHOWTEXT
%FSAMP CAIR RHOAIR
% asignar valor por defecto, por si se dejan por defecto sin modificar
open_or_closed_model= 'open';
%FSAMP=96000;
%CAIR = 344;
%RHOAIR = 1.21;
SHOWTEXT = 4;
specorder = 2;
difforder = 1;
int_or_ext_model = 'int';
EDcalcmethod = 'n';
% FSAMP CAIR RHOAIR 
save info open_or_closed_model SHOWTEXT specorder difforder int_or_ext_model EDcalcmethod '-append'



% --- Outputs from this function are returned to the command line.
function varargout = simulador5_env_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;



% function fsamp_Callback(hObject, eventdata, handles)
% 
% FSAMP=str2num(get(hObject,'String'));
% save info FSAMP '-append'

% --- Executes during object creation, after setting all properties.
% function fsamp_CreateFcn(hObject, eventdata, handles)
% 
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end



% function cair_Callback(hObject, eventdata, handles)
% 
% CAIR=str2num(get(hObject,'String'));
% save info CAIR '-append'

% --- Executes during object creation, after setting all properties.
% function cair_CreateFcn(hObject, eventdata, handles)
% 
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end



% function rhoair_Callback(hObject, eventdata, handles)
% 
% RHOAIR=str2num(get(hObject,'String'));
% save info RHOAIR '-append'

% --- Executes during object creation, after setting all properties.
% function rhoair_CreateFcn(hObject, eventdata, handles)
% 
% if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
%     set(hObject,'BackgroundColor','white');
% end


% --- Executes on selection change in showtext.
function showtext_Callback(hObject, eventdata, handles)

SHOWTEXT = get(hObject,'Value');
save info SHOWTEXT '-append'
        
% --- Executes during object creation, after setting all properties.
function showtext_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in atras.
function atras_Callback(hObject, eventdata, handles)

simulador4_env;
delete(handles.output);  %cierra la ventana anterior

% --- Executes on button press in ayuda.
function ayuda_Callback(hObject, eventdata, handles)

ayuda_simulador5_env;


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


% --- Executes on button press in continuar.
function continuar_Callback(hObject, eventdata, handles)

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
calcirs = 0;            % If you want to run the second calculation step (construct IRs), set the value 1,
                        % otherwise 0.
calctfs = 0;            %añadido por mi
% sources y receivers en sources_receivers.mat
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
save info directsound elemsize nedgesubs calcpaths calcirs calctfs skipcorners Rstart '-append'
simulador6_env;
delete(handles.output);  %cierra la ventana anterior

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

% --- Executes when selected object is changed in diffr.
function diffr_SelectionChangeFcn(hObject, eventdata, handles)

if hObject == handles.svensson
     EDcalcmethod ='n';
 elseif hObject == handles.vanderkooy
      EDcalcmethod ='v';
 elseif hObject == handles.kirchoff
    EDcalcmethod ='k';
 end
save info EDcalcmethod '-append'
