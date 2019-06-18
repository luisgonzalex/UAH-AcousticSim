function varargout = simulador4_env(varargin)
% SIMULADOR4_ENV MATLAB code for simulador4_env.fig
%  Generate  sources_receivers.mat with the number of sources and receivers
%  and the coordinates:
% nsources: number of sources
% nreceivers: number of receivers
% sources: matrix [nsources,3] with the coordinate x, y, z
% receivers: matrix [nreceivers,3] with the coordinate x, y, z

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @simulador4_env_OpeningFcn, ...
                   'gui_OutputFcn',  @simulador4_env_OutputFcn, ...
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


% --- Executes just before simulador4_env is made visible.
function simulador4_env_OpeningFcn(hObject, eventdata, handles, varargin)


% Choose default command line output for simulador4_env
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);



% --- Outputs from this function are returned to the command line.
function varargout = simulador4_env_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;



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


% --- Executes on button press in atras.
function atras_Callback(hObject, eventdata, handles)

simulador3_env;
delete(handles.output);  %cierra la ventana anterior

% --- Executes on button press in continuar.
function continuar_Callback(hObject, eventdata, handles)

nsources=handles.nsources;
nreceivers=handles.nreceivers;
sources=handles.sources;
receivers=handles.receivers;
save sources_receivers nsources nreceivers sources receivers
simulador5_env;
delete(handles.output);  %cierra la ventana anterior


function coor_source_Callback(hObject, eventdata, handles)

sources=str2num(get(hObject,'String'));
nsources=handles.nsources;
if length(sources) > nsources
    for i=1:nsources
        coor_src1(i,:)=sources(i,:);
        handles.sources=coor_src1;
    end
end
if length(sources) < nsources
    nsources=length(sources);
end
if length(sources) == nsources
   handles.sources=sources; 
end
% si en lugar de tres coordenadas se meten 2?? rellenar con 0s....(sin
% hacer)
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
if length(receivers) > nreceivers
    for i=1:nreceivers
        coor_receiver1(i,:)=receivers(i,:);
        handles.receivers=coor_receiver1;
    end
end
if length(receivers) < nreceivers
    nreceivers=length(receivers);
end
if length(receivers) == nreceivers
   handles.receivers=receivers; 
end
% si en lugar de tres coordenadas se meten 2?? rellenar con 0s....(sin
% hacer)
handles.nreceivers=nreceivers;
guidata(hObject,handles);



% --- Executes during object creation, after setting all properties.
function coor_receiver_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in representation_mio.
function representation_mio_Callback(hObject, eventdata, handles)
nsources=handles.nsources;
nreceivers=handles.nreceivers;
sources=handles.sources;
receivers=handles.receivers;
%close Figure 1
figure(1)
representation_room;
representation_src_rcv(nsources,nreceivers,sources,receivers);
