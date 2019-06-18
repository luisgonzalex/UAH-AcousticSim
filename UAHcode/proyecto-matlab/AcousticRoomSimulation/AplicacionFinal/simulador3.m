function varargout = simulador3(varargin)
% SIMULADOR3 MATLAB code for simulador3.fig

% Plot paths of simulation with differents options

% Use functions:
        % EDB2plotmodel()
        % EDB2extrnums()
        % representation_specularpaths()
        % representation_diffpaths()
        % representation_onepath()
        % representation_paths()
        
        % resultados.txt: resultados de la simulación generados en simulador2 en la llamada a EDB2main
        % sources_receivers.mat
        % info.mat
        % .mat generados en simulación
        
% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @simulador3_OpeningFcn, ...
                   'gui_OutputFcn',  @simulador3_OutputFcn, ...
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


% --- Executes just before simulador3 is made visible.
function simulador3_OpeningFcn(hObject, eventdata, handles, varargin)

% INICIALIZACIÓN VARIABLES:

% Opciones del modelo:
handles.normal=0;
handles.num_planos=0;
handles.num_edges=0;
handles.fuentes=0;
handles.receptores=0;
handles.corners=0;

% Representar fuente/receptor:
handles.plot_source=1;
handles.plot_receive=1;

% Representar rayos: 
handles.fuente_imagen=0;
u=0;
handles.u=u;
v=0;
handles.v=v;
w=0;
handles.w=w;
z=0;
handles.z=z;
todos_rayos=1;
un_rayo=0;
onlydifrac=0;
onlyrefle=0;
rayorder=0;
handles.todos_rayos=todos_rayos;
handles.un_rayo=un_rayo;
handles.onlydifrac=onlydifrac;
handles.onlyrefle=onlyrefle;
handles.rayorder=rayorder;

rayo_num=1;
handles.rayo_num=rayo_num;
order_repres=1;
handles.order_repres=order_repres;
order_onlyspec=1;
handles.order_onlyspec=order_onlyspec;
order_diff=1;
handles.order_diff=order_diff;


% Genera lista de receptores y fuentes según el número introducido en simulador3_env.m:
load('sources_receivers.mat');
for i=1:nsources
    list_sources{i}=i;
end
for i=1:nreceivers
    list_receivers{i}=i;
end
set(handles.plot_src,'String',list_sources);
set(handles.plot_rcv,'String',list_receivers);

num_rcv=1;
num_src=1;
handles.num_rcv=num_rcv;
handles.num_src=num_src;

separacion_point_edges=1;
handles.separacion_point_edges=separacion_point_edges;


% Choose default command line output for simulador3
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);


% --- Outputs from this function are returned to the command line.
function varargout = simulador3_OutputFcn(hObject, eventdata, handles) 

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in normales.
function normales_Callback(hObject, eventdata, handles)
val_normal=get(hObject,'Value');
if val_normal==1
    normal=4;
elseif val_normal==0
    normal=0;
end
handles.normal=normal;
guidata(hObject,handles);

% --- Executes on button press in num_planos.
function num_planos_Callback(hObject, eventdata, handles)
val_num_planos=get(hObject,'Value');
if val_num_planos==1
    num_planos=8;
elseif val_num_planos==0
    num_planos=0;
end
handles.num_planos=num_planos;
guidata(hObject,handles);

% --- Executes on button press in num_edges.
function num_edges_Callback(hObject, eventdata, handles)
val_num_edges=get(hObject,'Value');
if val_num_edges==1
    num_edges=16;
elseif val_num_edges==0
    num_edges=0;
end
handles.num_edges=num_edges;
guidata(hObject,handles);

% --- Executes on button press in fuentes.
function fuentes_Callback(hObject, eventdata, handles)
val_fuentes=get(hObject,'Value');
if val_fuentes==1
    fuentes=1;
elseif val_fuentes==0
    fuentes=0;
end
handles.fuentes=fuentes;
guidata(hObject,handles);

% --- Executes on button press in receptores.
function receptores_Callback(hObject, eventdata, handles)
val_receptores=get(hObject,'Value');
if val_receptores==1
    receptores=2;
elseif val_receptores==0
    receptores=0;
end
handles.receptores=receptores;
guidata(hObject,handles);

% --- Executes on button press in esquinas.
function esquinas_Callback(hObject, eventdata, handles)
val_corners=get(hObject,'Value');
if val_corners==1
    corners=32;
elseif val_corners==0
  corners=0;
end
handles.corners=corners;
guidata(hObject,handles);

% --- Executes on button press in borrar.
function borrar_Callback(hObject, eventdata, handles)
plotoptions=0;
load('info.mat')
eddatafile=[Filepath Filestem '_eddata.mat'];
axis(handles.axes1);
plotoptions3=EDB2extrnums('-37.5 30');
plotoptions2=[];
EDB2plotmodel(eddatafile,plotoptions,plotoptions2,plotoptions3)
axis equal;
% Borra campo que indica que rayos se han representado:
set(handles.rayosdibujados,'String','');
% Los checks no se desactivan...
guidata(hObject,handles);

% --- Executes on button press in resultados_simulacion.
function resultados_simulacion_Callback(hObject, eventdata, handles)
% RESULTADOS_SIMULACION EN .TXT:
 !notepad resultados.txt

% --- Executes on button press in dibujar_model.
function dibujar_model_Callback(hObject, eventdata, handles)
normal=handles.normal;
num_planos=handles.num_planos;
num_edges=handles.num_edges;
fuentes=handles.fuentes;
receptores=handles.receptores;
corners=handles.corners;

plotoptions=normal+num_planos+num_edges+fuentes+receptores+corners;

load('info.mat')
eddatafile=[Filepath Filestem '_eddata.mat'];

axis(handles.axes1);
plotoptions3=EDB2extrnums('-37.5 30');
num_src=handles.num_src;
num_rcv=handles.num_rcv;
size_src=size(num_src);
size_rcv=size(num_rcv);
if size_src(1) < size_rcv(1)
    dif=size_rcv(1)-size_src(1);
    num_src=[num_src zeros]';
elseif size_src(1) > size_rcv(1)
    dif=size_src(1)-size_rcv(1);
    num_rcv=[num_rcv zeros]';
end
plotoptions2=[num_src num_rcv];
EDB2plotmodel(eddatafile,plotoptions,plotoptions2,plotoptions3)
axis equal;
guidata(hObject,handles);

% --- Executes on selection change in num_rayo.
function num_rayo_Callback(hObject, eventdata, handles)
rayo_num=get(hObject,'Value');
handles.rayo_num=rayo_num;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function num_rayo_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in plot_src.
function plot_src_Callback(hObject, eventdata, handles)
plot_source=get(hObject,'Value');
u=0;
handles.u=u;
handles.plot_source=plot_source;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function plot_src_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in plot_rcv.
function plot_rcv_Callback(hObject, eventdata, handles)
plot_receive=get(hObject,'Value');
 u=0;
 handles.u=u;
handles.plot_receive=plot_receive;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function plot_rcv_CreateFcn(hObject, eventdata, handles)
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

function num_rcv_Callback(hObject, eventdata, handles)
num_rcv=str2num(get(hObject,'String'))';
handles.num_rcv=num_rcv;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function num_rcv_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function num_src_Callback(hObject, eventdata, handles)
num_src=str2num(get(hObject,'String'))';
handles.num_src=num_src;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function num_src_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in dibujar_paths.
function dibujar_paths_Callback(hObject, eventdata, handles)
plot_receive=handles.plot_receive;
plot_source=handles.plot_source;
fuente_imagen=handles.fuente_imagen;
todos_rayos=handles.todos_rayos;
un_rayo=handles.un_rayo;
onlydifrac=handles.onlydifrac;
onlyrefle=handles.onlyrefle;
rayorder=handles.rayorder;
rayo_num=handles.rayo_num;
order_repres=handles.order_repres;
order_onlyspec=handles.order_onlyspec;
order_diff=handles.order_diff;
separacion_point_edges=handles.separacion_point_edges;

load('info.mat')
eddatafile=[Filepath Filestem '_eddata.mat'];
reflpathsfile=[Filepath Filestem '_' num2str(plot_source) '_' num2str(plot_receive) '_edpaths.mat'];

if todos_rayos == 1  % Opción Dibujar todos los rayos
    api=0;
    n_rays=representation_paths(eddatafile,reflpathsfile,fuente_imagen,plot_source,plot_receive,separacion_point_edges,api);
    axis equal;
    set(handles.txt_rays,'String','El numero total de rayos representados es:');
    set(handles.rayosdibujados,'String',n_rays);
elseif un_rayo == 1 % Opción Representar rayo número
    representation_onepath(eddatafile,reflpathsfile,fuente_imagen,rayo_num,plot_source,plot_receive,separacion_point_edges);
    axis equal;
    set(handles.txt_rays,'String','Los rayos que ya ha representado son:');
    str=get(handles.rayosdibujados,'String');
    if str ~= 0
        set(handles.rayosdibujados,'String',[str ' ' num2str(rayo_num)]);
    else
     set(handles.rayosdibujados,'String',num2str(rayo_num));
    end
elseif rayorder == 1 % Opción Representar rayos de orden
    rep=representation_specularpaths(eddatafile,reflpathsfile,order_repres,fuente_imagen);
    dibujar=representation_diffpaths(eddatafile,reflpathsfile,order_repres,fuente_imagen,plot_source,plot_receive,separacion_point_edges);
    axis equal;
    dibujar=[rep' dibujar];
    set(handles.txt_rays,'String','Los rayos representados con el orden seleccionado son:');
    for i=1:length(dibujar)
        str=get(handles.rayosdibujados,'String');
        if str == 0
            set(handles.rayosdibujados,'String', num2str(dibujar(i)));
        else
            set(handles.rayosdibujados,'String',[str ' ' num2str(dibujar(i))]);
        end
    end
elseif onlydifrac == 1 % Opción Difracciones de orden
    dibujar=representation_diffpaths(eddatafile,reflpathsfile,order_diff,fuente_imagen,plot_source,plot_receive,separacion_point_edges);
    axis equal;
    set(handles.txt_rays,'String','Los rayos que tienen algún orden de difracción es:');
    for i=1:length(dibujar)
        str=get(handles.rayosdibujados,'String');
        if str == 0
            set(handles.rayosdibujados,'String', num2str(dibujar(i)));
        else
            set(handles.rayosdibujados,'String',[str ' ' num2str(dibujar(i))]);
        end
    end
elseif onlyrefle == 1 % Opción Reflexiones de orden
    rep=representation_specularpaths(eddatafile,reflpathsfile,order_onlyspec,fuente_imagen);
    axis equal;
    set(handles.txt_rays,'String','Los rayos representados son:');
    rep=rep';
    for i=1:length(rep)
        str=get(handles.rayosdibujados,'String');
        if str == 0
            set(handles.rayosdibujados,'String', num2str(rep(i)));
        else
            set(handles.rayosdibujados,'String',[str ' ' num2str(rep(i))]);
        end
    end
end

guidata(hObject,handles);


% --- Executes on button press in fuente_imagen.
function fuente_imagen_Callback(hObject, eventdata, handles)
fuente_imagen=get(hObject,'Value');
handles.fuente_imagen=fuente_imagen;
guidata(hObject,handles);


% --- Executes when selected object is changed in typerayos.
function typerayos_SelectionChangeFcn(hObject, eventdata, handles)
u=handles.u;
v=handles.v;
w=handles.w;
z=handles.z;
if hObject == handles.all_ray  % Opción Dibujar todos los rayos
     todos_rayos=1;
     un_rayo=0;
     onlydifrac=0;
     onlyrefle=0;
     rayorder=0;
elseif hObject == handles.one_ray  % Opción Representar rayo número
    todos_rayos=0;
    un_rayo=1;
    onlydifrac=0;
    onlyrefle=0;
    rayorder=0;
    load('info.mat')
    plot_receive=num2str(handles.plot_receive);
    plot_source=num2str(handles.plot_source);
    reflpathsfile=[Filepath Filestem '_' plot_source '_' plot_receive '_edpaths.mat'];
    load(reflpathsfile);
    num_rayos=size(reflpaths);
    if u==0  % Desplegable imprime los números posibles si no han sido dibujados antes
        for i=1:num_rayos(1)
            list_rayos{i}=i;
        end
        set(handles.num_rayo,'String',list_rayos);
        u=1;
    end
 elseif hObject == handles.ordenrayos % Opción Representar rayos de orden
     rayorder=1;
     onlyrefle=0;
     todos_rayos=0;
     un_rayo=0;
     onlydifrac=0;
     
    load('info.mat')
    plot_receive=num2str(handles.plot_receive);
    plot_source=num2str(handles.plot_source);
    reflpathsfile=[Filepath Filestem '_' plot_source '_' plot_receive '_edpaths.mat'];
    load(reflpathsfile);
    max_orden=length(pathtypevec(1,:));
    if v==0
        for i=1:max_orden
            list_rayos{i}=i;
        end
        set(handles.repres_orden,'String',list_rayos);
        v=1;
    end
    
 elseif hObject == handles.soloreflexiones % Opción Reflexiones de orden
     onlyrefle=1;
     rayorder=0;
     todos_rayos=0;
     un_rayo=0;
     onlydifrac=0;
     
    load('info.mat')
    plot_receive=num2str(handles.plot_receive);
    plot_source=num2str(handles.plot_source);
    reflpathsfile=[Filepath Filestem '_' plot_source '_' plot_receive '_edpaths.mat'];
    load(reflpathsfile);
    max_orden=length(pathtypevec(1,:));
    if w==0
        for i=1:max_orden
            list_rayos{i}=i;
        end
        set(handles.ordensoloreflex,'String',list_rayos);
        w=1;
    end
 elseif hObject == handles.solodifracc  % Opción Difracciones de orden
     onlydifrac=1;
     onlyrefle=0;
     rayorder=0;
     todos_rayos=0;
     un_rayo=0;
     load('info.mat')
    plot_receive=num2str(handles.plot_receive);
    plot_source=num2str(handles.plot_source);
    reflpathsfile=[Filepath Filestem '_' plot_source '_' plot_receive '_edpaths.mat'];
    load(reflpathsfile);
    max_orden=length(pathtypevec(1,:));
    if z==0
        for i=1:max_orden
            list_rayos{i}=i;
        end
        set(handles.ordendiffr,'String',list_rayos);
        w=1;
    end
end
 
handles.u=u;
handles.v=v;
handles.w=w;
handles.z=z;
handles.todos_rayos=todos_rayos;
handles.un_rayo=un_rayo;
handles.onlydifrac=onlydifrac;
handles.onlyrefle=onlyrefle;
handles.rayorder=rayorder;
guidata(hObject,handles);


% --- Executes on button press in autoescalado.
function autoescalado_Callback(hObject, eventdata, handles)
autoescalado=get(hObject,'Value');
if autoescalado == 1
    axis normal;
elseif autoescalado == 0
    axis equal;
end

% --- Executes on selection change in ordensoloreflex.
function ordensoloreflex_Callback(hObject, eventdata, handles)
order_onlyspec=get(hObject,'Value');
handles.order_onlyspec=order_onlyspec;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function ordensoloreflex_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ordendiffr.
function ordendiffr_Callback(hObject, eventdata, handles)
order_diff=get(hObject,'Value');
handles.order_diff=order_diff;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function ordendiffr_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes during object creation, after setting all properties.
function repres_orden_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

% --- Executes on selection change in repres_orden.
function repres_orden_Callback(hObject, eventdata, handles)
order_repres=get(hObject,'Value');
handles.order_repres=order_repres;
guidata(hObject,handles);

% --- Executes during object creation, after setting all properties.
function rayosdibujados_CreateFcn(hObject, eventdata, handles)


% --- Executes during object creation, after setting all properties.
function txt_rays_CreateFcn(hObject, eventdata, handles)



function sep_dist_edges_Callback(hObject, eventdata, handles)
separacion_point_edges=str2num(get(hObject,'String'))';
handles.separacion_point_edges=separacion_point_edges;
guidata(hObject,handles);


% --- Executes during object creation, after setting all properties.
function sep_dist_edges_CreateFcn(hObject, eventdata, handles)

if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
