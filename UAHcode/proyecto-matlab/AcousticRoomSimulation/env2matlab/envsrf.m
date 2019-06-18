% Convert .env and .srf TO "struct datos" MATLAB (data.mat) 
% version 7 (use fget)
clear all

%% open file .env
[archivo,ruta] = uigetfile('*.env','Please select the .env file');
if ~isstr(archivo)
    return
end
[temp1,nombre,extension] = fileparts(archivo);
envfile = [ruta,archivo];

%% read .env
fid = fopen(envfile,'r');
if fid == -1
    error(['ERROR: The file could not be opened.'])
end

tline = fgets(fid);

while ischar(tline)
   tline = fgets(fid);
   iden=strtok(tline);
   disp(iden)
   %numsurfaces
   if strcmpi(iden,'numsurfaces')
        pos = findstr(tline,'=');
        if ~isempty(pos)
          aux_surf=tline(pos+1:end); 
            ind = find( ~isspace(aux_surf) ); % posicion indice sin espacio   
            if isempty(ind)
             numsurfaces = [];        
            else
            numsurfaces = str2double(aux_surf( ind(1):ind(end) ));
            end
          clear aux_surf
        else
            error(['ERROR: Not found field: numsurfaces.'])
        end
        %content
   elseif strcmpi(iden,'content')
       pos = findstr(tline,'=');
       aux_content=tline(pos+1:end);
       ind = find( ~isspace(aux_content) ); % posicion indice sin espacio   
            if isempty(ind)
             content = [];        
            else
           content = aux_content( ind(1):ind(end) );
            end
            clear aux_content
       %version
   elseif strcmpi(iden,'version')
       pos = findstr(tline,'=');
       aux_version=tline(pos+1:end); 
       ind = find( ~isspace(aux_version) ); % posicion indice sin espacio   
            if isempty(ind)
             version = [];        
            else
            version = str2double(aux_version( ind(1):ind(end) ));
            end
        clear aux_version
     %comment
   elseif strcmpi(iden,'comment')
       pos = findstr(tline,'=');
       aux_comment=tline(pos+1:end);
       ind = find( ~isspace(aux_comment) ); % posicion indice sin espacio   
            if isempty(ind)
             comment = [];        
            else
           comment = aux_comment( ind(1):ind(end) );
            end
            clear aux_comment
      %dirEnvironments
   elseif strcmpi(iden,'direnvironments')
        pos = findstr(tline,'=');
       dirEnvironments=tline(pos+1:end);
       ind = find( ~isspace(dirEnvironments) ); % posicion indice sin espacio   
            if isempty(ind)
             dirEnvironments = [];        
            else
            dirEnvironments = dirEnvironments( ind(1):ind(end) );
            end
       %% .srf:   �sacarlo fuera del while cambiando la condicion del while?
   elseif strcmpi(iden,'surface0')
       pos = findstr(tline,'=');
       aux_s=tline(pos+1:end);
       ind = find( ~isspace(aux_s) ); % posicion indice sin espacio 
       if isempty(ind)
             aux_s0 = [];        
       else
            aux_s0 = aux_s( ind(1):ind(end) );
       end
            clear aux_s
        break
   end
            
end


surfaces.surface0.name=aux_s0;
ruta_srf0=[ruta,surfaces.surface0.name];
eval([' fid0 = fopen( ruta_srf0 ,''r'');']);
    if fid0 == -1
        error(['ERROR: The file ',aux_s0,' could not be opened.'])
    end
tline0 = fgets(fid0);

while ischar(tline0)
   tline0 = fgets(fid0);
   iden0=strtok(tline0);
   % numvertexes
   if strcmpi(iden0,'numvertexes')
        pos = findstr(tline0,'=');
        if ~isempty(pos)
          aux_vert=tline0(pos+1:end); 
            ind = find( ~isspace(aux_vert) ); % posicion indice sin espacio   
            if isempty(ind)
             numvertexes = [];        
            else
            numvertexes = str2double(aux_vert( ind(1):ind(end) ));
            end
          surfaces.surface0.numvertexes=numvertexes;
          clear aux_vert
        else
            error('ERROR: Not found field: numvertexes.')
        end
   elseif strcmpi(iden0,'vertex0')
       pos = findstr(tline0,'=');
       aux_v=tline0(pos+1:end);
       ind = find( ~isspace(aux_v) ); % posicion indice sin espacio 
       if isempty(ind)
             aux_v0 = [];        
       else
            aux_v0 = aux_v( ind(1):ind(end) );
       end
            clear aux_v
        break
   end
end

surfaces.surface0.vertex0=sscanf(aux_v0,'%d');

   for j=1:numvertexes-1
        tline0 = fgets(fid0);
        iden0=strtok(tline0);
        aux0=['vertex' num2str(j)];
     if strcmpi(iden0,aux0)
       pos = findstr(tline0,'=');
       aux_v=tline0(pos+1:end);
       ind = find( ~isspace(aux_v) ); % posicion indice sin espacio 
       if isempty(ind)
             aux_v = [];        
       else
            aux_v = aux_v( ind(1):ind(end) );
       end
    aux_v = sscanf(aux_v,'%d');
    eval(['surfaces.surface0.' iden0 '= aux_v;']);
     end
   end

   fclose(fid0);
%% resto de .srf
 for k = 1:numsurfaces-1
     tline = fgets(fid);
     iden=strtok(tline);
     aux=['surface' num2str(k)];
     if strcmpi(iden,aux)
       pos = findstr(tline,'=');
       aux_s=tline(pos+1:end);
       ind = find( ~isspace(aux_s) ); % posicion indice sin espacio 
       if isempty(ind)
             aux_s = [];        
       else
            aux_s = aux_s( ind(1):ind(end) );
       end
       
    eval(['surfaces.' iden '.name= aux_s;']);
    eval(['ruta_srfx=[ruta,surfaces.' iden '.name];'])
    eval([' fidx = fopen( ruta_srfx ,''r'');']);
    if fidx == -1
        error(['ERROR: The file ',aux_s,' could not be opened.'])
    end
tlinex = fgets(fidx);
while ischar(tlinex)
   tlinex = fgets(fidx);
   idenx=strtok(tlinex);
   % numvertexes
   if strcmpi(idenx,'numvertexes')
        pos = findstr(tlinex,'=');
        if ~isempty(pos)
          aux_vert=tlinex(pos+1:end); 
            ind = find( ~isspace(aux_vert) ); % posicion indice sin espacio   
            if isempty(ind)
             numvertexes = [];        
            else
            numvertexes = str2double(aux_vert( ind(1):ind(end) ));
            end
          clear aux_vert
         eval(['surfaces.' aux '.numvertexes=numvertexes;']);
        else
            error('ERROR: Not found field: numvertexes.')
        end
        
   elseif strcmpi(idenx,'vertex0')
       pos = findstr(tlinex,'=');
       aux_v=tlinex(pos+1:end);
       ind = find( ~isspace(aux_v) ); % posicion indice sin espacio 
       if isempty(ind)
             aux_vx = [];        
       else
            aux_vx = aux_v( ind(1):ind(end) );
       end
            clear aux_v
        break
   end
end

eval(['surfaces.' iden '.vertex0 = sscanf(aux_vx,''%d'');']);


eval(['aux_vertexes=surfaces.' aux '.numvertexes;']);
   for p=1:aux_vertexes-1
        tlinexx = fgets(fidx);
        idenxx=strtok(tlinexx);
        auxxx=['vertex' num2str(p)];
     if strcmpi(idenxx,auxxx)
       pos = findstr(tlinexx,'=');
       aux_v=tlinexx(pos+1:end);
       ind = find( ~isspace(aux_v) ); % posicion indice sin espacio 
       if isempty(ind)
             aux_v = [];        
       else
            aux_v = aux_v( ind(1):ind(end) );
       end
    aux_v = sscanf(aux_v,'%d');
   eval(['surfaces.' iden '.' idenxx '= aux_v;']);
     end
   end

   fclose(fidx);
     end
 end

fclose(fid);

save data numsurfaces comment content version dirEnvironments surfaces
