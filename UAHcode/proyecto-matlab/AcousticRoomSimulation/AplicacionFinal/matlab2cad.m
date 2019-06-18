function matlab2cad()
% Convert struct MATLAB (data_struct.mat) TO .CAD

load('data_struct.mat')

corn=1;
corners1(1,:)=[0,0,0];
aux = 0;

%% PARTE 1: lee estructura
% planecorners: contiene matriz donde filas definen cada uno de los planos, y los números definen el número de corner
% corners1: coordenadas de cada uno de los corners que han aparecido,aunque estén repetidos
% plane_true [número de corners (repetidos tambien),numero de planos]:coloca un 1 en el lugar

for i=0:numsurfaces-1
    eval(['numvert=surfaces.surface' num2str(i) '.numvertexes;']);
    for j=0:numvert-1
        eval(['aux_c = surfaces.surface' num2str(i) '.vertex' num2str(j) ';']);
        aux_c=aux_c';
        corners1(corn,:)=aux_c;
        planecorners(i+1,j+1)=corn;
        plane_true(j+aux+1,i+1) = 1;
        corn=corn+1;
    end
    aux = aux + numvert;
end


ncorners1=size(corners1);
ncorners1=ncorners1(1);
corners_clean=corners1;

u=1;
% elimino corners repetidos:
for p=1:ncorners1
    for r=1:ncorners1
        repeat=isequal(corners_clean(p,:),corners_clean(r,:));
        if p==r
            repeat=0;
        end
        if repeat==1
            ncorners_clean(u)=r;
            u=u+1;
        end
    end
    
    if ncorners_clean~=0
        %sutituye numero de corner en definicion de plano:
        for t=1:length(ncorners_clean)
            if ncorners_clean(t) < p
                continue
            end
            plane_true(p,:)=plane_true(p,:)+plane_true(ncorners_clean(t),:);
            plane_true(ncorners_clean(t),:)=plane_true(ncorners_clean(t),:)-plane_true(ncorners_clean(t),:);
        end
        ncorners_clean=0;
        u=1;
    end
end

% creo matriz corners:
j=1;
true = plane_true;
h=1;
for y=1:ncorners1
    a=sum(plane_true(y,:));
    if a ~= 0
        corners(j,:)=corners1(y,:);
        j=j+1;
    end
    if a == 0
        f0(h)=y;
        h=h+1;
    end
end
%elimino filas todo cero de plane_true y guardo en true
true(f0,:)=[];

% creo vector cornernumbers:
for i=1:length(corners)
    cornernumbers(i)=i;
end
cornernumbers=cornernumbers';

[tamp1,tamp2]=size(planecorners);
% creo matriz planecorners:
for t=1:numsurfaces
    e=find(true(:,t)==1);
    if length(e) < tamp2
        e = [e' zeros(1,tamp2-length(e))];
    end
    planecorners(t,:)=e;
    clear e;
end


clear a aux aux_c corn corners1 corners_clean e f0 h i j ncorners1 ncorners_clean plane_true r repeat t true u y
% END PARTE 1:
% cornernumbers: vector [número de corners (sin repetir),1]
% corners: coordenadas de cada uno de los cornernumbers
% planecorners: matriz [número de planos, número de vertices] con los corners que forman cada plano

%% PARTE 2: coloca corners en la definición de planecorners para que sean consecutivos

for ii=1:numsurfaces
    plano=planecorners(ii,:);
    pos_nozero=find(plano ~= 0);
    plano=[plano(pos_nozero) plano(pos_nozero(1:(length(pos_nozero)/2)))];
    for jj = 1:length(pos_nozero)
        co1numb = plano(jj);
        co2numb = plano(jj+1);
        co3numb = plano(jj+2);
        vec1 = corners(co1numb,:) - corners(co2numb,:);
        vec2 = corners(co3numb,:) - corners(co2numb,:);
        nvec = EDB2cross(vec1.',vec2.').';
        nveclen(jj) = norm(nvec);
        if nveclen(jj) > 0
            nvectorlist(jj,:) = nvec./nveclen(jj);
        end
    end
    clear plano co1numb co2numb co3numb vec1 vec2 nvec nveclen
    
    % si nvectorlist tiene toda la fila con el mismo signo es valido si no, hay
    % que cambiar el orden de planecorners!!!
    pos = find(nvectorlist ~= 0);
    normpos=nvectorlist(pos)>0;
    normneg=nvectorlist(pos)<0;
    if(sum(normpos)==0 | sum(normneg)==0)   % esta bien el orden si uno de ellos es 0
        
    else  % cambiar orden donde cambie el valor en solo un lugar
        if (nvectorlist(pos(1))~=nvectorlist(pos(2)))
            aux=planecorners(ii,2);
            planecorners(ii,2)=planecorners(ii,1);
            planecorners(ii,1)=aux;
            continue;   % para que no haga dos veces si ocurre dos veces y no evalue los elseif
        elseif (nvectorlist(pos(2))~=nvectorlist(pos(3)))
            aux=planecorners(ii,3);
            planecorners(ii,3)=planecorners(ii,2);
            planecorners(ii,2)=aux;
            continue;
        elseif (nvectorlist(pos(3))~=nvectorlist(pos(4)))
            aux=planecorners(ii,4);
            planecorners(ii,4)=planecorners(ii,3);
            planecorners(ii,3)=aux;
        end
    end
end

clear aux ii jj pos nvectorlist normpos normneg

%% PARTE 3: calcula el número de obstáculos, incluyendo al recinto
for p=1: length(corners)
    [posx,posy]=find(planecorners== p);
    planes2vert(p,:)=posx';
end

o=1;
obst_aux=planes2vert(1,:);
[valor,ind_unico]=unique(obst_aux);
aux_vertex(o) = 0;
for r=2:length(corners)
    % añadir planos
    valor_repetidos = setdiff( planes2vert(r,:),valor);
    if(length(  valor_repetidos ) == 3)
        aux_vertex(o+1) = length(valor);
        obst(o,:) = valor(aux_vertex(o)+1:aux_vertex(o+1));
        o = o + 1;
    end
    obst_aux = [obst_aux planes2vert(r,:)];
    %Comparar que no se repitan
    [valor,ind_unico]=unique(obst_aux);
end
obst(o,:) = valor(aux_vertex(o)+1:length(valor));
% fallo sin solucionar: si un obstaculo tiene mas o menos de 6 planos

%% PARTE 4: producto cruzado para ver si normal apunta hacia dentro del recinto
%definir limites recinto
numsurfaces_ini=1;
numsurfaces_fin=0;
for g=1:o
    
    plane=planecorners(obst(g,1),:);
    order_corners=corners(plane,:);
    for p=2:length(obst(g,:))
        plane=planecorners(obst(g,p),:);
        order_corners=[order_corners; corners(plane,:)];
    end
    
    minx=min(order_corners(:,1));
    maxx=max(order_corners(:,1));
    miny=min(order_corners(:,2));
    maxy=max(order_corners(:,2));
    minz=min(order_corners(:,3));
    maxz=max(order_corners(:,3));
    
    numsurfaces_fin=numsurfaces_fin+length(obst(g,:));
    if g <= 1
        for j=numsurfaces_ini:numsurfaces_fin
            vec1 = corners(planecorners(j,3),:) - corners(planecorners(j,2),:);
            vec2 = corners(planecorners(j,1),:) - corners(planecorners(j,2),:);
            nvec = EDB2cross(vec1.',vec2.').';
            %compruebo coordenada y:
            if nvec(1,2) ~= 0
                coor=2;
                if corners(planecorners(j,2),coor) == miny
                    if nvec(1,2) > 0
                        continue;
                    else
                        turn=1;
                    end
                elseif  corners(planecorners(j,2),coor) == maxy
                    if nvec(1,2) < 0
                        continue;
                    else
                        turn=1;
                    end
                else   % si no coincide con el max ni el min es que está dentro y da igual como se defina
                    continue;
                    warning('ha entrado al else');
                end
            end
            
            %compruebo coordenada x:
            if nvec(1,1) ~= 0
                coor=1;
                if corners(planecorners(j,2),coor) == minx
                    if nvec(1,1) > 0
                        continue;
                    else
                        turn=1;
                    end
                elseif  corners(planecorners(j,2),coor) == maxx
                    if nvec(1,1) < 0
                        continue;
                    else
                        turn=1;
                    end
                else   % si no coincide con el max ni el min es que está dentro y da igual como se defina
                    continue;
                end
            end
            
            %compruebo coordenada z:
            if nvec(1,3) ~= 0
                coor=3;
                if corners(planecorners(j,2),coor) == minz
                    if nvec(1,3) > 0
                        continue;
                    else
                        turn=1;
                    end
                elseif  corners(planecorners(j,2),coor) == maxz
                    if nvec(1,3) < 0
                        continue;
                    else
                        turn=1;
                    end
                else   % si no coincide con el max ni el min es que está dentro y da igual como se defina
                    continue;
                end
            end
            %se comprueba si está activado turn:
            if turn==1  % se cambia sentido de definición de corners en cada plano
                planecorners(j,:) = fliplr(planecorners(j,:));
                turn=0;
            end
        end
        numsurfaces_ini=numsurfaces_fin+1;
        clear order_corners
        
    elseif  g > 1
        for j=numsurfaces_ini:numsurfaces_fin
            vec1 = corners(planecorners(j,3),:) - corners(planecorners(j,2),:);
            vec2 = corners(planecorners(j,1),:) - corners(planecorners(j,2),:);
            nvec = EDB2cross(vec1.',vec2.').';
            
            %compruebo coordenada y:
            if nvec(1,2) ~= 0
                coor=2;
                if corners(planecorners(j,2),coor) == miny
                    if nvec(1,2) < 0
                        continue;
                    else
                        turn=1;
                    end
                elseif  corners(planecorners(j,2),coor) == maxy
                    if nvec(1,2) > 0
                        continue;
                    else
                        turn=1;
                    end
                else   % si no coincide con el max ni el min es que está dentro y da igual como se defina
                    continue;
                    warning('ha entrado al else');
                end
            end
            
            %compruebo coordenada x:
            if nvec(1,1) ~= 0
                coor=1;
                if corners(planecorners(j,2),coor) == minx
                    if nvec(1,1) < 0
                        continue;
                    else
                        turn=1;
                    end
                elseif  corners(planecorners(j,2),coor) == maxx
                    if nvec(1,1) > 0
                        continue;
                    else
                        turn=1;
                    end
                else   % si no coincide con el max ni el min es que está dentro y da igual como se defina
                    continue;
                end
            end
            
            %compruebo coordenada z:
            if nvec(1,3) ~= 0
                coor=3;
                if corners(planecorners(j,2),coor) == minz
                    if nvec(1,3) < 0
                        continue;
                    else
                        turn=1;
                    end
                elseif  corners(planecorners(j,2),coor) == maxz
                    if nvec(1,3) > 0
                        continue;
                    else
                        turn=1;
                    end
                else   % si no coincide con el max ni el min es que está dentro y da igual como se defina
                    continue;
                end
            end
            %se comprueba si está activado turn:
            if turn==1  % se cambia sentido de definición de corners en cada plano
                planecorners(j,:) = fliplr(planecorners(j,:));
                turn=0;
            end
        end
        numsurfaces_ini=numsurfaces_fin+1;
        clear order_corners
    end
end
%% PARTE 5: escribe archivo archivo_CAD.cad
escale = 1e-3;  % convierte mm de .env y .srf a metros!!
corners = corners*escale;
fid=fopen('environment.cad','w');
fprintf(fid,'CAD-file (version CATT-Acoustic v6)\n\n');
fprintf(fid,'%s\n',comment);
fprintf(fid,'%s\n\n',content);
fprintf(fid,'%%CORNERS\n\n');
for k = 1:length(cornernumbers)
    fprintf(fid,'\t %d \t %f \t %f \t %f \n', k, corners(k,1), corners(k,2), corners(k,3) );
end
fprintf(fid,'\n%%PLANES\n\n');
for k = 1:length(planecorners)
    fprintf(fid,'\t %d / /RIGID\n', k );
    fprintf(fid,'\t %d %d %d %d\n\n', planecorners(k,1), planecorners(k,2), planecorners(k,3), planecorners(k,4) );
end
%fprintf(fid,'%%SOURCES\n\n0 OMNI:SD0\n   0.0000000   0.0000000   0.0001000\n   0.0000000   0.0000000   0.0000000\n  85.0  88.0  91.0  94.0  97.0  100.0\n\n%%RECEIVERS\n\n   1   0.0000000   0.0000000   10.0000000\n');
fprintf(fid,'\n%%EOF');
fclose(fid);


