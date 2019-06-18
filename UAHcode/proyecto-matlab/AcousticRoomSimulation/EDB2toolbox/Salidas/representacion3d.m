close all
clear all


load('data.mat');
load('mioorigin_1_1_edpaths.mat');
specextradata=full(specextradata);
%% DIBUJA RECINTO: (sin mesa)


for i=0:numsurfaces-2  % -2 para que no dibuje la mesa!
    eval(['aux_surface=surfaces.surface' num2str(i) ';']);
    eval(['aux_nvertexes=aux_surface.numvertexes;']);
    for j=0:aux_nvertexes-1
      eval(['vertex = aux_surface.vertex' num2str(j) ';']);  
      x(j+1)=vertex(1);
      y(j+1)=vertex(2);
      z(j+1)=vertex(3);
    end
    x=x*1e-3;
    y=y*1e-3;
    z=z*1e-3;
    
    index=zeros(1,5);
    index(1,:)=[1 2 3 4 1];
    plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)))
    hold all
end

%% dibuja fuente y receptor
plot3(S(1,1),S(1,2),S(1,3),'ro')  %dibuja la fuente
plot3(R(1,1),R(1,2),R(1,3),'go')  %dibuja receptor

%% dibuja rayos:
for p=1:length(pathtypevec) %para saber cuando hay m�s de un punto de reflejo (orden mayor que 1)
    if pathtypevec(p,2)~=0
        d=p;
        break;
    end
end

for j=1:length(specextradata)
    if pathtypevec(j) == 102 %rayo directo
       x=[specextradata(j,1) specextradata(j,4)];
       y=[specextradata(j,2) specextradata(j,5)];
       z=[specextradata(j,3) specextradata(j,6)];
       index=zeros(1,2);
       index(1,:)=[1 2];
       plot3(x(index(1,:)),y(index(1,:)),z(index(1,:))); 
    elseif pathtypevec(j) == 115 %reflexion especular
        if j < d
         x=[R(1,1) specextradata(j,4)];
       y=[R(1,2) specextradata(j,5)];
       z=[R(1,3) specextradata(j,6)];
       index=zeros(1,2);
       index(1,:)=[1 2];
       plot3(x(index(1,:)),y(index(1,:)),z(index(1,:))); 
       
       x=[S(1,1) specextradata(j,4)];
       y=[S(1,2) specextradata(j,5)];
       z=[S(1,3) specextradata(j,6)];
       index=zeros(1,2);
       index(1,:)=[1 2];
       plot3(x(index(1,:)),y(index(1,:)),z(index(1,:))); 
        else  %mayor orden de difracci�n
          x=[R(1,1) specextradata(j,7) specextradata(j,4) S(1,1)];
       y=[R(1,2) specextradata(j,8) specextradata(j,5) S(1,2)];
       z=[R(1,3) specextradata(j,9) specextradata(j,6) S(1,3)];
       index=zeros(1,4);
       index(1,:)=[1 2 3 4];
       plot3(x(index(1,:)),y(index(1,:)),z(index(1,:))); 
       
        end
    end  
end

      