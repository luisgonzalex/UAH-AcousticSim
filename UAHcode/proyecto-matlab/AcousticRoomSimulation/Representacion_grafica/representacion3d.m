close all
clear all
load('data.mat')

for i=0:numsurfaces-1
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
