function representation_room(extension)
% Representa el recinto y obstáculos del modelo a simular

if extension == 1 % .ENV
    
    load('data_struct.mat');

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

        index=zeros(1,aux_nvertexes+1);
        index(1,:)=[1:aux_nvertexes 1];%   1 2 3 4 1];
        plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),'LineWidth',2)
        hold on
    end

elseif extension == 2 % .CAD
    
    load('info.mat');
    desiredname = [Filepath,Filestem,'_cadgeo'];
    cadgeofile = EDB2readcad(CADfile,desiredname);
    eval(['load ',cadgeofile]);
    [nplanes,slask] = size(planecorners);

   for ii = 1:nplanes   
        v=size(planecorners(1,:));
        for j=0:v(2)-1
          x(j+1)=corners(planecorners(ii,j+1),1);
          y(j+1)=corners(planecorners(ii,j+1),2);
          z(j+1)=corners(planecorners(ii,j+1),3);
        end

        index=zeros(1,v(2)+1);
        index(1,:)=[1:v(2) 1];
        plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),'LineWidth',2)
        hold on
   end
end


      