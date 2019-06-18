function representation_onepath(eddatafile,reflpathsfile,fuente_imagen,num_rayo,plot_source,plot_receive,separacion_point_edges)

[eddatafilepath,eddatafile,fileext] = fileparts(eddatafile);
eddatafile = [eddatafile,fileext];
Filestem = EDB2strpend(eddatafile,'_eddata');

[reflpathsfilepath,reflpathsfile,fileext] = fileparts(reflpathsfile);
reflpathsfile = [reflpathsfile,fileext];

eval(['load ',eddatafilepath,filesep,eddatafile]);
eval(['load ',reflpathsfilepath,filesep,reflpathsfile]);

colors=['y' 'm' 'c' 'r' 'g' 'b' 'k'];
c=randi(7,1,1);


pathtype = pathtypevec(num_rayo,:);
reflpath = reflpaths(num_rayo,:);
order=length(find(pathtype)~=0);
path_data=full(specextradata(num_rayo,:));
exist100=0;

if firstdiffrow ~=0
difraccion = mainlistguide(firstdiffrow,2);
else
    difraccion = length(pathtypevec) + 1 ;
end

if num_rayo < difraccion % aqui se evaluan los rayos directos y especulares puros
    if fuente_imagen == 1
            if pathtype(1) == 115
                x = path_data((1*3)-2);
                y = path_data((1*3)-1);
                z = path_data((1*3));    
                plot3(x,y,z,'rp');
            end
    end
    for h=1:order
        if pathtype(h) == 102 %rayo directo
            x=[path_data(1) path_data(4)];
            y=[path_data(2) path_data(5)];
            z=[path_data(3) path_data(6)];
            index=zeros(1,2);
            index(1,:)=[1 2];
            plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
        elseif pathtype(h)== 115  %reflexión especular
            if h+1 > order  

                x=[R(1,1) path_data(4+(h-1)*3)];
                y=[R(1,2) path_data(5+(h-1)*3)];
                z=[R(1,3) path_data(6+(h-1)*3)];
                index=zeros(1,2);
                index(1,:)=[1 2];
                plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                
                if h == 1
                    x=[S(1,1) path_data(4+(h-1)*3)];
                    y=[S(1,2) path_data(5+(h-1)*3)];
                    z=[S(1,3) path_data(6+(h-1)*3)];
                    index=zeros(1,2);
                    index(1,:)=[1 2];
                    plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                    
                else
                    x=[path_data(4+(h-2)*3) path_data(4+(h-1)*3)];
                    y=[path_data(5+(h-2)*3) path_data(5+(h-1)*3)];
                    z=[path_data(6+(h-2)*3) path_data(6+(h-1)*3)];
                    index=zeros(1,2);
                    index(1,:)=[1 2];
                    plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                   
                end
            else
                if h==1
                    x=[S(1,1) path_data(4+(h-1)*3)];
                    y=[S(1,2) path_data(5+(h-1)*3)];
                    z=[S(1,3) path_data(6+(h-1)*3)];
                    index=zeros(1,2);
                    index(1,:)=[1 2];
                    plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                    if order >= 3   %dibuje rayos intermedios
                        x=[path_data(4+(h-1)*3) path_data(4+(h)*3)];
                        y=[path_data(5+(h-1)*3) path_data(5+(h)*3)];
                        z=[path_data(6+(h-1)*3) path_data(6+(h)*3)];
                        index=zeros(1,2);
                        index(1,:)=[1 2];
                        plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                    end
                else
                    x=[path_data(4+(h-1)*3) path_data(4+(h)*3)];
                    y=[path_data(5+(h-1)*3) path_data(5+(h)*3)];
                    z=[path_data(6+(h-1)*3) path_data(6+(h)*3)];
                    index=zeros(1,2);
                    index(1,:)=[1 2];
                    plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));     
                end
            end
        end
end

elseif num_rayo >= difraccion  % en cuanto hay algún orden de difraccion

    if pathtypevec(num_rayo,1) == 115 %primero va a haber una reflexion  y luego la difraccion
    diffedge=reflpaths(num_rayo,2);
    %linea con mayor grosor:
        corners_edge1=edgecorners(diffedge,1);
        corners_edge2=edgecorners(diffedge,2);
        iv = [corners_edge1;corners_edge2];
        h = plot3(corners(iv,1),corners(iv,2),corners(iv,3));
        set(h,'LineWidth',2)
        
    diffplane=reflpaths(num_rayo,1);
    midpoint = mean(corners(edgecorners(diffedge,1:2),:)); %punto medio del edge
long_edge = corners(edgecorners(diffedge,1),:) - corners(edgecorners(diffedge,2),:);
      %  pos_long_edge=find(long_edge ~= 0);
      
      vector1=corners(edgecorners(diffedge,1),:)-midpoint;
        vector2=corners(edgecorners(diffedge,2),:)-midpoint;
        vector1=vector1/norm(vector1);
        vector2=vector2/norm(vector2);  
      
      
%         midpoint1_1coor=midpoint(1,pos_long_edge)+separacion_point_edges; % se añade un metro al punto medio
%         midpoint1_2coor=midpoint(1,pos_long_edge)-separacion_point_edges;
%         midpoint1_1=midpoint;
%         midpoint1_2=midpoint;
%         midpoint1_1(pos_long_edge)=midpoint1_1coor;
%         midpoint1_2(pos_long_edge)=midpoint1_2coor;
%         
midpoint1_1coor=midpoint+vector1.*separacion_point_edges;    
midpoint1_1=midpoint;
            midpoint1_2=midpoint;


        r=2;
         while sum(abs(long_edge/2) > abs(midpoint-midpoint1_1coor ))
        
             midpoint1_1coor=midpoint1_1(r-1,:)+vector1.*separacion_point_edges; % se añade un metro al punto medio
               midpoint1_2coor=midpoint1_2(r-1,:)+vector2.*separacion_point_edges;
               
               midpoint1_1(r,:)=midpoint1_1coor;
               midpoint1_2(r,:)=midpoint1_2coor;
%             midpoint1_1coor=midpoint1_1(r-1,pos_long_edge)+separacion_point_edges; % se añade un metro al punto medio
%             midpoint1_2coor=midpoint1_2(r-1,pos_long_edge)-separacion_point_edges;
%             midpoint1_1(r,:)=midpoint1_1(r-1,:);
%             midpoint1_2(r,:)=midpoint1_2(r-1,:);
%             midpoint1_1(r,pos_long_edge)=midpoint1_1coor;
%             midpoint1_2(r,pos_long_edge)=midpoint1_2coor;
            r=r+1;
        end
         midpoint1_1(r-1,:)=[];
        midpoint1_2(r-1,:)=[];
        %midpoint=[midpoint;midpoint1_1;midpoint1_2];
        midpoint=[midpoint;midpoint1_1(2:end,:);midpoint1_2(2:end,:)];
        
     corner2=corners(planecorners(diffplane,2),:);
     corner1=corners(planecorners(diffplane,1),:);
     corner3=corners(planecorners(diffplane,3),:);
     plane = [corner2   corner1-corner2   corner3-corner2];
     
     IS=specextradata(num_rayo,1:3);
     
     if fuente_imagen == 1
         plot3(IS(1),IS(2),IS(3),'rp')
     end
     
    for j=1:length(midpoint(:,1))     
     line = [midpoint(j,:)  IS-midpoint(j,:)];
    point = intersectLinePlane(line, plane, 1e-14);
        x=[S(1,1) point(1) midpoint(j,1) R(1,1)];
         y=[S(1,2) point(2) midpoint(j,2) R(1,2)];
         z=[S(1,3) point(3) midpoint(j,3) R(1,3)];
         index=zeros(1,4);
         index(1,:)=[1 2 3 4];
         plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
          hold on;
    end


elseif pathtypevec(num_rayo,1) == 100 %primero va a haber una difraccion  

   if pathtypevec(num_rayo,2) == 115
    diffedge=reflpaths(num_rayo,1);
    %linea con mayor grosor:
            corners_edge1=edgecorners(diffedge,1);
            corners_edge2=edgecorners(diffedge,2);
            iv = [corners_edge1;corners_edge2];
            h = plot3(corners(iv,1),corners(iv,2),corners(iv,3));
            set(h,'LineWidth',2)
            
    diffplane=reflpaths(num_rayo,2);
    midpoint = mean(corners(edgecorners(diffedge,1:2),:)); 
long_edge = corners(edgecorners(diffedge,1),:) - corners(edgecorners(diffedge,2),:);
           % pos_long_edge=find(long_edge ~= 0);
         
           vector1=corners(edgecorners(diffedge,1),:)-midpoint;
        vector2=corners(edgecorners(diffedge,2),:)-midpoint;
        vector1=vector1/norm(vector1);
        vector2=vector2/norm(vector2);
        
        
        midpoint1_1coor=midpoint+vector1.*separacion_point_edges;
        
        midpoint1_1=midpoint;
            midpoint1_2=midpoint;
        
%             midpoint1_1coor=midpoint(1,pos_long_edge)+separacion_point_edges; % se añade un metro al punto medio
%             midpoint1_2coor=midpoint(1,pos_long_edge)-separacion_point_edges;
%             midpoint1_1=midpoint;
%             midpoint1_2=midpoint;
%             midpoint1_1(pos_long_edge)=midpoint1_1coor;
%             midpoint1_2(pos_long_edge)=midpoint1_2coor;
            
            r=2;
            while sum(abs(long_edge/2) > abs(midpoint-midpoint1_1coor ))        
                 midpoint1_1coor=midpoint1_1(r-1,:)+vector1.*separacion_point_edges; % se añade un metro al punto medio
               midpoint1_2coor=midpoint1_2(r-1,:)+vector2.*separacion_point_edges;
               
               midpoint1_1(r,:)=midpoint1_1coor;
               midpoint1_2(r,:)=midpoint1_2coor;
%                 midpoint1_1coor=midpoint1_1(r-1,pos_long_edge)+separacion_point_edges; % se añade un metro al punto medio
%                 midpoint1_2coor=midpoint1_2(r-1,pos_long_edge)-separacion_point_edges;
%                 midpoint1_1(r,:)=midpoint1_1(r-1,:);
%                 midpoint1_2(r,:)=midpoint1_2(r-1,:);
%                 midpoint1_1(r,pos_long_edge)=midpoint1_1coor;
%                 midpoint1_2(r,pos_long_edge)=midpoint1_2coor;
                r=r+1;
            end
           midpoint1_1(r-1,:)=[];
            midpoint1_2(r-1,:)=[];

            %midpoint=[midpoint;midpoint1_1;midpoint1_2];
              midpoint=[midpoint;midpoint1_1(2:end,:);midpoint1_2(2:end,:)];
              
     corner2=corners(planecorners(diffplane,2),:);
     corner1=corners(planecorners(diffplane,1),:);
     corner3=corners(planecorners(diffplane,3),:);
     plane = [corner2   corner1-corner2   corner3-corner2];
     
     RS=specextradata(num_rayo,4:6);
     if fuente_imagen == 1
         plot3(RS(1),RS(2),RS(3),'rp')
     end
     for j=1:length(midpoint(:,1))
     line = [midpoint(j,:)  RS-midpoint(j,:)];
    point = intersectLinePlane(line, plane, 1e-14);
        x=[S(1,1) midpoint(j,1) point(1) R(1,1)];
                 y=[S(1,2) midpoint(j,2) point(2) R(1,2)];
                 z=[S(1,3) midpoint(j,3) point(3) R(1,3)];
         index=zeros(1,4);
         index(1,:)=[1 2 3 4];
         plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
        
          hold on;
     end
   elseif pathtypevec(num_rayo,2) == 100
       diffedge1=reflpaths(num_rayo,1);
       corners_edge1=edgecorners(diffedge1,1);
            corners_edge2=edgecorners(diffedge1,2);
            iv = [corners_edge1;corners_edge2];
             h = plot3(corners(iv,1),corners(iv,2),corners(iv,3));
             set(h,'LineWidth',2)
             
       diffedge2=reflpaths(num_rayo,2);
       corners_edge1=edgecorners(diffedge2,1);
     corners_edge2=edgecorners(diffedge2,2);
         iv = [corners_edge1;corners_edge2];
         h = plot3(corners(iv,1),corners(iv,2),corners(iv,3));
         set(h,'LineWidth',2)
         midpoint_edge1 = mean(corners(edgecorners(diffedge1,1:2),:)); 
        midpoint_edge2 = mean(corners(edgecorners(diffedge2,1:2),:)); 
        long_edge1 = corners(edgecorners(diffedge1,1),:) - corners(edgecorners(diffedge1,2),:);
        long_edge2 = corners(edgecorners(diffedge2,1),:) - corners(edgecorners(diffedge2,2),:);
%         pos_long_edge1=find(long_edge1 ~= 0);
%         pos_long_edge2=find(long_edge2 ~= 0);

vector1_edge1=corners(edgecorners(diffedge1,1),:)-midpoint_edge1;
        vector2_edge1=corners(edgecorners(diffedge1,2),:)-midpoint_edge1;
        vector1_edge1=vector1_edge1/norm(vector1_edge1);
        vector2_edge1=vector2_edge1/norm(vector2_edge1);
        
       vector1_edge2=corners(edgecorners(diffedge2,1),:)-midpoint_edge2;
        vector2_edge2=corners(edgecorners(diffedge2,2),:)-midpoint_edge2;
        vector1_edge2=vector1_edge2/norm(vector1_edge2);
        vector2_edge2=vector2_edge2/norm(vector2_edge2);

%         midpoint_edge1_1_1coor=midpoint_edge1(1,pos_long_edge1) + separacion_point_edges; 
%         midpoint_edge1_1_2coor=midpoint_edge1(1,pos_long_edge1) - separacion_point_edges;
%         midpoint_edge1_1_1=midpoint_edge1;
%         midpoint_edge1_1_2=midpoint_edge1;
%         midpoint_edge1_1_1(pos_long_edge1)=midpoint_edge1_1_1coor;
%         midpoint_edge1_1_2(pos_long_edge1)=midpoint_edge1_1_2coor;
%         
%         midpoint_edge2_1_1coor=midpoint_edge2(1,pos_long_edge2) + separacion_point_edges; 
%         midpoint_edge2_1_2coor=midpoint_edge2(1,pos_long_edge2) - separacion_point_edges;
%         midpoint_edge2_1_1=midpoint_edge2;
%         midpoint_edge2_1_2=midpoint_edge2;
%         midpoint_edge2_1_1(pos_long_edge2)=midpoint_edge2_1_1coor;
%         midpoint_edge2_1_2(pos_long_edge2)=midpoint_edge2_1_2coor;
    
 midpoint_edge1_1_1coor=midpoint_edge1+vector1_edge1.*separacion_point_edges; 
    midpoint_edge2_1_1coor=midpoint_edge2+vector1_edge2.*separacion_point_edges; 
    
            midpoint_edge1_1_1=midpoint_edge1;
            midpoint_edge1_1_2=midpoint_edge1;
            
            midpoint_edge2_1_1=midpoint_edge2;
            midpoint_edge2_1_2=midpoint_edge2;

        r=2;
      while sum(abs(long_edge1/2) > abs(midpoint_edge1-midpoint_edge1_1_1coor )) && sum(abs(long_edge2/2) > abs(midpoint_edge2-midpoint_edge2_1_1coor ))
        midpoint_edge1_1_1coor=midpoint_edge1_1_1(r-1,:)+vector1_edge1.*separacion_point_edges; % se añade un metro al punto medio
               midpoint_edge1_1_2coor=midpoint_edge1_1_2(r-1,:)+vector2_edge1.*separacion_point_edges;
               
               midpoint_edge2_1_1coor=midpoint_edge2_1_1(r-1,:)+vector1_edge2.*separacion_point_edges; % se añade un metro al punto medio
               midpoint_edge2_1_2coor=midpoint_edge2_1_2(r-1,:)+vector2_edge2.*separacion_point_edges;
               
               midpoint_edge1_1_1(r,:)=midpoint_edge1_1_1coor;
               midpoint_edge1_1_2(r,:)=midpoint_edge1_1_2coor;
               
               midpoint_edge2_1_1(r,:)=midpoint_edge2_1_1coor;
               midpoint_edge2_1_2(r,:)=midpoint_edge2_1_2coor;
          
%             midpoint_edge1_1_1coor=midpoint_edge1_1_1(r-1,pos_long_edge1)+separacion_point_edges; 
%             midpoint_edge1_1_2coor=midpoint_edge1_1_2(r-1,pos_long_edge1)-separacion_point_edges;
%             midpoint_edge1_1_1(r,:)=midpoint_edge1_1_1(r-1,:);
%             midpoint_edge1_1_2(r,:)=midpoint_edge1_1_2(r-1,:);
%             midpoint_edge1_1_1(r,pos_long_edge1)=midpoint_edge1_1_1coor;
%             midpoint_edge1_1_2(r,pos_long_edge1)=midpoint_edge1_1_2coor;
%             
%             midpoint_edge2_1_1coor=midpoint_edge2_1_1(r-1,pos_long_edge2)+separacion_point_edges; 
%             midpoint_edge2_1_2coor=midpoint_edge2_1_2(r-1,pos_long_edge2)-separacion_point_edges;
%             midpoint_edge2_1_1(r,:)=midpoint_edge2_1_1(r-1,:);
%             midpoint_edge2_1_2(r,:)=midpoint_edge2_1_2(r-1,:);
%             midpoint_edge2_1_1(r,pos_long_edge2)=midpoint_edge2_1_1coor;
%             midpoint_edge2_1_2(r,pos_long_edge2)=midpoint_edge2_1_2coor;
%             
            r=r+1;
        end
        midpoint_edge1_1_1(r-1,:)=[];
        midpoint_edge1_1_2(r-1,:)=[];
     %   midpoint_edge1=[midpoint_edge1;midpoint_edge1_1_1;midpoint_edge1_1_2];
             midpoint_edge1=[midpoint_edge1;midpoint_edge1_1_1(2:end,:);midpoint_edge1_1_2(2:end,:)];
             
        midpoint_edge2_1_1(r-1,:)=[];
        midpoint_edge2_1_2(r-1,:)=[];
      %  midpoint_edge2=[midpoint_edge2;midpoint_edge2_1_1;midpoint_edge2_1_2];
      midpoint_edge2=[midpoint_edge2;midpoint_edge2_1_1(2:end,:);midpoint_edge2_1_2(2:end,:)];
    %%%%%%%%%%%
        for j=1:length(midpoint_edge1(:,1))
         x=[S(1,1) midpoint_edge1(j,1) midpoint_edge2(j,1) R(1,1)];
         y=[S(1,2) midpoint_edge1(j,2) midpoint_edge2(j,2) R(1,2)];
         z=[S(1,3) midpoint_edge1(j,3) midpoint_edge2(j,3) R(1,3)];
         index=zeros(1,4);
         index(1,:)=[1 2 3 4];
         plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c)); 
    
        
          hold on;
        end
       
       
   end
    end
    
    
end
    


