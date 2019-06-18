function dibujar=representation_diffpaths(eddatafile,reflpathsfile,order_diff,fuente_imagen,plot_source,plot_receive,separacion_point_edges)

[eddatafilepath,eddatafile,fileext] = fileparts(eddatafile);
eddatafile = [eddatafile,fileext];
Filestem = EDB2strpend(eddatafile,'_eddata');

[reflpathsfilepath,reflpathsfile,fileext] = fileparts(reflpathsfile);
reflpathsfile = [reflpathsfile,fileext];

eval(['load ',eddatafilepath,filesep,eddatafile])
eval(['load ',reflpathsfilepath,filesep,reflpathsfile])

if firstdiffrow ~=0
    
    colors=['y' 'm' 'c' 'r' 'g' 'b' 'k'];
    c=1;
    
    
    [tipos,nada]=size(mainlistguide);
    
    
    dibujar=0;
    
    for d=firstdiffrow:tipos
        
        comb=mainlistguidepattern(d,:);
        
        first_diff=mainlistguide(d,2);
        lasttype=mainlistguide(d,3);
        
        rep_diff=first_diff:lasttype;
        
        if strcmp(comb(1,1),'s') && order_diff == 1 %primero va a haber una reflexion  y luego la difraccion
            if dibujar ~= 0
                dibujar=[dibujar rep_diff];
            else
                dibujar=rep_diff;
            end
            for i=1:length(rep_diff)
                diffedge=reflpaths(rep_diff(i),2);
                
                %linea con mayor grosor:
                corners_edge1=edgecorners(diffedge,1);
                corners_edge2=edgecorners(diffedge,2);
                iv = [corners_edge1;corners_edge2];
                h = plot3(corners(iv,1),corners(iv,2),corners(iv,3));
                set(h,'LineWidth',2)
                
                diffplane=reflpaths(rep_diff(i),1);
                
                midpoint = mean(corners(edgecorners(diffedge,1:2),:)); %punto medio del edge
                long_edge = corners(edgecorners(diffedge,1),:) - corners(edgecorners(diffedge,2),:);
                
                
                vector1=corners(edgecorners(diffedge,1),:)-midpoint;
                vector2=corners(edgecorners(diffedge,2),:)-midpoint;
                vector1=vector1/norm(vector1);
                vector2=vector2/norm(vector2);
                
                midpoint1_1coor=midpoint+vector1.*separacion_point_edges;
                midpoint1_1=midpoint;
                midpoint1_2=midpoint;
                
                
                r=2;
                while sum(abs(long_edge/2) > abs(midpoint-midpoint1_1coor ))
                    midpoint1_1coor=midpoint1_1(r-1,:)+vector1.*separacion_point_edges; % se añade un metro al punto medio
                    midpoint1_2coor=midpoint1_2(r-1,:)+vector2.*separacion_point_edges;
                    
                    midpoint1_1(r,:)=midpoint1_1coor;
                    midpoint1_2(r,:)=midpoint1_2coor;
                    
                    
                    r=r+1;
                end
                midpoint1_1(r-1,:)=[];
                midpoint1_2(r-1,:)=[];
                
                midpoint=[midpoint;midpoint1_1(2:end,:);midpoint1_2(2:end,:)];
                
                corner2=corners(planecorners(diffplane,2),:);
                corner1=corners(planecorners(diffplane,1),:);
                corner3=corners(planecorners(diffplane,3),:);
                plane = [corner2   corner1-corner2   corner3-corner2];
                
                IS=specextradata(rep_diff(i),1:3);
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
                    
                end
                hold on;
                if c<length(colors)
                    c=c+1;
                else
                    c=1;
                end
            end
            
            
        elseif strcmp(comb(1,1),'d') %primero va a haber una difraccion  
            if order_diff == 2  && strcmp(comb(1,2),'d')  % difraccion de orden 2 (la maxima posible)
                if dibujar ~= 0
                    dibujar=[dibujar rep_diff];
                else
                    dibujar=rep_diff;
                end
                for i=1:length(rep_diff)
                    diffedge1=reflpaths(rep_diff(i),1);
                    corners_edge1=edgecorners(diffedge1,1);
                    corners_edge2=edgecorners(diffedge1,2);
                    iv = [corners_edge1;corners_edge2];
                    h = plot3(corners(iv,1),corners(iv,2),corners(iv,3));
                    set(h,'LineWidth',2)
                    
                    diffedge2=reflpaths(rep_diff(i),2);
                    
                    corners_edge1=edgecorners(diffedge2,1);
                    corners_edge2=edgecorners(diffedge2,2);
                    iv = [corners_edge1;corners_edge2];
                    h = plot3(corners(iv,1),corners(iv,2),corners(iv,3));
                    set(h,'LineWidth',2)
                    
                    
                    midpoint_edge1 = mean(corners(edgecorners(diffedge1,1:2),:));
                    midpoint_edge2 = mean(corners(edgecorners(diffedge2,1:2),:));
                    long_edge1 = corners(edgecorners(diffedge1,1),:) - corners(edgecorners(diffedge1,2),:);
                    long_edge2 = corners(edgecorners(diffedge2,1),:) - corners(edgecorners(diffedge2,2),:);
                    
                    
                    vector1_edge1=corners(edgecorners(diffedge1,1),:)-midpoint_edge1;
                    vector2_edge1=corners(edgecorners(diffedge1,2),:)-midpoint_edge1;
                    vector1_edge1=vector1_edge1/norm(vector1_edge1);
                    vector2_edge1=vector2_edge1/norm(vector2_edge1);
                    
                    vector1_edge2=corners(edgecorners(diffedge2,1),:)-midpoint_edge2;
                    vector2_edge2=corners(edgecorners(diffedge2,2),:)-midpoint_edge2;
                    vector1_edge2=vector1_edge2/norm(vector1_edge2);
                    vector2_edge2=vector2_edge2/norm(vector2_edge2);
                    
                    
                    
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
                        
                        
                        
                        r=r+1;
                    end
                    midpoint_edge1_1_1(r-1,:)=[];
                    midpoint_edge1_1_2(r-1,:)=[];
                    
                    midpoint_edge1=[midpoint_edge1;midpoint_edge1_1_1(2:end,:);midpoint_edge1_1_2(2:end,:)];
                    
                    midpoint_edge2_1_1(r-1,:)=[];
                    midpoint_edge2_1_2(r-1,:)=[];
                    
                    midpoint_edge2=[midpoint_edge2;midpoint_edge2_1_1(2:end,:);midpoint_edge2_1_2(2:end,:)];
                    
                    for j=1:length(midpoint_edge1(:,1))
                        x=[S(1,1) midpoint_edge1(j,1) midpoint_edge2(j,1) R(1,1)];
                        y=[S(1,2) midpoint_edge1(j,2) midpoint_edge2(j,2) R(1,2)];
                        z=[S(1,3) midpoint_edge1(j,3) midpoint_edge2(j,3) R(1,3)];
                        index=zeros(1,4);
                        index(1,:)=[1 2 3 4];
                        plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                        if c<length(colors)
                            c=c+1;
                        else
                            c=1;
                        end
                    end
                    hold on;
                    if c<length(colors)
                        c=c+1;
                    else
                        c=1;
                    end
                end
                
            elseif order_diff == 1 && strcmp(comb(1,2),'s')
                if dibujar ~= 0
                    dibujar=[dibujar rep_diff];
                else
                    dibujar=rep_diff;
                end
                for i=1:length(rep_diff)
                    diffedge=reflpaths(rep_diff(i),1);
                    
                    corners_edge1=edgecorners(diffedge,1);
                    corners_edge2=edgecorners(diffedge,2);
                    iv = [corners_edge1;corners_edge2];
                    h = plot3(corners(iv,1),corners(iv,2),corners(iv,3));
                    set(h,'LineWidth',2)
                    
                    diffplane=reflpaths(rep_diff(i),2);
                    midpoint = mean(corners(edgecorners(diffedge,1:2),:)); %punto medio del edge
                    long_edge = corners(edgecorners(diffedge,1),:) - corners(edgecorners(diffedge,2),:);
                    
                    
                    vector1=corners(edgecorners(diffedge,1),:)-midpoint;
                    vector2=corners(edgecorners(diffedge,2),:)-midpoint;
                    vector1=vector1/norm(vector1);
                    vector2=vector2/norm(vector2);
                    
                    
                    
                    midpoint1_1coor=midpoint+vector1.*separacion_point_edges; % se añade un metro al punto medio
                    midpoint1_1=midpoint;
                    midpoint1_2=midpoint;
                    
                    r=2;
                    
                    while sum(abs(long_edge/2) > abs(midpoint-midpoint1_1coor ))
                        midpoint1_1coor=midpoint1_1(r-1,:)+vector1.*separacion_point_edges; % se añade un metro al punto medio
                        midpoint1_2coor=midpoint1_2(r-1,:)+vector2.*separacion_point_edges;
                        
                        midpoint1_1(r,:)=midpoint1_1coor;
                        midpoint1_2(r,:)=midpoint1_2coor;
                        
                        r=r+1;
                    end
                    midpoint1_1(r-1,:)=[];
                    midpoint1_2(r-1,:)=[];
                    
                    
                    midpoint=[midpoint;midpoint1_1(2:end,:);midpoint1_2(2:end,:)];
                    
                    corner2=corners(planecorners(diffplane,2),:);
                    corner1=corners(planecorners(diffplane,1),:);
                    corner3=corners(planecorners(diffplane,3),:);
                    plane = [corner2   corner1-corner2   corner3-corner2];
                    
                    RS=specextradata(rep_diff(i),4:6);
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
                        
                    end
                    hold on;
                    if c<length(colors)
                        c=c+1;
                    else
                        c=1;
                    end
                end
                
                
            end
            
        end
    end
else
    dibujar=0;
end



