function [n_rays,raytype,paths]=representation_paths(eddatafile,reflpathsfile,fuente_imagen,plot_source,plot_receive,separacion_point_edges,api)

[eddatafilepath,eddatafile,fileext] = fileparts(eddatafile);
eddatafile = [eddatafile,fileext];
Filestem = EDB2strpend(eddatafile,'_eddata');

[reflpathsfilepath,reflpathsfile,fileext] = fileparts(reflpathsfile);
reflpathsfile = [reflpathsfile,fileext];

eval(['load ',eddatafilepath,filesep,eddatafile])
eval(['load ',reflpathsfilepath,filesep,reflpathsfile])

if api == 1
    ind_paths=1;
    path1=0;
else
    raytype=0;
    paths=0;
end

colors=['y' 'm' 'c' 'r' 'g' 'b' 'k'];
c=1;

n_rays=length(pathtypevec);

if firstdiffrow ~= 0
    difraccion = mainlistguide(firstdiffrow,2);
else
    difraccion=n_rays+1;
end

%specextradata=full(specextradata);

for p=1:difraccion-1
    order=length(find(pathtypevec(p,:))~=0);
    for h=1:order
        if fuente_imagen == 1
            if pathtypevec(p,h) == 115
                x = specextradata(p,(h*3)-2);
                y = specextradata(p,(h*3)-1);
                z = specextradata(p,(h*3));
                plot3(x,y,z,'rp');
                hold on;
            end
        end
        if pathtypevec(p,h) == 102 %rayo directo
            x=[specextradata(p,1) specextradata(p,4)];
            y=[specextradata(p,2) specextradata(p,5)];
            z=[specextradata(p,3) specextradata(p,6)];
            index=zeros(1,2);
            index(1,:)=[1 2];
            if api == 1
                path1=[full(x(1)) full(y(1)) full(z(1)) full(x(2)) full(y(2)) full(z(2))];
            else
                plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                hold on;
            end
        elseif pathtypevec(p,h)== 115  %reflexión especular
            if h+1 > order
                x=[R(1,1) specextradata(p,4+(h-1)*3)];
                y=[R(1,2) specextradata(p,5+(h-1)*3)];
                z=[R(1,3) specextradata(p,6+(h-1)*3)];
                index=zeros(1,2);
                index(1,:)=[1 2];
                if api == 1
                    if path1==0
                        if order ~= 1
                            path1=[full(x(2)) full(y(2)) full(z(2))  R(1,1) R(1,2) R(1,3)];
                        end
                    else
                        if order == 2
                            path1=[path1 full(x(2)) full(y(2)) full(z(2)) R(1,1) R(1,2) R(1,3)];
                        else
                            path1=[path1 R(1,1) R(1,2) R(1,3)];
                        end
                    end
                else
                    plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                    hold on;
                end
                if h == 1
                    x=[S(1,1) specextradata(p,4+(h-1)*3)];
                    y=[S(1,2) specextradata(p,5+(h-1)*3)];
                    z=[S(1,3) specextradata(p,6+(h-1)*3)];
                    index=zeros(1,2);
                    index(1,:)=[1 2];
                    if api == 1
                        if path1==0
                            path1=[S(1,1) S(1,2) S(1,3) full(x(2)) full(y(2)) full(z(2)) R(1,1) R(1,2) R(1,3)];
                        else
                            
                        end
                    else
                        plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                        hold on;
                    end
                else
                    x=[specextradata(p,4+(h-2)*3) specextradata(p,4+(h-1)*3)];
                    y=[specextradata(p,5+(h-2)*3) specextradata(p,5+(h-1)*3)];
                    z=[specextradata(p,6+(h-2)*3) specextradata(p,6+(h-1)*3)];
                    index=zeros(1,2);
                    index(1,:)=[1 2];
                    if api == 1
                    else
                        plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                        hold on;
                    end
                end
            else
                if h==1
                    x=[S(1,1) specextradata(p,4+(h-1)*3)];
                    y=[S(1,2) specextradata(p,5+(h-1)*3)];
                    z=[S(1,3) specextradata(p,6+(h-1)*3)];
                    index=zeros(1,2);
                    index(1,:)=[1 2];
                    if api == 1
                        if path1==0
                            path1=[full(x(1)) full(y(1)) full(z(1)) full(x(2)) full(y(2)) full(z(2))];
                        else
                            path1=[path1 full(x(1)) full(y(1)) full(z(1)) full(x(2)) full(y(2)) full(z(2))];
                        end
                    else
                        plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                        hold on;
                    end
                    if order >= 3   %dibuje rayos intermedios
                        x=[specextradata(p,4+(h-1)*3) specextradata(p,4+(h)*3)];
                        y=[specextradata(p,5+(h-1)*3) specextradata(p,5+(h)*3)];
                        z=[specextradata(p,6+(h-1)*3) specextradata(p,6+(h)*3)];
                        index=zeros(1,2);
                        index(1,:)=[1 2];
                        if api == 1
                            if path1==0
                                path1=[ full(x(2)) full(y(2)) full(z(2))];
                            else
                                path1=[path1 full(x(2)) full(y(2)) full(z(2))];
                            end
                        else
                            plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                        end
                    end
                else
                    x=[specextradata(p,4+(h-1)*3) specextradata(p,4+(h)*3)];
                    y=[specextradata(p,5+(h-1)*3) specextradata(p,5+(h)*3)];
                    z=[specextradata(p,6+(h-1)*3) specextradata(p,6+(h)*3)];
                    index=zeros(1,2);
                    index(1,:)=[1 2];
                    if api == 1
                        if path1==0
                            path1=[full(x(2)) full(y(2)) full(z(2))];
                        else
                            path1=[path1 full(x(2)) full(y(2)) full(z(2))];
                        end
                    else
                        plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                        hold on;
                    end
                end
            end
        end
    end
    
    c=c+1;
    if c > 7
        c=1;
    end
    
    if api==1
        paths(ind_paths,:)={path1};
        ind_paths=ind_paths+1;
        path1=0;
    end
end

raytype=pathtypevec(1:difraccion-1,:);
ind_raytype=difraccion;
ind_pathtypevec=difraccion;


% DIFRACCIONES:

if firstdiffrow ~= 0
    [tipos,nada]=size(mainlistguide);
    
    for d=firstdiffrow:tipos
        
        comb=mainlistguidepattern(d,:);
        
        first_diff=mainlistguide(d,2);
        lasttype=mainlistguide(d,3);
        rep_diff=first_diff:lasttype;
        
        if strcmp(comb(1,1),'s') %primero va a haber una reflexion  y luego la difraccion
            for i=1:length(rep_diff)
                
                diffedge=reflpaths(rep_diff(i),2);
                
                corners_edge1=edgecorners(diffedge,1);
                corners_edge2=edgecorners(diffedge,2);
                iv = [corners_edge1;corners_edge2];
                if api == 0
                    h = plot3(corners(iv,1),corners(iv,2),corners(iv,3));
                    set(h,'LineWidth',2)
                end
                
                diffplane=reflpaths(rep_diff(i),1);
                midpoint = mean(corners(edgecorners(diffedge,1:2),:));
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
                    if api == 1
                        
                        path1=[x(1) y(1) z(1) x(2) y(2) z(2) x(3) y(3) z(3) x(4) y(4) z(4)];
                        
                        paths(ind_paths,:)={path1};
                        path1=0;
                        ind_paths=ind_paths+1;
                        raytype(ind_raytype,:)=pathtypevec(ind_pathtypevec,:);
                        ind_raytype=ind_raytype+1;
                        
                    else
                        index=zeros(1,4);
                        index(1,:)=[1 2 3 4];
                        plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                    end
                    
                end
                if api==1
                    ind_pathtypevec=ind_pathtypevec+1;
                end
                if c<length(colors)
                    c=c+1;
                else
                    c=1;
                end
                hold on;
            end
            
        elseif strcmp(comb(1,1),'d') %primero va a haber una difraccion  y luego reflexion
            if strcmp(comb(1,2),'s')
                for i=1:length(rep_diff)
                    diffedge=reflpaths(rep_diff(i),1);
                    
                    corners_edge1=edgecorners(diffedge,1);
                    corners_edge2=edgecorners(diffedge,2);
                    iv = [corners_edge1;corners_edge2];
                    if api == 0
                        h = plot3(corners(iv,1),corners(iv,2),corners(iv,3));
                        
                        set(h,'LineWidth',2)
                    end
                    diffplane=reflpaths(rep_diff(i),2);
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
                        
                        midpoint1_1coor=midpoint1_1(r-1,:)+vector1.*separacion_point_edges;
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
                        if api == 1
                            
                            path1=[x(1) y(1) z(1) x(2) y(2) z(2) x(3) y(3) z(3) x(4) y(4) z(4)];
                            
                            paths(ind_paths,:)={path1};
                            path1=0;
                            ind_paths=ind_paths+1;
                            raytype(ind_raytype,:)=pathtypevec(ind_pathtypevec,:);
                            ind_raytype=ind_raytype+1;
                        else
                            index=zeros(1,4);
                            index(1,:)=[1 2 3 4];
                            plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                        end
                    end
                    if api==1
                        ind_pathtypevec=ind_pathtypevec+1;
                    end
                    hold on;
                    if c<length(colors)
                        c=c+1;
                    else
                        c=1;
                    end
                end
                
            elseif strcmp(comb(1,2),'d')
                for i=1:length(rep_diff)
                    diffedge1=reflpaths(rep_diff(i),1);
                    
                    corners_edge1=edgecorners(diffedge1,1);
                    corners_edge2=edgecorners(diffedge1,2);
                    iv = [corners_edge1;corners_edge2];
                    h = plot3(corners(iv,1),corners(iv,2),corners(iv,3));
                    if api == 0
                        set(h,'LineWidth',2)
                    end
                    
                    diffedge2=reflpaths(rep_diff(i),2);
                    
                    corners_edge1=edgecorners(diffedge2,1);
                    corners_edge2=edgecorners(diffedge2,2);
                    iv = [corners_edge1;corners_edge2];
                    if api == 0
                        h = plot3(corners(iv,1),corners(iv,2),corners(iv,3));
                        set(h,'LineWidth',2)
                    end
                    
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
                        
                        midpoint_edge1_1_1coor=midpoint_edge1_1_1(r-1,:)+vector1_edge1.*separacion_point_edges;
                        midpoint_edge1_1_2coor=midpoint_edge1_1_2(r-1,:)+vector2_edge1.*separacion_point_edges;
                        
                        midpoint_edge2_1_1coor=midpoint_edge2_1_1(r-1,:)+vector1_edge2.*separacion_point_edges;
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
                        if api == 1
                            
                            path1=[x(1) y(1) z(1) x(2) y(2) z(2) x(3) y(3) z(3) x(4) y(4) z(4)];
                            
                            paths(ind_paths,:)={path1};
                            path1=0;
                            ind_paths=ind_paths+1;
                            raytype(ind_raytype,:)=pathtypevec(ind_pathtypevec,:);
                            ind_raytype=ind_raytype+1;
                        else
                            index=zeros(1,4);
                            index(1,:)=[1 2 3 4];
                            plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                            hold on;
                        end
                        
                    end
                    if api==1
                        ind_pathtypevec=ind_pathtypevec+1;
                    end
                    if c<length(colors)
                        c=c+1;
                    else
                        c=1;
                    end
                end
            else    % difracción de un orden (100 0)
                for i=1:length(rep_diff)
                    diffedge1=reflpaths(rep_diff(i),1);
                    
                    corners_edge1=edgecorners(diffedge1,1);
                    corners_edge2=edgecorners(diffedge1,2);
                    iv = [corners_edge1;corners_edge2];
                    h = plot3(corners(iv,1),corners(iv,2),corners(iv,3));
                    if api == 0
                        set(h,'LineWidth',2)
                    end
                    
              
                    midpoint_edge1 = mean(corners(edgecorners(diffedge1,1:2),:));
                    long_edge1 = corners(edgecorners(diffedge1,1),:) - corners(edgecorners(diffedge1,2),:);

                    
                    
                    vector1_edge1=corners(edgecorners(diffedge1,1),:)-midpoint_edge1;
                    vector2_edge1=corners(edgecorners(diffedge1,2),:)-midpoint_edge1;
                    vector1_edge1=vector1_edge1/norm(vector1_edge1);
                    vector2_edge1=vector2_edge1/norm(vector2_edge1);
                    
                    midpoint_edge1_1_1coor=midpoint_edge1+vector1_edge1.*separacion_point_edges;

                    
                    midpoint_edge1_1_1=midpoint_edge1;
                    midpoint_edge1_1_2=midpoint_edge1;
                    

                    
                    r=2;
                    while sum(abs(long_edge1/2) > abs(midpoint_edge1-midpoint_edge1_1_1coor )) 
                        
                        midpoint_edge1_1_1coor=midpoint_edge1_1_1(r-1,:)+vector1_edge1.*separacion_point_edges;
                        midpoint_edge1_1_2coor=midpoint_edge1_1_2(r-1,:)+vector2_edge1.*separacion_point_edges;
                        
                        midpoint_edge1_1_1(r,:)=midpoint_edge1_1_1coor;
                        midpoint_edge1_1_2(r,:)=midpoint_edge1_1_2coor;

                        r=r+1;
                    end
                    midpoint_edge1_1_1(r-1,:)=[];
                    midpoint_edge1_1_2(r-1,:)=[];
                    
                    midpoint_edge1=[midpoint_edge1;midpoint_edge1_1_1(2:end,:);midpoint_edge1_1_2(2:end,:)];
  
                    for j=1:length(midpoint_edge1(:,1))
                        x=[S(1,1) midpoint_edge1(j,1) R(1,1)];
                        y=[S(1,2) midpoint_edge1(j,2) R(1,2)];
                        z=[S(1,3) midpoint_edge1(j,3) R(1,3)];
                        if api == 1
                            
                            path1=[x(1) y(1) z(1) x(2) y(2) z(2) x(3) y(3) z(3)];
                            
                            paths(ind_paths,:)={path1};
                            path1=0;
                            ind_paths=ind_paths+1;
                            raytype(ind_raytype,:)=pathtypevec(ind_pathtypevec,:);
                            ind_raytype=ind_raytype+1;
                        else
                            index=zeros(1,3);
                            index(1,:)=[1 2 3];
                            plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                            hold on;
                        end
                        
                    end
                    if api==1
                        ind_pathtypevec=ind_pathtypevec+1;
                    end
                    if c<length(colors)
                        c=c+1;
                    else
                        c=1;
                    end
                end
                
            end
            
        end
    end
end






