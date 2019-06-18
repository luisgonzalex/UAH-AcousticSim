function rep=representation_specularpaths(eddatafile,reflpathsfile,order,fuente_imagen)
[eddatafilepath,eddatafile,fileext] = fileparts(eddatafile);
eddatafile = [eddatafile,fileext];
Filestem = EDB2strpend(eddatafile,'_eddata');

[reflpathsfilepath,reflpathsfile,fileext] = fileparts(reflpathsfile);
reflpathsfile = [reflpathsfile,fileext];

eval(['load ',eddatafilepath,filesep,eddatafile])
eval(['load ',reflpathsfilepath,filesep,reflpathsfile])

colors=['y' 'm' 'c' 'r' 'g' 'b' 'k'];
c=1;

if directsoundrow == 1
    first=2;
else
    first=1;
end


[ncombs,maxorder]=size(pathtypevec);

if firstdiffrow ~= 0
    ultima_refl=mainlistguide(firstdiffrow-1,3);
    reflexiones=first:ultima_refl;
elseif firstdiffrow == 0
    reflexiones=first:ncombs;
end

if order < maxorder
    rep=find(pathtypevec(reflexiones,order+1)==0);
    busc=find(pathtypevec(reflexiones,order)==0);
    if busc ~= 0   % para orden 2 sirve para 3??¿
        repet=[rep' busc'];
        repet=repet';
        [unicos,ind_unicos]=unique(repet);
        ind_repetidos = setdiff(1:length(repet),ind_unicos);
        repet=repet(ind_repetidos,1);
        rep=rep(length(repet)+1:end);
    end
elseif order==maxorder
    rep=find(pathtypevec(reflexiones,order)==115);
end

%specextradata=full(specextradata);

if fuente_imagen == 1
    for i=1:length(rep)
        x = specextradata(rep(i),1);
        y = specextradata(rep(i),2);
        z = specextradata(rep(i),3);
        plot3(x,y,z,'rp')
    end
end

for p=1:length(rep)
    % order=length(find(pathtypevec(p,:))~=0);
    for h=1:order
        if h+1 > order
            x=[R(1,1) specextradata(rep(p),4+(h-1)*3)];
            y=[R(1,2) specextradata(rep(p),5+(h-1)*3)];
            z=[R(1,3) specextradata(rep(p),6+(h-1)*3)];
            index=zeros(1,2);
            index(1,:)=[1 2];
            plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
            hold on;
            if h == 1
                x=[S(1,1) specextradata(rep(p),4+(h-1)*3)];
                y=[S(1,2) specextradata(rep(p),5+(h-1)*3)];
                z=[S(1,3) specextradata(rep(p),6+(h-1)*3)];
                index=zeros(1,2);
                index(1,:)=[1 2];
                plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                hold on;
            else
                x=[specextradata(rep(p),4+(h-2)*3) specextradata(rep(p),4+(h-1)*3)];
                y=[specextradata(rep(p),5+(h-2)*3) specextradata(rep(p),5+(h-1)*3)];
                z=[specextradata(rep(p),6+(h-2)*3) specextradata(rep(p),6+(h-1)*3)];
                index=zeros(1,2);
                index(1,:)=[1 2];
                plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                hold on;
            end
        else
            if h==1
                x=[S(1,1) specextradata(rep(p),4+(h-1)*3)];
                y=[S(1,2) specextradata(rep(p),5+(h-1)*3)];
                z=[S(1,3) specextradata(rep(p),6+(h-1)*3)];
                index=zeros(1,2);
                index(1,:)=[1 2];
                plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                hold on;
                if order >= 3   %dibuje rayos intermedios
                    x=[specextradata(rep(p),4+(h-1)*3) specextradata(rep(p),4+(h)*3)];
                    y=[specextradata(rep(p),5+(h-1)*3) specextradata(rep(p),5+(h)*3)];
                    z=[specextradata(rep(p),6+(h-1)*3) specextradata(rep(p),6+(h)*3)];
                    index=zeros(1,2);
                    index(1,:)=[1 2];
                    plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                end
            else
                x=[specextradata(rep(p),4+(h-1)*3) specextradata(rep(p),4+(h)*3)];
                y=[specextradata(rep(p),5+(h-1)*3) specextradata(rep(p),5+(h)*3)];
                z=[specextradata(rep(p),6+(h-1)*3) specextradata(rep(p),6+(h)*3)];
                index=zeros(1,2);
                index(1,:)=[1 2];
                plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)),colors(c));
                hold on;
                
            end
            
            
        end
    end
    if c<7
        c=c+1;
    else
        c=1;
    end
end