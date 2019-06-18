%dibuja BENCH_Lsp_Kessel.cad y rayo directo y reflejado del ejemplo

close all
clear all
load('prueba3_1_1_edpaths.mat')


pathtypevec_char=[char(pathtypevec(:,1)) char(pathtypevec(:,2))];
specextradata_f=full(specextradata);
%CORNERS

%    1   -0.2000000   -0.4400000   -0.3200000
%    2   0.2000000   -0.4400000   -0.3200000
%    3   0.2000000   0.2000000   -0.3200000
%    4   -0.2000000   0.2000000   -0.3200000
%    5   -0.2000000   -0.4400000   0.0000000
%    6   0.2000000   -0.4400000   0.0000000
%    7   0.2000000   0.2000000   0.0000000
%    8   -0.2000000   0.2000000   0.0000000
c= [-0.2000000   -0.4400000   -0.3200000; 0.2000000   -0.4400000   -0.3200000;0.2000000   0.2000000   -0.3200000;-0.2000000   0.2000000   -0.3200000; -0.2000000   -0.4400000   0.0000000;0.2000000   -0.4400000   0.0000000;0.2000000   0.2000000   0.0000000;-0.2000000   0.2000000   0.0000000];

%PLANES

%   1 / /RIGID
p=[   1  4  3  2  ;5  6  7  8 ;1 2 6 5 ;3 4 8 7;2 3 7 6;1 5 8 4]
% 
%   2 / /RIGID
%   5  6  7  8  
% 
%   3 / /RIGID
%   1 2 6 5 
% 
%   4 / /RIGID
%   3 4 8 7
% 
%   5 / /RIGID
%   2 3 7 6
%   
%   6 / /RIGID
%   1 5 8 4
numsurfaces=6;


for i=0:numsurfaces-1
     for j=0:3 
       x(j+1)=c(p(i+1,j+1),1);
       y(j+1)=c(p(i+1,j+1),2);
      z(j+1)=c(p(i+1,j+1),3);
     end

    index=zeros(1,5);
    index(1,:)=[1 2 3 4 1];
    plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)))
    hold all
end
plot3(S(1,1),S(1,2),S(1,3),'ro')  %dibuja la fuente
plot3(R(1,1),R(1,2),R(1,3),'go')  %dibuja receptor

%dibuja rayo directo:
image_source_direct=[0 0 1.000e-4];
image_receiver_direct=[0 0 10];
x=[image_source_direct(1) image_receiver_direct(1)];
y=[image_source_direct(2) image_receiver_direct(2)];
z=[image_source_direct(3) image_receiver_direct(3)];
index=zeros(1,2);
index(1,:)=[1 2];
plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)));

%dibuja rayo reflejado¿?:
image_source_specular=[0 0 -1.000e-04];
plot3(image_source_specular(1,1),image_source_specular(1,2),image_source_specular(1,3),'gx')
image_receiver_specular=[0 0 0];
x=[S(1) image_receiver_specular(1) R(1)];
y=[S(1) image_receiver_specular(2) R(2)];
z=[S(1) image_receiver_specular(3) R(3)];
index=zeros(1,3);
index(1,:)=[1 2 3];
plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)));

%plano en que refleja: al estar el receptor fuera del recinto...
plot3(c(5,1),c(5,2),c(5,3),'bx');
plot3(c(6,1),c(6,2),c(6,3),'bx');
plot3(c(7,1),c(7,2),c(7,3),'bx');
plot3(c(8,1),c(8,2),c(8,3),'bx');