close all
load('data.mat')

x=[surfaces.surface0.vertex0(1),surfaces.surface0.vertex1(1),surfaces.surface0.vertex2(1),surfaces.surface0.vertex3(1)];
y=[surfaces.surface0.vertex0(2),surfaces.surface0.vertex1(2),surfaces.surface0.vertex2(2),surfaces.surface0.vertex3(2)];
z=[surfaces.surface0.vertex0(3),surfaces.surface0.vertex1(3),surfaces.surface0.vertex2(3),surfaces.surface0.vertex3(3)];
x=x*1e-3;
y=y*1e-3;
z=z*1e-3;

index=zeros(1,5);
index(1,:)=[1 2 3 4 1];
plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)))


%dibuja una linea:

x=[x(1) 0];
y=[y(1) 0];
z=[z(1) 0];
index=zeros(1,2);
index(1,:)=[1 2];
figure
plot3(x(index(1,:)),y(index(1,:)),z(index(1,:)));