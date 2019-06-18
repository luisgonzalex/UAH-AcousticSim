function representation_src_rcv(nsources,nreceivers,sources,receivers)
% representa fuentes y/o receptores en el recinto representado con representation_room.m
% fuentes con circulo en rojo
% receptores con circulo en verde

for i=1:nsources
    plot3(sources(i,1),sources(i,2),sources(i,3),'ro','LineWidth',2)  %dibuja la fuente
    hold on
end
for i=1:nreceivers
    plot3(receivers(i,1),receivers(i,2),receivers(i,3),'go','LineWidth',2)  %dibuja receptor
    hold on
end


