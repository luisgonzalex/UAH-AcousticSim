function [tiporayo,trayectorias]=getpaths(sources,receivers,infos)


%%sources: matriz [nsources,3] con las coordenadas de las fuentes
%%receivers: matriz [nreceivers,3] con las coordenas de los receptores
%%infos: estructura creada con InitSimulation()

api=1;
infos.sources=sources;
infos.receivers=receivers;

EDB2main(api,infos);


nsources=length(sources(:,1));
nreceivers=length(receivers(:,1));

eddatafile=[infos.Filepath infos.Filestem '_eddata.mat'];

for i=1:nsources
    for j=1:nreceivers
        reflpathsfile=[infos.Filepath infos.Filestem '_' num2str(i) '_' num2str(j) '_edpaths.mat'];
        [n_rays,raytype,paths]=representation_paths(eddatafile,reflpathsfile,0,i,j,infos.separacion_point_edges,api);
        eval(['tiporayo.source' num2str(i) '_receiver' num2str(j) '=raytype;']);
        eval(['trayectorias.source' num2str(i) '_receiver' num2str(j) '=paths;']);
    end
end



