Add any files (in subdirs if you want) here.


--env2matlab
--------envsrf.m : convierte .env y .srf a una estructura de datos de MATLAB y la guarda en data.mat
--------matlab2cad_v3.m : Convierte data.mat al formato usado en .CAD. Para ello dentro tiene la primera parte en variables, una segunda en la que se hace el producto cruzado para modificar orden de definición de corners en cada plano, y una tercera parte donde se escribe un archivo .CAD  
--------EDB2cross.m : función de EDT que calcula el producto cruzado (lo usa matlab2cad_v3.m)
--------data.mat: valores de IdiapRoom
--Representacion_grafica
--------BoxPlot3.m : ejemplo para dibujar gráficos 3D
--------representacion3d_1pared.m : cogiendo data.mat dibuja una única pared
--------representacion3d.m : cogiendo data.mat dibuja en 3D la sala
--------ejemplo_dibujoEDT.m : cogiendo pruebaX_1_1_edpaths.mat dibuja recinto del ejemplo de EDT, fuente, receptor,rayo directo, rayo reflejado