% Edge diffraction toolbox
% Version EDB2 24 April 2015
%
% From version EDB1, all functions started with EDB1 rather than EDB. The B1, B2 etc will indicate
% version number.
%
% Main programs
%  EDB2main       	  - Calculates the specular and edge diffraction IRs and saves them in files.
%
% Read CAD file and geometry pre-processing
%  EDB2readsetup       - Runs the EDsetup m-file and saves all the settings in a mat-file.
%  EDB2readcad         - Reads a file of type .CAD and saves all the geometry data in a mat-file. 
%  EDB2readac          - EDB2readac - Reads a file of type .AC (made by e.g. Invis AC3D) and saves all the geometry data in a mat-file. 
%  EDB2ed2geox         - Calculates 2nd- and higher-order edge-related geom. parameters.
%  EDB2edgeox          - Calculates some plane- and edge-related geometrical parameters.
%  EDB2SorRgeo         - Calculates some source- or receiver-related geometrical parameters.
%
% Find sound paths
%  EDB2findpathsISESx  - Finds all possible paths that include direct sound, specular, diffraction.
%  EDB2diff2ISES       - Gives list of paths that includes a 2nd-order diff. combination.
%  EDB2diffISESx       - Gives list of paths that includes a 1st-order diff. combination.
%  EDB2diffNISES       - Gives list of paths that includes an N-order diff. combination.
%  EDB2directsound     - Checks if the direct sound is valid.
%  EDB2findISEStree    - Constructs a list ("an ISES tree") of all possible specular-diffraction combinations.
%  EDB2findis          - Returns the image source coordinates.
%  EDB2checkobstrpaths - Checks obstruction for a list of paths, or check if specular reflections are valid.
%  EDB2checkobstr_pointtoedge      - Checks obstruction for a list of point-to-edge paths.
%  EDB2checkobstr_edgetoedge       - Checks obstruction for a list of edge-to-edge paths.
%  EDB2chkISvisible    - Checks if paths from a set of IS to a set of R pass through their refl. planes.
%  EDB2chkISvisiblex   - Checks if paths from a set of IS to a set of R pass through their refl. planes. Extended version.
%  EDB2getedgepoints   - Calculates a number of edge coordinates.
%  EDB2infrontofplane  - Checks if a point is in front of, in plane with, or behind a list of planes. 
%  EDB2poinpla         - Detects if one or more points are inside a number of finite planes. 
%  EDB2poinplax        - Detects if one or more points are inside a number of finite planes. Extended version.
%  EDB2speculISES      - Finds the valid specular reflections by checking visibility and obstruction.
%
% Create impulse responses
%  EDB2makeirs         - Constructs IRs from a list of paths in an edpathsfile.
%  EDB2betaoverml      - Integrand function which is called for num. int. of ED IR.
%  EDB2irfromslotlist  - Creates an IR from a list of sample numbers and amplitudes
%  EDB2wedge1st_int    - Computes the 1st-order diffraction impulse response.
%  EDB2wedge2nd 		  - Computes the second-order diffraction impulse responses.
%  EDB2wedgeN 		  - Computes the Nth-order diffraction impulse responses.
%  EDB2quadstep        - Recursive numerical integration routine for EDB2wedge1st_int.
%
% Create transfer functions
%  EDB2maketfs         - Constructs TFs from a list of paths in an edpathsfile.
%  EDB2betaoverml_fd   - Integrand function which is called for num. int. of ED TF, first order.
%  EDB2betaoverml2_fd  - Integrand function which is called for num. int. of ED TF, second order.
%  EDB2wedge1st_fd     - Computes the 1st-order diffraction transfer function.
%  EDB2wedge2nd_fd 	   - Computes the second-order diffraction transfer function.
%  EDB2wedge3rd_intcorefd   - Integrand function for third-order ED TF. Such integration is
%						 extremely slow and therefore third-order diffraction is shut off in
%					     the toolbox.
%
% Various plot functions
%  EDB2makemovie       - Makes a movie based on IRs.
%  EDB2makesnapshot    - Makes a snapshot of a sound field based on IRs.
%  EDB2plotmodel       - Plots a model which is given in an eddatafile.
%  EDB2plotstackedIRs  - Plots a number of IRs in a stacked way.
%  EDB2plotstackedTFs  - Plots a number of TFs in a stacked way.
%  EDB2plotpath 		  - Plots a model which is given in an eddatafile and one or more of the paths in an edreflpathsfile.
%
% Various string functions
%  EDB2extrnums        - Extracts the numerical values in a text-string, separated by a given delimiter.
%  EDB2fistr           - Reads a text file, line by line, until a given textstring is found.
%  EDB2myvalueinput    - Prints a text on the screen and waits for a numerical to be input.
%  EDB2strpend         - Removes a specified ending from a filename.
%
% Various numerical functions
%  EDB2calcdist        - Gives distances between sets of points in 3D space.
%  EDB2compress7       - 
%  EDB2coordtrans1     - Transforms one set of cartesian coordinates to edge-related cylindrical coordinates.
%  EDB2coordtrans2     - Transforms two sets of cartesian coordinates to edge-related cylindrical coordinates.
%  EDB2creindexmatrix  - Creates a matrix with index numbers.
%  EDB2cross           - Stripped down version of Matlab's built in function cross.
%
% Function for editing paths
%  EDB2editpathlist    - Does some semi-automatic editing of an edpathsfile.
