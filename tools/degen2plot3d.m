% Program to convert camber geometry in *.m degenGeom file 
% generated by OpenVSP to Plot3D format

function [Xright,Yright,Zright] = degen2plot3d(degenGeomFile)

% Remove .m extension if present in degenGeomFile
if (degenGeomFile(end-1:end) == '.m')
  commandName = degenGeomFile(1:end-2);
else
  commandName = degenGeomFile;
end

% Run degenGeom matlab file to get dataset
run(commandName);

% Check no. of geometries
nGeo = size(degenGeom,2);

% Extract camber surface coordinates of right wing
Xright = degenGeom(1).plate.x;
Yright = degenGeom(1).plate.y;
Zright = degenGeom(1).plate.zCamber;

% Flip order of X coordinates of right wing
Xright = flip(Xright,2);
Zright = flip(Zright,2);

% if (nGeo == 2)
%   % Extract camber surface coordinates of left wing
%   Xleft = degenGeom(2).plate.x;
%   Yleft = degenGeom(2).plate.y;
%   Zleft = degenGeom(2).plate.zCamber;
%   
%   % Flip order of X coordinates of left wing
%   Xleft = flip(Xleft,2);
%   Zleft = flip(Zleft,2);
%   
%   % Flip order of Y coordinates of left wing
%   Yleft = flip(Yleft,1);
%   Xleft = flip(Xleft,1);
%   Zleft = flip(Zleft,1);
% end

% Mirror right wing to get left wing

nx = size(Xright,1);
ny = size(Xright,2);
nz = 1;

% Write to PLOT3D format
fileID = fopen([commandName '.xyz'],'w');
fprintf(fileID,'%u %u %u\n',nx,ny,nz)
fprintf(fileID,'%15.7f',Xright)
fprintf(fileID,'%15.7f',Yright)
fprintf(fileID,'%15.7f',Zright)
fclose(fileID);

return;
