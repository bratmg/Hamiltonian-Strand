%% ##################################################################
%
% FindNodeAndConnectivity.m
%
% Code to extract the node points and triangle connectivity from
% a STL file
% ###################################################################
close all;clear all;clc
% ===================================================================

% execute C code that strips the text from the STL file
system('./execute.sh');

% load the normal and vertex positions
data = load('./dataPoints.dat');

% total numbers of triangular faces
numTriangle = data(1,1);

% create arrays for normals and vertices
normal = zeros(numTriangle,3);
vertex = zeros(3*numTriangle,3);
data   = data(2:end,:); % ignore first line (contains numTriangle)

dataCount = 0;
vCount    = 0;
for i = 1:numTriangle
   
   % load the normal and vertex count
   normal(i,:)        = data(dataCount+1,:);
   vertex(vCount+1,:) = data(dataCount+2,:);
   vertex(vCount+2,:) = data(dataCount+3,:);
   vertex(vCount+3,:) = data(dataCount+4,:);
   
   % update counters
   vCount    = vCount + 3;
   dataCount = dataCount + 4;
   
end

% ===================================================================
% Find repeating nodes
% ===================================================================

% create empty nodes and connectivity arrays
coord = [];
conn  = [];

coordCount = 0;
for i = 1:vCount
   
   % obtain row and column index in conn matrix
   connRow = ceil(i/3);
   connCol = mod(i,3);
   if(connCol == 0); connCol = 3; end
   
   % isolate the current vertex
   tempVertex = vertex(i,:);
   
   if(isempty(coord)==0) % coord array is non-empty
      
      % Loop through all the elements in coord
      lenCoord  = size(coord,1);
      flagExist = 0;
      for j = 1:lenCoord
         
         % check if tempVertex matches a coord array
         if(   coord(j,1)==tempVertex(1) && ...
               coord(j,2)==tempVertex(2) && ...
               coord(j,3)==tempVertex(3))
            
            conn(connRow,connCol) = j; % nodeID
            flagExist             = 1;
            
            % if the vertex is new. add it to the coord list
            % and conn list
         end
         
         
      end
      
      % If the coord does not exist, add vertex to list of coord
      if(flagExist == 0)
         coord(coordCount + 1,:) = tempVertex;
         conn(connRow,connCol)   = coordCount + 1;
         coordCount              = coordCount + 1;
      end
      
   else % coord array is empty
      coord(coordCount + 1,:) = tempVertex;
      conn(connRow,connCol)   = coordCount + 1;
      coordCount              = coordCount + 1;
      
   end
   
   
end

% ===================================================================
% plotting and write to file
% ===================================================================
close all

figure(1)
trimesh(conn,coord(:,1),coord(:,2),coord(:,3))
axis('equal')

fid = fopen('coordRobin.dat','w+');
for i = 1:coordCount
   fprintf(fid,'%f %f %f \n',coord(i,1),coord(i,2),coord(i,3));
end
fclose(fid);

fid = fopen('connRobin.dat','w+');
connCount = size(conn,1);
for i = 1:connCount
   fprintf(fid,'%d %d %d \n',conn(i,1),conn(i,2),conn(i,3));
end
fclose(fid);



% ###################################################################
% END OF FILE
% ###################################################################