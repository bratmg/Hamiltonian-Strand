%% ##################################################################
% Code to try the tessellation of a sphere
% ###################################################################
clear all;clc
% ===================================================================
maxLevel       = 0;

% ===================================================================
% Create the original icosahedron
% ===================================================================
theta  = 26.56505117707799*pi/180.0;
ctheta = cos(theta);
stheta = sin(theta);

% vertex count
vCount = 0;

% lower vertex
vCount = vCount + 1;
vertex(vCount,:) = [ 0 0 -1 ];

% lower pentagon
phi = pi/5;
for i = 1:5
   vCount           = vCount + 1;
   vertex(vCount,:) = [ ctheta*cos(phi) ctheta*sin(phi) -stheta];
   phi              = phi + 2*pi/5;
end

% upper pentagon
phi = 0;
for i = 1:5
   vCount           = vCount + 1;
   vertex(vCount,:) = [ ctheta*cos(phi) ctheta*sin(phi) +stheta];
   phi              = phi + 2*pi/5;
end

% upper vertex
vCount = vCount + 1;
vertex(vCount,:) = [ 0 0 1 ];

% ===================================================================
% Create a list of triangles IDs and the corresponding node IDs
% Node ID is the index number in vertex array
% ===================================================================
triList = [...
    1  3  2;
    1  4  3;
    1  5  4;
    1  6  5;
    1  2  6;
    2  3  8;
    3  4  9;
    4  5 10;
    5  6 11;
    6  2  7;
    2  8  7;
    3  9  8;
    4 10  9;
    5 11 10;
    6  7 11;
    7  8 12;
    8  9 12;
    9 10 12;
   10 11 12;
   11  7 12
   ];
triListCount = length(triList);
% ===================================================================
% Add another level of triangular subdivision
% ===================================================================
for n = 1:maxLevel

   numTriangle  = length(triList); % number of triangles
   eListCount   = 0;
   eList        = []; % empty edge list (contains the two nodes and the middle node)
   
   triList2     = [];
   triListCount = 0;

   for i = 1:numTriangle
      node1 = triList(i,1);
      node2 = triList(i,2);
      node3 = triList(i,3);

      
      % check if node1 and node2 exist   
      nodes                = [node1 node2];
      [node4,eList,vertex] = FindNodeInEList(nodes,eList,vertex);

      % check if node2 and node3 exist   
      nodes                = [node2 node3];
      [node5,eList,vertex] = FindNodeInEList(nodes,eList,vertex);

      % check if node1 and node3 exist   
      nodes                = [node1 node3];
      [node6,eList,vertex] = FindNodeInEList(nodes,eList,vertex);

      % create the triList array
      triList2(triListCount+1,:) = [node1 node4 node6]; 
      triList2(triListCount+2,:) = [node4 node2 node5];
      triList2(triListCount+3,:) = [node6 node5 node3];
      triList2(triListCount+4,:) = [node6 node4 node5];

      triListCount = triListCount + 4;
      
   end

   triList = triList2;

end
% ===================================================================
% plotting and write to file
% ===================================================================
%close all

figure
trimesh(triList,vertex(:,1),vertex(:,2),vertex(:,3))
axis('square')

lenVertex = length(vertex);
fid = fopen('coordSphere.dat','w+');
fprintf(fid,'%d\n',lenVertex);
for i = 1:lenVertex
   fprintf(fid,'%f %f %f \n',vertex(i,1),vertex(i,2),vertex(i,3));
end
fclose(fid);

fid = fopen('connSphere.dat','w+');
fprintf(fid,'%d\n',triListCount);
for i = 1:triListCount
   fprintf(fid,'%d %d %d \n',triList(i,1),triList(i,2),triList(i,3));
end
fclose(fid);

% ###################################################################
% END OF FILE
% ###################################################################
