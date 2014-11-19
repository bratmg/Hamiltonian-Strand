%% ##################################################################
% PrepareLoops.m
% Written by Beatrice Roget
%
% Note: There is a domain specific condition in Tedges (0.61 and 1.1)
%       QEVec1 and QEVec2 - should it be upto 35 instead of 30
%
% NOTE: There is a special correction for spherical surface which
%       pushes the points onto the sphere surface
%       Search for 'sphCorr'
% ###################################################################
close all;clear all;clc
% ===================================================================
% Read in  coord and connectivity files
% ===================================================================
disp('Reading coord and connectivity files...');

pts = load('../robin/coordRobin.dat');
%pts = load('naca1/naca_coord.dat');
%pts = load('../SphereTessellation/coordSphere.dat');
%pts = load('../TestGrid/coordSquareGrid.dat');
nn  = length(pts(:,1)); % number of nodes
%pts = [pts zeros(nn,1)];
x   = pts(:,1);
y   = pts(:,2);
z   = pts(:,3);

% The three numbers in a given line correspond to the
% vertices of a triangle. Numbers are referenced to
% the lines in naca_coord.dat.
tri = load('../robin/connRobin.dat');
%tri = load('naca1/naca_conn.dat');
%tri = load('../SphereTessellation/connSphere.dat');
%tri = load('../TestGrid/connSquareGrid.dat');
nt  = length(tri(:,1)); % number of triangles

disp('==================================');
disp(' ENSURE THAT iSURFACE IS ACCURATE ');
disp('==================================');

figure
hold on
h=trimesh(tri,x,y,z);
set(h,'EdgeColor','r');
axis equal

% ===================================================================
% Invert connectivity list
% ===================================================================
disp('Inverting connectivity list...');

% reshaped as row-wise. Because of the transpose used on tri.
List1  = reshape(tri(1:nt,1:3)',3*nt,1);
n1     = nt; % total number of triangles
n2     = nn; % total number of nodes
Index1 = (1:3:3*n1+1);

% Given the starting and ending index of a triangle, the nodes of that
% triangle are known. The following function, however, outputs data
% such that, the triangles containing a given node are now listed
[List2,Index2] = ReverseList(List1,3*nt,n1,n2,Index1);

% ===================================================================
% Go through vertices and find loops
% (cell loops and vertex loops)
% ===================================================================
ikc = 0; ikv = 0;

Icloop = ones(nn,1);
Ivloop = ones(nn,1);
disp('Finding loops (of cells and vertices)...');

% Loop over total number of nodes
for i=1:nn
   
   k   = 99;            % not sure why 99 (guess just a large number)
   k1  = k;
   k2  = k+3;
   i1  = Index2(i);     % start Index2 for a given node (vertex)
   i2  = Index2(i+1)-1; % end Index2 for a given node (based on number of connected edges)
   nc  = i2-i1+1;       % number of cells (edges/connections)
   
   % Creates the list of cell connectivities
   % e.g.: triangles connected to node i
   %        128 100   1
   %        157 128   1
   %        157   1   2
   CC  = tri(List2(i1:i2),:);
   
   % Contains the memory for the OTHER two node IDs for a given triangle
   CC2 = zeros(nc,2);
   
   % Loop for each of the connected edges (cells)
   for j=1:nc
      
      % if the 1st point of three is the vertex itself
      % the other two are 2 and 3
      if (CC(j,1)==i)
         CC2(j,1) = CC(j,2);
         CC2(j,2) = CC(j,3);
         continue
      end
      
      % if the 2nd point of three is the vertex itself
      % the other two are 3 and 1
      if (CC(j,2)==i)
         CC2(j,1) = CC(j,3);
         CC2(j,2) = CC(j,1);
         continue
      end
      
      % Else condition
      CC2(j,1) = CC(j,1);
      CC2(j,2) = CC(j,2);
   end % end j = 1:nc
   
   % If the number of connections is only 1
   % Haven't yet encountered this condition
   if (nc==1)
      ncl = 1;
      nv  = 2;
      CLoops(ikc+1:ikc+ncl) = List2(i1);
      VLoops(ikv+1:ikv+nv)  = CC2(1,1:2);
   else % if more that 1 connection to a node (more common)
      
      % cloop = Cell loop
      cloop1    = [];        % begin as empty loop
      cloop2(1) = List2(i1); % the first element is the node itself
      kloop1    = 0;         % counter for number of elements in cloop1
      kloop2    = 1;         % counter for number of elements in cloop2
      
      % V1 and V2 vertex hardcoded to 1st row of CC2
      V1        = CC2(1,1);
      V2        = CC2(1,2);
      EV(k+1)   = V1;  % end vertices for the loop (k was 99)
      EV(k+2)   = V2;  % The last two end vertices are CC2(1,1) and CC2(1,2)
      iact      = ones(1,nc);%ones(1:nc);
      iact(1)   = 0;
      nact      = nc-1; % number of active edge connections
      icont     = 1;    % flag for continuation
      
      % while continuation flag is 1
      % This is a fancy algorithm to build the vertices from a list of
      % CC values. The chain builds to left if V1 and to the right if V2
      % For e.g., stop at i == 907 (look for ppt in the folder explaining
      % this algorithm)
      while (icont==1)
         % loop from 2 to remaining number of edges (cells)
         for j=2:nc
            if (iact(j)==1)
               % loop for other two vertices
               for kk=1:2
                  
                  % The (3-kk) basically switches to the other column
                  % of the same row in CC2. If kk = 1, switches to 2,
                  % is kk = 2, switches to 1
                  if (CC2(j,kk) == V1)
                     V1               = CC2(j,3-kk);   % switch V1 to other column in same row
                     EV(k1)           = V1;            % fill EV from element 99 towards 1
                     k1               = k1-1;          % decrement k1
                     iact(j)          = 0;             % change activation for that nc columns
                     nact             = nact-1;        % decrement nact by 1
                     cloop1(kloop1+1) = List2(i1+j-1); % append the node ID into cloop1
                     kloop1           = kloop1+1;      % increment kloop1 counter
                     break;
                  elseif (CC2(j,kk) == V2)
                     V2               = CC2(j,3-kk);
                     EV(k2)           = V2;            % fill EV from element 102 towards Inf
                     k2               = k2+1;          % increment k2
                     iact(j)          = 0;             % change activation for that nc columns
                     nact             = nact-1;        % decrement nact by 1
                     cloop2(kloop2+1) = List2(i1+j-1); % append the node ID into cloop2
                     kloop2           = kloop2+1;      % increment kloop2 counter
                     break;
                  end
                  
               end % for loop kk
            end % if iact
            
            % Termination condition for icont based of it there are
            % any remaining active edges left (nact)
            if (nact==0)
               icont=0;
               break;
            end
            
         end % for loop (j)
      end % while condition (icont)
      
      % Cell loops
      ncl                   = kloop1 + kloop2; % number of cells
      CLoops(ikc+1:ikc+ncl) = [cloop1(kloop1:-1:1) cloop2(1:kloop2)];
      
      % Vertex loops
      nv                    = k2-k1-1; % number of vertices
      VLoops(ikv+1:ikv+nv)  = EV(k1+1:k2-1); % select the appropriate loops
      
   end
   
   ikc         = ikc+ncl;       % increment total cell counter
   ikv         = ikv+nv;        % increment total vertex (node) counter
   Icloop(i+1) = Icloop(i)+ncl; % index of cell loop
   Ivloop(i+1) = Ivloop(i)+nv;  % index of vertex loop
   
end

% ===================================================================
% Assign a color to each vertex so that neighbors
% are all of different color
% This is the implementation of the Greedy Colouring Algorithm
% ===================================================================

disp('Assigning a color to each vertex ...');

colr   = zeros(nn,1);
maxcol = 1;
clr    = ('brgmkcykbrgmcykbrgmcykbrgmcykbrgmcykbrgmcykbrgmcykbrgmcykbrgmcyk');

% loop over total number of nodes
for i=1:nn
   
   % vertex loop centered around node i
   loop    = VLoops(Ivloop(i):Ivloop(i+1)-1);
   lenLoop = length(loop);
   okcol   = (1:maxcol);
   for k=1:lenLoop
      if (colr(loop(k))>0)
         okcol(colr(loop(k)))=maxcol+1;
      end
   end
   colr(i)=min(okcol);
   maxcol=max(maxcol,colr(i));
end

% Empty colour index array
ColIndx=[];

for i=1:maxcol
   % Identify the nodes of a particular colour
   eval(sprintf('indx=find(colr==%d);',i));
   
   % number of nodes of a particular colour
   ncol(i) = length(indx);
   
   ColIndx = [ColIndx;indx]; %#ok<AGROW>
end

% ===================================================================
% Extract edges from triangulation
% ===================================================================
disp('Finding triangles edges ...');

% Find the edges of the triangles
% neT - total number of edges
% edge array (contains the two nodes, triangle the edge belongs to and
% if the edge is repeated)
[neT,Tedges] = FindEdges(3,tri',nt);

% ===================================================================
% Invert edge list to create list of edges per triangle
% ===================================================================
disp('Reversing list of edges...');

Index1(1) = 1;
k         = 0; % runningcounter for List1

% Loop over all the identified edges
for i=1:neT
   
   Index1(i+1) = Index1(i);
   
   % The edges belong to SOME triangle
   % (not sure when the condition will fail)
   if (Tedges(3,i)>0)
      k           = k+1;
      List1(k)    = Tedges(3,i);
      Index1(i+1) = Index1(i+1)+1;
   end
   
   % if the edge is shared by two triangles
   if (Tedges(4,i)>0)
      k           = k+1;
      List1(k)    = Tedges(4,i);
      Index1(i+1) = Index1(i+1)+1;
   end
end

n1 = neT; % number of edges
n2 = nt;  % number of triangles

% Given an edge, using the starting and ending index, triangles the
% edges belong to are known. The following function, however, outputs data
% such that, the edge IDs of each triangle are now listed
[Elist2,Eindex2] = ReverseList(List1,k,n1,n2,Index1);

Elist2 = reshape(Elist2,3,nt);

% ===================================================================
% Create new vertices on each edge, and create quad edges
% ===================================================================
disp('Creating new vertices and quad edges...');

aa     = 0.5;
bb     = 0.5;
cc     = 0.5;
k1     = nn; % total number of nodes
k2     = 0;

% extend pts array to include original data pts and edges
pts    = [pts;zeros(3*neT,3)];

% arrays includes provision for ___
Qedges = zeros(4*neT+18*nt,6);

% Loop over total number of edges to generate 1/4, 1/2 and
% 3/4 points along each edge and create basic connectivity info
for i=1:neT
   i1 = Tedges(1,i); % one of the end nodes of the edge
   i2 = Tedges(2,i); % other end node of the edge
   A  = pts(i1,:);   % (x,y,z) coordinate of i1
   B  = pts(i2,:);   % (x,y,z) coordinate of i2
   
   % Is this a domain specific condition. Must check.
   %if ((Tedges(4,i)==0)&&(abs(A(2))<0.61)&&(abs(A(1))<1.1))
   %    isurface=1;
   %else
   isurface=0;
   %end
   
   %    if (isurface==1)
   %       pts(i1,2) = sign(A(2))*naca(A(1));
   %       pts(i2,2) = sign(B(2))*naca(B(1));
   %    end
   %    A = pts(i1,:);
   %    B = pts(i2,:);
   %
   % ================================================================
   %      0.25   0.25    0.25   0.25
   %    |------|------||------|------|
   %    A     A1      A2      A3    B
   % ================================================================
   
   A1             = A+aa*(B-A)/2; % quarter position from A
   A2             = (B+A)/2;      % half way between A and B
   A3             = B+aa*(A-B)/2; % quarter position from B
   
   pts(k1+1,:)    = A1; % append the three newly formed
   pts(k1+2,:)    = A2; % points along the edge A-B
   pts(k1+3,:)    = A3; % to the pts array after 'nn' rows
   
   % The following indices are a one-to-one map to the indices
   % of the edges and nodes to the pts array
   % Qedges(k2+1) [ A1 - A  ]
   % Qedges(k2+2) [ A2 - A1 ]
   % Qedges(k2+3) [ A2 - A3 ]
   % Qedges(k2+4) [ A3 - B  ]
   Qedges(k2+1,1) = k1+1;
   Qedges(k2+1,2) = i1;
   Qedges(k2+2,1) = k1+2;
   Qedges(k2+2,2) = k1+1;
   Qedges(k2+3,1) = k1+2;
   Qedges(k2+3,2) = k1+3;
   Qedges(k2+4,1) = k1+3;
   Qedges(k2+4,2) = i2;
   
   % make sure airfoil shape is accurate
   % domain specific section of the code
   if (isurface==1)
      pts(k1+1,2)=sign(A1(2))*naca(A1(1));
      pts(k1+2,2)=sign(A2(2))*naca(A2(1));
      pts(k1+3,2)=sign(A3(2))*naca(A3(1));
   end
   
   k1 = k1+3; % increment k1 + 3 (3 spaces used for A1, A2, and A3)
   k2 = k2+4; % increment k2 + 4 (4 used in Qedge for connectivity info)
end

% total number of elements in Q (hopefully lesser than 4*neT+18*nt )
neQ = k2;

% ===================================================================
% Go through triangles and create quad cells
% ===================================================================
disp('Creating quad cells ...');

pts   = [pts;zeros(7*nt,3)]; % append the pts array
Qconn = zeros(12*nt,4);      % quad-connectivity
k1    = nn+3*neT;            % reset counter pointer to after data points and edges
k2    = 0;
k3    = 4*neT;

% loop over total number of triangles
for i=1:nt
   v  = tri(i,:);  % nodes of a given triangle
   iA = v(1);      % node ID of point A
   iB = v(2);      % node ID of point B
   iC = v(3);      % node ID of point C
   
   A  = pts(iA,:); % [x,y,z] coordinate of point A
   B  = pts(iB,:); % [x,y,z] coordinate of point B
   C  = pts(iC,:); % [x,y,z] coordinate of point C
   
   e  = Elist2(:,i); % lists the edge IDs of a given triangle
   
   % Find the coordinate and indices of points on triangle edges
   for j=1:3 % loop through each of the three edges
      edg    = Tedges(:,e(j));  % list Cols 1-4 of the particular edge
      istart = nn + 3*(e(j)-1); % edge index in pts array
      is2    = 4*(e(j)-1);      % counter for Qedges ( set of 4 )
      
      % Qedges: (refer hardcopy documentation)
      %  The 3rd and 4th column refers to the index in Qconn array which
      %    contains the particular edge
      %  The 5th and 6th column refers to the edge number within the said
      %    loop - therefore runs from 1-4 (part of a quad cell)
      % If the edge has the 1st node as node A
      if (edg(1) == iA)
         
         
         % If the edge has the 2nd node as node B (i.e., A-B)
         if (edg(2) == iB)
            
            
            % Scan to the relevant edge counter in pts arrays
            % to retrive A1, A2, and A3 (D, E, and F, respectively)
            iD = istart+1; D = pts(iD,:);
            iE = istart+2; E = pts(iE,:);
            iF = istart+3; F = pts(iF,:);
            
            Qedges(is2+1,4) = k2+1;
            Qedges(is2+2,4) = k2+2;
            Qedges(is2+3,3) = k2+3;
            Qedges(is2+4,3) = k2+4;
            Qedges(is2+1,6) = 1;
            Qedges(is2+2,6) = 1;
            Qedges(is2+3,5) = 1;
            Qedges(is2+4,5) = 1;
            
         elseif (edg(2) == iC) % If the edge has the 2nd node as node C (i.e., A-C)
            
            
            % Scan to the relevant edge counter in pts arrays
            % to retrive A1, A2, and A3 (D, E, and F, respectively)
            iL = istart+1; L = pts(iL,:);
            iK = istart+2; K = pts(iK,:);
            iJ = istart+3; J = pts(iJ,:);
            
            Qedges(is2+1,3) = k2+1;
            Qedges(is2+2,3) = k2+9;
            Qedges(is2+3,4) = k2+8;
            Qedges(is2+4,4) = k2+7;
            Qedges(is2+1,5) = 4;
            Qedges(is2+2,5) = 2;
            Qedges(is2+3,6) = 2;
            Qedges(is2+4,6) = 3;
         end
      end
      
      % If the edge has the 1st node as node B
      if (edg(1)==iB)
         
         % If the edge has the 2nd node as node A (i.e., B-A)
         if (edg(2)==iA)
            
            % Scan to the relevant edge counter in pts arrays
            % to retrive A1, A2, and A3 (D, E, and F, respectively)
            iF = istart+1; F = pts(iF,:);
            iE = istart+2; E = pts(iE,:);
            iD = istart+3; D = pts(iD,:);
            
            Qedges(is2+1,3) = k2+4;
            Qedges(is2+2,3) = k2+3;
            Qedges(is2+3,4) = k2+2;
            Qedges(is2+4,4) = k2+1;
            Qedges(is2+1,5) = 1;
            Qedges(is2+2,5) = 1;
            Qedges(is2+3,6) = 1;
            Qedges(is2+4,6) = 1;
            
         elseif (edg(2) == iC) % If the edge has the 2nd node as node B (i.e., B-C)
            
            iG = istart+1; G = pts(iG,:);
            iH = istart+2; H = pts(iH,:);
            iI = istart+3; I = pts(iI,:);
            
            Qedges(is2+1,4) = k2+4;
            Qedges(is2+2,4) = k2+5;
            Qedges(is2+3,3) = k2+6;
            Qedges(is2+4,3) = k2+7;
            Qedges(is2+1,6) = 2;
            Qedges(is2+2,6) = 2;
            Qedges(is2+3,5) = 2;
            Qedges(is2+4,5) = 2;
         end
      end
      
      % If the edge has the 1st node as node C
      if (edg(1)==iC)
         
         % If the edge has the 2nd node as node A (i.e., C-A)
         if (edg(2)==iA)
            
            iJ = istart+1; J = pts(iJ,:);
            iK = istart+2; K = pts(iK,:);
            iL = istart+3; L = pts(iL,:);
            
            Qedges(is2+1,4) = k2+7;
            Qedges(is2+2,4) = k2+8;
            Qedges(is2+3,3) = k2+9;
            Qedges(is2+4,3) = k2+1;
            Qedges(is2+1,6) = 3;
            Qedges(is2+2,6) = 2;
            Qedges(is2+3,5) = 2;
            Qedges(is2+4,5) = 4;
            
         elseif (edg(2)==iB) % If the edge has the 2nd node as node B (i.e., C-B)
            
            iI = istart+1; I = pts(iI,:);
            iH = istart+2; H = pts(iH,:);
            iG = istart+3; G = pts(iG,:);
            
            Qedges(is2+1,3) = k2+7;
            Qedges(is2+2,3) = k2+6;
            Qedges(is2+3,4) = k2+5;
            Qedges(is2+4,4) = k2+4;
            Qedges(is2+1,5) = 2;
            Qedges(is2+2,5) = 2;
            Qedges(is2+3,6) = 2;
            Qedges(is2+4,6) = 2;
         end
      end
      
      
   end
   
   
   
   % Define interior points and add to vertex list
   O = (A+B+C)/3;  iO = k1+3;
   M = O+bb*(A-O); iM = k1+1;
   N = O+cc*(E-O); iN = k1+2;
   P = O+cc*(K-O); iP = k1+4;
   Q = O+bb*(B-O); iQ = k1+5;
   R = O+cc*(H-O); iR = k1+6;
   S = O+bb*(C-O); iS = k1+7;
   
   % Add to pts list after the entries of the
   % points and the edge points (A1, A2, and A3)
   pts(k1+1:k1+7,:) = [M;N;O;P;Q;R;S];
   
   % internal loops within each triangle
   % Set of 12 counter-clockwise loops
   % (Refer to fig - hardcopy)
   Qconn(k2+1,1:4)   = [iA iD iM iL];
   Qconn(k2+2,1:4)   = [iD iE iN iM];
   Qconn(k2+3,1:4)   = [iE iF iQ iN];
   Qconn(k2+4,1:4)   = [iF iB iG iQ];
   Qconn(k2+5,1:4)   = [iQ iG iH iR];
   Qconn(k2+6,1:4)   = [iR iH iI iS];
   Qconn(k2+7,1:4)   = [iS iI iC iJ];
   Qconn(k2+8,1:4)   = [iS iJ iK iP];
   Qconn(k2+9,1:4)   = [iP iK iL iM];
   Qconn(k2+10,1:4)  = [iM iN iO iP];
   Qconn(k2+11,1:4)  = [iN iQ iR iO];
   Qconn(k2+12,1:4)  = [iO iR iS iP];
   
   % Append internal edges to Qedges
   % The first 12 loops contains edges of the original triangle
   Qedges(k3+1,1:6)  = [iM iD k2+2  k2+1  4 2];
   Qedges(k3+2,1:6)  = [iN iE k2+3  k2+2  4 2];
   Qedges(k3+3,1:6)  = [iQ iF k2+4  k2+3  4 2];
   Qedges(k3+4,1:6)  = [iQ iG k2+5  k2+4  1 3];
   Qedges(k3+5,1:6)  = [iR iH k2+6  k2+5  1 3];
   Qedges(k3+6,1:6)  = [iS iI k2+7  k2+6  1 3];
   Qedges(k3+7,1:6)  = [iS iJ k2+8  k2+7  1 4];
   Qedges(k3+8,1:6)  = [iP iK k2+9  k2+8  1 3];
   Qedges(k3+9,1:6)  = [iM iL k2+1  k2+9  3 3];
   Qedges(k3+10,1:6) = [iP iM k2+10 k2+9  4 4];
   Qedges(k3+11,1:6) = [iO iN k2+11 k2+10 4 2];
   Qedges(k3+12,1:6) = [iR iQ k2+5  k2+11 4 2];
   
   % The following loops are interal loops to the triangle
   Qedges(k3+13,1:6) = [iN iQ k2+11 k2+3  1 3];
   Qedges(k3+14,1:6) = [iO iR k2+12 k2+11 1 3];
   Qedges(k3+15,1:6) = [iP iS k2+8  k2+12 4 3];
   Qedges(k3+16,1:6) = [iR iS k2+12 k2+6  2 4];
   Qedges(k3+17,1:6) = [iO iP k2+10 k2+12 3 4];
   Qedges(k3+18,1:6) = [iN iM k2+2  k2+10 3 1];
   
   % Updates for the running counters
   k1 = k1 +  7; % Qedges
   k2 = k2 + 12; % Qconn
   k3 = k3 + 18; % Qedges
   
end

neQ1 = neQ;
neQ  = neQ+18*nt;

% ===================================================================
% Find loops of quad edges
% Qloops contain the indices of the edges they form.
% Can be linked to Qedge and therefore Qconn
% ===================================================================
disp('Finding loops of quad edges ...');

% Empty initial arrays for the Q (quad) loops
Q1loops = [];
Q2loops = [];
IQloop  = ones(nn+1,1);

% Running counters for the Q loops
kk1 = 0;
kk2 = 0;

% Loop through the total number of nodes
for ii=1:nn
   
   i      = ColIndx(ii);                       % colour index of the loop
   ncl    = Icloop(i+1)-Icloop(i);             % total number of cells in the loop
   CL     = CLoops(Icloop(i):Icloop(i)+ncl-1); % triangle IDs in the loop
   QEvec1 = zeros(30,1); % should this be 35 as there can be 7 triangles
   QEvec2 = zeros(30,1);
   
   % Loop through the cell (triangles) in the loop
   % REFER DOCUMENTATION FOR THIS SECTION OF LOOP BUILDING
   for j=1:ncl
      
      
      cell = CL(j);          % particular cell ID
      v    = tri(cell,:);    % node IDs comprising of the cell
      iA   = v(1);           % node A ID
      iB   = v(2);           % node B ID
      iC   = v(3);           % node C ID
      e    = Elist2(:,cell); % list of edges for the given cell
      
      % if the node A ID is the colour index of the cell (triangle) it
      % belongs to ___
      if (iA==i)
         
         % Loop through the edges of the triangle
         for k=1:3
            
            edg = Tedges(:,e(k)); % 1-4 edge array for the particle edge
            is2 = 4*(e(k)-1);     % Index of the edge in Qedge array
            
            % Edge A-B
            if (edg(1)==iA)
               if (edg(2)==iB)
                  QEvec1(5*(j-1)+1) = is2 + 3;
                  QEvec2(5*(j-1)+1) = is2 + 4;
               end
            end
            
            % Edge B-A
            if (edg(1)==iB)
               if (edg(2)==iA)
                  QEvec1(5*(j-1)+1) = is2 + 2;
                  QEvec2(5*(j-1)+1) = is2 + 1;
               end
            end
            
            % Edge C-A
            if (edg(1)==iC)
               if (edg(2)==iA)
                  QEvec1(5*(j-1)+5) = is2 + 2;
                  QEvec2(5*(j-1)+5) = is2 + 1;
               end
            end
            
            % Edge A-C
            if (edg(1)==iA)
               if (edg(2)==iC)
                  QEvec1(5*(j-1)+5) = is2 + 3;
                  QEvec2(5*(j-1)+5) = is2 + 4;
               end
            end
         end
         
         
         QEvec1(5*(j-1)+(2:4)) = neQ1 + 18*(cell-1) + (13:15);
         QEvec2(5*(j-1)+(2:4)) = neQ1 + 18*(cell-1) + (4:6);
         
      elseif (iB==i)
         
         for k=1:3
            
            edg = Tedges(:,e(k));
            is2 = 4*(e(k)-1);
            
            if (edg(1)==iB)
               if (edg(2)==iC)
                  QEvec1(5*(j-1)+1) = is2 + 3;
                  QEvec2(5*(j-1)+1) = is2 + 4;
               end
            end
            
            if (edg(1)==iC)
               if (edg(2)==iB)
                  QEvec1(5*(j-1)+1) = is2 + 2;
                  QEvec2(5*(j-1)+1) = is2 + 1;
               end
            end
            
            if (edg(1)==iA)
               if (edg(2)==iB)
                  QEvec1(5*(j-1)+5) = is2 + 2;
                  QEvec2(5*(j-1)+5) = is2 + 1;
               end
            end
            
            if (edg(1)==iB)
               if (edg(2)==iA)
                  QEvec1(5*(j-1)+5) = is2 + 3;
                  QEvec2(5*(j-1)+5) = is2 + 4;
               end
            end
         end
         
         QEvec1(5*(j-1)+(2:4)) = neQ1 + 18*(cell-1) + (16:18);
         QEvec2(5*(j-1)+(2:4)) = neQ1 + 18*(cell-1) + (7:9);
         
      elseif (iC==i)
         for k=1:3
            
            edg = Tedges(:,e(k));
            is2 = 4*(e(k)-1);
            
            if (edg(1)==iC)
               if (edg(2)==iA)
                  QEvec1(5*(j-1)+1) = is2 + 3;
                  QEvec2(5*(j-1)+1) = is2 + 4;
               end
            end
            
            if (edg(1)==iA)
               if (edg(2)==iC)
                  QEvec1(5*(j-1)+1) = is2 + 2;
                  QEvec2(5*(j-1)+1) = is2 + 1;
               end
            end
            
            if (edg(1)==iB)
               if (edg(2)==iC)
                  QEvec1(5*(j-1)+5) = is2 + 2;
                  QEvec2(5*(j-1)+5) = is2 + 1;
               end
            end
            
            if (edg(1)==iC)
               if (edg(2)==iB)
                  QEvec1(5*(j-1)+5) = is2 + 3;
                  QEvec2(5*(j-1)+5) = is2 + 4;
               end
            end
         end
         
         QEvec1(5*(j-1)+(2:4)) = neQ1 + 18*(cell-1) + (10:12);
         QEvec2(5*(j-1)+(2:4)) = neQ1 + 18*(cell-1) + (1:3);
         
      else
         disp('error');
      end
   end
   
   % ================================================================
   % Build full quad edges loop (without rep.)
   % QEVec1 and QEVec2 are local arrays for a given node
   % Q1 and Q2 loops are accumulative loops for all nodes, which are
   % indexed based on IQloop. Single index for both loops are
   % sufficient as Q1 and Q2 are of the same length
   % ================================================================
   
   Q1loop(kk1+1) = QEvec1(1); % kk1 - running counter for edges of a loop 1
   
   % loop through number of triangles for a given node
   % Build Inner loop
   for j=1:ncl
      Q1loop(kk1+2:kk1+5) = QEvec1(5*(j-1)+(2:5));
      kk1                 = kk1 + 4;
   end
   kk1 = kk1 + 1;
   
   Q2loop(kk2+1)=QEvec2(1); % kk2 - running counter for edges of a loop 2
   
   % loop through number of triangles for a given node
   % Build Outer loop
   for j=1:ncl
      Q2loop(kk2+2:kk2+5) = QEvec2(5*(j-1)+(2:5));
      kk2                 = kk2+4;
   end
   kk2 = kk2 + 1;
   
   % start and end index of Q1 and Q2 loop for a given node
   IQloop(ii+1) = IQloop(ii) + ncl*4 + 1;
   
end

nloops1 = length(Q1loop(1,:)); % number of Q1 loops/ Q2 loops
nloops  = 2*nloops1;           % Total number of Q1 and Q2 loops
Qloops  = zeros(nloops,1);
IQloops = ones(2*nn+1,1);

k  = 0;
kk = 1;

% Run through all the nodes
% Qloops are ordered as
% [Qloop1 (node1) | Qloop2 (node1) | Qloop1 (node2) | Qloop2 (node2) ... ]
for i=1:nn
   i1  = IQloop(i);     % start index for the loop
   i2  = IQloop(i+1)-1; % end index for the loop
   nel = i2-i1+1;       % total number of edges for the loop
   
   % Ordering of Q1loop
   Qloops(k+1:k+nel) = Q1loop(i1:i2);
   k                 = k+nel;
   
   % Ordering of Q2loop
   Qloops(k+1:k+nel) = Q2loop(i1:i2);
   k                 = k+nel;
   
   % Ordering of index for inner loop
   IQloops(kk+1)     = IQloops(kk)+nel;
   kk                = kk+1;
   
   % Ordering of index for outer loop (same as inner loop)
   IQloops(kk+1)     = IQloops(kk)+nel;
   kk                = kk+1;
   
end


% ===================================================================
% Creating strand meshes
% ==================================================================
numStrandPoints = 1;
strandTemplate  = linspace(0,1,numStrandPoints);
numStrands      = length(Qloops);
strandNormal    = zeros(numStrands,3);
strandPts       = zeros(numStrands*numStrandPoints,3);

pts   = [pts;zeros(numStrands*(numStrandPoints-1),3)]; % append the pts array
k1    = nn+3*neT;            % reset counter pointer to after data points and edges
k2    = 0;
k3    = 4*neT;

% initialize counters
kk = 0;
ik = 0;
k  = 0;
jk = 0;
for i = 1:nn % loop through the nodes
   for ii = 1:2 % loop for inner and outer
      
      ik    = ik+1;
      nel   = IQloops(ik+1)-IQloops(ik);
      
      for j=1:nel
         iA = Qedges(Qloops(kk+1),1);
         iB = Qedges(Qloops(kk+1),2);
         
         % mid point of edge (WRONG FOR ANYTHING BUT A SPHERE!)
         midPoint            = (pts(iA,:)+pts(iB,:))/2;
         norm                = sqrt(sum(midPoint.*midPoint));
         strandNormal(k+1,:) = midPoint./norm; % unit normal vector
         
         
         for jj = 1:numStrandPoints
            strandPoints(jk+1,:) = midPoint;% + strandNormal(k+1,:).*strandTemplate(jj);
            jk = jk+1;
         end
         
         k  = k+1;
         kk = kk+1;
         
      end
      
   end
end
%%
% ===================================================================
% Outputs and data writing
% ===================================================================
disp('plotting quad loops...');

figure
hold on
kk = 0;
ik = 0;
for i=1:nn
   for ii=1:2 % inner and outer loop
      if(clr(colr(ColIndx(i)))=='b')
         ik    = ik+1;
         nel   = IQloops(ik+1)-IQloops(ik);
         
         kkstart = kk;
         for jj = 1:numStrandPoints
            
            Ploop = zeros(100,3);
            k     = 0;
            kk    = kkstart;
            for j=1:nel
               
               iA = Qedges(Qloops(kk+1),1);
               iB = Qedges(Qloops(kk+1),2);
               
               midPoint = 0.5*(pts(iA,:)+pts(iB,:)); % mid point of the edge
               normal   = midPoint./sqrt(sum(midPoint.*midPoint)); % unit normal vector
               
               % extruding to form "strand" loops
               Ploop(k+1,:) = midPoint;% + normal.*strandTemplate(jj); 
               
               
               k  = k+1;
               kk = kk+1;
            end
            Ploop = Ploop(1:k,:);
            h     = plot3(Ploop(:,1),Ploop(:,2),Ploop(:,3));
            set(h,'Color',clr(colr(ColIndx(i))),'LineStyle','-','Marker','.');
         end
      end
      
   end
end
axis equal
%plot3(strandPoints(:,1),strandPoints(:,2),strandPoints(:,3),'k.','MarkerSize',5)
%plot3(strandPoints(1:50,1),strandPoints(1:50,2),strandPoints(1:50,3),'k.-')


% ===================================================================
% Writing out Tecplot data
% ===================================================================
disp('Writing out Tecplot data...');
% ===================================================================
fid = fopen('quads.dat','w');
fprintf(fid,'%s\n','VARIABLES="X","Y","Z"');
eval(sprintf('fprintf(fid,''%%s\\n'',''ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'');',length(pts(:,1)),12*nt));
lenPts = length(pts(:,1));
for i=1:lenPts
   fprintf(fid,'%g %g %g\n',pts(i,:));
end

lenQconn = length(Qconn(:,1));
for i=1:lenQconn
   fprintf(fid,' %d %d %d %d \n',Qconn(i,1),Qconn(i,2),Qconn(i,3),Qconn(i,4));
end
fclose(fid);
stop
% ===================================================================
disp('Writing out general data...');
% I assume that the -1's in the output into file is to accomodate the
% 1 index of Matlab and 0 index of C/C++
% ===================================================================
fid=fopen('QuadData/coord.dat','w');
fprintf(fid,'%d\n',length(pts(:,1)));

for i=1:lenPts
   fprintf(fid,'%g %g %g\n',pts(i,:));
end
fclose(fid);

fid=fopen('QuadData/conn.dat','w');
fprintf(fid,'%d\n',length(Qconn(:,1)));
for i=1:lenQconn
   fprintf(fid,' %d %d %d %d \n',Qconn(i,1)-1,Qconn(i,2)-1,Qconn(i,3)-1,Qconn(i,4)-1);
end
fclose(fid);

fid=fopen('QuadData/qedges.dat','w');
fprintf(fid,'%d\n',neQ);
for i=1:neQ
   fprintf(fid,' %d %d %d %d %d %d \n',Qedges(i,1:6)-1);
end
fclose(fid);

fid=fopen('QuadData/qloops.dat','w');
fprintf(fid,'%d\n',nloops);
for i=1:nloops
   fprintf(fid,' %d\n',Qloops(i)-1);
end
fclose(fid);
fid=fopen('QuadData/iqloops.dat','w');
fprintf(fid,'%d\n',2*nn+1);
for i=1:2*nn+1
   fprintf(fid,' %d\n',IQloops(i)-1);
end
fclose(fid);

fid=fopen('QuadData/ncolors.dat','w');
fprintf(fid,'%d\n',maxcol);
for i=1:maxcol
   fprintf(fid,' %d\n',ncol(i));
end
fclose(fid);

k=0;
for i=1:maxcol
   disp(sprintf('Preparing Tecplot File for loops %d ...',i)); %#ok<DSPS>
   eval(sprintf('fid=fopen(''loops%d.dat'',''w'');',i));
   P=[];
   conn=[];
   kk=0;
   for j=1:2*ncol(i) % go through all loops of color i
      k=k+1;
      i1=IQloops(k);
      i2=IQloops(k+1)-1;
      for ii=i1:i2
         ie=Qloops(ii,1);
         iA=Qedges(ie,1);
         iB=Qedges(ie,2);
         kk=kk+1;
         P=[P;(pts(iA,:)+pts(iB,:))*0.5]; %#ok<AGROW>
         if (ii<i2)
            conn=[conn;kk kk+1]; %#ok<AGROW>
         end
      end
   end
   
   lenP    = length(P(:,1));
   lenConn = length(conn(:,1));
   
   fprintf(fid,'%s\n','VARIABLES="X","Y","Z"');
   eval(sprintf('fprintf(fid,''%%s\\n'',''ZONE DATAPACKING=POINT, NODES=%d, ELEMENTS=%d, ZONETYPE=FELINESEG'');',...
      lenP,lenConn));
   
   for j=1:lenP
      fprintf(fid,'%g %g %g\n',P(j,:));
   end
   for j=1:lenConn
      fprintf(fid,'%d %d \n',conn(j,1),conn(j,2));
   end
   fclose(fid);
   
end
% ###################################################################
% END OF FILE
% ###################################################################
