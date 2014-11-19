%% ############################
% PrepareLoops.m
% Written by Dr. Jayanarayanan Sitaraman
% ##############################
close all
clear all

%-----------------------------------------
% Read in  coord and connectivity files
%------------------------------------------
disp('Reading coord and connectivity files...');
iSphere = 0;
%pts = load('../SphereTessellation/coordSphere.dat');
pts = load('../robin/coordRobin.dat');
nn=length(pts(:,1));

x=pts(:,1);
y=pts(:,2);
z=pts(:,3);

% node index list of each triangles
%tri = load('../SphereTessellation/connSphere.dat');
tri = load('../robin/connRobin.dat');

nt=length(tri(:,1)); % number of triangles


[nn nt]
dsd

figure
hold on
h=trimesh(tri,x,y,z);
set(h,'EdgeColor','k');
axis equal


lenPts = length(pts(:,1));
lenTri = length(tri);
fid = fopen('originalMesh.dat','w');
fprintf(fid,'%s\n','VARIABLES="X","Y","Z"');
eval(sprintf('fprintf(fid,''%%s\\n'',''ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE'');',lenPts,lenTri));

for i=1:lenPts
    fprintf(fid,'%g %g %g\n',pts(i,:));
end

for i=1:lenTri
    fprintf(fid,'%d %d %d\n',tri(i,:));
end
fclose(fid);


%-----------------------------------------
% Invert connectivity list
%-----------------------------------------
disp('Inverting connectivity list...');
% reshaped colume -> row-wise
List1=reshape(tri(1:nt,1:3)',3*nt,1);

n1=nt; % number of traiangles
n2=nn; % number of nodes
Index1=[1:3:3*n1+1];

% List2 is a list of cell index that consists specific node.
[List2,Index2]=ReverseList(List1,3*nt,n1,n2,Index1);

%-----------------------------------------
% Go through vertices and find loops
% (cell loops and vertex loops)
%-----------------------------------------
ikc=0;ikv=0;
Icloop=ones(nn,1);
Ivloop=ones(nn,1);
disp('Finding loops (of cells and vertices)...');

% Loop over total number of nodes
for i=1:nn
   
   k=99;
   k1=k;
   k2=k+3;
   i1=Index2(i);
   i2=Index2(i+1)-1;
   nc=i2-i1+1;
   
   % CC: the node list of the cell that surrounding node.
   CC=tri(List2(i1:i2),:);
   
   CC2=zeros(nc,2);
   
   % Store Other two nodes index expept the center node.
   for j=1:nc
      if (CC(j,1)==i)
         CC2(j,1)=CC(j,2);
         CC2(j,2)=CC(j,3);
         continue
      end
      if (CC(j,2)==i)
         CC2(j,1)=CC(j,3);
         CC2(j,2)=CC(j,1);
         continue
      end
      CC2(j,1)=CC(j,1);
      CC2(j,2)=CC(j,2);
   end
   
   % find vortex loop (VLoops): find vortex index that surroundings specific node starting from node1
   % find cell loop (CLoops): find cell index that surroundings specific node starting from node 1
   if (nc==1) % if the number of neighboring cell is only 1
      ncl=1;
      nv=2;
      CLoops(ikc+1:ikc+ncl)=List2(i1);
      VLoops(ikv+1:ikv+nv)=CC2(1,1:2);
   else % more than 1 neighboring cells
      cloop1=[];
      cloop2(1)=List2(i1);
      kloop1=0;
      kloop2=1;
      V1=CC2(1,1);
      V2=CC2(1,2);
      EV(k+1)=V1;  % end vertices for the loop
      EV(k+2)=V2;
      iact=ones(1,nc);
      iact(1)=0;
      nact=nc-1;
      icont=1;
      
      while (icont==1)
         for j=2:nc
            if (iact(j)==1)
               for kk=1:2
                  if (CC2(j,kk)==V1)
                     V1=CC2(j,3-kk);
                     EV(k1)=V1; k1=k1-1;
                     iact(j)=0;
                     nact=nact-1;
                     cloop1(kloop1+1)=List2(i1+j-1);
                     kloop1=kloop1+1;
                     break;
                  elseif (CC2(j,kk)==V2)
                     V2=CC2(j,3-kk);
                     EV(k2)=V2; k2=k2+1;
                     iact(j)=0;
                     nact=nact-1;
                     cloop2(kloop2+1)=List2(i1+j-1);
                     kloop2=kloop2+1;
                     break;
                  end %if
               end %for
            end %if
            
            % Termination condition (remaining active edge)
            if (nact==0)
               icont=0;
               break;
            end
            
         end % for loop (j)
      end % while condition (icont)
      
      % Cell loops
      ncl=kloop1+kloop2;
      CLoops(ikc+1:ikc+ncl)=[cloop1(kloop1:-1:1) cloop2(1:kloop2)];
      
      % Vertex loops
      nv=k2-k1-1;
      VLoops(ikv+1:ikv+nv)=EV(k1+1:k2-1);
   end
   
   
   % ikc is total number of surrounding cells
   % ikv is total number of surrouding nodes
   % Icloop and Ivloop may be for check process
   % Icloop: add the number of surrounding cells at node number location
   % Ivloop: add the number of surrounding nodes at node number location
   ikc=ikc+ncl;
   ikv=ikv+nv;
   Icloop(i+1)=Icloop(i)+ncl;
   Ivloop(i+1)=Ivloop(i)+nv;
   
end

%------------------------------------------------
% Assign a color to each vertex
%------------------------------------------------
% (so that neighbors are all of different color)

disp('Assigning a color to each vertex ...');
colr=zeros(nn,1);
maxcol=1;
clr=['brgmkcykbrgmcykbrgmcykbrgmcykbrgmcykbrgmcykbrgmcykbrgmcykbrgmcyk'];

% colr : each node has color that is defined as numbers (1-6)
% colr of each node is defined by selecting minimum okcol
for i=1:nn
   
   % vertex loop centered around node i
   loop=VLoops(Ivloop(i):Ivloop(i+1)-1);
   okcol=[1:maxcol];
   
   
   for k=1:length(loop)
      if (colr(loop(k))>0)
         okcol(colr(loop(k)))=maxcol+1;
      end
   end
   colr(i)=min(okcol);
   maxcol=max(maxcol,colr(i));
   
end

% ncol : the total number of nodes at each color
% ColIndx : arrange the node number according to the same color number
ColIndx=[];
for i=1:maxcol
   eval(sprintf('indx=find(colr==%d);',i));
   ncol(i)=length(indx);
   ColIndx=[ColIndx;indx];
end


%--------------------------------------
% Extract edges from triangulation
%--------------------------------------
% find node number that consists triangle cell
% and find neighbor cell index for each face.
% neT : Total number of edges
disp('Finding triangles edges ...');
[neT,Tedges]=FindEdges(3,tri',nt);


%-------------------------------------------------------
% Invert edge list to create list of edges per triangle
%-------------------------------------------------------
disp('Reversing list of edges...');
Index1(1)=1;
k=0;
for i=1:neT
   Index1(i+1)=Index1(i);
   
   % The edges belong to SOME triangle
   if (Tedges(3,i)>0)
      k=k+1;
      List1(k)=Tedges(3,i);
      Index1(i+1)=Index1(i+1)+1;
   end
   
   % if the edge is shared by two triangles
   if (Tedges(4,i)>0)
      k=k+1;
      List1(k)=Tedges(4,i);
      Index1(i+1)=Index1(i+1)+1;
   end
end

n1=neT; %number of edges
n2=nt; % number of triagnels

%Edge indexes that consists triangles
[Elist2,Eindex2]=ReverseList(List1,k,n1,n2,Index1);
Elist2=reshape(Elist2,3,nt);

%--------------------------------------------------------
% Create new vertices on each edge, and create quad edges
%--------------------------------------------------------
disp('Creating new vertices and quad edges...');
aa=0.5;
bb=0.5;
cc=0.5;
k1=nn;
k2=0;
pts=[pts;zeros(3*neT,3)];

Qedges=zeros(4*neT+18*nt,6);

% neT is the number of faces
% nt is the number of cells
% A, B is coordinate of nodes
for i=1:neT
   i1=Tedges(1,i);
   i2=Tedges(2,i);
   A=pts(i1,:);
   B=pts(i2,:);
   
   % ys code
   % to define surface face?(isurface)
   %   if ((Tedges(4,i)==0) &&(abs(A(2))<0.61)&&(abs(A(1))<1.1))
   %     isurface=1;
   %   else
   %     isurface=0;
   %   end
   
   % ys code
   % divide each face by 3 element, A1 A2 A3
   % and the points coordinate is saved after nnode.
   % Qedge is connectivity of nodes that consists same face.
   A1=A+aa*(B-A)/2;
   A2=(B+A)/2;
   A3=B+aa*(A-B)/2;
   
   if(iSphere==1)
      A1 = A1./sqrt(sum(A1.*A1));
      A2 = A2./sqrt(sum(A2.*A2));
      A3 = A3./sqrt(sum(A3.*A3));
   end
   
   % ================================================================
   %      0.25   0.25    0.25   0.25
   %    |------|------||------|------|
   %    A     A1      A2      A3    B
   % ================================================================
   
   
   pts(k1+1,:)=A1;
   pts(k1+2,:)=A2;
   pts(k1+3,:)=A3;
   
   
   % Qedges(k2+1) [ A1 - A  ]
   % Qedges(k2+2) [ A2 - A1 ]
   % Qedges(k2+3) [ A2 - A3 ]
   % Qedges(k2+4) [ A3 - B  ]
   Qedges(k2+1,1)=k1+1;
   Qedges(k2+1,2)=i1;
   Qedges(k2+2,1)=k1+2;
   Qedges(k2+2,2)=k1+1;
   Qedges(k2+3,1)=k1+2;
   Qedges(k2+3,2)=k1+3;
   Qedges(k2+4,1)=k1+3;
   Qedges(k2+4,2)=i2;
   
   % ys code
   % k1 is total number of node
   % k2 is total number of face
   k1=k1+3;
   k2=k2+4;
   
end

%total number of elements in Q
neQ=k2;

% Go through triangles and create quad cells
%-------------------------------------------------------
disp('Creating quad cells ...');
pts=[pts;zeros(7*nt,3)];
Qconn=zeros(12*nt,4);
k1=nn+3*neT;
k2=0;
k3=4*neT;

for i=1:nt
   v=tri(i,:);
   iA=v(1);iB=v(2);iC=v(3);
   A=pts(iA,:);B=pts(iB,:);C=pts(iC,:);
   e=Elist2(:,i);
   
   % Find coord and indices of points on triangle edges
   for j=1:3
      edg=Tedges(:,e(j));
      istart=nn+3*(e(j)-1);
      is2=4*(e(j)-1);
      
      % for the case of edg(1) is same as 1st node
      if (edg(1)==iA)
         if (edg(2)==iB)
            
            iD=istart+1;D=pts(iD,:);
            iE=istart+2;E=pts(iE,:);
            iF=istart+3;F=pts(iF,:);
            
            Qedges(is2+1,4)=k2+1;
            Qedges(is2+2,4)=k2+2;
            Qedges(is2+3,3)=k2+3;
            Qedges(is2+4,3)=k2+4;
            Qedges(is2+1,6)=1;
            Qedges(is2+2,6)=1;
            Qedges(is2+3,5)=1;
            Qedges(is2+4,5)=1;
            
         elseif (edg(2)==iC)
            iL=istart+1;L=pts(iL,:);
            iK=istart+2;K=pts(iK,:);
            iJ=istart+3;J=pts(iJ,:);
            Qedges(is2+1,3)=k2+1;
            Qedges(is2+2,3)=k2+9;
            Qedges(is2+3,4)=k2+8;
            Qedges(is2+4,4)=k2+7;
            Qedges(is2+1,5)=4;
            Qedges(is2+2,5)=2;
            Qedges(is2+3,6)=2;
            Qedges(is2+4,6)=3;
         end
      end
      
      % for the case of edg(1) is same as 2nd node
      if (edg(1)==iB)
         if (edg(2)==iA)
            iF=istart+1;F=pts(iF,:);
            iE=istart+2;E=pts(iE,:);
            iD=istart+3;D=pts(iD,:);
            Qedges(is2+1,3)=k2+4;
            Qedges(is2+2,3)=k2+3;
            Qedges(is2+3,4)=k2+2;
            Qedges(is2+4,4)=k2+1;
            Qedges(is2+1,5)=1;
            Qedges(is2+2,5)=1;
            Qedges(is2+3,6)=1;
            Qedges(is2+4,6)=1;
         elseif (edg(2)==iC)
            iG=istart+1;G=pts(iG,:);
            iH=istart+2;H=pts(iH,:);
            iI=istart+3;I=pts(iI,:);
            Qedges(is2+1,4)=k2+4;
            Qedges(is2+2,4)=k2+5;
            Qedges(is2+3,3)=k2+6;
            Qedges(is2+4,3)=k2+7;
            Qedges(is2+1,6)=2;
            Qedges(is2+2,6)=2;
            Qedges(is2+3,5)=2;
            Qedges(is2+4,5)=2;
         end
      end
      
      % for the case of edg(1) is same as 3rd node
      if (edg(1)==iC)
         if (edg(2)==iA)
            iJ=istart+1;J=pts(iJ,:);
            iK=istart+2;K=pts(iK,:);
            iL=istart+3;L=pts(iL,:);
            Qedges(is2+1,4)=k2+7;
            Qedges(is2+2,4)=k2+8;
            Qedges(is2+3,3)=k2+9;
            Qedges(is2+4,3)=k2+1;
            Qedges(is2+1,6)=3;
            Qedges(is2+2,6)=2;
            Qedges(is2+3,5)=2;
            Qedges(is2+4,5)=4;
         elseif (edg(2)==iB)
            iI=istart+1;I=pts(iI,:);
            iH=istart+2;H=pts(iH,:);
            iG=istart+3;G=pts(iG,:);
            Qedges(is2+1,3)=k2+7;
            Qedges(is2+2,3)=k2+6;
            Qedges(is2+3,4)=k2+5;
            Qedges(is2+4,4)=k2+4;
            Qedges(is2+1,5)=2;
            Qedges(is2+2,5)=2;
            Qedges(is2+3,6)=2;
            Qedges(is2+4,6)=2;
         end
      end
      
   end
   
   % Define interior points and add to vertex list
   
   % D,E,F,G,H,I,J,K,L is the node on the face.
   % M,N,O,P,Q,R,S is the node in the face.
   % Qconn : node to cell connectivity of quad cells.
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   % Qedges : node to (left cell and face index) or right cell and face index)
   % node1, node2, R_Cell_index, L_Cell_index, face_index_R_Cell, face_index_L_Cell
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   
   O=(A+B+C)/3;iO=k1+3;
   
   
   M=O+bb*(A-O);iM=k1+1;
   N=O+cc*(E-O);iN=k1+2;
   P=O+cc*(K-O);iP=k1+4;
   Q=O+bb*(B-O);iQ=k1+5;
   R=O+cc*(H-O);iR=k1+6;
   S=O+bb*(C-O);iS=k1+7;
   
   % Push onto surface of sphere
   if(iSphere==1)
      M = M./sqrt(sum(M.*M));
      N = N./sqrt(sum(N.*N));
      O = O./sqrt(sum(O.*O));
      P = P./sqrt(sum(P.*P));
      Q = Q./sqrt(sum(Q.*Q));
      R = R./sqrt(sum(R.*R));
      S = S./sqrt(sum(S.*S));
      
   end
   
   pts(k1+1:k1+7,:)=[M;N;O;P;Q;R;S];
   Qconn(k2+1,1:4)=[iA iD iM iL];
   Qconn(k2+2,1:4)=[iD iE iN iM];
   Qconn(k2+3,1:4)=[iE iF iQ iN];
   Qconn(k2+4,1:4)=[iF iB iG iQ];
   Qconn(k2+5,1:4)=[iQ iG iH iR];
   Qconn(k2+6,1:4)=[iR iH iI iS];
   Qconn(k2+7,1:4)=[iS iI iC iJ];
   Qconn(k2+8,1:4)=[iS iJ iK iP];
   Qconn(k2+9,1:4)=[iP iK iL iM];
   Qconn(k2+10,1:4)=[iM iN iO iP];
   Qconn(k2+11,1:4)=[iN iQ iR iO];
   Qconn(k2+12,1:4)=[iO iR iS iP];
   Qedges(k3+1,1:6)=[iM iD k2+2 k2+1 4 2];
   Qedges(k3+2,1:6)=[iN iE k2+3 k2+2 4 2];
   Qedges(k3+3,1:6)=[iQ iF k2+4 k2+3 4 2];
   Qedges(k3+4,1:6)=[iQ iG k2+5 k2+4 1 3];
   Qedges(k3+5,1:6)=[iR iH k2+6 k2+5 1 3];
   Qedges(k3+6,1:6)=[iS iI k2+7 k2+6 1 3];
   Qedges(k3+7,1:6)=[iS iJ k2+8 k2+7 1 4];
   Qedges(k3+8,1:6)=[iP iK k2+9 k2+8 1 3];
   Qedges(k3+9,1:6)=[iM iL k2+1 k2+9 3 3];
   Qedges(k3+10,1:6)=[iP iM k2+10 k2+9 4 4];
   Qedges(k3+11,1:6)=[iO iN k2+11 k2+10 4 2];
   Qedges(k3+12,1:6)=[iR iQ k2+5 k2+11 4 2];
   Qedges(k3+13,1:6)=[iN iQ k2+11 k2+3 1 3];
   Qedges(k3+14,1:6)=[iO iR k2+12 k2+11 1 3];
   Qedges(k3+15,1:6)=[iP iS k2+8 k2+12 4 3];
   Qedges(k3+16,1:6)=[iR iS k2+12 k2+6 2 4];
   Qedges(k3+17,1:6)=[iO iP k2+10 k2+12 3 4];
   Qedges(k3+18,1:6)=[iN iM k2+2 k2+10 3 1];
   k1=k1+7;
   k2=k2+12;
   k3=k3+18;
   
end
neQ1=neQ;
neQ=neQ+18*nt;



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%find normal vectors at each face for STRAND GRID
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
nlayer = 10;
init_length = 0.04;
stretch = 1.0;

nvec_cell = zeros(nt,3);
for i=1:nt
   
   nn1=tri(i,1);
   nn2=tri(i,2);
   nn3=tri(i,3);
   
   xx1=pts(nn1,1);
   yy1=pts(nn1,2);
   zz1=pts(nn1,3);
   
   xx2=pts(nn2,1);
   yy2=pts(nn2,2);
   zz2=pts(nn2,3);
   
   xx3=pts(nn3,1);
   yy3=pts(nn3,2);
   zz3=pts(nn3,3);
   
   dx2 = xx2 - xx1;
   dx3 = xx3 - xx1;
   
   dy2 = yy2 - yy1;
   dy3 = yy3 - yy1;
   
   dz2 = zz2 - zz1;
   dz3 = zz3 - zz1;
   
   nvec_cell(i,1) = dy2*dz3-dy3*dz2;
   nvec_cell(i,2) = dz2*dx3-dx2*dz3;
   nvec_cell(i,3) = dx2*dy3-dy2*dx3;
   
   nvec_cell(i,:) = nvec_cell(i,:)/(nvec_cell(i,1)^2 + nvec_cell(i,2)^2 + nvec_cell(i,3)^2)^0.5;
   
end

nvec_node= zeros(nn,3);
ave_nvec = zeros(nn,3);

%calculate normal vector at each node of tri elements
for i=1:nn
   ncl = Icloop(i+1)-Icloop(i);
   CL=CLoops(Icloop(i):Icloop(i)+ncl-1);
   
   for j=1:ncl
      cell = CL(j);
      ave_nvec(i,:) = ave_nvec(i,:) + nvec_cell(cell,:);
   end
   ave_nvec(i,:) = ave_nvec(i,:)/ncl;
   nvec_node(i,:) = ave_nvec(i,:)/(ave_nvec(i,1)^2+ave_nvec(i,2)^2+ave_nvec(i,3)^2)^0.5;
end


%Normal vector smoothing and calculating node points of tri elements for
%layers
for k=1:20 % iteration of avaeraing. (user's constant)
   for i=1:nn
      %  find neighbor vectors
      loop=VLoops(Ivloop(i):Ivloop(i+1)-2); %last node index is duplicate.(all node is interior nodes)
      nvert = length(loop);
      
      %  averaging neighbor vector
      iflag = 0;
      nvec = (nvec_node(i,1)^2 + nvec_node(i,2)^2 + nvec_node(i,3)^2)^0.5;
      for j=1:nvert
         crt_vec(j,:) = (nvec_node(i,:)+nvec_node(j,:))/2;
         crt = (crt_vec(j,1)^2+crt_vec(j,2)^2+crt_vec(j,3)^2)^0.5;
         
         if (abs(crt - nvec)<=0.3) iflag = 1;
         end
      end
      
      if(iflag==0)
         
         for j=1:nvert
            vert = loop(j);
            nvec_node(i,:) = nvec_node(i,:) + nvec_node(vert,:);
         end
         
         nvec_node(i,:) = nvec_node(i,:)/(nvert+1);
         nvec_node(i,:) = nvec_node(i,:)/(nvec_node(i,1)^2+nvec_node(i,2)^2+nvec_node(i,3)^2)^0.5;
         
      end
      
   end
end

% pts2 : node points are stored at "pts2" except surface nodes.
% First 2nd layer tri nodes -> 2nd layer quads nodes
% 3rd layer tri nodes -> 3rd layer quads nodes......
% tri nodes of each layer is calculated using normal vectors.
leng=0;
leng_strand = init_length;
for i=1:nlayer
   for j=1:nn
      pts2(j+leng,:)= pts(j,:) + i*nvec_node(j,:)*leng_strand;
   end
   pts2 = [pts2;zeros(3*neT+7*nt,3)];
   leng = length(pts2);
   leng_strand = leng_strand*stretch;
end

x=pts2(:,1);
y=pts2(:,2);
z=pts2(:,3);


%--------------------------------------------------------
% REPEAT MAKING QUAD CELLS IN TRI CELLS
%--------------------------------------------------------
disp('Creating new vertices and quad edges...');
aa=0.5;
bb=0.5;
cc=0.5;
k1=nn;
k2=0;

for j=1:nlayer
   
   k1=nn;
   k2=0;
   for i=1:neT
      i1=Tedges(1,i);
      i2=Tedges(2,i);
      A=pts2(i1+(j-1)*(nn+3*neT+7*nt),:);
      B=pts2(i2+(j-1)*(nn+3*neT+7*nt),:);
      
      % divide each face by 3 element, A1 A2 A3
      % and the points coordinate is saved after nnode.
      % Qedge is connectivity of nodes that consists same face.
      A1=A+aa*(B-A)/2;
      A2=(B+A)/2;
      A3=B+aa*(A-B)/2;
      
      % ================================================================
      %      0.25   0.25    0.25   0.25
      %    |------|------||------|------|
      %    A     A1      A2      A3    B
      % ================================================================
      
      pts2(k1+1+(j-1)*(nn+3*neT+7*nt),:)=A1;
      pts2(k1+2+(j-1)*(nn+3*neT+7*nt),:)=A2;
      pts2(k1+3+(j-1)*(nn+3*neT+7*nt),:)=A3;
      
      k1=k1+3;
      k2=k2+4;
      
   end
end

%total number of elements in Q
neQ=k2;

%-------------------------------------------------------
% Go through triangles and create quad cells
%-------------------------------------------------------
disp('Creating quad cells ...');
%pts2=[pts2;zeros(7*nt,3)];
k1=nn+3*neT;
k2=0;
k3=4*neT;

for k=1:nlayer
   
   k1=nn+3*neT;
   k2=0;
   k3=4*neT;
   for i=1:nt
      v=tri(i,:);
      iA=v(1);iB=v(2);iC=v(3);
      A=pts2(iA+(k-1)*(nn+3*neT+7*nt),:);B=pts2(iB+(k-1)*(nn+3*neT+7*nt),:);C=pts2(iC+(k-1)*(nn+3*neT+7*nt),:);
      e=Elist2(:,i);
      
      % Find coord and indices of points on triangle edges
      for j=1:3
         edg=Tedges(:,e(j));
         istart=nn+3*(e(j)-1);
         is2=4*(e(j)-1);
         
         % ys code
         % Qedges :
         
         % for the case of edg(1) is same as 1st node
         if (edg(1)==iA)
            if (edg(2)==iB)
               
               iD=istart+1;D=pts2(iD+(k-1)*(nn+3*neT+7*nt),:);
               iE=istart+2;E=pts2(iE+(k-1)*(nn+3*neT+7*nt),:);
               iF=istart+3;F=pts2(iF+(k-1)*(nn+3*neT+7*nt),:);
               
               
               
            elseif (edg(2)==iC)
               iL=istart+1;L=pts2(iL+(k-1)*(nn+3*neT+7*nt),:);
               iK=istart+2;K=pts2(iK+(k-1)*(nn+3*neT+7*nt),:);
               iJ=istart+3;J=pts2(iJ+(k-1)*(nn+3*neT+7*nt),:);
               
            end
         end
         
         % for the case of edg(1) is same as 2nd node
         if (edg(1)==iB)
            if (edg(2)==iA)
               iF=istart+1;F=pts2(iF+(k-1)*(nn+3*neT+7*nt),:);
               iE=istart+2;E=pts2(iE+(k-1)*(nn+3*neT+7*nt),:);
               iD=istart+3;D=pts2(iD+(k-1)*(nn+3*neT+7*nt),:);
               
            elseif (edg(2)==iC)
               iG=istart+1;G=pts2(iG+(k-1)*(nn+3*neT+7*nt),:);
               iH=istart+2;H=pts2(iH+(k-1)*(nn+3*neT+7*nt),:);
               iI=istart+3;I=pts2(iI+(k-1)*(nn+3*neT+7*nt),:);
            end
         end
         
         % for the case of edg(1) is same as 3rd node
         if (edg(1)==iC)
            if (edg(2)==iA)
               iJ=istart+1;J=pts2(iJ+(k-1)*(nn+3*neT+7*nt),:);
               iK=istart+2;K=pts2(iK+(k-1)*(nn+3*neT+7*nt),:);
               iL=istart+3;L=pts2(iL+(k-1)*(nn+3*neT+7*nt),:);
            elseif (edg(2)==iB)
               iI=istart+1;I=pts2(iI+(k-1)*(nn+3*neT+7*nt),:);
               iH=istart+2;H=pts2(iH+(k-1)*(nn+3*neT+7*nt),:);
               iG=istart+3;G=pts2(iG+(k-1)*(nn+3*neT+7*nt),:);
            end
         end
         
      end
      
      % Define interior points and add to vertex list
      
      % ys code
      % D,E,F,G,H,I,J,K,L is the node on the face.
      % M,N,O,P,Q,R,S is the node in the face.
      % Qconn : node to cell connectivity of quad cells.
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      % Qedges : node to (left cell and face index) or right cell and face index)
      % node1, node2, R_Cell_index, L_Cell_index, face_index_R_Cell, face_index_L_Cell
      %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
      
      O=(A+B+C)/3;iO=k1+3;
      
      M=O+bb*(A-O);iM=k1+1;
      N=O+cc*(E-O);iN=k1+2;
      P=O+cc*(K-O);iP=k1+4;
      Q=O+bb*(B-O);iQ=k1+5;
      R=O+cc*(H-O);iR=k1+6;
      S=O+bb*(C-O);iS=k1+7;
      
      pts2(k1+1+(k-1)*(nn+3*neT+7*nt):k1+7+(k-1)*(nn+3*neT+7*nt),:)=[M;N;O;P;Q;R;S];
      
      k1=k1+7;
      k2=k2+12;
      k3=k3+18;
      
   end
   
end

% REPEAT END
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%-------------------------------------------
% Find loops of quad edges
%-------------------------------------------
disp('Finding loops of quad edges ...');
Q1loops=[];
Q2loops=[];
IQloop=ones(nn+1,1);
kk1=0;
kk2=0;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% only for the initial surface mesh layer because same connectivity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ncl : number of surrounding cells at certain node
for ii=1:nn
   i=ColIndx(ii);
   ncl=Icloop(i+1)-Icloop(i);
   CL=CLoops(Icloop(i):Icloop(i)+ncl-1);
   QEvec1=zeros(30,1);
   QEvec2=zeros(30,1);
   
   % ys code
   % Elist2 : cell to face connectivity from cell number 1
   for j=1:ncl
      cell=CL(j);
      v=tri(cell,:);iA=v(1);iB=v(2);iC=v(3);
      e=Elist2(:,cell);
      
      
      if (iA==i) % ys code
         for k=1:3
            edg=Tedges(:,e(k));
            is2=4*(e(k)-1);
            
            if (edg(1)==iA)
               if (edg(2)==iB)
                  QEvec1(5*(j-1)+1)=is2+3;
                  QEvec2(5*(j-1)+1)=is2+4;
               end
            end
            
            if (edg(1)==iB)
               if (edg(2)==iA)
                  QEvec1(5*(j-1)+1)=is2+2;
                  QEvec2(5*(j-1)+1)=is2+1;
               end
            end
            
            if (edg(1)==iC)
               if (edg(2)==iA)
                  QEvec1(5*(j-1)+5)=is2+2;
                  QEvec2(5*(j-1)+5)=is2+1;
               end
            end
            
            if (edg(1)==iA)
               if (edg(2)==iC)
                  QEvec1(5*(j-1)+5)=is2+3;
                  QEvec2(5*(j-1)+5)=is2+4;
               end
            end
            
         end    % ys code for k=1:3
         
         QEvec1(5*(j-1)+[2:4])=neQ1+18*(cell-1)+[13:15];
         QEvec2(5*(j-1)+[2:4])=neQ1+18*(cell-1)+[4:6];
         
         
      elseif (iB==i) %ys code
         
         
         for k=1:3
            edg=Tedges(:,e(k));
            is2=4*(e(k)-1);
            
            if (edg(1)==iB)
               if (edg(2)==iC)
                  QEvec1(5*(j-1)+1)=is2+3;
                  QEvec2(5*(j-1)+1)=is2+4;
               end
            end
            
            if (edg(1)==iC)
               if (edg(2)==iB)
                  QEvec1(5*(j-1)+1)=is2+2;
                  QEvec2(5*(j-1)+1)=is2+1;
               end
            end
            
            if (edg(1)==iA)
               if (edg(2)==iB)
                  QEvec1(5*(j-1)+5)=is2+2;
                  QEvec2(5*(j-1)+5)=is2+1;
               end
            end
            
            if (edg(1)==iB)
               if (edg(2)==iA)
                  QEvec1(5*(j-1)+5)=is2+3;
                  QEvec2(5*(j-1)+5)=is2+4;
               end
            end
            
            
         end
         
         
         QEvec1(5*(j-1)+[2:4])=neQ1+18*(cell-1)+[16:18];
         QEvec2(5*(j-1)+[2:4])=neQ1+18*(cell-1)+[7:9];
         
      elseif (iC==i) % ys code
         
         for k=1:3
            
            edg=Tedges(:,e(k));
            is2=4*(e(k)-1);
            
            if (edg(1)==iC)
               if (edg(2)==iA)
                  QEvec1(5*(j-1)+1)=is2+3;
                  QEvec2(5*(j-1)+1)=is2+4;
               end
            end
            
            if (edg(1)==iA)
               if (edg(2)==iC)
                  QEvec1(5*(j-1)+1)=is2+2;
                  QEvec2(5*(j-1)+1)=is2+1;
               end
            end
            
            if (edg(1)==iB)
               if (edg(2)==iC)
                  QEvec1(5*(j-1)+5)=is2+2;
                  QEvec2(5*(j-1)+5)=is2+1;
               end
            end
            
            if (edg(1)==iC)
               if (edg(2)==iB)
                  QEvec1(5*(j-1)+5)=is2+3;
                  QEvec2(5*(j-1)+5)=is2+4;
               end
            end
         end % ys code for k=1:3
         
         QEvec1(5*(j-1)+[2:4])=neQ1+18*(cell-1)+[10:12];
         QEvec2(5*(j-1)+[2:4])=neQ1+18*(cell-1)+[1:3];
         
      else % ys code
         disp('error');
      end
      
   end % ys code for j=1:ncl
   
   
   
   % Build full quad edges loop (without rep.)
   
   % ys code
   % Q1loop is same as QEvec1 except reputation numbers.
   Q1loop(kk1+1)=QEvec1(1);
   for j=1:ncl
      Q1loop(kk1+2:kk1+5)=QEvec1(5*(j-1)+[2:5]);
      kk1=kk1+4;
   end
   kk1=kk1+1;
   
   
   Q2loop(kk2+1)=QEvec2(1);
   for j=1:ncl
      Q2loop(kk2+2:kk2+5)=QEvec2(5*(j-1)+[2:5]);
      kk2=kk2+4;
   end
   kk2=kk2+1;
   
   
   % ys code: IQloop is number of vector for each nodes case.
   IQloop(ii+1)=IQloop(ii)+ncl*4+1;
end

nloops1=length(Q1loop(1,:));
nloops=2*nloops1;
Qloops=zeros(nloops,1);
IQloops=ones(2*nn+1,1);
k=0;
kk=1;


% Qloops = (Q1loop + Q2loop) for surrouding loop of node #1 +
% (Q1loop + Q2loop) for surrounding loop of node #2 + ....
for i=1:nn
   i1=IQloop(i);
   i2=IQloop(i+1)-1;
   nel=i2-i1+1;
   
   % ordering of Q1loop
   Qloops(k+1:k+nel)=Q1loop(i1:i2);
   k=k+nel;
   
   % ordering of Q2loop
   Qloops(k+1:k+nel)=Q2loop(i1:i2);
   k=k+nel;
   
   % ordering of index for inner loop
   IQloops(kk+1)=IQloops(kk)+nel;
   kk=kk+1;
   
   % ordering of index for outer loop (same as inner loop)
   IQloops(kk+1)=IQloops(kk)+nel;
   kk=kk+1;
end

disp('plotting quad loops...');
figure
hold on
kk=0;
ik=0;

for i=1:nn
   for ii=1:2
      
      if(clr(colr(ColIndx(i)))=='b')
         ik=ik+1;
         nel=IQloops(ik+1)-IQloops(ik);
         Ploop=zeros(100,3);
         k=0;
         
         for j=1:nel
            
            % get node index that consists faces of qloop
            % iA, iB
            
            iA=Qedges(Qloops(kk+1),1);
            iB=Qedges(Qloops(kk+1),2);
            
            %%
            Ploop(k+1,:)=(pts(iA,:)+pts(iB,:))/2;
            k=k+1;
            kk=kk+1;
            
            %Ploop(k+1,:) = Ploop(k+1,:)./sqrt(sum(Ploop(k+1,:).*Ploop(k+1,:)));
            
         end
         Ploop=Ploop(1:k,:);
         h=plot3(Ploop(:,1),Ploop(:,2),Ploop(:,3));
         
         set(h,'Color',clr(colr(ColIndx(i))),'LineStyle','-','Marker','.');
      end
   end
end
axis equal


disp('Writing out Tecplot data...quads.dat');
fid=fopen('quads.dat','w');
fprintf(fid,'%s\n','VARIABLES="X","Y","Z"');
eval(sprintf('fprintf(fid,''%%s\\n'',''ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'');',length(pts(:,1)),12*nt));
for i=1:length(pts(:,1))
   fprintf(fid,'%g %g %g\n',pts(i,:));
end
for i=1:length(Qconn(:,1))
   fprintf(fid,' %d %d %d %d \n',Qconn(i,1),Qconn(i,2),Qconn(i,3),Qconn(i,4));
end
fclose(fid);

% write strand grids
fid=fopen('Strandgrid.dat','w');
% strand point of each layer

for j=1:length(Qconn(:,1))
   fprintf(fid,'%s\n','VARIABLES="X","Y","Z"');
   fprintf(fid,'%s\n','ZONE');
   
   %averaing for surface nodes
   xxx = 0;
   yyy = 0;
   zzz = 0;
   for i=1:4
      xxx = xxx + pts(Qconn(j,i),1);
      yyy = yyy + pts(Qconn(j,i),2);
      zzz = zzz + pts(Qconn(j,i),3);
   end
   
   aaa(j,1) = xxx/4;
   aaa(j,2) = yyy/4;
   aaa(j,3) = zzz/4;
   
   strand_grid = aaa(j,:);
   fprintf(fid,'%d %d %d\n', strand_grid);
   
   %averaging for layer nodes
   for k=1:nlayer
      xxx = 0;
      yyy = 0;
      zzz = 0;
      
      for i=1:4
         xxx = xxx + pts2(Qconn(j,i)+(k-1)*(nn+3*neT+7*nt),1);
         yyy = yyy + pts2(Qconn(j,i)+(k-1)*(nn+3*neT+7*nt),2);
         zzz = zzz + pts2(Qconn(j,i)+(k-1)*(nn+3*neT+7*nt),3);
      end
      
      aaa(j,1) = xxx/4;
      aaa(j,2) = yyy/4;
      aaa(j,3) = zzz/4;
      
      strand_grid = aaa(j,:);
      fprintf(fid,'%d %d %d\n', strand_grid);
   end
   
end
fclose(fid);
disp('Complete strand grid ...');


%%%%  To write node points of layers %%%
disp('Writing out Tecplot data...quads_layer');
fid=fopen('quads_layer.dat','w');
fprintf(fid,'%s\n','VARIABLES="X","Y","Z"');
eval(sprintf('fprintf(fid,''%%s\\n'',''ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL'');',length(pts(:,1))*nlayer,12*nt*nlayer));
for i=1:length(pts2(:,1))
   fprintf(fid,'%g %g %g\n',pts2(i,:));
end

for i=1:length(Qconn(:,1))
   fprintf(fid,' %d %d %d %d \n',Qconn(i,1),Qconn(i,2),Qconn(i,3),Qconn(i,4));
end

for j=1:nlayer-1
   for i=1:length(Qconn(:,1))
      fprintf(fid,' %d %d %d %d \n',Qconn(i,1)+j*(nn+3*neT+7*nt),Qconn(i,2)+j*(nn+3*neT+7*nt),Qconn(i,3)+j*(nn+3*neT+7*nt),Qconn(i,4)+j*(nn+3*neT+7*nt));
   end
end

fclose(fid);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

disp('Writing out general data...');
fid=fopen('QuadData/coord.dat','w');
fprintf(fid,'%d\n',length(pts(:,1)));

for i=1:length(pts(:,1))
   fprintf(fid,'%g %g %g\n',pts(i,:));
end
fclose(fid);
fid=fopen('QuadData/conn.dat','w');
fprintf(fid,'%d\n',length(Qconn(:,1)));
for i=1:length(Qconn(:,1))
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
strandPoints =[];

for i=1:maxcol
   disp(sprintf('Preparing Tecplot File for loops %d ...',i));
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
         P=[P;(pts(iA,:)+pts(iB,:))*0.5];
         if (ii<i2)
            conn=[conn;kk kk+1];
         end
      end
      
      
   end
   
   
   fprintf(fid,'%s\n','VARIABLES="X","Y","Z"');
   eval(sprintf('fprintf(fid,''%%s\\n'',''ZONE DATAPACKING=POINT, NODES=%d, ELEMENTS=%d, ZONETYPE=FELINESEG'');',...
      length(P(:,1)),length(conn(:,1))));
   
   for j=1:length(P(:,1))
      fprintf(fid,'%g %g %g\n',P(j,:));
   end
   for j=1:length(conn(:,1))
      fprintf(fid,'%d %d \n',conn(j,1),conn(j,2));
   end
   fclose(fid);
   
end





k=0;
%% Repeat making HAM loop for layers
for i=1:maxcol
   disp(sprintf('Preparing Tecplot File for loops %d ...',i));
   eval(sprintf('fid=fopen(''loops_layer%d.dat'',''w'');',i));
   P=[];
   conn=[];
   kk=0;
   for j=1:2*ncol(i) % go through all loops of color i
      k=k+1;
      i1=IQloops(k);
      i2=IQloops(k+1)-1;
      
      for jj=1:nlayer
         for ii=i1:i2
            ie=Qloops(ii,1);
            iA=Qedges(ie,1);
            iB=Qedges(ie,2);
            kk=kk+1;
            
            P=[P;(pts2(iA+(jj-1)*(nn+3*neT+7*nt),:)+pts2(iB+(jj-1)*(nn+3*neT+7*nt),:))*0.5];
            
            if (ii<i2)
               conn=[conn;kk kk+1];
            end
         end
      end
   end
   
   
   fprintf(fid,'%s\n','VARIABLES="X","Y","Z"');
   eval(sprintf('fprintf(fid,''%%s\\n'',''ZONE DATAPACKING=POINT, NODES=%d, ELEMENTS=%d, ZONETYPE=FELINESEG'');',...
      length(P(:,1)),length(conn(:,1))));
   
   for j=1:length(P(:,1))
      fprintf(fid,'%g %g %g\n',P(j,:));
   end
   for j=1:length(conn(:,1))
      fprintf(fid,'%d %d \n',conn(j,1),conn(j,2));
   end
   fclose(fid);
   
end





