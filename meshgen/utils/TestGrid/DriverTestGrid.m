%% ##################################################################
%
% Test code to generate a 2D simple triangular grid and write out
% the _coord.dat and _conn.dat file for use in the Hamilton Path Code
%
% Nov 7: Bharath
% ###################################################################
close all;clear all;clc
% ===================================================================

itype = 2; % 1 - right-angled, 2 - isoceles

if( itype ~= 1 && itype ~= 2)
   disp('Incorrect itype choice');break;
end

Nx   = 53;
Ny   = 53;
Ntot = Nx*Ny;

xmin = 0;
xmax = 10;
ymin = 0;
ymax = 10;

xPos = linspace(xmin,xmax,Nx);
yPos = linspace(ymin,ymax,Ny);

deltax = abs(xPos(2)-xPos(1));
deltay = abs(yPos(2)-yPos(1));

% initialization
if(itype==2)
   Nextra = ceil(0.5*Ny);
   Ntot   = Ntot + Nextra;
end
pointCoord = zeros(Ntot,3); % x, y, and z position

% Writing the data points to file
if(itype==1); factor = 0.0;end
if(itype==2); factor = 0.5;end
for i = 1:Nx
   for j = 1:Ny
      index               = i + (j-1)*Nx;
      if (i==1 && mod(j,2)==1 && itype==2)
         pointCoord(index,:) = [xPos(i)+0.5*deltax-factor*deltax*mod(j,2)  yPos(j) 0];
      else
         pointCoord(index,:) = [xPos(i)-factor*deltax*mod(j,2)  yPos(j) 0];
      end
   end % j loop
end % i loop

% add extra points at the end for the isolceles triangles
if(itype==2)
   index = Nx*Ny;
   for j = 1:Ny
      if (mod(j,2)==1)
         index = index+1;
         pointCoord(index,:) = [xmax ymin+(j-1)*deltay 0];
      end % if loop
   end % j loop
end

fid = fopen('coordSquareGrid.dat','w+');
fprintf(fid,'%d \n',Ntot);
for i = 1:Ntot
   fprintf(fid,'%f %f %f\n',pointCoord(i,1),pointCoord(i,2),pointCoord(i,3));
end
fclose(fid);

% writing the connectivity information for the triangles
numSquares  = (Nx-1)*(Ny-1);
numTriangle = 2*numSquares;

% initialization
connIndex = zeros(numTriangle,3);
count     = 0;

for i = 1:Nx-1
   for j = 1:Ny-1
      
      if(itype==1) % Right angled triangles
         count  = count + 1;
         index1 = (i  ) + (j-1)*Nx;
         index2 = (i  ) + (j  )*Nx;
         index3 = (i+1) + (j-1)*Nx;
         
      elseif(itype==2) % Equilateral triangles
         if(mod(j,2)==1)
            count  = count + 1;
            index1 = (i  ) + (j-1)*Nx;
            index2 = (i  ) + (j  )*Nx;
            index3 = (i+1) + (j-1)*Nx;
         else
            count  = count + 1;
            index1 = (i  ) + (j-1)*Nx;
            index2 = (i  ) + (j  )*Nx;
            index3 = (i+1) + (j  )*Nx;
         end
      end % itype
      
      connIndex(count,:) = [index1 index2 index3];
   end
end

for i = 2:Nx
   for j = 1:Ny-1
      
      if(itype == 1) % right angled triangles
         count  = count + 1;
         index1 = (i  ) + (j-1)*Nx;
         index2 = (i-1) + (j  )*Nx;
         index3 = (i  ) + (j  )*Nx;
         
      elseif(itype == 2) % equilateral triangles
         if(mod(j,2)==1)
            count  = count + 1;
            index1 = (i  ) + (j-1)*Nx;
            index2 = (i-1) + (j  )*Nx;
            index3 = (i  ) + (j  )*Nx;
         else
            count  = count + 1;
            index1 = (i  ) + (j-1)*Nx;
            index2 = (i-1) + (j-1)*Nx;
            index3 = (i  ) + (j  )*Nx;
         end
      end
      
      connIndex(count,:) = [index1 index2 index3];
   end
end

% add the connectivity for the right-hand edge right angled cells
if(itype==2)
   iExtra = 1;
   index = Nx*Ny;
   for j = 1:Ny-1
      if (mod(j,2)==1)
         count  = count + 1;
         index1 = Nx    + (j-1)*Nx;
         index2 = Nx    + (j  )*Nx;
         index3 = index + iExtra;
         iExtra = iExtra + 1;
         connIndex(count,:) = [index1 index2 index3];
      end
   end
   iExtra = 2;
   index = Nx*Ny;
   for j = 1:Ny-1
      if (mod(j,2)==0)
         count  = count + 1;
         index1 = Nx    + (j-1)*Nx;
         index2 = Nx    + (j  )*Nx;
         index3 = index + iExtra;
         iExtra = iExtra + 1;
         connIndex(count,:) = [index1 index2 index3];
      end
   end % j loop
end 

lenConn = length(connIndex);
fid = fopen('connSquareGrid.dat','w+');
fprintf(fid,'%d \n',lenConn);
for i = 1:lenConn
   fprintf(fid,'%d %d %d \n',connIndex(i,1),connIndex(i,2),connIndex(i,3));
end
fclose(fid);

figure(1)
axis('square')
trimesh(connIndex,pointCoord(:,1),pointCoord(:,2));
% 
% figure(2)
% axis('square')
% plot(pointCoord(:,1),pointCoord(:,2),'r.');

% ###################################################################
% END OF FILE
% ###################################################################