clc
clear all

nx = 27;
ny = 27;
%nz = 10;

xi = -0;
yi = -0;
%zi = 0;

xend = 10;
yend = 10;
%zend = 1;

dx = (xend-xi)/nx;
dy = (yend-yi)/ny;
%dz = (zend-zi)/nz;
%square_grid program

%for k=1:nz
for j=1:ny
   
   for i=1:nx
      x(i,j) = xi + dx*(i-1);
      y(i,j) = yi + dy*(j-1);
      %z(i,j) = zi + dz*(k-1);
   end
end
%end


for j=1:ny
   
   if rem(j,2)==0
      for i=2:nx
         x(i,j) = x(i,j) - dx*0.5;
         %z(i,j) = zi + dz*(k-1);
      end
      x(nx+1,j) = xend -dx;
      y(nx+1,j) = y(i,j);
   end
end


% coordinate
fid = fopen('squaregrid.dat','w');
fprintf(fid,'VARIABLES="X","Y"\n');
fprintf(fid,'ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FETRIANGLE\n',nx*ny+ny/2,(nx-1)*(ny-1)*2+ny-1);

for j=1:ny
   for i=1:nx
      fprintf(fid,'%e %e\n',x(i,j),y(i,j));
   end
end

for j=1:ny
   if rem(j,2) ==0
      fprintf(fid,'%e %e\n',x(nx+1,j),y(nx+1,j));
   end
end

iflag = 0;
for j=1:ny-1
   if rem(j,2) ~=0
      iflag = iflag+1;
      fprintf(fid,'%d %d %d\n',nx*(j-1)+nx,nx*ny+iflag,nx*j+nx);
   else
      fprintf(fid,'%d %d %d\n',nx*(j-1)+nx,nx*ny+iflag,nx*j+nx);
   end
end


for j=1:ny-1
   if rem(j,2)==0
      for i=1:nx-1
         fprintf(fid,'%d %d %d\n',i+nx*(j-1),(i+1)+nx*(j-1),(i)+nx*j);
         fprintf(fid,'%d %d %d\n',(i+1)+nx*(j-1),(i+1)+nx*j,i+nx*j);
      end
   else
      for i=1:nx-1
         fprintf(fid,'%d %d %d\n',i+nx*(j-1),(i+1)+nx*(j-1),(i+1)+nx*j);
         fprintf(fid,'%d %d %d\n',i+nx*(j-1),(i+1)+nx*j,i+nx*j);
      end
   end
end

fclose(fid);

%% write data

icoord = 0;
for j=1:ny
   for i=1:nx
      icoord       = icoord + 1;
      xpos(icoord) = x(i,j);
      ypos(icoord) = y(i,j);
%       fprintf(fid,'%e %e %e\n',x(i,j),y(i,j),0);
   end
end

for j=1:ny
   if rem(j,2) ==0
      icoord       = icoord + 1;
      xpos(icoord) = x(nx+1,j);
      ypos(icoord) = y(nx+1,j);
%       fprintf(fid,'%e %e %e\n',x(nx+1,j),y(nx+1,j),0);
   end
end

fid= fopen('squaregridcoord.dat','w+');
fprintf(fid,'%d\n',icoord);
for i = 1:icoord
   fprintf(fid,'%f %f %f\n',xpos(i),ypos(i),0);
end
fclose(fid);

%connectivity
iconn = 0;
conn  = [];


iflag = 0;
for j=1:ny-1
   if rem(j,2) ~=0
      iflag = iflag+1;
      
      iconn         = iconn+1;
      conn(iconn,:) = [nx*(j-1)+nx,nx*ny+iflag,nx*j+nx];
      
%       fprintf(fid,'%d %d %d\n',nx*(j-1)+nx,nx*ny+iflag,nx*j+nx);
   else
      iconn         = iconn+1;
      conn(iconn,:) = [nx*(j-1)+nx,nx*ny+iflag,nx*j+nx];
%       fprintf(fid,'%d %d %d\n',nx*(j-1)+nx,nx*ny+iflag,nx*j+nx);
   end
end

for j=1:ny-1
   if rem(j,2)==0
      for i=1:nx-1
         iconn         = iconn+1;
         conn(iconn,:) = [i+nx*(j-1),(i+1)+nx*(j-1),(i)+nx*j];
%          fprintf(fid,'%d %d %d\n',i+nx*(j-1),(i+1)+nx*(j-1),(i)+nx*j);
         
         iconn         = iconn+1;
         conn(iconn,:) = [(i+1)+nx*(j-1),(i+1)+nx*j,i+nx*j];
%          fprintf(fid,'%d %d %d\n',(i+1)+nx*(j-1),(i+1)+nx*j,i+nx*j);
      end
   else
      for i=1:nx-1
         iconn         = iconn+1;
         conn(iconn,:) = [i+nx*(j-1),(i+1)+nx*(j-1),(i+1)+nx*j];
%          fprintf(fid,'%d %d %d\n',i+nx*(j-1),(i+1)+nx*(j-1),(i+1)+nx*j);
         
         iconn         = iconn+1;
         conn(iconn,:) = [i+nx*(j-1),(i+1)+nx*j,i+nx*j];
%          fprintf(fid,'%d %d %d\n',i+nx*(j-1),(i+1)+nx*j,i+nx*j);
      end
   end
end

fid = fopen('squaregridconn.dat','w+');
fprintf(fid,'%d\n',iconn);
for i = 1:iconn
   fprintf(fid,'%d %d %d\n',conn(i,1),conn(i,2),conn(i,3));
end
fclose(fid);

figure(1)
trimesh(conn,xpos,ypos);

disp('program end');



