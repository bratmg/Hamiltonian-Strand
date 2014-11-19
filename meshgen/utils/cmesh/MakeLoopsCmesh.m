clear all;
close all;
x=load('x.dat');
y=load('y.dat');
nra=length(x(:,1));
nth=length(x(1,:));
n1=27;

% Create list of points
pts=[reshape(x',nra*nth,1) reshape(y',nra*nth,1)];

% Make Quads
Qconn=zeros((nra-1)*(nth-1),4);
k=0;
for j=1:nra-1
  for i=1:nth-1
    k=k+1;
    i1=(j-1)*nth+i;
    i2=i1+1;
    if (i1<=n1) ii1=nth+1-i1; else ii1=i1; end
    if (i2<=n1) ii2=nth+1-i2; else ii2=i2; end
    Qconn(k,:)=[ii1 ii2 i2+nth i1+nth];
  end
end

% Make edges
Qedges=zeros(nth*(nra-1),6);
k=0;
for j=1:nra-1
  for i=1:nth
    k=k+1;
    i1=(j-1)*nth+i;
    i2=i1+nth;
    c1=(j-1)*(nth-1)+i-1;
    c2=c1+1;
    e1=2;e2=4;
    if (i==1)
      c1=0;e1=0;
    elseif(i==nth)
      c2=0;e2=0;
    end
    if (i1<=n1) i1=nth+1-i1; end
    if (i2<=n1) i2=nth+1-i2; end
    Qedges(k,:)=[i1 i2 c1 c2 e1 e2]; 
  end
end

Qedges=[Qedges;zeros((n1-1)*(2*nra-1),6)];
for i=1:n1-1
  for j=nra-1:-1:1
    k=k+1;
    i1=j*nth+i+1;
    i2=i1-1;
    c1=(j-1)*(nth-1)+i;
    c2=c1+nth-1;
    e1=3;e2=1;
    if (j==nra-1) 
      c2=0;e2=0; 
    end
    Qedges(k,1:6)=[i1 i2 c1 c2 e1 e2];
  end
  k=k+1;
  Qedges(k,1:6)=[nth-i nth-i+1 nth-i i  1 1];
  for j=1:nra-1
    k=k+1;
    i1=nth*(j+1)-i;
    i2=i1+1;
    c1=(nth-1)*(j+1)-i+1;
    c2=c1-(nth-1);
    e1=1;e2=3;
    if (j==nra-1) 
      c1=0;e1=0; 
    end
    Qedges(k,1:6)=[i1 i2 c1 c2 e1 e2];
  end
end

na=nth-1-2*(n1-1);
Qedges=[Qedges;zeros(nra*na,6)];
for i=n1+1:n1+na
  for j=1:nra
    i1=nth*(j-1)+i-1;
    i2=i1+1;
    if (i1<=n1) i1=nth+1-i1; end
    if (i2<=n1) i2=nth+1-i2; end
    c1=(nth-1)*(j-1)+i-1;
    c2=c1-(nth-1);
    e1=1;
    e2=3;
    if (j==nra) 
      c1=0;e1=0; 
    elseif (j==1)
      c2=0;e2=0;
    end
    k=k+1;
    Qedges(k,1:6)=[i1 i2 c1 c2 e1 e2];
  end
end

% Make loops
nedges=length(Qedges(:,1));
Qloops=zeros(nedges,1);
IQloops=ones(n1-1+nra-1+na+1,1);
k=0;
for j=1:nra-1
  for i=1:nth
    k=k+1;
    Qloops(k)=k;
  end
  Qloops(k-nth+1:k)=Qloops(k:-1:k-nth+1);
  IQloops(j+1)=IQloops(j)+nth;
end
istart=nra;
for i=1:n1-1
  for j=1:2*nra-1
    k=k+1;
    Qloops(k)=k;
  end
  IQloops(istart+i)=IQloops(istart+i-1)+2*nra-1;
end
istart=istart+n1-1;
for i=1:na
  for j=1:nra
    k=k+1;
    Qloops(k)=k;
  end
  IQloops(istart+i)=IQloops(istart+i-1)+nra;
end

maxcol=2;
ncol(1)=nra-1;
ncol(2)=n1-1+na;
    
npts=length(pts(:,1));
nquads=length(Qconn(:,1));
nEloops=length(Qloops(:,1));
nloops=nra-1+n1-1+na;
neQ=length(Qedges(:,1));

%
disp('Writing out Tecplot data...');
fid=fopen('quads.dat','w');
fprintf(fid,'%s\n','VARIABLES="X","Y","Z"');
fprintf(fid,'ZONE N=%d, E=%d, DATAPACKING=POINT, ZONETYPE=FEQUADRILATERAL\n', ...
	     npts,nquads);

for i=1:npts
  fprintf(fid,'%1.15g %1.15g %1.15g \n',pts(i,1:2),0.0);
end
for i=1:nquads
  fprintf(fid,' %d %d %d %d \n',Qconn(i,1:4));
end
fclose(fid);

system('mkdir QuadData')

disp('Writing out general data...');
fid=fopen('QuadData/coord.dat','w');
fprintf(fid,'%d\n',npts);
for i=1:npts
fprintf(fid,'%1.15g %1.15g %1.15g\n',pts(i,1:2),0.0);
end
fclose(fid);
fid=fopen('QuadData/conn.dat','w');
fprintf(fid,'%d\n',nquads);
for i=1:nquads
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
fprintf(fid,'%d\n',nEloops);
for i=1:nEloops
  fprintf(fid,' %d\n',Qloops(i)-1);
end
fclose(fid);
fid=fopen('QuadData/iqloops.dat','w');
fprintf(fid,'%d\n',nloops+1);
for i=1:nloops+1
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
  disp(sprintf('Preparing Tecplot File for loops %d ...',i));
  eval(sprintf('fid=fopen(''loops%d.dat'',''w'');',i));
  P=[];
  conn=[];
  kk=0;
  iloop=0;
  while (iloop<ncol(i))
    k=k+1;
    iloop=iloop+1;
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
    fprintf(fid,'%1.15g %1.15g %1.15g \n',P(j,1:2),0.0);
  end
  for j=1:length(conn(:,1))
    fprintf(fid,'%d %d \n',conn(j,1),conn(j,2));
  end
  fclose(fid);

end

    

