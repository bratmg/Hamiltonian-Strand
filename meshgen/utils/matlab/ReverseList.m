%% ##################################################################
% ReverseList.m
%
% Written by Jayanarayanan Sitaraman
% Commented by Bharath Govindarajan
% ###################################################################
function [List2,Index2]=ReverseList(List1,N,n1,n2,Index1)

% Creates reverse list and associated index list
%
% List1 is a         list of elements of Set 2 per element of Set 1 (size N)
% List2 is a reverse list of elements of Set 1 per element of Set 2 (size N)
% n1 is the number of elements in Set 1
% n2 is the number of elements in Set 2
% Index1 is the starting index in List1 for each element of Set 1 (size n1+1)
% Index2 is the starting index in List2 for each element of Set 2 (size n2+1)

% n2 is number of data points in _coord.dat file
for i=1:n2
   Ninst(i)=0;
end

% N is (3*nt) - total vertex list
% the following loop builds a list of the number of
% times a given vertex appears in the _conn.dat file.
% Implies the number of edges connecting a given vertex
% e.g.: if Ninst(230) = 7, then vertex ID 230 is connected to 7 edges
for i=1:N
   Ninst(List1(i))=Ninst(List1(i))+1;
end



% First point in Index2 set to 1
Index2(1)=1;

% Loop from 2 to n2+1
% n2 is number of data points in _coord.dat file
% Index2(i) contains the Index2(i-1) + the number of
% edges connected to vertex(i-1)
for i=2:n2+1
   Index2(i)=Index2(i-1)+Ninst(i-1);
end

k2(1:n2) = 0;
k1       = 0;

% Loop over n1 (lines in _conn.dat file)
for i1=1:n1
   
   % The loop runs from the 1:3 or x:x+2
   for j=Index1(i1):Index1(i1+1)-1   
      k1     = k1+1;      % incremental counter
      
      % effectively goes through all the elements in _conn.dat file
      % n1*3 (3 = size of j loop)
      i2     = List1(k1); % store in i2 the value of 
      k2(i2) = k2(i2)+1;  % k2 performs the function similar to Ninst
      
      List2(Index2(i2)+k2(i2)-1) = i1;
      
   end
end

% ===================================================================
% E.g. List2 = [218,315,316||13,234,216||13,38,96||..]
%     Index2 = [ 1, 4, 7, ... ] (need not be offset by 3)
%
% The lines in _conn.dat file with vertex ID i can be found in 
% List2 by looking at indiced from Index(i) to Index(i+1)-1
%
% Therefore, Vertex ID 2 is contained in lines 13, 246, and 216
% ===================================================================

% ###################################################################
% END OF FILE
% ###################################################################
