%% ##################################################################
%
% InsertEdge.m
%
% Function to insert an edge
%
% Inputs:  nn    - total number of data points
%          me    - 
%          eloc  - two vertices of the edge
%          ci    - polygon ID (1,2,3,...)
%          ne1   - running counter for the number of edges
%          iptr1 - running counter for number of edges associated with a node
%          edge1 - running counter for edge array
%
% Outputs: ne    - running counter for the number of edges
%          iptr  - running counter for number of edges associated with a node
%          edge  - running counter for edge array
%
% Written by Jayanarayanan Sitaraman
% Comments added by Bharath Govindarajan
% ###################################################################
function [ne,iptr,edge]=InsertEdge(nn,me,eloc,ci,ne1,iptr1,edge1)

ne   = ne1;
iptr = iptr1;
edge = edge1;
e1   = eloc;  % temp variable that mirrors eloc

% Swap procedure to ensure e1(2) > e1(1)
if (e1(1)>e1(2)) 
   te    = e1(1);
   e1(1) = e1(2);
   e1(2) = te;
end   

ip = iptr(e1(1));


while (ip>0)

   e2 = edge(1:2,ip)';

   % Swap procedure to ensure e2(2) > e2(1)
   if (e2(1)>e2(2))
     te    = e2(1);
     e2(1) = e2(2);
     e2(2) = te;
   end

   % if both edges e1 and e2 are identical
   if (sum(abs(e1-e2))==0)
      
     edge(4,ip) = ci;
     return;
   end

   ip = edge(5,ip);
   
   
end

ne           = ne+1;
edge(1:2,ne) = eloc;
edge(3,ne)   = ci;
edge(5,ne)   = iptr(e1(1));
iptr(e1(1))  = ne;


return
% ###################################################################
% END OF FILE
% ###################################################################