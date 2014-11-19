%% ##################################################################
%
% FindEdges.m
%
% Function to find the edges of a polygon
%
% Inputs:  nv    - number of vertices
%          cc    - set of nv sided polygons
%          nc    - total number of polygons
%
% Outputs: ne    - number of edges
%          edges - list of edges
%
% Written by Jayanarayanan Sitaraman
% Comments added by Bharath Govindarajan
% ###################################################################
function [ne,edges]=FindEdges(nv,cc,nc)

me = nc*4; % maximum number of edges (I think just a safe initial upper bound)
nn = 0;

% These loops find 'nn' which turns out to be the total
% number of data points
for i=1:nc % loop over total number of polygons
  for m=1:nv % loop over each vertex
    nn=max(nn,cc(m,i));
  end
end

iptr  = zeros(nn,1);
etmp  = zeros(5,me);
iflag = zeros(nc,1);
ne    = 0; % set in InsertEdge function

% loop over total number of polygons
for i=1:nc
%  disp(sprintf('processing triangle %d out of %d (%g)',i,nc,i/nc*100));

   % loop over number of vertices per polygon
   for j=1:nv

      jp1     = mod(j,nv)+1; % sets the number [1,2,3,..,nv] (nv = 3 for triangle)
      % e.g., j = [1,2,3]; jp1 = [2,3,1]

      % two ends of an edge (3 for triangle)
      eloc(1) = cc(j,i);
      eloc(2) = cc(jp1,i);
      
      
      [ne,iptr,etmp] = InsertEdge(nn,me,eloc,i,ne,iptr,etmp);

   end
end

% first four enteries of every column
% Cols 1 and 2 are the node IDs of the edge
% Col 3 is the polygon/triangle ID the edge belongs to
% Col 4 is triangle ID with which the edge is shared (if any, 0 is none)
edges=etmp(1:4,1:ne); % 'ne' should be hopefully less than 'me'

% ###################################################################
% END OF FILE
% ###################################################################