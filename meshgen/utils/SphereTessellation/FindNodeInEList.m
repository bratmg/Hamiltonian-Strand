%% ##################################################################
% Function that identifies if a mid-node has been previously created
% and if not, creates a new one. Required to identify the connectivity
% of the triangles
% ###################################################################
function [nodeID,eListOut,vertexOut] = FindNodeInEList(nodes,eListIn,vertexIn)

   eListOut  = eListIn;
   vertexOut = vertexIn;

	nodeMin = min(nodes); % minimum of the two node values
   nodeMax = nodes(find(nodes~=nodeMin));

   % identify the index in eList where this particular edge exist
   index = [];
   temp1 = [];

   if(isempty(eListIn)==0)
      temp1 = find(eListIn(:,1)==nodeMin);
   end 

   % if the first node exists, try to find the second node
   if(isempty(temp1)==0)   
      temp2 = find(eListIn(temp1,2)==nodeMax);
      index = temp1(temp2);
   end      


   if(isempty(index)==1)
      % If empty index. Implies the edge does not yet exist. Therefore
      % add to eList and create nodeID for the mid node.
      % if empty index create node 4
      tempVertex = 0.5*(vertexIn(nodeMin,:) + vertexIn(nodeMax,:));
      dist       = sqrt(tempVertex(1)^2 + tempVertex(2)^2 + tempVertex(3)^2);

      % add the node to vertex list
      vCount              = size(vertexIn,1) + 1;
      vertexOut(vCount,:) = 1/dist.*tempVertex;

      % add this edge to eList. An indication that the mid node exists
      eListCount             = size(eListIn,1) + 1;
      eListOut(eListCount,:) = [nodeMin nodeMax vCount];

      % node ID for node4
      nodeID = vCount;
   else
      % both the nodes exist and therefore the middle node has been
      % created previously and the node ID is the third column
      % in eList(index) array
      nodeID = eListIn(index,3);
   end
% ###################################################################
% END OF FILE
% ###################################################################