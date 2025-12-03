function [INE, degrees, localIDs] = IENtoINE(IEN)
% Compute INE, node-elements incidence array (dual to IEN up to
% permutations within rows). 
%
% INPUT:
% IEN - #El-by-#NpE array, element-node incidence matrix. 
%
% OUTPUT:
% INE - #Nodes-by-#EpN array, node-element incidence matrix, row INE(i,:)
%   corresponds to ith node and contains IDs of elements that share it.
%   Since the number of incident elements varies over nodes, the unused
%   entries are filled with zeros.
% degrees - #Nodes-by-1 matrix containing degrees/valencies of nodes,
%   degrees(i) is the number of elements in the mesh node i is incident to.
% localIDs - #Nodes-by-#EpN array, localIDs(i,j) is the local ID of node i
%   within the element INE(i,j). Thus, INE(i,:) and localIDs(i,:) contain
%   respectively the row and column subscripts of each entry of node i
%   within IEN.  
%
% Examples:
% [INE, degrees, localIDs] = IENtoINE(IEN)
% IEN2 = IENtoINE(INE);
% norm(sort(IEN2,2)-sort(IEN,2)) % equals to zero

numNodes=max(IEN(:));
numNpE=size(IEN,2);

% Let's determine the maximum degree/valency of nodes
% to allocate enough storage for INE
degrees=zeros(numNodes,1);
for i=1:numNodes
    degrees(i)=nnz(IEN==i);
end

numEpN=max(degrees);
INE=zeros(numNodes,numEpN);
localIDs=zeros(numNodes,numEpN);
degrees=zeros(numNodes,1); % we use its entries as pointers, so reset it

for i=1:size(IEN,1) % Iterate over elements
    curNodes=nonzeros(IEN(i,:))';
    degrees(curNodes)=degrees(curNodes)+1;
    INE( sub2ind([numNodes,numEpN],curNodes',degrees(curNodes)) ) = i;
    % it was assumed that there are no repetitions in curNodes
    localIDs( sub2ind([numNodes,numEpN],curNodes',degrees(curNodes)) ) = find(IEN(i,:));
end
end %function