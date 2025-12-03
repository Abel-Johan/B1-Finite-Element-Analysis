function grad_n=recoveryPatch(u2,nodeCoords,IEN,elementType,BIEN,varargin)
% Recover gradients at nodes from given nodal field values using
% superconvergent patch recovery (SPR) technique.
%
% INPUT:
% u2 - #Nodes-by-#Eq matrix containing nodal field values (e.g.
%   displacements); 
% nodeCoords - #Nodes-by-#Dim matrix, coordinates of nodes.
% IEN - #El-by-#NpE array, element-node incidence matrix. 
% elementType - single string, the type of elements.
% BIEN - a 1-by-#BR cell array containing boundary incidence matrices.
% INE (optional) - #Nodes-by-#EpN array, node-element incidence matrix, row
%   INE(i,:) corresponds to ith node and contains IDs of elements that
%   share it.
% degrees (optional) - #Nodes-by-1 matrix containing degrees (valencies) of
%   nodes, i.e. the number of elements in the mesh each node belongs to.
%
% OUTPUT:
% grad_n - #Eq*#Dim-by-#Nodes matrix containing nodal values of gradients
%   recovered using the procedure, so that grad_n((i-1)*#Dim+j,m)
%   corresponds to du_i/dx_j evaluated at m-th node.
%
% COMMENTS:
%   INE and degrees should be provided together. If not provided, they are
%   computed using IENtoINE

%% Compute gradient values at Gauss Points
elDat=elementData(elementType);
numGP=elDat.numGPSPR;
[gradu,GPCoords]=recoveryGP(u2,nodeCoords,IEN,elementType,numGP);
grad_n=recoveryPatchGP(gradu,GPCoords,nodeCoords,IEN,elementType,BIEN,varargin);
end %function