function e_i=recoveryEvaluateStrains(u2,nodeCoords,IEN,elementType,el,ptNatCoords)
% Interpolate strains from nodal displacements at points given by natural
% coordinates.
%
% INPUT:
% u2 - #Nodes-by-#Eq matrix containing nodal field values (e.g.
%   displacements);
% nodeCoords - #Nodes-by-#Dim matrix, coordinates of nodes.
% IEN - #El-by-#NpE array, element-node incidence matrix. 
% elementType - single string, the type of elements.
% el - k-by-1 element IDs of the query points.
% ptNatCoords - k-by-#ElDim natural coordinates of the query points.
%
% OUTPUT:
% e_i - 3-by-k or 6-by-k matrix containing values of strains, column
%   e_i(:,m) corresponds to strains evaluated at m-th query point.

%% parsing input
numEq=size(u2,2);
e_i=recoveryEvaluateGradients(u2,nodeCoords,IEN,elementType,el,ptNatCoords);
if numEq==2
    A=[1 0 0 0;0 0 0 1;0 1 1 0];
elseif numEq==3
    A=[1 0 0 0 0 0 0 0 0;...
       0 0 0 0 1 0 0 0 0;...
       0 0 0 0 0 0 0 0 1;...
       0 0 0 0 0 1 0 1 0;...
       0 0 1 0 0 0 1 0 0;...
       0 1 0 1 0 0 0 0 0];
else
    error('2D or 3D elasticity problem expected');
end
e_i=A*e_i;
end%function