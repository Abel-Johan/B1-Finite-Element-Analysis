function gradu_i=recoveryEvaluateGradients(u2,nodeCoords,IEN,elementType,el,ptNatCoords)
% Interpolate gradients from nodal field values at points given by natural
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
% gradu_i - #Eq*#Dim-by-k matrix containing values of gradients, 
%   gradu2_i((i-1)*#Dim+j,m) equals to du_i/dx_j evaluated at m-th
%   query point.

%% 
numEq=size(u2,2);
numPts=size(ptNatCoords,1);
if exist(['shapeFunctions' elementType],'file')==2
    shapeFunctions=str2func(['shapeFunctions' elementType]);
else
    error(['Can''t find shape functions for specified element type: "' ...
        elementType '"' ]);
end
numDim=nargin(shapeFunctions);
%% subroutine
    function Vals=interpolateGradient(elementNodeCoords,nodalVals,curPt0)
        [~, naturalDerivatives]=shapeFunctions(curPt0);
        J=elementNodeCoords'*naturalDerivatives'; %Jacobian from natural to spatial coordinates
        spatialDerivatives=naturalDerivatives'/J; %same as naturalDerivatives'*inv(J)
        Vals=spatialDerivatives'*nodalVals;
        Vals=Vals';
        Vals=Vals(:)';
    end %subfunction

%% main 
gradu_i=zeros(numEq*numDim,numPts);
for i=1:numPts
    gradu_i(:,i)=interpolateGradient(nodeCoords(IEN(el(i),:),:),u2(IEN(el(i),:),:),ptNatCoords(i,:));
end %loop over query pts

end%function