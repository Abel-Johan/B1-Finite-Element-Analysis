function gradu_n=recoveryAvg(u2,nodeCoords,IEN,elementType)
% Recover gradients at nodes by averaging values obtained by interpolation
% within adjacent elements.
%
% INPUT:
% u2 - #Nodes-by-#Eq matrix containing nodal field values (e.g.
%   displacements);
% nodeCoords - #Nodes-by-#Dim matrix, coordinates of nodes.
% IEN - #El-by-#NpE array, element-node incidence matrix. 
% elementType - single string, the type of elements.
%
% OUTPUT:
% gradu_n - #Eq*#Dim-by-#Nodes matrix containing nodal values of gradients
%   recovered using the procedure, so that gradu_n((i-1)*#Dim+j,m)
%   equals to du_i/dx_j evaluated at m-th node.

%% 
nodeNatCoords=elementData(elementType,'nodeNatCoords');

numEq=size(u2,2);
numNodes=size(nodeCoords,1);
numDim=size(nodeCoords,2);
numEl=size(IEN,1);
%%
gradu_n=zeros(numDim*numEq,numNodes);
cnt=zeros(1,numNodes);
for i=1:numEl
    gradu_n(:,IEN(i,:))=gradu_n(:,IEN(i,:))+...
        recoveryEvaluateGradients(u2,nodeCoords,IEN,elementType,...
        repmat(i,size(nodeNatCoords,1),1),nodeNatCoords);
    cnt(IEN(i,:))=cnt(IEN(i,:))+1;
end
gradu_n=gradu_n./repmat(cnt,numDim*numEq,1);
end %function