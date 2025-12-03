function e_n=recoveryAvgEng(u2,nodeCoords,IEN,elementType)
% Recover strains at nodes by averaging values obtained by interpolation
% within adjacent elements.
%
% INPUT:
% u2 - #Nodes-by-#Eq matrix of displacements.
% nodeCoords - #Nodes-by-#Dim matrix, coordinates of nodes.
% IEN - #El-by-#NpE array, element-node incidence matrix. 
% elementType - single string, the type of elements.
%
% OUTPUT:
% e_n - 3-by-#Nodes or 6-by-#El*#GP matrix containing nodal values of
%   strains recovered using the procedure, so that column e_n(:,i)
%   equals to strains evaluated at i-th node.

%% 
nodeNatCoords=elementData(elementType,'nodeNatCoords');

numEq=size(u2,2);
numNodes=size(nodeCoords,1);
numDim=size(nodeCoords,2);
numEl=size(IEN,1);
%%
if numEq==2
    e_n=zeros(3,numNodes);
elseif numEq==3
    e_n=zeros(6,numNodes);
else
    error('2D or 3D elasticity problem expected');
end
cnt=zeros(1,numNodes);
for i=1:numEl
    e_n(:,IEN(i,:))=e_n(:,IEN(i,:))+...
        recoveryEvaluateStrains(u2,nodeCoords,IEN,elementType,...
        repmat(i,size(nodeNatCoords,1),1),nodeNatCoords);
    cnt(IEN(i,:))=cnt(IEN(i,:))+1;
end
e_n=e_n./repmat(cnt,size(e_n,1),1);
end %function