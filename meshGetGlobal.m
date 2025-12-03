function ptCoords=meshGetGlobal(nodeCoords,IEN,elementType,el,ptNatCoords)
% Map natural coordinates to global coordinates
%
% INPUT:
% nodeCoords - #Nodes-by-#Dim matrix, coordinates of nodes.
% IEN - #El-by-#NpE array, element-node incidence matrix. 
% elementType - single string, the type of elements.
% el - k-by-1 element IDs of the query points.
% ptNatCoords - k-by-#ElDim natural coordinates of the query points.
%
% OUTPUT:
% ptCoords - k-by-#Dim global coordinates of the query points.
%
% COMMENTS:
% Since isoparametric elements are used, the procedure is identical to
% recoveryEvaluateField();
%% 
numDim=size(nodeCoords,2);
numPts=size(ptNatCoords,1);

if exist(['shapeFunctions' elementType],'file')==2
    shapeFunctions=str2func(['shapeFunctions' elementType]);
else
    error(['Can''t find shape functions for specified element type: "' ...
        elementType '"' ]);
end

%% subroutine
    function curPt=getGlobal(elementNodeCoords,curPt0)
        shapeFunctionsVals=shapeFunctions(curPt0);
        curPt=shapeFunctionsVals*elementNodeCoords;
    end %subfunction

%% main 
ptCoords=zeros(numPts,numDim);
for i=1:numPts
    elementNodeCoords=nodeCoords(IEN(el(i),:),:);
    ptCoords(i,:)=getGlobal(elementNodeCoords,ptNatCoords(i,:));
end %loop over pts

end%function