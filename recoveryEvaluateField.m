function u2_i=recoveryEvaluateField(u2,IEN,elementType,el,ptNatCoords)
% Interpolate values at points given by natural coordinates from nodal
% values.
%
% INPUT:
% u2 - #Nodes-by-#Eq matrix containing nodal field values (e.g.
%   displacements);
% IEN - #El-by-#NpE array, element-node incidence matrix. 
% elementType - single string, the type of elements.
% el - k-by-1 element IDs of the query points.
% ptNatCoords - k-by-#ElDim natural coordinates of the query points.
%
% OUTPUT:
% u2_i - #Eq*#Dim-by-k matrix containing requested values, column
%   gradu2_i((i-1)*#Dim+j,m) corresponds to du_i/dx_j evaluated at m-th
%   query point.
%
% COMMENTS:
% Since isoparametric elements are used, the procedure is identical to
% meshGetGlobal();

%% 
numEq=size(u2,2);
numPts=size(ptNatCoords,1);

if exist(['shapeFunctions' elementType],'file')==2
    shapeFunctions=str2func(['shapeFunctions' elementType]);
else
    error(['Can''t find shape functions for specified element type: "' ...
        elementType '"' ]);
end

%% subroutine
    function nodalVals=interpolateValue(nodalVals,curPt0)
        shapeFunctionsVals=shapeFunctions(curPt0);
        nodalVals=shapeFunctionsVals*nodalVals;
    end %subfunction

%% main 
u2_i=zeros(numPts,numEq);
for i=1:numPts
    u2_i(i,:)=interpolateValue(u2(IEN(el(i),:),:),ptNatCoords(i,:));
end %loop over pts

end%function