function [gradu_GP,GPCoords,varargout]=recoveryGP(u2,nodeCoords,IEN,elementType,varargin)
% Recover gradients at Gauss points by interpolation of nodal field values.
%
% INPUT:
% u2 - #Nodes-by-#Eq matrix containing nodal field values (e.g.
%   displacements);
% nodeCoords - #Nodes-by-#Dim matrix, coordinates of nodes.
% IEN - #El-by-#NpE array, element-node incidence matrix. 
% elementType - single string, the type of elements.
% numGP (optional) - #GP, number of Gauss points. Default value - 1.
%
% OUTPUT:
% gradu_GP - #Eq*#Dim-by-#El*#GP matrix containing values of gradients
%   recovered using the procedure, so that gradu((i-1)*#Dim+j,(m-1)*#GP+l)
%   equals to du_i/dx_j evaluated at l-th Gauss point of m-th element. 
% GPCoords - #El*#GP-by-#Dim matrix containing coordinates of Gauss points,
%   so that GPCoords((m-1)*#GP+l,:) are the coordinates of l-th Gauss point
%   of m-th element.
% GPNatCoords - #El*#GP-by-#ElDim matrix containing natural coordinates of
%   Gauss points, so that GPNatCoords((m-1)*#GP+l,:) are the coordinates of
%   l-th Gauss point of m-th element. el - #El*#GP-by-1 matrix containing
% element IDs of corresponding Gauss
%   points, i.e. first #GP rows contain 1, followed by #GP rows containing
%   2 and so on.
%
% Examples:
% [gradu,GPCoords] = recoveryGP(u2,nodeCoords,IEN,elementType,numGP) 
% [gradu,GPCoords] = recoveryGP(u2,nodeCoords,IEN,elementType) 
% [gradu,GPCoords,GPNatCoords,el] = recoveryGP(...) 


numEl=size(IEN,1);

if ~isempty(varargin)
    numGP=varargin{1};
else 
    numGP=1;
end

[gaussPts, ~]=gaussPoints(elementType,numGP);
el=repmat([1:numEl],numGP,1);
el=el(:);
ptNatCoords=repmat(gaussPts,numEl,1);
gradu_GP=recoveryEvaluateGradients(u2,nodeCoords,IEN,elementType,el,ptNatCoords);
GPCoords=recoveryEvaluateField(nodeCoords,IEN,elementType,el,ptNatCoords);

if nargout>=3
    varargout{1}=ptNatCoords;
    varargout{2}=el;
end

end %function