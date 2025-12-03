function [e_GP,GPCoords,varargout]=recoveryGPEng(u2,nodeCoords,IEN,elementType,varargin)
% Recover strains at Gauss points by interpolation of displacements.
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
% e_GP - 3-by-#El*#GP or 6-by-#El*#GP  matrix containing values of
%   strains recovered using the procedure, so that e_GP(i,(m-1)*#GP+l)
%   equals to i-th strain evaluated at l-th Gauss point of m-th element. 
% GPCoords - #El*#GP-by-#Dim matrix containing coordinates of Gauss points,
%   so that GPCoords((m-1)*#GP+l,:) are the coordinates of l-th Gauss point
%   of m-th element.
% GPNatCoords - #El*#GP-by-#DimNat matrix containing natural coordinates of
%   Gauss points, so that GPNatCoords((m-1)*#GP+l,:) are the coordinates of
%   l-th Gauss point of m-th element.
% el - #El*#GP-by-1 matrix containing element IDs of corresponding Gauss
%   points, i.e. first #GP rows contain 1, followed by #GP rows containing
%   2 and so on.
%
% Examples:
% [e_GP,GPCoords]=recoveryGPEng(u2,nodeCoords,IEN,elementType,numGP) 
% [e_GP,GPCoords]=recoveryGPEng(u2,nodeCoords,IEN,elementType)
% [e_GP,GPCoords,GPNatCoords,el] = recoveryGPEng(...) 



numEl=size(IEN,1);

if ~isempty(varargin)
    numGP=varargin{1};
else 
    numGP=1;
end

[gaussPts, ~]=gaussPoints(elementType,numGP);
el=repmat([1:numEl],1,numGP);
el=el(:);
ptNatCoords=repmat(gaussPts,numEl,1);
e_GP=recoveryEvaluateStrains(u2,nodeCoords,IEN,elementType,el,ptNatCoords);
GPCoords=recoveryEvaluateField(nodeCoords,IEN,elementType,el,ptNatCoords);

if nargout>=3
    varargout{1}=ptNatCoords;
    varargout{2}=el;
end

end %function