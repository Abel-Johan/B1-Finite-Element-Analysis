function [K, argout] = formStiffnessMatrix_inc(nodeCoords, IEN, elementType, numGP, CMatrix, u, du, s)
% Compute stiffness matrix. Differs from formStiffnessMatrix() in that
%   displacement increment du is provided for evaluation of CMatrix.
%
% INPUT:
% nodeCoords - #Nodes-by-#Dim matrix, coordinates of nodes.
% IEN - #El-by-#NpE array, element-node incidence matrix. 
% elementType - single string, the type of elements.
% numGP - single scalar, number of Gauss points per element.
% CMatrix - a function handle returning a  #Eq*#Dim-by-#Eq*#Dim matrix
%   (grad(u)-->grad(u)*). The CMatrix input is x_ (1-by-#Dim), u_
%   (1-by-#Eq), gradu_ (#Dim-by-#Eq), s_ (1-by-k), du_ (1-by-#Eq), dgradu_
%   (#Dim-by-#Eq), - values of position, displacement, displacement
%   gradient, state variables, displacement increment and displacement
%   gradient increment at a Gauss points. x_, u_, gradu_, ddu_, dgradu_ are
%   interpolated; s_ is provided by the user. Function CMatrix may have a
%   second output argument, which is a 1-by-#AO cell array. This additional
%   output is then collected at Gauss points and returned at the global
%   level as additional output of formStiffnessMatrix().
% u - #Nodes-by-#Eq real array, nodal displacements provided for evaluation
%   of CMatrix.
% du - #Nodes-by-#Eq real array, nodal displacement increments provided for
%   evaluation of CMatrix.
% s - #El*#GP-by-k real array, state variables at GPs.
%   s((m-1)*#GP+l,:) is passed to CMatrix at l-th Gauss point
%   of m-th element. Default value - all zeros, k=1.
% 
% OUTPUT:
% K - (#Nodes*#Eq)-by-(#Nodes*#Eq) stiffness matrix
%
% COMMENTS:
% This rountine simply rearranges input and output parameters and invokes
% formStiffnessMatrix().
%
% Example:
%   K = formStiffnessMatrix_inc(nodeCoords, IEN, elementType, numGP, CMat, u, du, s)
%   K = formStiffnessMatrix_inc(nodeCoords, IEN, elementType, numGP, CMat, [], du, s)
%   K = formStiffnessMatrix_inc(nodeCoords, IEN, elementType, numGP, CMat, u, du, [])
%   [K, argout] = formStiffnessMatrix_inc(...)

numEl=size(IEN,1); %total number of elements in the mesh
numNodes=size(nodeCoords,1); %total number of nodes in the mesh 
numDim=size(nodeCoords,2); %the dimensionality of the underlying space
numEq=numel(u)/numNodes;

u_du=[u du];

    function [CMat, argout]=CMatrix2(x,udu_loc,gradudu_loc,s)
        u_loc=udu_loc(:,1:numEq);
        du_loc=udu_loc(numEq+1:numEq*2);
        gradu_loc=gradudu_loc(:,1:numEq*numDim);
        graddu_loc=gradudu_loc(numEq*numDim+1:numEq*numDim*2);
        [CMat, argout]=CMatrix(x,u_loc,gradu_loc,s,du_loc,graddu_loc);
    end

[K, argout] = formStiffnessMatrix(nodeCoords, IEN, elementType, numGP, @CMatrix2, ...
    u_du, s);

end %function