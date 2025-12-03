function [F, argout]= formInternalForce(nodeCoords, IEN, elementType, numGP, CMatrix,...
    u, varargin)
% Compute equivalent internal nodal force vector from given nodal 
% displacements.
%
% INPUT:
% nodeCoords - #Nodes-by-#Dim matrix, coordinates of nodes.
% IEN - #El-by-#NpE array, element-node incidence matrix. 
% elementType - single string, the type of elements.
% numGP - single scalar, number of Gauss points per element.
% CMatrix - a #Eq2*#Dim-by-#Eq1*#Dim matrix (grad(u1)-->grad(u2)*) or a
%   function handle returning such matrix. The CMatrix input is
%   x_ (1-by-#Dim), u_ (1-by-#Eq), gradu_ (#Dim-by-#Eq), s_ (1-by-k) -
%   values of position, displacement, displacement gradient and state
%   variables at a Gauss points. x_, u_, and gradu_ are interpolated; s_ is
%   provided by the user. 
%   Function CMatrix may have a second output argument, which is a 1-by-#AO
%   cell array. This additional output is then collected at Gauss points
%   and returned at the global level as additional output of
%   formStiffnessMatrix().
% u - #Nodes-by-#Eq1 matrix, nodal displacements provided for
%   evaluation of CMatrix. Default value - all zeros.
% s (optional) - #El*#GP-by-k matrix, state variables at GPs.
%   s((m-1)*#GP+l,:) is passed to CMatrix at l-th Gauss point
%   of m-th element. Default value - all zeros, k=1.
% 
% OUTPUT:
% F - #Nodes*#Eq2-by-1 force vector
%
% Examples:
%   F = formInternalForce2(nodeCoords, IEN, elementType, numGP, CMat, u)
%   F = formInternalForce2(nodeCoords, IEN, elementType, numGP, CMat, u, s)
%   [F, argout] = formInternalForce2(...)

% Implemented as a shortcut for formInternalForce2()

numEl=size(IEN,1); %total number of elements in the mesh
numNodes=size(nodeCoords,1); %total number of nodes in the mesh 
numDim=size(nodeCoords,2); %the dimensionality of the underlying space

switch numDim % check geometric dimensionality
    case {2,3}
    otherwise
        error('This function only works for 2D or 3D problems')
end
numEq1=size(u,2);

if isa(CMatrix,'function_handle')
    isnonconstant=1;
    % check if additional output is requested and generated
    if (nargout>1) && ( nargout(CMatrix)>1 || nargout(CMatrix)<0 )
        additionalOutput=1;
        % need to know dimensions of additional output
        [C_sample, args_loc]=...
            CMatrix(nodeCoords(1,:),u(1,:),zeros(1,numEq1*numDim),s(1,:));
        numAO=numel(args_loc);
        sizeAO=zeros(numAO,2);
        argout=cell(1,numAO); % initialize global storage for additional output
        for i=1:numAO
            sizeAO(i,:)=size(args_loc{i});
            argout{i}=zeros(numEl*numGP*sizeAO(i,1),sizeAO(i,2));
        end
    elseif (nargout>1)
        error('Additional output is requested, but CMatrix does not seem to return it');
    else
        C_sample=...
            CMatrix(nodeCoords(1,:),u(1,:),zeros(1,numEq1*numDim),s(1,:));
    end
elseif isa(BetaMatrix,'function_handle')
    isnonconstant=1;
    C_sample=CMatrix; % Number of scalar equations
    CMatrix=@(x,u,gradu,s)(CMatrix);
else
    isnonconstant=0;
    C_sample=CMatrix; % Number of scalar equations
    if (s_provided==1)
        warning('Arguments were provided for evaluation of constant stiffness matrix (possible error)');
    end
end
numEq2=size(C_sample,1)/numDim; % Number of scalar equations in output

BetaMatrix=zeros(numEq2*numDim,numEq1);

[F, argout]= formInternalForce2(nodeCoords, IEN, elementType, numGP, CMatrix,...
    BetaMatrix, u, varargin);

end %function