function M = formMassMatrix(nodeCoords, IEN, elementType, numGP, MMatrix,varargin)
% Compute mass matrix.
%
% INPUT:
% nodeCoords - #Nodes-by-#Dim matrix, coordinates of nodes.
% IEN - #El-by-#NpE array, element-node incidence matrix. 
% elementType - single string, the type of elements.
% numGP - single scalar, number of Gauss points per element.
% MMatrix - a #Eq-by-#Eq matrix (du/dt-->du/dt*) or a
%   function handle returning such matrix. The function input is
%   x_ (1-by-#Dim), u_ (1-by-#Eq), gradu_ (#Dim-by-#Eq), s_ (1-by-k) -
%   values of position, displacement, displacement gradient and state
%   variables at a Gauss points. x_, u_, and gradu_ are interpolated; s_ is
%   provided by the user. 
% u (optional) - #Nodes-by-#Eq matrix, nodal displacements provided for
%   evaluation of MMatrix. Default value - all zeros.
% s (optional) - #El*#GP-by-k matrix, state variables at GPs.
%   s((m-1)*#GP+l,:) is passed to MMatrix at l-th Gauss point
%   of m-th element. Default value - all zeros, k=1.
% 
% OUTPUT:
% M - (#Nodes*#Eq)-by-(#Nodes*#Eq) mass matrix
%
% Examples:
%   M = formMassMatrix(nodeCoords, IEN, elementType, numGP, MMat)
%   M = formMassMatrix(nodeCoords, IEN, elementType, numGP, MMat, u)
%   M = formMassMatrix(nodeCoords, IEN, elementType, numGP, MMat, [], s)
%   M = formMassMatrix(nodeCoords, IEN, elementType, numGP, MMat, u, s)

numEl=size(IEN,1); %total number of elements in the mesh
numNodes=size(nodeCoords,1); %total number of nodes in the mesh 
numDim=size(nodeCoords,2); %the dimensionality of the underlying space
%% ==== Parsing and validating input ==== %
switch numDim % check geometric dimensionality
    case {2,3}
    otherwise
        error('This function only works for 2D or 3D problems')
end

if numel(varargin)==2 % parse u, s if provided
    u_provided=1; u=varargin{1};
    s_provided=1; s=varargin{2};
elseif numel(varargin)==1
    u_provided=1; u=varargin{1};
    s_provided=0; s=[];
elseif numel(varargin)==0
    u_provided=0; u=[];
    s_provided=0; s=[];
else 
    error('Too many input argument');
end
if isempty(u)
    u=zeros(numNodes,1);
end
if isempty(s)
    s=zeros(numEl*numGP,1);
end

if isa(MMatrix,'function_handle')
    isnonconstant=1;
    M_sample=...
        MMatrix(nodeCoords(1,:),u(1,:),zeros(1,numEqU*numDim),s(1,:));
else
    isnonconstant=0;
    M_sample=MMatrix;
    if (u_provided==1)||(s_provided==1)
        warning('Arguments were provided for evaluation of constant stiffness matrix (possible error)');
    end
end

numEq=length(M_sample); % Number of scalar equations
numGDoF=numNodes*numEq; %total number of displacement degrees of freedom

if exist(['shapeFunctions' elementType],'file')==2
    shapeFunctions=str2func(['shapeFunctions' elementType]);
    % Assign a function handle to shapeFunctions
else
    error(['Can''t find shape functions for specified element type: "' ...
        elementType '"' ]);
end

%% ==== The following subroutine assembles element local mass matrix ====
    function MLocal = formLocalMass(elementNodeCoords, shapeFunctions,...
            gaussPts, gaussWeights,MMat)
    % forms a mass matrix for a given element
        numNpE=size(elementNodeCoords,1); %number of nodes in this element
        numLocalDoF=numNpE*numEq; %number of element degrees of freedom
        MLocal=zeros(numLocalDoF);
        H=zeros(numEq,numLocalDoF);
        for l=1:numel(gaussWeights)
            [shapeFunctionsVals, naturalDerivatives]=shapeFunctions(gaussPts(l,1),gaussPts(l,2));
            J=elementNodeCoords'*naturalDerivatives'; %Jacobian from natural to spatial coordinates
            for ieq=1:numEq
                H(ieq,ieq:numEq:end)=shapeFunctionsVals;
            end
            MLocal=MLocal+H'*MMat*H*sqrt(det(J'*J))*gaussWeights(l);
        end
    end %formLocalMass

    function MLocal = formLocalMassVar(elementNodeCoords, shapeFunctions,...
            gaussPts, gaussWeights,MMat,elementNodeVals,elementStVars)
    % forms a mass matrix for a given element
        numNpE=size(elementNodeCoords,1); %number of nodes in this element
        numLocalDoF=numNpE*numEq; %number of element degrees of freedom
        MLocal=zeros(numLocalDoF);
        H=zeros(numEq,numLocalDoF);
        for l=1:numel(gaussWeights)
            [shapeFunctionsVals, naturalDerivatives]=shapeFunctions(gaussPts(l,1),gaussPts(l,2));
            J=elementNodeCoords'*naturalDerivatives'; %Jacobian from natural to spatial coordinates
            for ieq=1:numEq
                H(ieq,ieq:numEq:end)=shapeFunctionsVals;
            end
            
            elementGPCoords=shapeFunctionsVals*elementNodeCoords;
            elementGPVals=shapeFunctionsVals*elementNodeVals;
            graduVals=spatialDerivatives'*elementNodeVals;
            MMatVals=MMat(elementGPCoords,elementGPVals,graduVals,elementStVars(l,:));
            MLocal=MLocal+H'*MMatVals*H*sqrt(det(J'*J))*gaussWeights(l);
        end
    end %formLocalMass


%% ==== Assembly of global mass matrix ====
[gaussPts, gaussWeights]=gaussPoints(elementType,numGP);
M=zeros(numGDoF,numGDoF); %initialize global mass matrix
if ~isnonconstant %linear problem
    for m=1:numEl %loop over elements
        elementNodesIDs=IEN(m,:); %global IDs of current element's nodes
        elementNodeCoords=nodeCoords(elementNodesIDs,:); %current element's nodes
        Mlocal = formLocalMass(elementNodeCoords, shapeFunctions,...
            gaussPts, gaussWeights,MMatrix);
        elementDoFIDs=node2DoFs(elementNodesIDs',numEq);
        M(elementDoFIDs,elementDoFIDs)=M(elementDoFIDs,elementDoFIDs)+Mlocal;
    end %loop over elements
else %nonlinear problem
    for m=1:numEl %loop over elements
        elementNodesIDs=IEN(m,:); %global IDs of current element's nodes
        elementNodeCoords=nodeCoords(elementNodesIDs,:); %current element's nodes
        elementNodeVals=u(elementNodesIDs,:); %displacement values at current el's nodes
        elementStVars=s( (numGP*(m-1)+1):(numGP*m) , :); %state variables at current el's Gauss points
        Mlocal = formLocalMassVar(elementNodeCoords, shapeFunctions,...
            gaussPts, gaussWeights,MMatrix,elementNodeVals,elementStVars);
        elementDoFIDs=node2DoFs(elementNodesIDs',numEq);
        M(elementDoFIDs,elementDoFIDs)=M(elementDoFIDs,elementDoFIDs)+Mlocal;
    end %loop over elements
end
end %function