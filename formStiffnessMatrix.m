function [K, argout] = formStiffnessMatrix(nodeCoords, IEN, elementType, numGP, CMatrix, varargin)
% Compute stiffness matrix.
%
% INPUT:
% nodeCoords - #Nodes-by-#Dim matrix, coordinates of nodes.
% IEN - #El-by-#NpE array, element-node incidence matrix. 
% elementType - single string, the type of elements.
% numGP - single scalar, number of Gauss points per element.
% CMatrix - a #Eq*#Dim-by-#Eq*#Dim matrix (grad(u)-->grad(u)*) or a
%   function handle returning such matrix. The CMatrix input is
%   x_ (1-by-#Dim), u_ (1-by-#Eq), gradu_ (#Dim-by-#Eq), s_ (1-by-k) -
%   values of position, displacement, displacement gradient and state
%   variables at a Gauss points. x_, u_, and gradu_ are interpolated; s_ is
%   provided by the user. 
%   Function CMatrix may have a second output argument, which is a 1-by-#AO
%   cell array. This additional output is then collected at Gauss points
%   and returned at the global level as additional output of
%   formStiffnessMatrix().
% u (optional) - #Nodes-by-#Eq matrix, nodal displacements provided for
%   evaluation of CMatrix. Default value - all zeros.
% s (optional) - #El*#GP-by-k matrix, state variables at GPs.
%   s((m-1)*#GP+l,:) is passed to CMatrix at l-th Gauss point
%   of m-th element. Default value - all zeros, k=1.
% 
% OUTPUT:
% K - (#Nodes*#Eq)-by-(#Nodes*#Eq) stiffness matrix
% args - 1-by-#AO cell array, where #AO is the number of additional outputs of
%   CMatrix(). args{i} is a #El*#GP*n1(i)-by-n2(i) array, where [n1(i)
%   n2(i)] is the size of i-th additional output of CMatrix().
%
% COMMENTS:
% This routine uses the gradH-matrix. See formStiffnessMatrixEng for a
% similar procedure using the B-matrix.
%
% Examples:
%   K = formStiffnessMatrix(nodeCoords, IEN, elementType, numGP, CMat)
%   K = formStiffnessMatrix(nodeCoords, IEN, elementType, numGP, CMat, u)
%   K = formStiffnessMatrix(nodeCoords, IEN, elementType, numGP, CMat, [], s)
%   K = formStiffnessMatrix(nodeCoords, IEN, elementType, numGP, CMat, u, s)
%   [K, argout] = formStiffnessMatrix(...)


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
numEqU=size(u,2);

additionalOutput=0;
if isa(CMatrix,'function_handle')
    isnonconstant=1;
    % check if additional output is requested and generated
    if (nargout>1) && ( nargout(CMatrix)>1 || nargout(CMatrix)<0 )
        additionalOutput=1;
        % need to know dimensions of additional output
        [C_sample, args_loc]=...
            CMatrix(nodeCoords(1,:),u(1,:),zeros(1,numEqU*numDim),s(1,:));
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
            CMatrix(nodeCoords(1,:),u(1,:),zeros(1,numEqU*numDim),s(1,:));
    end
else
    isnonconstant=0;
    C_sample=CMatrix; % Number of scalar equations
    if (u_provided==1)||(s_provided==1)
        warning('Arguments were provided for evaluation of constant stiffness matrix (possible error)');
    end
end

numEq=length(C_sample)/numDim; % Number of scalar equations
numGDoF=numNodes*numEq; %total number of displacement degrees of freedom

if abs(round(numEq)-numEq)>1e-6 %check if numEq is integer
    error('CMatrix dimension is inconsistent');
end

if exist(['shapeFunctions' elementType],'file')==2
    shapeFunctions=str2func(['shapeFunctions' elementType]);
    % Assign a function handle to shapeFunctions
else
    error(['Can''t find shape functions for specified element type: "' ...
        elementType '"' ]);
end
 
%% ==== Local stiffness matrix assembly ==
    function KLocal = formLocalStiffness(elementNodeCoords, shapeFunctions,...
            gaussPts, gaussWeights,CMat)
    % forms a stiffness matrix for a given element
        numNpE=size(elementNodeCoords,1); %number of nodes in this element
        numLocalDoF=numNpE*numEq; %number of element degrees of freedom
        KLocal=zeros(numLocalDoF);
        gradH=zeros(numDim*numEq,numLocalDoF);
        for l=1:numel(gaussWeights) %loop over Gauss points
            [~, naturalDerivatives]=shapeFunctions(gaussPts(l,:));
            %Jacobian from natural to spatial coordinates
            J=elementNodeCoords'*naturalDerivatives'; 
            spatialDerivatives=naturalDerivatives'/J;
            if numDim==2
                for ieq=1:numEq
                    gradH((ieq-1)*numDim+1,ieq:numEq:end)=spatialDerivatives(:,1)';
                    gradH((ieq-1)*numDim+2,ieq:numEq:end)=spatialDerivatives(:,2)';
                end
            elseif numDim==3
                for ieq=1:numEq
                    gradH((ieq-1)*numDim+1,ieq:numEq:end)=spatialDerivatives(:,1)';
                    gradH((ieq-1)*numDim+2,ieq:numEq:end)=spatialDerivatives(:,2)';
                    gradH((ieq-1)*numDim+3,ieq:numEq:end)=spatialDerivatives(:,3)';
                end
            else
                error('2D or 3D problem expected');
            end
            
            KLocal=KLocal+gradH'*CMat*gradH*sqrt(det(J'*J))*gaussWeights(l);
                
        end %loop over Gauss points
    end %formLocalStiffness

    function [KLocal, args_el] = formLocalStiffnessVar(elementNodeCoords, shapeFunctions,...
            gaussPts, gaussWeights,CMat,elementNodeVals,elementStVars)
    % forms a stiffness matrix for a given element
        numNpE=size(elementNodeCoords,1); %number of nodes in this element
        numLocalDoF=numNpE*numEq; %number of element degrees of freedom
        KLocal=zeros(numLocalDoF);
        gradH=zeros(numDim*numEq,numLocalDoF);
        H=zeros(numEq,numLocalDoF);
        if additionalOutput==1
            args_el=cell(1,numAO);
            for i_=1:numAO
                args_el{i_}=zeros(numGP*sizeAO(i_,1),sizeAO(i_,2));
            end
        end
        for l=1:numel(gaussWeights) %loop over Gauss points
            [shapeFunctionsVals, naturalDerivatives]=shapeFunctions(gaussPts(l,:));
            J=elementNodeCoords'*naturalDerivatives'; %Jacobian from natural to spatial coordinates
            spatialDerivatives=naturalDerivatives'/J;
            if numDim==2
                for ieq=1:numEq
                    gradH((ieq-1)*numDim+1,ieq:numEq:end)=spatialDerivatives(:,1)';
                    gradH((ieq-1)*numDim+2,ieq:numEq:end)=spatialDerivatives(:,2)';
                end
            elseif numDim==3
                for ieq=1:numEq
                    gradH((ieq-1)*numDim+1,ieq:numEq:end)=spatialDerivatives(:,1)';
                    gradH((ieq-1)*numDim+2,ieq:numEq:end)=spatialDerivatives(:,2)';
                    gradH((ieq-1)*numDim+3,ieq:numEq:end)=spatialDerivatives(:,3)';
                end
            else
                error('2D or 3D problem expected');
            end
            for ieq=1:numEq
                H(ieq,ieq:numEq:end)=shapeFunctionsVals;
            end
            
            elementGPCoords=shapeFunctionsVals*elementNodeCoords;
            elementGPVals=shapeFunctionsVals*elementNodeVals;
            graduVals=spatialDerivatives'*elementNodeVals;
            if additionalOutput==1
                [CMatVals, args_loc]=...
                    CMat(elementGPCoords,elementGPVals,graduVals,elementStVars(l,:));
                for i_=1:numAO
                    args_el{i_}(sizeAO(i_,1)*(l-1)+1:sizeAO(i_,1)*l,:)=...
                        args_loc{i_};
                end
            else
                CMatVals=...
                    CMat(elementGPCoords,elementGPVals,graduVals,elementStVars(l,:));
            end
            KLocal=KLocal+gradH'*CMatVals*gradH*sqrt(det(J'*J))*gaussWeights(l);
        end %loop over Gauss points
    end %formLocalStiffness

%% ==== Assembly of global stiffness matrix ====
[gaussPts, gaussWeights]=gaussPoints(elementType,numGP);
K=zeros(numGDoF,numGDoF); % initialize global stiffness matrix
if ~isnonconstant %linear problem, constant CMatrix
    for m=1:numEl %loop over elements
        elementNodesIDs=IEN(m,:); %global IDs of current element's nodes
        elementNodeCoords=nodeCoords(elementNodesIDs,:); %current element's nodes
        Klocal = formLocalStiffness(elementNodeCoords, shapeFunctions,...
            gaussPts, gaussWeights,CMatrix); %current element's stiffness matrix
        elementDoFIDs=node2DoFs(elementNodesIDs',numEq);
        K(elementDoFIDs,elementDoFIDs)=K(elementDoFIDs,elementDoFIDs)+Klocal;
    end %loop over elements
else %nonlinear problem or varying CMatrix
    for m=1:numEl %loop over elements
        elementNodesIDs=IEN(m,:); %global IDs of current element's nodes
        elementNodeCoords=nodeCoords(elementNodesIDs,:); %current element's nodes
        elementNodeVals=u(elementNodesIDs,:); %displacement values at current el's nodes
        elementStVars=s( (numGP*(m-1)+1):(numGP*m) , :); %state variables at current el's Gauss points
        %get current element's stiffness matrix [and additional output]
        if additionalOutput==1
            [Klocal, args_el]=formLocalStiffnessVar(elementNodeCoords, shapeFunctions,...
                gaussPts, gaussWeights,CMatrix,elementNodeVals,elementStVars);
            for i=1:numAO
                argout{i}( (numGP*sizeAO(i,1)*(m-1)+1) : numGP*sizeAO(i,1)*m , : )=...
                    args_el{i};
            end
        else
            Klocal = formLocalStiffnessVar(elementNodeCoords, shapeFunctions,...
                gaussPts, gaussWeights,CMatrix,elementNodeVals,elementStVars);
        end
        elementDoFIDs=node2DoFs(elementNodesIDs',numEq);
        K(elementDoFIDs,elementDoFIDs)=K(elementDoFIDs,elementDoFIDs)+Klocal;
    end %loop over elements
end
end %function