function [F, argout]= formInternalForce2(nodeCoords, IEN, elementType, numGP, CMatrix,...
    BetaMatrix, u, varargin)
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
% BetaMatrix - a #Eq2*#Dim-by-#Eq1 matrix (u1-->grad(u2)) or a
%   function handle returning such matrix. The input parameters of
%   BetaMatrix are similar to those of CMatrix. Additional output is not
%   collected.
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
%   F = formInternalForce2(nodeCoords, IEN, elementType, numGP, CMat, BMat, u)
%   F = formInternalForce2(nodeCoords, IEN, elementType, numGP, CMat, BMat, u, s)
%   [F, argout] = formInternalForce2(...)

numEl=size(IEN,1); %total number of elements in the mesh
numNodes=size(nodeCoords,1); %total number of nodes in the mesh 
numDim=size(nodeCoords,2); %the dimensionality of the underlying space
%% ==== Parsing and validating input ==== %
switch numDim % check geometric dimensionality
    case {2,3}
    otherwise
        error('This function only works for 2D or 3D problems')
end

u_provided=1;
if numel(varargin)==1
    s_provided=1; s=varargin{2};
else %numel(varargin)==0
    s_provided=0; s=[];
end
if isempty(s)
    s=zeros(numEl*numGP,1);
end

numEq1=size(u,2); % Number of scalar equations in input
additionalOutput=0;
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
    if ~isa(BetaMatrix,'function_handle')
        BetaMatrix=@(x,u,gradu,s)(BetaMatrix);
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

if isa(BetaMatrix,'function_handle')
    B_sample=...
        BetaMatrix(nodeCoords(1,:),u(1,:),zeros(1,numEq1*numDim),s(1,:));
else
    B_sample=BetaMatrix;
end

numEq2=size(C_sample,1)/numDim; % Number of scalar equations in output
numGDoF2=numNodes*numEq2; %total number of displacement degrees of freedom

if (numEq1~=size(B_sample,2))||(numEq2~=size(B_sample,1)/numDim)
    error('Dimensions of CMatrix and BetaMatrix are inconsistent');
end

if (abs(round(numEq1)-numEq1)>1e-6)||(abs(round(numEq2)-numEq2)>1e-6)
    error('CMatrix dimension is inconsistent');
end

if exist(['shapeFunctions' elementType],'file')==2
    shapeFunctions=str2func(['shapeFunctions' elementType]);
    % Assign a function handle to shapeFunctions
else
    error(['Can''t find shape functions for specified element type: "' ...
        elementType '"' ]);
end
%% ==== The following subroutines assemble element local stiffness matrix ==
    function FLocal = formLocalForce(elementNodeCoords, shapeFunctions2d,...
            gaussPts, gaussWeights,CMat,BetaMat,elementNodeVals)
    % forms a stiffness matrix for a given element
        numNpE=size(elementNodeCoords,1); %number of nodes in this element
        numLocalDoF1=numNpE*numEq1; %number of element degrees of freedom (input)
        numLocalDoF2=numNpE*numEq2; %number of element degrees of freedom (output)
        FLocal=zeros(numLocalDoF2,1);
        gradH1=zeros(numDim*numEq1,numLocalDoF1);
        gradH2=zeros(numDim*numEq2,numLocalDoF2);
        H1=zeros(numEq1,numLocalDoF1);
        H2=zeros(numEq2,numLocalDoF2);
        for l=1:numel(gaussWeights) %loop over Gauss points
            [shapeFunctionsVals, naturalDerivatives]=shapeFunctions2d(gaussPts(l,1),gaussPts(l,2));
            J=elementNodeCoords'*naturalDerivatives'; %Jacobian from natural to spatial coordinates
            spatialDerivatives=naturalDerivatives'/J; %same as naturalDerivatives'*inv(J)
            if numDim==2
                for ieq=1:numEq1
                    gradH1((ieq-1)*numDim+1,ieq:numEq1:end)=spatialDerivatives(:,1)';
                    gradH1((ieq-1)*numDim+2,ieq:numEq1:end)=spatialDerivatives(:,2)';
                end
                for ieq=1:numEq2
                    gradH2((ieq-1)*numDim+1,ieq:numEq2:end)=spatialDerivatives(:,1)';
                    gradH2((ieq-1)*numDim+2,ieq:numEq2:end)=spatialDerivatives(:,2)';
                end
            elseif numDim==3
                for ieq=1:numEq1
                    gradH1((ieq-1)*numDim+1,ieq:numEq1:end)=spatialDerivatives(:,1)';
                    gradH1((ieq-1)*numDim+2,ieq:numEq1:end)=spatialDerivatives(:,2)';
                    gradH1((ieq-1)*numDim+2,ieq:numEq1:end)=spatialDerivatives(:,2)';
                end
                for ieq=1:numEq2
                    gradH2((ieq-1)*numDim+1,ieq:numEq2:end)=spatialDerivatives(:,1)';
                    gradH2((ieq-1)*numDim+2,ieq:numEq2:end)=spatialDerivatives(:,2)';
                    gradH2((ieq-1)*numDim+2,ieq:numEq2:end)=spatialDerivatives(:,2)';
                end
            else
                error('2D or 3D problem expected');
            end
            for ieq=1:numEq1
                H1(ieq,ieq:numEq1:end)=shapeFunctionsVals;
            end
            for ieq=1:numEq2
                H2(ieq,ieq:numEq2:end)=shapeFunctionsVals;
            end

            u_loc=elementNodeVals';
            u_loc=u_loc(:);
            FLocal=FLocal+gradH2'*(CMat*gradH1 + BetaMat*H1)*u_loc*sqrt(det(J'*J))*gaussWeights(l);
                
        end %loop over Gauss points
    end %formLocalStiffness

    function [FLocal, args_el] = formLocalForceVar(elementNodeCoords, shapeFunctions2d,...
            gaussPts, gaussWeights,CMat,BetaMat,elementNodeVals,elementStVars)
    % forms a stiffness matrix for a given element
        numNpE=size(elementNodeCoords,1); %number of nodes in this element
        numLocalDoF1=numNpE*numEq1; %number of element degrees of freedom (input)
        numLocalDoF2=numNpE*numEq2; %number of element degrees of freedom (output)
        FLocal=zeros(numLocalDoF2,1);
        gradH1=zeros(numDim*numEq1,numLocalDoF1);
        gradH2=zeros(numDim*numEq2,numLocalDoF2);
        H1=zeros(numEq1,numLocalDoF1);
        H2=zeros(numEq2,numLocalDoF2);
        if additionalOutput==1
            args_el=cell(1,numAO);
            for i_=1:numAO
                args_el{i_}=zeros(numGP*sizeAO(i_,1),sizeAO(i_,2));
            end
        end
        for l=1:numel(gaussWeights) %loop over Gauss points
            [shapeFunctionsVals, naturalDerivatives]=shapeFunctions2d(gaussPts(l,1),gaussPts(l,2));
            J=elementNodeCoords'*naturalDerivatives'; %Jacobian from natural to spatial coordinates
            spatialDerivatives=naturalDerivatives'/J; %same as naturalDerivatives'*inv(J)
            if numDim==2
                for ieq=1:numEq1
                    gradH1((ieq-1)*numDim+1,ieq:numEq1:end)=spatialDerivatives(:,1)';
                    gradH1((ieq-1)*numDim+2,ieq:numEq1:end)=spatialDerivatives(:,2)';
                end
                for ieq=1:numEq2
                    gradH2((ieq-1)*numDim+1,ieq:numEq2:end)=spatialDerivatives(:,1)';
                    gradH2((ieq-1)*numDim+2,ieq:numEq2:end)=spatialDerivatives(:,2)';
                end
            elseif numDim==3
                for ieq=1:numEq1
                    gradH1((ieq-1)*numDim+1,ieq:numEq1:end)=spatialDerivatives(:,1)';
                    gradH1((ieq-1)*numDim+2,ieq:numEq1:end)=spatialDerivatives(:,2)';
                    gradH1((ieq-1)*numDim+2,ieq:numEq1:end)=spatialDerivatives(:,2)';
                end
                for ieq=1:numEq2
                    gradH2((ieq-1)*numDim+1,ieq:numEq2:end)=spatialDerivatives(:,1)';
                    gradH2((ieq-1)*numDim+2,ieq:numEq2:end)=spatialDerivatives(:,2)';
                    gradH2((ieq-1)*numDim+2,ieq:numEq2:end)=spatialDerivatives(:,2)';
                end
            else
                error('2D or 3D problem expected');
            end
            for ieq=1:numEq1
                H1(ieq,ieq:numEq1:end)=shapeFunctionsVals;
            end
            for ieq=1:numEq2
                H2(ieq,ieq:numEq2:end)=shapeFunctionsVals;
            end
            
            elementGPCoords=shapeFunctionsVals*elementNodeCoords;
            elementGPVals=shapeFunctionsVals*elementNodeVals;
            graduVals=spatialDerivatives'*elementNodeVals;
            CMatVals=CMat(elementGPCoords,elementGPVals,graduVals,elementStVars(l,:));
            BetaMatVals=BetaMat(elementGPCoords,elementGPVals,graduVals,elementStVars(l,:));
            
            u_loc=elementNodeVals';
            u_loc=u_loc(:);
            FLocal=FLocal+...
                gradH2'*(CMatVals*gradH1 + BetaMatVals*H1)*u_loc...
                *sqrt(det(J'*J))*gaussWeights(l);
        end %loop over Gauss points
    end %formLocalStiffness

%% ==== Assembly of global stiffness matrix ====
[gaussPts, gaussWeights]=gaussPoints(elementType,numGP);
F=zeros(numGDoF2,1); % initialize global stiffness matrix

if ~isnonconstant %linear problem, constant CMatrix
    for m=1:numEl %loop over elements
        elementNodesIDs=IEN(m,:); %global IDs for current element's nodes
        elementNodeCoords=nodeCoords(elementNodesIDs,:); %current element's nodes
        elementNodeVals=u(elementNodesIDs,:); %displacement values at current el's nodes
        Flocal = formLocalForce(elementNodeCoords, shapeFunctions,...
            gaussPts, gaussWeights,CMatrix,BetaMatrix, elementNodeVals);
        elementDoFIDs2=node2DoFs(elementNodesIDs',numEq2);
        F(elementDoFIDs2)=F(elementDoFIDs2)+Flocal;
    end 
else %nonlinear problem or varying CMatrix
    for m=1:numEl %loop over elements
        elementNodesIDs=IEN(m,:); %global IDs for current element's nodes
        elementNodeCoords=nodeCoords(elementNodesIDs,:); %current element's nodes
        elementNodeVals=u(elementNodesIDs,:); %displacement values at current el's nodes
        elementStVars=s( (numGP*(m-1)+1):(numGP*m) , :); %state variables at current el's Gauss points
        if additionalOutput==1
            [Flocal, args_el] = formLocalForceVar(elementNodeCoords, shapeFunctions,...
                gaussPts, gaussWeights,CMatrix,BetaMatrix, elementNodeVals, elementStVars);
            for i=1:numAO
                argout{i}( (numGP*sizeAO(i,1)*(m-1)+1) : numGP*sizeAO(i,1)*m , : )=...
                    args_el{i};
            end
        else
            Flocal = formLocalForceVar(elementNodeCoords, shapeFunctions,...
                gaussPts, gaussWeights,CMatrix,BetaMatrix, elementNodeVals, elementStVars);
        end
        
        elementDoFIDs2=node2DoFs(elementNodesIDs',numEq2);
        F(elementDoFIDs2)=F(elementDoFIDs2)+Flocal;
    end %loop over elements
end
   

end %function