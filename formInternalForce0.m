function F = formInternalForce0(nodeCoords, IEN, elementType, s)
% Compute equivalent internal nodal force vector from stress values
% given at Gauss points. 
%
% INPUT:
% nodeCoords - #Nodes-by-#Dim matrix, coordinates of nodes.
% IEN - #El-by-#NpE array, element-node incidence matrix. 
% elementType - single string, the type of elements.
% s - #El*#GP-by-#Eq*#Dim matrix, stress values at Gauss points. Each row
%   s((m-1)*#GP+l,:) corresponds to the stress tensor of at l-th Gauss point
%   of m-th element.
% 
% OUTPUT:
% F - #Nodes*#Eq-by-1 force vector
%
% Examples:
%   F = formInternalForce0(nodeCoords, IEN, elementType, s)

numEl=size(IEN,1); %total number of elements in the mesh
numNodes=size(nodeCoords,1); %total number of nodes in the mesh 
numDim=size(nodeCoords,2); %the dimensionality of the underlying space
numEq=size(s,2)/numDim;
numGDoF=numNodes*numEq; %total number of degrees of freedom
numGP=size(s,1)/numEl;
%% ==== Parsing and validating input ==== %
switch numDim % check geometric dimensionality
    case {2,3}
    otherwise
        error('This function only works for 2D or 3D problems')
end

if exist(['shapeFunctions' elementType],'file')==2
    shapeFunctions=str2func(['shapeFunctions' elementType]);
    % Assign a function handle to shapeFunctions
else
    error(['Can''t find shape functions for specified element type: "' ...
        elementType '"' ]);
end

%% ==== The following subroutines assemble element local force vector ==
    function FLocal = formLocalForce(elementNodeCoords, shapeFunctions,...
            gaussPts, gaussWeights, elementStressGP)
    % forms a stiffness matrix for a given element
        numNpE=size(elementNodeCoords,1); %number of nodes in this element
        numLocalDoF=numNpE*numEq; %number of element degrees of freedom (input)
        FLocal=zeros(numLocalDoF,1);
        gradH=zeros(numDim*numEq,numLocalDoF);
        for l=1:numel(gaussWeights) %loop over Gauss points
            [~, naturalDerivatives]=shapeFunctions(gaussPts(l,:));
            J=elementNodeCoords'*naturalDerivatives'; %Jacobian from natural to spatial coordinates
            spatialDerivatives=naturalDerivatives'/J; %same as naturalDerivatives'*inv(J)
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

            FLocal=FLocal+gradH'*elementStressGP(l,:)'*sqrt(det(J'*J))*gaussWeights(l);
                
        end %loop over Gauss points
    end %formLocalStiffness

%% ==== Assembly of global force vector ====
[gaussPts, gaussWeights]=gaussPoints(elementType,numGP);
F=zeros(numGDoF,1); % initialize global stiffness matrix

for m=1:numEl %loop over elements
    elementNodesIDs=IEN(m,:); %global IDs for current element's nodes
    elementNodeCoords=nodeCoords(elementNodesIDs,:); %current element's nodes
    elementGPStress=s(numGP*(m-1)+1:numGP*m,:);
    Flocal = formLocalForce(elementNodeCoords, shapeFunctions,...
        gaussPts, gaussWeights, elementGPStress);
    elementDoFIDs=node2DoFs(elementNodesIDs',numEq);
    F(elementDoFIDs)=F(elementDoFIDs)+Flocal;
end %loop over elements

end %function