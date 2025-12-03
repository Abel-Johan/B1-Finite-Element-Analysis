function F = formBodyForceVector(nodeCoords, IEN, elementType, numGP, bodyForce)
% Compute body force vector.
% INPUT:
% nodeCoords - #Nodes-by-#Dim matrix, coordinates of nodes.
% IEN - #El-by-#NpE array, element-node incidence matrix. 
% elementType - single string, the type of elements.
% numGP - single scalar, number of Gauss points per element.
% bodyForce - #Eq-by-1 matrix or a function handle returning such matrix.
%   The function input is 1-by-#Dim position vector.
%
% OUTPUT:
% F - (#nodes*#Eq)-by-1 body force vector
%
% Examples:
%   F = formBodyForceVector(nodeCoords, IEN, elementType, numGP, [0 -9.8*8.05])
%   F = formBodyForceVector(nodeCoords, IEN, elementType, numGP, bodyForce)

numEl=size(IEN,1); %total number of elements in the mesh
numNodes=size(nodeCoords,1); %total number of nodes in the mesh 
numEq=numel(bodyForce(nodeCoords(1,:))); % Number of scalar equations
numDoF=numNodes*numEq; %total number of degrees of freedom

%% ==== Parsing and validating input ==== %
if exist(['shapeFunctions' elementType],'file')==2
    shapeFunctions=str2func(['shapeFunctions' elementType]);
    % Assign a function handle to shapeFunctions
else
    error(['Can''t find shape functions for specified element type: "' ...
        elementType '"' ]);
end

%% ==== The following subroutine assembles element local body force vector==
    function FbLocal = formLocalBodyForce(elementNodeCoords, shapeFunctions,...
            gaussPts, gaussWeights, bodyForce)
        % forms a body force vector a given element
        numNpE=size(elementNodeCoords,1); %number of nodes in this element
        numLocalDoF=numNpE*numEq; %number of element degrees of freedom
        FbLocal=zeros(numLocalDoF,1);
        H=zeros(numEq,numLocalDoF);
        for l=1:numel(gaussWeights) %loop over Gauss points
            [shapeFunctionsVals, naturalDerivatives]=shapeFunctions(gaussPts(l,:));
            bodyForceVals=bodyForce(shapeFunctionsVals*elementNodeCoords);
            % body force evaluated at the Gauss point in global coordinates
            J=elementNodeCoords'*naturalDerivatives'; %Jacobian from natural to spatial coordinates
            for ieq=1:numEq
                H(ieq,ieq:numEq:end)=shapeFunctionsVals;
            end
            FbLocal=FbLocal+H'*bodyForceVals*sqrt(det(J'*J))*gaussWeights(l);
        end %loop over Gauss points
    end %formLocalBodyForce

%% ==== Assembly of global body force vector ====
[gaussPts, gaussWeights]=gaussPoints(elementType,numGP);
F=zeros(numDoF,1); %initialize global force vector
for m=1:numEl %loop over elements
    elementNodesIDs=IEN(m,:); %global IDs of current element's nodes
    elementNodeCoords=nodeCoords(elementNodesIDs,:); %current element's nodes
    Flocal = formLocalBodyForce(elementNodeCoords, shapeFunctions,...
        gaussPts, gaussWeights, bodyForce);
    elementDoFIDs=node2DoFs(elementNodesIDs',numEq);
    F(elementDoFIDs)=F(elementDoFIDs)+Flocal;
end %loop over elements
end %function