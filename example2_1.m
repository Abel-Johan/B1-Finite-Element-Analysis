% Steady state mechanical equivilibrium and steady state heat equation 
% on a rectangular domain (linear plane stress problem).
%
% Elastic and thermal parts are united in a single elliptic equation (N=3)
% and are solved at once. Note that the resulting stiffness matrix is not
% symmetric.

%% ===== Problem Parameters ===== %
domain=[0 0 2 2];
% Material properties
youngsModulus = 193e6;
poissonsRatio = 0.253;
properties=elasticProperties('youngsModulus',youngsModulus,'poissonsRatio',poissonsRatio);
CMatrix_E=properties.CPlaneStressFull;
% Heat conductivity tensor 
CondCoeff=16.2;
CMatrix_T=CondCoeff*[1 0; 0 1];
% Combined CMatrix
CMatrix=blkdiag(CMatrix_E,CMatrix_T);
% Thermal expansion
TEcoeff=9.7e-6;
TEtensor=TEcoeff*[1 0; 0 1]; %thermal expansion
TEtensortr=TEtensor';
BetaMatrix=CMatrix*[zeros(6,2) [-TEtensortr(:); zeros(2,1)]]; % u-->gradu map

%% ===== FEM Solver Parameters ===== %
elementType='2dQ2';
elementType1d='1dQ2';
numGP=9;
numGP1d=2;
elementSize=.4; %element diameter
% misc
numDim=2; % A 2D problem
numEq=3; % Number of scalar equations: 2 for elastic part, 1 for thermal

%% ===== Specifying Body Force and Boundary conditions =====
bodyForce2d=@(x,y)([0*x 0*x 0*x]'); % distributed body force and heat source
bodyForce2d=@(x)(bodyForce2d(x(:,1),x(:,1)));

boundaryTraction2d=@(x,y)([0*x 0*x 300*y]'); %boundary traction and heat flux
boundaryTraction2d=@(x)(boundaryTraction2d(x(:,1),x(:,2)));

boundaryDisplacement2d=@(x,y)[0*x 0*x 000*(x<=1e-6)]'; %boundary displacement and temperature
boundaryDisplacement2d=@(x)(boundaryDisplacement2d(x(:,1),x(:,2)));
isDirichlet=[0 0 0; 0 0 0;0 0 0;1 1 1];

% Mesh generation
[nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect2d(domain, elementType, elementSize);
BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);

% misc
numEl=size(IEN,1); %total number of elements in the mesh
numNodes=size(nodeCoords,1); %total number of nodes in the mesh 
numGDoF=numNodes*numEq; %total number of displacement/temp degrees of freedom

% Drawing the original mesh
hold off;
drawElements(nodeCoords,IEN,elementType)
hold on;
view(2);

% Matrix, vector assembly
K = formStiffnessMatrix2(nodeCoords, IEN, elementType, numGP, CMatrix, BetaMatrix);
Fb = formBodyForceVector(nodeCoords, IEN, elementType, numGP, bodyForce2d);

%% Applying Boundary Conditions
[u_prescribed, Fs, prescribedDoF, freeDoF]=...
    formBC(nodeCoords,BIEN,elementType1d,numGP1d,...
    boundaryTraction2d,boundaryDisplacement2d,isDirichlet);

%%
F=Fs+Fb;
KK=K(freeDoF,freeDoF);
if ~isempty(prescribedDoF)
    FF=F(freeDoF)-K(freeDoF,prescribedDoF)*u_prescribed(prescribedDoF);
else
    FF=F(freeDoF);
end

% Solving
u=zeros(numGDoF,1);
u(freeDoF)=KK\FF;
u(prescribedDoF)=u_prescribed(prescribedDoF);

%% Recovery and visualisation
u2=reshape(u,[numEq numel(u)/numEq])';
u2_E=u2(:,[1 2]); % nodal mechanical displacements
u2_T=u2(:,3); % nodal temperature
[gradu]=recoveryPatch(u2_E(:,[1 2]),nodeCoords,IEN,elementType,BIEN);
s2=CMatrix_E*gradu;

% Drawing the results
hold on;
factor=1e2;
figure(1); clf;
subplot(2,1,1);
drawElementsInterp(nodeCoords+u2_E*factor,IEN,elementType,u2_T);
title('T, Temperature distribution');
view(2)
subplot(2,1,2);
drawElementsInterp(nodeCoords+u2_E*factor,IEN,elementType,s2(1,:)+s2(1,:));
title('I_1=\sigma_{xx} + \sigma_{yy} stress invariant');
view(2)