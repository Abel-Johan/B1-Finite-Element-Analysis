% Steady state mechanical equilibrium and steady-state heat equation 
% on a rectangular domain (linear plane stress problem).
%
% Elastic and thermal parts are represented by two equations, which are
% solved separately. 

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
% Thermal expansion
TEcoeff=9.7e-6;
TEtensor=TEcoeff*[1 0; 0 1]; %thermal expansion
TEtensortr=TEtensor';
BetaMatrix=CMatrix_E*TEtensortr(:); % th-->gradu map

%% ===== FEM Solver Parameters ===== %
elementType='2dQ1';
elementType1d='1dQ1';
numGP=4;
numGP1d=2;
elementSize=.4; %element diameter
% misc
numDim=2; % A 2D problem
numEq_E=2; % Number of scalar equations: 2 for elastic part, 1 for thermal
numEq_T=1; % Number of scalar equations: 2 for elastic part, 1 for thermal

%% ===== Specifying Body Force and Boundary conditions =====
bodyForce2d_E=@(x,y)([0*x 0*x]'); % distributed body force
bodyForce2d_E=@(x)(bodyForce2d_E(x(:,1),x(:,1)));
bodyForce2d_T=@(x,y)([0*x]'); % distributed heat source
bodyForce2d_T=@(x)(bodyForce2d_T(x(:,1),x(:,1)));

boundaryTraction2d_E=@(x,y)([0*x 0*x]'); %boundary traction 
boundaryTraction2d_E=@(x)(boundaryTraction2d_E(x(:,1),x(:,2)));
boundaryTraction2d_T=@(x,y)([300*y]'); %boundary heat flux
boundaryTraction2d_T=@(x)(boundaryTraction2d_T(x(:,1),x(:,2)));

boundaryDisplacement2d_E=@(x,y)[0*x 0*x]'; %boundary displacement and temperature
boundaryDisplacement2d_E=@(x)(boundaryDisplacement2d_E(x(:,1),x(:,2)));
boundaryDisplacement2d_T=@(x,y)[000*(x<=1e-6)]'; %boundary displacement and temperature
boundaryDisplacement2d_T=@(x)(boundaryDisplacement2d_T(x(:,1),x(:,2)));
isDirichlet_E=[0;0;0;1];
isDirichlet_T=[0;0;0;1];

%% Mesh generation
[nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect2d(domain, elementType, elementSize);
BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);

% misc
numEl=size(IEN,1); %total number of elements in the mesh
numNodes=size(nodeCoords,1); %total number of nodes in the mesh 
numGDoF_E=numNodes*numEq_E; %total number of displacement/temp degrees of freedom
numGDoF_T=numNodes*numEq_T; %total number of displacement/temp degrees of freedom

% Drawing the original mesh
hold off;
drawElements(nodeCoords,IEN,elementType)
hold on;
view(2);

%% Matrix, vector assembly
K_E = formStiffnessMatrix(nodeCoords, IEN, elementType, numGP, CMatrix_E);
Fb_E = formBodyForceVector(nodeCoords, IEN, elementType, numGP, bodyForce2d_E);

K_T = formStiffnessMatrix(nodeCoords, IEN, elementType, numGP, CMatrix_T);
Fb_T = formBodyForceVector(nodeCoords, IEN, elementType, numGP, bodyForce2d_T);

%% Applying Boundary Conditions
[u_prescribed_E, Fs_E, prescribedDoF_E, freeDoF_E]=...
    formBC(nodeCoords,BIEN,elementType1d,numGP1d,...
    boundaryTraction2d_E,boundaryDisplacement2d_E,isDirichlet_E);
        
[u_prescribed_T, Fs_T, prescribedDoF_T, freeDoF_T]=...
    formBC(nodeCoords,BIEN,elementType1d,numGP1d,...
    boundaryTraction2d_T,boundaryDisplacement2d_T,isDirichlet_T);
        
%% Solution
% Thermal part
F_T=Fs_T+Fb_T;
KK_T=K_T(freeDoF_T,freeDoF_T);
if ~isempty(prescribedDoF_T)
    FF_T=F_T(freeDoF_T)-K_T(freeDoF_T,prescribedDoF_T)*u_prescribed_T(prescribedDoF_T);
else
    FF_T=F_T(freeDoF_T);
end

u_T=zeros(numGDoF_T,1);
u_T(freeDoF_T)=KK_T\FF_T;
u_T(prescribedDoF_T)=u_prescribed_T(prescribedDoF_T);

% Elastic part
F_TE = formInternalForce2(nodeCoords, IEN, elementType, numGP, zeros(4,2), BetaMatrix, u_T);

F_E=Fs_E+Fb_E+F_TE;
KK_E=K_E(freeDoF_E,freeDoF_E);
if ~isempty(prescribedDoF_E)
    FF_E=F_E(freeDoF_E)-K_E(freeDoF_E,prescribedDoF_E)*u_prescribed_E(prescribedDoF_E);
else
    FF_E=F_E(freeDoF_E);
end

u_E=zeros(numGDoF_E,1);
u_E(freeDoF_E)=KK_E\FF_E;
u_E(prescribedDoF_E)=u_prescribed_E(prescribedDoF_E);

%% Recovery and visualisation
u2_E=reshape(u_E,[numEq_E numel(u_E)/numEq_E])';
[gradu]=recoveryPatch(u2_E,nodeCoords,IEN,elementType,BIEN);
s2=CMatrix_E*gradu;

% Drawing the results
hold on;
factor=1e2;
figure(1); clf;
subplot(2,1,1);
drawElementsInterp(nodeCoords+u2_E*factor,IEN,elementType,u_T);
title('T, Temperature distribution');
view(2)
subplot(2,1,2);
drawElementsInterp(nodeCoords+u2_E*factor,IEN,elementType,s2(1,:)+s2(1,:));
title('I_1=\sigma_{xx} + \sigma_{yy} stress invariant');
view(2)