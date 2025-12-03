% Steady state heat diffusion equation on a rectangular domain
% BC: Dirichlet on the left side (T=100), Neumann elsewhere (Flux=100*y on
% the right side, Flux=0 on the top and bottom sides).
% Distributed heat source q=5000 in rectangular region 0.5<=x<=0.7, y<=0.2.

%% ===== Problem Parameters ===== %
domain=[0 0 2 2];
% Heat diffusion tensor 
CMatrix=[1 0; 0 1];

%% ===== FEM Solver Parameters ===== %
elementType='2dQ1';
numGP=4; % to suit both linear and quadratic elements
elementType1d='1dQ1';
numGP1d=2; % to suit both linear and quadratic elements
elementSize=.4; %element diameter
% misc
numDim=2; % A 2D problem
numEq=1; % Number of scalar equations 

%% ===== Specifying Body Force and Boundary conditions =====
bodyForce2d=@(x,y)([5e3*((.5<=x)&(x<=.7)&(y<=.2))]'); %distributed heat source
bodyForce2d=@(x)(bodyForce2d(x(:,1),x(:,2)));

boundaryTraction2d=@(x,y)([(100*y)*(x>=2-1e-6)]'); %boundary flux
boundaryTraction2d=@(x)(boundaryTraction2d(x(:,1),x(:,2)));

boundaryDisplacement2d=@(x,y)[100*(x<=1e-6)]'; %boundary temperature
boundaryDisplacement2d=@(x)(boundaryDisplacement2d(x(:,1),x(:,2)));
isDirichlet=[0;0;0;1];

%% Mesh generation
[nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect2d(domain, elementType, elementSize);
boundaryIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);

% misc
numEl=size(IEN,1); %total number of elements in the mesh
numNodes=size(nodeCoords,1); %total number of nodes in the mesh 
numGDoF=numNodes*numEq; %total number of displacement degrees of freedom

% Drawing the original mesh
hold off;
drawElements(nodeCoords,IEN,elementType)
hold on;
view(2);

% Matrix, vector assembly
K = formStiffnessMatrix(nodeCoords, IEN, elementType, numGP, CMatrix);
Fb = formBodyForceVector(nodeCoords, IEN, elementType, numGP, bodyForce2d);

%% Applying Boundary Conditions
[u_prescribed, Fs, prescribedDoF, freeDoF]=...
    formBC(nodeCoords,boundaryIEN,elementType1d,numGP1d,...
    boundaryTraction2d,boundaryDisplacement2d,isDirichlet);

%% Solution
F=Fs+Fb;
KK=K(freeDoF,freeDoF);
if ~isempty(prescribedDoF)
    FF=F(freeDoF)-K(freeDoF,prescribedDoF)*u_prescribed(prescribedDoF);
else
    FF=F(freeDoF);
end

u=zeros(numGDoF,1);
u(freeDoF)=KK\FF;
u(prescribedDoF)=u_prescribed(prescribedDoF);

%% Visualisation
hold on;
factor=1e1;
figure(1); clf;
drawElementsInterp(nodeCoords,IEN,elementType,u);
title('T, Temperature distribution');
view(2)