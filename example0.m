% Minimal example - plane stress equilibrium problem
%
% We consider a plane stress problem for a rectangular domain.
% Zero-displacement boundary condition is applied to the left boundary,
% traction is applied to the right boundary, top and bottom boundaries are
% load-free. 

clear; close all;
%% Mesh generation
elementType='2dQ1'; %define element type
elementType1d='1dQ1'; %define boundary element type - must be consistent
[nodeCoords, IEN, boundaryElementIDs, boundaryNodeLocalID]=...
    meshRect2d([0 0 2 2],elementType,[10 10]); %creating mesh 10-by-10 el.
%create boundary incidence arrays
BIEN=IENtoBIEN(IEN, boundaryElementIDs, boundaryNodeLocalID);

%% Matrix and vector assembly 
% Boundary conditions
%define function handles for BC
bndDisplacement=@(x)(zeros(size(x))'); %zero displacement
bndTraction0=@(x)(zeros(size(x))'); %zero traction
bndTraction1=@(x)(repmat([1; 1],size(x,1)));  %non-zero traction
%make cell arrays of function handles - boundary condition cell arrays must
%match the structure of BIEN 
bndTractions={bndTraction0,bndTraction1,bndTraction0,bndTraction0};
bndDisplacements=repmat({bndDisplacement},4,1); 
isDirichlet=[0; 0; 0; 1]; %define which boundary regions have Dirichlet BC

numGP1d=1; %number of Gauss Points for 1d boundary elements
% assemble boundary load vector for Neumann BC
% and evaluate displacements for Dirichlet BC
[u_prescribed, Fs, prescribedDoF, freeDoF]=...
    formBC(nodeCoords,BIEN,elementType1d,numGP1d,...
    bndTractions,bndDisplacements,isDirichlet);

% Body force
numGP=1; %number of Gauss Points for 2d elements
%define function handle for body force
bodyForce=@(x)(repmat([0; -1],size(x,1))); 
% assemble body force vector
Fb = formBodyForceVector(nodeCoords, IEN, elementType, numGP, bodyForce);

% Stiffness matrix
%define the matrix of elasticities
CMatrix=elasticProperties('youngsModulus',193e6,'poissonsRatio',0.253,...
    'CPlaneStressEng'); %elasticity tensor (Voigt notation)
%assemble stiffness matrix
K = formStiffnessMatrixEng(nodeCoords, IEN, elementType, numGP, CMatrix);

%% Solution 
F=Fb+Fs; %total load vector
%define equivalent load vector
FF=F(freeDoF)-K(freeDoF,prescribedDoF)*u_prescribed(prescribedDoF);
%define equivalent stiffness matrix
KK=K(freeDoF,freeDoF);
u=zeros(size(u_prescribed)); %initialise vector of displacements
%solve linear equations
u(freeDoF)=KK\FF;

%% Stress recovery 
u2=reshape(u,[2 numel(u)/2])'; 
%recover strains at element centroids
[strain, GPCoords]=recoveryGPEng(u2,nodeCoords,IEN,elementType,1);
%evaluate stresses at centroids
s2=CMatrix*strain;

%% Visualisation 
figure(1);clf;
drawElements(nodeCoords,IEN,elementType); %draw initial undeformed mesh 
hold on;
drawNodes(nodeCoords); %draw nodes
drawNodes(GPCoords,'x'); %draw Gauss points

figure(2);clf;
drawElements(nodeCoords,IEN,elementType); %draw initial undeformed mesh 
hold on;
factor=1e7; %scaling factor to make small deformations visible
%draw deformed mesh on top with alpha=.3
drawElements(nodeCoords+u2*factor,IEN,elementType,s2(1,:)',.3);
drawNodes(nodeCoords,BIEN{4},{'ks','filled'}); %draw pinned nodes
title('\sigma_x');