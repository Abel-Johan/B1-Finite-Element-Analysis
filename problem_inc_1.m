%%%%% Original file for Problem 1 (incremental)
clear; close all;
%% Mesh generation
elementType='2dQ1'; %define element type
domain=[0 0 4 1]; %[x0 y0 x1 y1]
numElements=[15 10]; %[numElX numElY]
[nodeCoords, IEN, boundaryElementIDs, boundaryNodeLocalID]=...
    meshRect2d(domain,elementType,numElements); %generate the mesh

%% Parameters
%elasticity tensor
CMatrix=elasticProperties('youngsModulus',193e6,'poissonsRatio',0.253,'CPlaneStressEng');

%% Matrix and vector assembly 
numNodes=size(nodeCoords,1); %total number of nodes
numDoFs=numNodes*2; %number of global degrees of freedom
% Assemble stiffness matrix
numGP=4; %number of Gauss points used in quadrature
K = formStiffnessMatrixEng(nodeCoords, IEN, elementType, numGP, CMatrix);
% Body force
Fb = zeros(numDoFs,1);

% Boundary conditions
leftNodes=find(nodeCoords(:,1)==domain(1)); 
rightNodes=find(nodeCoords(:,1)==domain(3));
leftXDoF=(leftNodes-1)*2+1;
leftYDoF=(leftNodes-1)*2+2;
rightXDoF=(rightNodes-1)*2+1;
rightYDoF=(rightNodes-1)*2+2;

prescribedDoF=[leftXDoF; leftYDoF]; %list of prescribed global DoFs
freeDoF=setdiff(1:numDoFs,prescribedDoF); %list of free DoFs
u_prescribed=zeros(numDoFs,1); 
% Surface traction
Fs = zeros(numDoFs,1); %define global force vector for surface traction

%% Solution 
F=Fb+Fs; %total load vector
%define the free part of load vector
FF=F(freeDoF)-K(freeDoF,prescribedDoF)*u_prescribed(prescribedDoF);
%define the free part of stiffness matrix
KK=K(freeDoF,freeDoF);
%solve linear equations
u=zeros(numDoFs,1);
u(freeDoF)=KK\FF;
u(prescribedDoF)=u_prescribed(prescribedDoF);

%% Stress recovery 
u2=reshape(u,[2 numel(u)/2])'; %reshape s.t. Ux=u2(:,1), Uy=u2(:,2)
%recover strains at centroids elements
[strain, GPCoords]=recoveryGPEng(u2,nodeCoords,IEN,elementType,1);
%evaluate stresses at the centroids
s2=CMatrix*strain;

%% Visualisation 
figure(1);clf; 
%draw initial undeformed mesh
drawElements(nodeCoords,IEN,elementType,0,0);
hold on;
factor=1e7; %scaling factor to amplify small deformations
%draw deformed mesh
drawElements(nodeCoords+u2*factor,IEN,elementType,s2(1,:)',.7);
title('\sigma_x')