function [d1, d2]=example6(y0,Tr)
% Run example6_main.m
% Deflection of a bimetal
% Steady state mechanical equilibrium and equilibrium heat equation on a
% rectangular domain (linear plane stress problem). The only aspect of
% thermomechanical coupling considered is heat expansion. 
%
% INPUT:
% y0 - a scalar 0<=y0<=2 defining the interface location between steel and brass.
% Tr - scalar, heat flux on the right edge of the domain.
%
% OUTPUT
% d1 - displacement of upper-right corner
% d2 - displacement of lower-right corner

if (nargin==0)&&(nargout==0)
    example6_main;
    return;
end

doplot=1;
factor=1e2;
%% ===== Problem Parameters ===== %
if nargin==0
    y0=1;
end
%
domain=[0 0 10 2];
%% ==== Material properties ==== 
youngsModulus1 = 193e6;
youngsModulus2 = 112e6;
poissonsRatio1 = 0.253;
poissonsRatio2 = 0.331;
properties1=elasticProperties('youngsModulus',youngsModulus1,'poissonsRatio',poissonsRatio1);
properties2=elasticProperties('youngsModulus',youngsModulus2,'poissonsRatio',poissonsRatio2);
CMatrix_E1=properties1.CPlaneStressFull;
CMatrix_E2=properties2.CPlaneStressFull;
    % this function defines position-dependent elasticity modulus
    function CMat=CMatrix_E(x,~,~,~)
        
        if nargin==0
            CMat=CMatrix_E1;
        else
            CMat=((x(2)>=y0)*CMatrix_E1+(x(2)<y0)*CMatrix_E2);
        end
    end

% Heat conductivity tensor
CondCoeff1=16.2;
CondCoeff2=121;
CMatrix_T1=CondCoeff1*[1 0; 0 1];
CMatrix_T2=CondCoeff2*[1 0; 0 1];
    % this function defines position-dependent heat conductivity tensor
    function CMat=CMatrix_T(x,~,~)
        if nargin==0
            CMat=CMatrix_T1;
        else
            CMat=((x(2)>=y0)*CMatrix_T1+(x(2)<y0)*CMatrix_T2);
        end
    end

% Thermal expansion
TEcoeff1=9.7e-6;
TEcoeff2=19.9e-6;
TEtensor1=TEcoeff1*[1 0; 0 1]; %thermal expansion
TEtensor2=TEcoeff2*[1 0; 0 1]; %thermal expansion
% this function defines position-dependent thermal expansion tensor
    function TEten=TEtensortr(x)
        if nargin==0
            TEten=TEtensor1;
        else
            TEten=((x(2)>=y0)*TEtensor1+(x(2)<y0)*TEtensor2);
        end
        TEten=TEten(:);
    end
BetaMatrix=@(x,u,du,s)-CMatrix_E(x,u,du,s)*TEtensortr(x); % th-->gradu map


%% ===== FEM Solver Parameters ===== %
elementType='2dQ1';
elementType1d='1dQ1';
numGP=4;
numGP1d=2;
elementSize=[21 21]; %number of elements
% elementSize=[15 15]; %number of elements
% misc
numEq_E=2; % Number of scalar equations: 2 for elastic part, 1 for thermal
numEq_T=1; % Number of scalar equations: 2 for elastic part, 1 for thermal

%% ===== Specifying Body Force and Boundary conditions =====
bodyForce2d_E=@(x,y)([0*x 0*x]'); % distributed body force
bodyForce2d_E=@(x)(bodyForce2d_E(x(:,1),x(:,1)));
bodyForce2d_T=@(x,y)([0*x]'); % distributed heat source
bodyForce2d_T=@(x)(bodyForce2d_T(x(:,1),x(:,1)));

boundaryTraction2d_E=@(x,y)([0*x 0*x]'); %boundary traction
boundaryTraction2d_E=@(x)(boundaryTraction2d_E(x(:,1),x(:,2)));
boundaryTraction2d_T=@(x,y)([00*y]'); %boundary heat flux
boundaryTraction2d_T=@(x)(boundaryTraction2d_T(x(:,1),x(:,2)));

boundaryDisplacement2d_E=@(x,y)[0*x 0*x]'; %boundary displacement
boundaryDisplacement2d_E=@(x)(boundaryDisplacement2d_E(x(:,1),x(:,2)));
boundaryDisplacement2d_T=@(x,y)[00*(x<=1e-6+domain(1))+Tr*(x>=-1e-6+domain(3))]'; %boundary temperature
boundaryDisplacement2d_T=@(x)(boundaryDisplacement2d_T(x(:,1),x(:,2)));
isDirichlet_E=[0 0; 0 0; 0 0; 1 1];
pinnedDoF=2; %y-DoF of the first node.
isDirichlet_T=[0;1;0;1];

%% Mesh generation
[nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect2d(domain, elementType, elementSize);
boundaryIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);

% misc
numEl=size(IEN,1); %#ok<NASGU> %total number of elements in the mesh
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
    formBC(nodeCoords,boundaryIEN,elementType1d,numGP1d,...
    boundaryTraction2d_E,boundaryDisplacement2d_E,isDirichlet_E,pinnedDoF);


[u_prescribed_T, Fs_T, prescribedDoF_T, freeDoF_T]=...
    formBC(nodeCoords,boundaryIEN,elementType1d,numGP1d,...
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
Ft_E = formInternalForce2(nodeCoords, IEN, elementType, numGP, zeros(4,2), BetaMatrix, u_T);

F_E=Fs_E+Fb_E-Ft_E;
KK_E=K_E(freeDoF_E,freeDoF_E);
if ~isempty(prescribedDoF_E)
    FF_E=F_E(freeDoF_E)-K_E(freeDoF_E,prescribedDoF_E)*u_prescribed_E(prescribedDoF_E);
else
    FF_E=F_E(freeDoF_E);
end

u_E=zeros(numGDoF_E,1);
u_E(freeDoF_E)=KK_E\FF_E;
u_E(prescribedDoF_E)=u_prescribed_E(prescribedDoF_E);

u2_E=reshape(u_E,[numEq_E numel(u_E)/numEq_E])';
[gradu_E_GP,~]=recoveryGP(u2_E,nodeCoords,IEN,elementType,1);

% Plotting the results
if doplot
    figure(1);clf;
    subplot(2,1,1);
    drawElementsInterp(nodeCoords,IEN,elementType,u_T);
    title({'T, temperature distribution','undeformed mesh','dotted - interface bw 2 materials'});
    line(domain([1,3]),[y0 y0],'Linestyle','--','linewidth',2)
    subplot(2,1,2);
    drawElements(nodeCoords+u2_E*factor,IEN,elementType,gradu_E_GP(2,:)+gradu_E_GP(3,:));
    title({'u_y+v_x, engineering shear strain', ['displacement factor=' num2str(factor)], ['y0=' num2str(y0)]});
    view(2);
    drawnow;
end

% upper-right and lower-right corners' node IDs
uc=find((nodeCoords(:,1)==domain(3))&(nodeCoords(:,2)==domain(4)));
lc=find((nodeCoords(:,1)==domain(3))&(nodeCoords(:,2)==domain(2)));

d1=u2_E(uc,:); %#ok<FNDSB>
d2=u2_E(lc,:); %#ok<FNDSB>
end%file