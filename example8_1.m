function [ts,u,v,a]=example8_1(elementSize,tfinal)
% Oscilation of cantilever beam - implicit Newmark's method
% Run example8_main
% Plane stress, rectangular domain
% We consider a cantilever beam which is being pulled and then released at
% the free end. 

if (nargin==0)&&(nargout==0)
    example8_main;
    return;
end
doplot=1;
displayevery=2; % Every k-th timestep to plot
factor=1e5;
%% ===== Problem Parameters ===== %
domain=[0 0; 10 3]; %rectangular domain
thickness=1; %thickness of a beam

% Material properties
youngsModulus = 193e6;
poissonsRatio = 0.253;
massDensity=7.999e3;
properties=elasticProperties('youngsModulus',youngsModulus,'poissonsRatio',poissonsRatio);
CMatrix=properties.CPlaneStressFull;
MMatrix=[1 0; 0 1]*massDensity*thickness;
numEq=2; % Number of scalar equations (elasticity)

%% ===== FEM Solver Parameters ===== %
elementType='2dP1';
gamma=1/2; % Newmark parameter
beta=1/4; % Newmark parameter
dt = 3e-02; % Fixed timestep
numGP=1;
elData=elementData(elementType);
bndElementType=elData.bndElementType;
bndElData=elementData(bndElementType);
numGPbnd=bndElData.numGPFull;

%% ===== BOUNDARY CONDITIONS DEFINITION =====
traction1=[0 2]; % pull-stage distributed force applied to the free end
traction2=[0 0]; % release-stage distributed force applied to the free end

bodyForce2d=@(x)([0*x(:,1) 0*x(:,1)]');
boundaryTraction2d_1=@(x)([0*x(:,1) 0*x(:,1)]');
boundaryTraction2d_2_1=@(x)...% pull
    ([(traction1(1)+0*x(:,1)) (traction1(2)+0*x(:,1))]');
boundaryTraction2d_2_2=@(x)...% release
    ([(traction2(1)+0*x(:,1)) (traction2(2)+0*x(:,1))]');
boundaryTraction2d_3=@(x)([0*x(:,1) 0*x(:,1)]');
boundaryTraction2d_4=@(x)([0*x(:,1) 0*x(:,1)]');

boundaryTractions2d_1={boundaryTraction2d_1,boundaryTraction2d_2_1,...
    boundaryTraction2d_3,boundaryTraction2d_4}; % pull
boundaryTractions2d_2={boundaryTraction2d_1,boundaryTraction2d_2_2,...
    boundaryTraction2d_3,boundaryTraction2d_4}; % release
boundaryDisplacement2d=@(x,y)([0*x(:,1) 0*x(:,1)]'); % x,y - k-by-1, output must be numEq-by-1 or numEq-by-k
isDirichlet=[0; 0; 0; 1]; % bottom, right, top, left boundaries of 2D domain

%% ===== Mesh generation .... solution =====
% ===== Mesh generation =====
[nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect2d(domain, elementType,elementSize);
BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);

% id of the lower-right corner
hlID=find(nodeCoords(:,1)==domain(2,1)&nodeCoords(:,2)==domain(1,2)); %#ok<NASGU>

% misc
numEl=size(IEN,1); %#ok<NASGU> %total number of elements in the mesh
numNodes=size(nodeCoords,1); %total number of nodes in the mesh
numGDoF=numNodes*numEq; %total number of global degrees of freedom (elasticity)

% ===== Matrix, vector assembly =====
K = formStiffnessMatrix(nodeCoords, IEN, elementType, numGP, CMatrix);
Fb = formBodyForceVector(nodeCoords, IEN, elementType, numGP, bodyForce2d);
M = formMassMatrix(nodeCoords, IEN, elementType, numGP, MMatrix);

% ===== Part 1 - Pull
% Applying Boundary Conditions
[u_prescribed, Fs_1, prescribedDoF, freeDoF]=...
    formBC(nodeCoords,BIEN,bndElementType,numGPbnd,...
    boundaryTractions2d_1,boundaryDisplacement2d,isDirichlet);
a_prescribed=u_prescribed;
v_prescribed=u_prescribed;

F=Fb+Fs_1;

% Defining the free parts
KK=K(freeDoF,freeDoF);
FF=F(freeDoF)...
    -K(freeDoF,prescribedDoF)*u_prescribed(prescribedDoF)...
    -M(freeDoF,prescribedDoF)*a_prescribed(prescribedDoF);

% Solving for displacements
u=zeros(numGDoF,1);
u(freeDoF)=KK\FF;
u(prescribedDoF)=u_prescribed(prescribedDoF);

% Part 2 - Release
% Applying Boundary Conditions
[u_prescribed, Fs_2, prescribedDoF, freeDoF]=...
    formBC(nodeCoords,BIEN,bndElementType,numGPbnd,...
    boundaryTractions2d_2,boundaryDisplacement2d,isDirichlet);
[prescribedDoF,freeDoF]=preventRigidMotion(nodeCoords,prescribedDoF,freeDoF);

% Initial conditions
u0=u;
v0=zeros(numGDoF,1);
a0=zeros(numGDoF,1);

F=Fb+Fs_2;

M1=diag(sum(M)); % Lumping mass matrix

ts=0:dt:tfinal;


[output]=solveNewmarkLa(beta,gamma,ts,freeDoF,prescribedDoF,...
    M1,[],K,F,...
    u0,v0,a0,...
    u_prescribed,v_prescribed,a_prescribed,...
    'increment',@outputFunction);
u=output{1};
v=output{2};
a=output{3};

% function passed to solveNewmarkLa() used for output and plotting
    function res=outputFunction(i,t,u,v,a,fDoF,pDoF,M,C,K,F,u0,v0,a0,up,vp,ap) %#ok<INUSL,INUSD>
        if doplot && mod(i-1,displayevery)==0
            % Plotting the results
            u2=reshape(u,[numEq numel(u)/numEq])';
            [gradu]=recoveryAvg(u2,nodeCoords,IEN,elementType);
            s2=CMatrix*gradu;
            
            figure(1);clf;
            drawElementsInterp(nodeCoords+u2*factor,IEN,elementType,s2(1,:));
            axis((eye(4)+.2*[-1 -1 0 0; 1 1 0 0; 0 0 -1 -1; 0 0 1 1])*domain(:));
            title('\sigma_{xx} stress component');
            view(2)
            drawnow;
        end
        res={u,v,a};
    end

end


