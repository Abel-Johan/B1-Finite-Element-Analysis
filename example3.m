% Steady state mechanical equilibrium and dynamics of heat diffusion on a
% rectangular domain (linear plane stress problem). The only aspect of
% thermomechanical coupling considered is heat expansion. 

doplot=1; % if we want plots
displayevery=100; % Every k-th timestep to plot
%% ===== Problem Parameters ===== %
domain=[0 0 2 2];
% Material properties
thickness=0.1; 
youngsModulus = 193e6;
poissonsRatio = 0.253;
properties=elasticProperties('youngsModulus',youngsModulus,'poissonsRatio',poissonsRatio);
CMatrix_E=properties.CPlaneStressFull;
% Heat conductivity tensor 
CondCoeff=16.2;
CMatrix_T=CondCoeff*[1 0; 0 1];
% Heat capacity 
massDensity=7.999e3;
specificCapacity=5e2;
Capacity=thickness*massDensity*specificCapacity;
% Combined CMatrix
CMatrix=blkdiag(CMatrix_E,CMatrix_T);
% Thermal expansion
TEcoeff=9.7e-6;
TEtensor=TEcoeff*[1 0; 0 1]; %thermal expansion tensor
TEtensortr=TEtensor';
BetaMatrix=CMatrix_E*TEtensortr(:); % th-->gradu map

%% ===== FEM Solver and Presentation Parameters ===== %
elementType='2dQ1';
elementType1d='1dQ1';
numGP=4;
numGP1d=2;
elementSize=.1; %element diameter
% misc
numDim=2; % A 2D problem
numEq_E=2; % Number of scalar equations: 2 for elastic part, 1 for thermal
numEq_T=1; % Number of scalar equations: 2 for elastic part, 1 for thermal
% time-stepping parameters
alpha=1; % generalized Trapezoidal method parameter
finalTime=1e5;
numTimesteps=800;
dt=finalTime/numTimesteps;
numTimesteps=finalTime/dt;
factor=1e2; %displacement display amplification factor


%% ===== Specifying Body Force and Boundary and Initial conditions =====
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

init_u0_2d_T=@(x,y)([000*ones(size(x))]); % x,y - column-vectors of same size
init_u0_2d_T=@(x)(init_u0_2d_T(x(:,1),x(:,2))); % x,y - column-vectors of same size

%% Mesh generation
[nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect2d(domain, elementType, elementSize);
BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);

% misc
numEl=size(IEN,1); %total number of elements in the mesh
numNodes=size(nodeCoords,1); %total number of nodes in the mesh 
numGDoF_E=numNodes*numEq_E; %total number of displacement/temp degrees of freedom
numGDoF_T=numNodes*numEq_T; %total number of displacement/temp degrees of freedom

% Drawing the original mesh
if doplot
    figure(1); clf;
    drawElements(nodeCoords,IEN,elementType)
    hold on;
    view(2);
end

%% Matrix, vector assembly
K_E = formStiffnessMatrix(nodeCoords, IEN, elementType, numGP, CMatrix_E);
Fb_E = formBodyForceVector(nodeCoords, IEN, elementType, numGP, bodyForce2d_E);

K_T = formStiffnessMatrix(nodeCoords, IEN, elementType, numGP, CMatrix_T);
Fb_T = formBodyForceVector(nodeCoords, IEN, elementType, numGP, bodyForce2d_T);
C_T = formMassMatrix(nodeCoords, IEN, elementType, numGP, Capacity);
u0_T=init_u0_2d_T(nodeCoords); %initial temperature distribution

%% Applying Boundary Conditions
[u_prescribed_E, Fs_E, prescribedDoF_E, freeDoF_E]=...
    formBC(nodeCoords,BIEN,elementType1d,numGP1d,...
    boundaryTraction2d_E,boundaryDisplacement2d_E,isDirichlet_E);
        
[u_prescribed_T, Fs_T, prescribedDoF_T, freeDoF_T]=...
    formBC(nodeCoords,BIEN,elementType1d,numGP1d,...
    boundaryTraction2d_T,boundaryDisplacement2d_T,isDirichlet_T);
v_prescribed_T=zeros(size(u_prescribed_T)); %velocity needs to be defined for the heat equation
u0_T(prescribedDoF_T)=u_prescribed_T(prescribedDoF_T);

%% Solution
% Thermal part
F_T=Fs_T+Fb_T;
KK_T=K_T(freeDoF_T,freeDoF_T);
CC_T=C_T(freeDoF_T,freeDoF_T);
if ~isempty(prescribedDoF_T)
    FF_T=F_T(freeDoF_T)-K_T(freeDoF_T,prescribedDoF_T)*u_prescribed_T(prescribedDoF_T);
else
    FF_T=F_T(freeDoF_T);
end

u_T=u0_T; %initial temperature and rate
v_T=zeros(numGDoF_T,1);
v_T(freeDoF_T)=CC_T\(FF_T+KK_T*u_T(freeDoF_T));

u_E=zeros(numGDoF_E,1);
u_E(prescribedDoF_E)=u_prescribed_E(prescribedDoF_E);
for i=1:numTimesteps
    % make a time-step
    u_old_T=u_T(freeDoF_T); 
    v_old_T=v_T(freeDoF_T);
    % effective right-hand side
    rhs=FF_T+CC_T*(u_old_T+(1-alpha)*dt*v_old_T)/(alpha*dt);
    % effective linear operator
    lhs=(CC_T+alpha*dt*KK_T)/(alpha*dt);
    u_T(freeDoF_T)=lhs\rhs; %solution
    v_T(freeDoF_T)=(u_T(freeDoF_T)-u_old_T-(1-alpha)*dt*v_old_T)/alpha/dt;
    
    % compute elastic componenets and plot results at selected timesteps
    if mod(i,displayevery)==1
        % Elastic part
        Ft_E = formInternalForce2(nodeCoords, IEN, elementType, numGP, ...
            zeros(4,2), BetaMatrix, u_T);
        F_E=Fs_E+Fb_E+Ft_E;
        KK_E=K_E(freeDoF_E,freeDoF_E);
        if ~isempty(prescribedDoF_E)
            FF_E=F_E(freeDoF_E)-K_E(freeDoF_E,prescribedDoF_E)*u_prescribed_E(prescribedDoF_E);
        else
            FF_E=F_E(freeDoF_E);
        end
        u_E(freeDoF_E)=KK_E\FF_E;
        
        [gradu_E,GPCoords]=recoveryGP(u_E,nodeCoords,IEN,elementType,1);
        u2_E=reshape(u_E,[numEq_E numel(u_E)/numEq_E])';
        if doplot
            % Plotting the results
            hold on;
            figure(1); clf;
            drawElementsInterp(nodeCoords+u2_E*factor,IEN,elementType,u_T);
            title({'Temperature distribution',['Deformed mesh, t=' num2str(dt*i)]});
            view(2)
            drawnow;
        end
    end
end


