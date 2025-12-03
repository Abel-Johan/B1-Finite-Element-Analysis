function example4_1
% Mechanical equilibrium of a linearly elastic plate.
% Here we solve equation and compare numerical results to the analytical
% exact value of stresses and displacements.
% For the derivation of the analyical solution see
% "plane_strain_analytic_benchmark.nb". 

% figure(1);clf(1);figure(2);clf(2);
close all;
do_plot=1;
%% ===== Problem Parameters ===== %
domain=[-2 -2 -1 -1]; %rectangular domain
% Material properties
youngsModulus = 193e6;
poissonsRatio = 0.253;
elProp=elasticProperties('youngsModulus',youngsModulus,'poissonsRatio',poissonsRatio);
CMatrix_E=elProp.CPlaneStrainFull;

%% ===== FEM Solver Parameters ===== %

%uncomment any of the below lines
% elementType='2dP1'; elementType1d='1dQ1'; numGP=1; numGP1d=2; 
%  elementType='2dP2'; elementType1d='1dQ2'; numGP=4; numGP1d=2;
% elementType='2dQ1'; elementType1d='1dQ1'; numGP=4; numGP1d=2;
% elementType='2dQ2r'; elementType1d='1dQ2'; numGP=4; numGP1d=2;
elementType='2dQ2r'; elementType1d='1dQ2'; numGP=9; numGP1d=2;
% elementType='2dQ2'; elementType1d='1dQ2'; numGP=9; numGP1d=2;

numEq_E=2; % Number of scalar equations (elasticity)
%% ===== Analytical solution ===== %
mu=youngsModulus/(2+2*poissonsRatio);mu=elProp.mu; %#ok<*NASGU>
nu=poissonsRatio;nu=elProp.poissonsRatio;

s11_=@(x,y)(1/9).*((-18)+(-1).*exp(1).^(3.*x).*cos(3.*y)+cos(3.*x).*cosh(3.* ...
    y));
s12_=@(x,y)(1/9).*(exp(1).^(3.*x).*sin(3.*y)+sin(3.*x).*sinh(3.*y));
s22_=@(x,y)(1/9).*((-18)+exp(1).^(3.*x).*cos(3.*y)+(-1).*cos(3.*x).*cosh(3.* ...
    y));
c1=(1/54).*exp(1).^(-6).*mu.^(-1).*((-1).*exp(1).^3.*sin(3)+(-2).* ...
    sin(6)+exp(1).^6.*((-54)+108.*nu+cos(3).*sinh(3)+(-1).*cos(6).* ...
    sinh(6)+3.*sin(6).*sinh(6)));
u_=@(x,y)(1/108).*exp(1).^(-6).*mu.^(-1).*(108.*exp(1).^6.*(((-1)+2.*nu).*( ...
    2+x)+c1.*mu.*(2+y))+2.*cos(6)+(-2).*exp(1).^(6+3.*x).*cos(3.*y)+( ...
    19+9.*y+(-1).*exp(1).^12.*(5+3.*y)).*sin(6)+2.*exp(1).^6.*cosh(3.* ...
    y).*sin(3.*x));
v_=@(x,y)(1/108).*exp(1).^(-6).*mu.^(-1).*((-108).*exp(1).^6.*(c1.*mu.*(2+ ...
    x)+(-1).*((-1)+2.*nu).*(2+y))+cos(6)+(-1).*(16+9.*x).*sin(6)+exp( ...
    1).^12.*((-1).*cos(6)+3.*(2+x).*sin(6))+2.*exp(1).^(6+3.*x).*sin( ...
    3.*y)+(-2).*exp(1).^6.*cos(3.*x).*sinh(3.*y));
ux_=@(x,y)(1/18).*mu.^(-1).*((-18)+36.*nu+(-1).*exp(1).^(3.*x).*cos(3.*y)+ ...
    cos(3.*x).*cosh(3.*y));
uy_=@(x,y)(1/36).*exp(1).^(-6).*mu.^(-1).*(36.*c1.*exp(1).^6.*mu+(-1).*((-3) ...
    +exp(1).^12).*sin(6)+2.*exp(1).^(6+3.*x).*sin(3.*y)+2.*exp(1).^6.* ...
    sin(3.*x).*sinh(3.*y));
vx_=@(x,y)(1/18).*exp(1).^(-6).*mu.^(-1).*((-1).*sin(6)+exp(1).^(6+3.*x).* ...
    sin(3.*y)+exp(1).^6.*((-18).*c1.*mu+sin(6).*sinh(6))+exp(1).^6.* ...
    sin(3.*x).*sinh(3.*y));
vy_=@(x,y)(1/18).*mu.^(-1).*((-18)+36.*nu+exp(1).^(3.*x).*cos(3.*y)+(-1).* ...
    cos(3.*x).*cosh(3.*y));

%% ===== Specifying Body Force and Boundary conditions =====
bodyForce2d_E=@(x,y)([0*x 0*x]'); % Elastic part
bodyForce2d_E=@(x)(bodyForce2d_E(x(:,1),x(:,2)));

boundaryTraction2d_1=@(x,y)([-s12_(x,y) -s22_(x,y)]');
boundaryTraction2d_1=@(x)boundaryTraction2d_1(x(:,1),x(:,2));

boundaryTraction2d_2=@(x,y)([s11_(x,y) s12_(x,y)]');
boundaryTraction2d_2=@(x)boundaryTraction2d_2(x(:,1),x(:,2));

boundaryTraction2d_3=@(x,y)([s12_(x,y) s22_(x,y)]');
boundaryTraction2d_3=@(x)boundaryTraction2d_3(x(:,1),x(:,2));

boundaryTraction2d_4=@(x,y)([-s11_(x,y) -s12_(x,y)]');
boundaryTraction2d_4=@(x)boundaryTraction2d_4(x(:,1),x(:,2));

boundaryTraction2d_E={boundaryTraction2d_1,boundaryTraction2d_2,boundaryTraction2d_3,boundaryTraction2d_4};
isDirichlet_E=[0 0; 0 0; 0 0; 0 0]; % <-- uncomment this for Neumann
% isDirichlet_E=[1 1; 1 1; 1 1; 1 1]; % <-- or this for Dirichlet

boundaryDisplacement2d_E=@(x,y)( [u_(x,y) v_(x,y)]' ); % x,y - k-by-1, output must be numEq-by-1 or numEq-by-k
boundaryDisplacement2d_E=@(x)(boundaryDisplacement2d_E(x(:,1),x(:,2)));

%% ============ ========== PART 1 ========== ============
    function [numEl,h,err_s_GP,err_s_C,err_s_N,err_s_SPR,err_u_N,t]=go(elementSize)
        tic;
        % ===== Mesh generation =====
        [nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect2d(domain, elementType,elementSize);
        BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);
        [INE,degrees]=IENtoINE(IEN);
        
        % misc
        numEl=size(IEN,1); %total number of elements in the mesh
        numNodes=size(nodeCoords,1); %total number of nodes in the mesh
        numGDoF_E=numNodes*numEq_E; %total number of global degrees of freedom (elasticity)
        
        % ===== Matrix, vector assembly =====
        K_E = formStiffnessMatrix(nodeCoords, IEN, elementType, numGP, CMatrix_E);
        Fb_E = formBodyForceVector(nodeCoords, IEN, elementType, numGP, bodyForce2d_E);
        
        % Applying Boundary Conditions
        [u_prescribed_E, Fs_E, prescribedDoF_E, freeDoF_E]=...
            formBC(nodeCoords,BIEN,elementType1d,numGP1d,...
            boundaryTraction2d_E,boundaryDisplacement2d_E,isDirichlet_E);
        
        % Preventing rigid body motion
        [prescribedDoF_E,freeDoF_E]=preventRigidMotion(nodeCoords,prescribedDoF_E,freeDoF_E);
        
        % total force
        F_E=Fb_E+Fs_E;
        
        % Defining the free parts of stiffness
        KK_E=K_E(freeDoF_E,freeDoF_E);
        FF_E=F_E(freeDoF_E)-K_E(freeDoF_E,prescribedDoF_E)*u_prescribed_E(prescribedDoF_E);
        
        % ===== Solution =====
        u_E=zeros(numGDoF_E,1);
        u_E(prescribedDoF_E)=u_prescribed_E(prescribedDoF_E);
        u_E(freeDoF_E)=KK_E\FF_E;
        
        u2_E=reshape(u_E,[numEq_E numel(u_E)/numEq_E])';
        u2_anal_N=[u_(nodeCoords(:,1),nodeCoords(:,2)) v_(nodeCoords(:,1),nodeCoords(:,2))];
        % ===== Stress recovery =====
        [gradu_GP, GPCoords]=recoveryGP(u2_E,nodeCoords,IEN,elementType,numGP);
        [gradu_C, CCoords]=recoveryGP(u2_E,nodeCoords,IEN,elementType,1);
        [gradu_N]=recoveryAvg(u2_E,nodeCoords,IEN,elementType);
        [gradu_SPR]=recoveryPatch(u2_E,nodeCoords,IEN,elementType,BIEN,INE,degrees);
        
        s2_GP=CMatrix_E*gradu_GP;
        s2_C=CMatrix_E*gradu_C;
        s2_N=CMatrix_E*gradu_N;
        s2_SPR=CMatrix_E*gradu_SPR;
        s2_anal_GP=[s11_(GPCoords(:,1),GPCoords(:,2))'; s12_(GPCoords(:,1),GPCoords(:,2))';...
            s12_(GPCoords(:,1),GPCoords(:,2))'; s22_(GPCoords(:,1),GPCoords(:,2))'];
        s2_anal_C=[s11_(CCoords(:,1),CCoords(:,2))'; s12_(CCoords(:,1),CCoords(:,2))';...
            s12_(CCoords(:,1),CCoords(:,2))'; s22_(CCoords(:,1),CCoords(:,2))'];
        s2_anal_N=[s11_(nodeCoords(:,1),nodeCoords(:,2))'; s12_(nodeCoords(:,1),nodeCoords(:,2))';...
            s12_(nodeCoords(:,1),nodeCoords(:,2))'; s22_(nodeCoords(:,1),nodeCoords(:,2))'];
        
        % ===== Errors =====
        %This is an approximation of the L2-norm on uniform grids.
        err_s_GP=norm((s2_GP-s2_anal_GP)/sqrt(numel(s2_GP)));
        err_s_C=norm((s2_C-s2_anal_C)/sqrt(numel(s2_C)));
        err_s_N=norm((s2_N-s2_anal_N)/sqrt(numel(s2_N)));
        err_s_SPR=norm((s2_SPR-s2_anal_N)/sqrt(numel(s2_SPR)));
        err_u_N=norm((u2_E-u2_anal_N)/sqrt(numel(u2_E)));
        
        curElCoords=nodeCoords(IEN(1,:),:);
        h=norm(curElCoords(1,:)-curElCoords(2,:));
        
        t=toc;
        % ===== Plots ===== %
        if do_plot==1
            factor=1e6;
            figure(3); clf;
            drawElementsInterp(nodeCoords+u2_E*factor,IEN,elementType,s2_SPR(4,:),.5*ones(size(IEN,1),1));
            hold on;
            view(2); title('yy-component of stress')
            drawnow;
        end
    end

%% Trying different meshes
elementSize={[6,6],[15,15],[24,24],[32,32],[48,48],[64,64]};
elementSize={[6,6],[15,15],[24,24],[32,32]};

hs=zeros(1,numel(elementSize));
numEls=hs;
err_ss_GP=hs;
err_ss_C=hs;
err_ss_N=hs;
err_ss_SPR=hs;
err_us=hs;
t=hs;

for i=1:numel(elementSize)
    [numEls(i),hs(i),err_ss_GP(i),err_ss_C(i),err_ss_N(i),err_ss_SPR(i),err_us(i),t(i)]=go(elementSize{i});
end

numEls,hs,err_ss_GP,err_ss_C,err_ss_N,err_ss_SPR,err_us %#ok<NOPRT>

%% plotting convergence

xdata=log(hs);
ydatau=log(err_us);
ydatas_GP=log(err_ss_GP);
ydatas_C=log(err_ss_C);
ydatas_N=log(err_ss_N);
ydatas_SPR=log(err_ss_SPR);

slopeu=diff(ydatau)./diff(xdata);
slopes_GP=diff(ydatas_GP)./diff(xdata);
slopes_C=diff(ydatas_C)./diff(xdata);
slopes_N=diff(ydatas_N)./diff(xdata);
slopes_SPR=diff(ydatas_SPR)./diff(xdata);

figure(1);set(gcf,'Position',[432.2000   84.2000  882.4000  681.6000]);
subplot(2,1,1);hold off;
loglog(hs,err_us,'--*k'); hold on;
grid on;
xlabel('h, max element diameter')
ylabel('Displacement error')
axis equal
title({['Element type = ', elementType ',  numGP = ' num2str(numGP)],...
    ['slope= ' num2str(fliplr(slopeu))]});
axis equal

subplot(2,1,2);hold off;
loglog(hs,err_ss_GP,'--*r');hold on
loglog(hs,err_ss_C,'--*b');
loglog(hs,err_ss_N,'--*k');
loglog(hs,err_ss_SPR,'--*g');
grid on;
xlabel('h, max element diameter')
ylabel('Stress error')
title({['slope@GP= ' num2str(fliplr(slopes_GP))],...
    ['slope@C= ' num2str(fliplr(slopes_C))],...
    ['slope@N= ' num2str(fliplr(slopes_N))],...
    ['slope@N(SPR)= ' num2str(fliplr(slopes_SPR))]});
legend('at GP','at Centroids','at Nodes','at Nodes (SPR)')
disp(['Respective elapsed time'])
disp(num2str(t))
disp(['Total time:' num2str(sum(t))])
end
