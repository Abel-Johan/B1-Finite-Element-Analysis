function example5_1
% Mechanical equilibrium - plate with a hole
% Here we generate mesh using MATLAB pdetool, solve equation and compare
% numerical results to the analytical solution

% Domain
a=.5; %radius of the hole
L=1; %size of the plate

elementType='2dP1';
elementType1d='1dQ1';
numGP=4;
numGP1d=2;

numEq_E=2; % Number of scalar equations (elasticity)

%% Material properties
youngsModulus = 193e6;
poissonsRatio = 0.253;
properties=elasticProperties('youngsModulus',youngsModulus,'poissonsRatio',poissonsRatio);
CMatrix_E=properties.CPlaneStrainFull;
%% ===== Analytical solution ===== %
mu=properties.mu;
nu=properties.poissonsRatio;

s11_=@(x,y)(1/2).*(x.^2+y.^2).^(-4).*(2.*(x.^2+y.^2).^4+3.*a.^4.*(x.^4+(-6).* ...
    x.^2.*y.^2+y.^4)+a.^2.*((-5).*x.^6+7.*x.^4.*y.^2+13.*x.^2.*y.^4+ ...
    y.^6));
s12_=@(x,y)a.^2.*x.*y.*(x.^2+y.^2).^(-4).*((-5).*x.^4+(-2).*x.^2.*y.^2+3.* ...
    y.^4+6.*a.^2.*(x+(-1).*y).*(x+y));
s22_=@(x,y)(1/2).*a.^2.*(x.^2+y.^2).^(-4).*((-3).*a.^2.*(x.^4+(-6).*x.^2.* ...
    y.^2+y.^4)+(x.^2+y.^2).*(x.^4+(-12).*x.^2.*y.^2+3.*y.^4));
u_=@(x,y)(1/4).*mu.^(-1).*(x.^2+y.^2).^(-3/2).*((-2).*((-1)+nu).*x.*(x.^2+ ...
    y.^2).^(1/2).*(2.*a.^2+x.^2+y.^2)+a.^2.*((-1).*a.^2+x.^2+y.^2).* ...
    cos(3.*atan2(y,x)));
v_=@(x,y)(1/4).*mu.^(-1).*(x.^2+y.^2).^(-3/2).*((-2).*y.*(x.^2+y.^2).^(1/2) ...
    .*(a.^2.*(1+(-2).*nu)+nu.*(x.^2+y.^2))+a.^2.*((-1).*a.^2+x.^2+ ...
    y.^2).*sin(3.*atan2(y,x)));
ux_=@(x,y)(1/4).*mu.^(-1).*(x.^2+y.^2).^(-5/2).*(2.*((-1)+nu).*(x.^2+y.^2) ...
    .^(1/2).*(2.*a.^2.*(x+(-1).*y).*(x+y)+(-1).*(x.^2+y.^2).^2)+a.^2.* ...
    ((-1).*x.*((-3).*a.^2+x.^2+y.^2).*cos(3.*atan2(y,x))+3.*y.*((-1).* ...
    a.^2+x.^2+y.^2).*sin(3.*atan2(y,x)))); %#ok<NASGU>
uy_=@(x,y)(1/2).*a.^2.*mu.^(-1).*x.*y.*(x.^2+y.^2).^(-4).*(6.*a.^2.*(x+(-1) ...
    .*y).*(x+y)+(x.^2+y.^2).*(((-9)+4.*nu).*x.^2+((-1)+4.*nu).*y.^2)); ...
     %#ok<NASGU>
vx_=@(x,y)(1/2).*a.^2.*mu.^(-1).*x.*y.*(x.^2+y.^2).^(-4).*(6.*a.^2.*(x+(-1) ...
    .*y).*(x+y)+(-1).*(x.^2+y.^2).*((1+4.*nu).*x.^2+((-7)+4.*nu).* ...
    y.^2)); %#ok<NASGU>
vy_=@(x,y)(1/4).*mu.^(-1).*(x.^2+y.^2).^(-4).*((-2).*nu.*(x.^2+y.^2).^4+(-3) ...
    .*a.^4.*(x.^4+(-6).*x.^2.*y.^2+y.^4)+a.^2.*(x.^2+y.^2).*((1+4.*nu) ...
    .*x.^4+(-12).*x.^2.*y.^2+(3+(-4).*nu).*y.^4)); %#ok<NASGU>

%% ===== Specifying Body Force and Boundary conditions =====
% 1 2 3 4 5 = right, top, bottom, left, inner hole
bodyForce2d_E=@(x,y)([0*x 0*x]'); % Elastic part
bodyForce2d_E=@(x)(bodyForce2d_E(x(:,1),x(:,2)));

boundaryTraction2d_1=@(x,y)([s11_(x,y) s12_(x,y)]');
boundaryTraction2d_1=@(x)boundaryTraction2d_1(x(:,1),x(:,2));

boundaryTraction2d_2=@(x,y)([s12_(x,y) s22_(x,y)]');
boundaryTraction2d_2=@(x)boundaryTraction2d_2(x(:,1),x(:,2));

boundaryTraction2d_3=@(x,y)([-s12_(x,y) -s22_(x,y)]');
boundaryTraction2d_3=@(x)boundaryTraction2d_3(x(:,1),x(:,2));

boundaryTraction2d_4=@(x,y)([-s11_(x,y) -s12_(x,y)]');
boundaryTraction2d_4=@(x)boundaryTraction2d_4(x(:,1),x(:,2));

boundaryTraction2d_5=@(x,y)([0*x 0*x]');
boundaryTraction2d_5=@(x)boundaryTraction2d_5(x(:,1),x(:,2));

boundaryTraction2d_E={boundaryTraction2d_1,boundaryTraction2d_2,boundaryTraction2d_3,boundaryTraction2d_4,boundaryTraction2d_5};
isDirichlet_E=[0 0; 0 0; 0 0; 0 0; 0 0]; % <-- uncomment this for Neumann
isDirichlet_E=[1 1; 1 1; 1 1; 1 1; 1 1]; % <-- or this for Dirichlet
% isDirichlet_E=[0 0; 0 0; 0 1; 1 0; 0 0]; % <-- or this for symmetry BC

boundaryDisplacement2d_E=@(x,y)( [u_(x,y) v_(x,y)]' ); % x,y - k-by-1, output must be numEq-by-1 or numEq-by-k
boundaryDisplacement2d_E=@(x)(boundaryDisplacement2d_E(x(:,1),x(:,2)));

%% ===== Creating geometry using Matlab PDEToolbox =====
csg_rect=[3 4 0 L L 0 0 0 L L]'; %define the rectangle
csg_circ=[1 0 0 a]'; %define the circle
gd=[csg_rect, [csg_circ; zeros(numel(csg_rect)-numel(csg_circ),1)]]; %combine
ns = char('rect1','circ1')'; %Give names to the shapes - namespace matrix
sf = 'rect1-circ1';% Specify the set formula
[dl,bt] = decsg(gd,sf,ns); %#ok<ASGLU> % produces decomposed geometry

%% ===== PART 1 =====
    function [numEl,err_s_GP,err_u]=go(hmax)
        % using PDEToolbox to generate triangular mesh
        [nodeCoords,e,t] = initmesh(dl,'Hmax',hmax);
        nodeCoords=nodeCoords'; % bring it to the form we use
        IEN=t(1:3,:)'; % bring it to the form we use
        
        INE=IENtoINE(IEN);
        [boundaryElementIDs, boundaryNodeLocalIDs]=edges2sublists(e,IEN,INE);
        
        % plot the mesh
        figure(1);clf;subplot(2,1,1);
        drawElements(nodeCoords,IEN,elementType);view(2); hold on;
        drawNodes(nodeCoords([e(1,:) e(2,:)],:)); hold off; colorbar;
        subplot(2,1,2); pdegplot(dl,'EdgeLabels','on'); axis tight; axis equal;colorbar;
        drawnow;
        
        numEl=size(IEN,1); %total number of elements in the mesh
        numNodes=size(nodeCoords,1); %total number of nodes in the mesh
        numGDoF_E=numNodes*numEq_E; %total number of global degrees of freedom (elasticity)
        
        % ===== Matrix, vector assembly =====
        K_E = formStiffnessMatrix(nodeCoords, IEN, elementType, numGP, CMatrix_E);
        Fb_E = formBodyForceVector(nodeCoords, IEN, elementType, numGP, bodyForce2d_E);
        
        boundaryIEN=IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);
        [u_prescribed_E, Fs_E, prescribedDoF_E, freeDoF_E]=...
            formBC(nodeCoords,boundaryIEN,elementType1d,numGP1d,...
            boundaryTraction2d_E,boundaryDisplacement2d_E,isDirichlet_E);
        
        % Preventing rigid body motion
        [prescribedDoF_E,freeDoF_E]=preventRigidMotion(nodeCoords,prescribedDoF_E,freeDoF_E);
        
        % total force
        F_E=Fb_E+Fs_E;
        
        % Defining the free parts of stiffness
        KK_E=K_E(freeDoF_E,freeDoF_E);
        if ~isempty(prescribedDoF_E)
            FF_E=F_E(freeDoF_E)-K_E(freeDoF_E,prescribedDoF_E)*u_prescribed_E(prescribedDoF_E);
        else
            FF_E=F_E;
            warning('Need to prescribe at least some DoFs');
        end
        
        % Solving for mechanical equilibrium
        u_E=zeros(numGDoF_E,1);
        u_E(prescribedDoF_E)=u_prescribed_E(prescribedDoF_E);
        u_E(freeDoF_E)=KK_E\FF_E;
        
        factor=1e8;
        u2_E=reshape(u_E,[numEq_E numel(u_E)/numEq_E])';
        u2_anal=[u_(nodeCoords(:,1),nodeCoords(:,2)) v_(nodeCoords(:,1),nodeCoords(:,2))];
        
        % ===== Stress recovery =====
        [gradu_GP, GPCoords]=recoveryGP(u2_E,nodeCoords,IEN,elementType,1);
        s2_GP=CMatrix_E*gradu_GP;
        s2_anal_GP=[s11_(GPCoords(:,1),GPCoords(:,2))'; s12_(GPCoords(:,1),GPCoords(:,2))';...
            s12_(GPCoords(:,1),GPCoords(:,2))'; s22_(GPCoords(:,1),GPCoords(:,2))'];
        
        % ===== Errors =====
        err_s_GP=norm((s2_GP-s2_anal_GP)/sqrt(numel(s2_GP)));
        err_u=norm((u2_E-u2_anal)/sqrt(numel(u2_E)));
        
        curElCoords=nodeCoords(IEN(1,:),:);
        h=norm(curElCoords(1,:)-curElCoords(2,:));%#ok<NASGU>

        % ====== Plots =====
        fprintf('\n Number of elements: %8i \n', numEl)
        
        figure(2); clf;
        subplot(2,1,1)
        drawElements(nodeCoords+u2_E*factor,IEN,elementType,s2_GP(1,:),.5*ones(size(IEN,1),1));
        view(2); title('xx-component of stress - numerical')
        drawnow;
        subplot(2,1,2)
        drawElements(nodeCoords+u2_E*factor,IEN,elementType,s2_anal_GP(1,:),.5*ones(size(IEN,1),1));
        view(2); title('xx-component of stress - analytical')
        drawnow;
    end %function go


%% ==========   ================== 
hmaxs=[.4 .3 .2 .15 .1 .06 .05 .04 .033];
numEls=zeros(size(hmaxs));
err_ss=numEls;
err_us=numEls;

for i=1:numel(hmaxs)
    [numEls(i),err_ss(i),err_us(i)]=go(hmaxs(i));
end

numEls,hmaxs,err_ss,err_us %#ok<NOPRT>

%% plotting convergence

xdata=log(hmaxs);
ydatau=log(err_us);
ydatas=log(err_ss);

slopeu=diff(ydatau)./diff(xdata) %#ok<NOPRT>
slopes=diff(ydatas)./diff(xdata) %#ok<NOPRT>

figure(3);
subplot(2,1,1);
loglog(hmaxs,err_us,'--*k');
grid on;
xlabel('h, max element diameter')
ylabel('Displacement error')
axis equal
title({['Element Type = ' elementType],...
    ['slope= ' num2str(fliplr(slopeu))]});

subplot(2,1,2);
loglog(hmaxs,err_ss,'--*k');
grid on;
xlabel('h, max element diameter')
ylabel('Stress error')
legend('at Centroids');
title(['slope= ' num2str(fliplr(slopes))]);
% axis equal

end