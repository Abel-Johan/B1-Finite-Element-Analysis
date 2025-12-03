function example7_2%%%% Mechanical equilibrium 3D and 2D - bending cantilever beam
% plane strain vs 3dQ1 elements

%% ===== Problem Parameters ===== %
domain=[0 0; 10 3]; %rectangular domain
thickness=1;
domain3d=[domain(1,:) 0;  domain(2,:) thickness]; %rectangular domain
% Material properties
youngsModulus = 193e6;
poissonsRatio = 0.253;
properties=elasticProperties('youngsModulus',youngsModulus,'poissonsRatio',poissonsRatio);
%% ====== BC
bodyForce2d_E=@(x,y)([0*x 0*x]'); % Elastic part
bodyForce2d_E=@(x)(bodyForce2d_E(x(:,1),x(:,2)));

boundaryTraction2d_1=@(x,y)([0*x 0*x]');
boundaryTraction2d_1=@(x)boundaryTraction2d_1(x(:,1),x(:,2));

boundaryTraction2d_2=@(x,y)([(500+0*x) (2+0*x)]');
boundaryTraction2d_2=@(x)boundaryTraction2d_2(x(:,1),x(:,2));

boundaryTraction2d_3=@(x,y)([0*x 0*x]');
boundaryTraction2d_3=@(x)boundaryTraction2d_3(x(:,1),x(:,2));

boundaryTraction2d_4=@(x,y)([0*x 0*x]');
boundaryTraction2d_4=@(x)boundaryTraction2d_4(x(:,1),x(:,2));

boundaryTraction2d_E={boundaryTraction2d_1,boundaryTraction2d_2,...
    boundaryTraction2d_3,boundaryTraction2d_4};

boundaryDisplacement2d_E=@(x,y)([0*x 0*x]'); % x,y - k-by-1, output must be numEq-by-1 or numEq-by-k
boundaryDisplacement2d_E=@(x)(boundaryDisplacement2d_E(x(:,1),x(:,2)));

%% ============ ========== PART 1 (3D) ========== ============

    function [u2_E,smpGrad,nodeCoords]=go3d(elementSize,smpPts)
        %% ===== Problem Parameters ===== %
        CMatrix_E=properties.CFull;
        %% ===== FEM Solver Parameters ===== %
        
        %uncomment any of the below lines
        elementType='3dQ1';
        elData=elementData(elementType);
        bndElementType=elData.bndElementType;
        numGP=elData.numGPFull;
        bndElData=elementData(bndElementType);
        numGPbnd=bndElData.numGPFull;
        
        numEq_E=3; % Number of scalar equations (elasticity)
        %% ===== Specifying Body Force and Boundary conditions =====
        bodyForce3d_E=@(x,y,z)([bodyForce2d_E([x y])' 0*z]'); % Elastic part
        bodyForce3d_E=@(x)(bodyForce3d_E(x(:,1),x(:,2),x(:,3)));
        
        boundaryTraction3d_1=@(x,y,z)([boundaryTraction2d_E{1}([x y])' 0*z]');
        boundaryTraction3d_1=@(x)boundaryTraction3d_1(x(:,1),x(:,2),x(:,3));
        
        boundaryTraction3d_2=@(x,y,z)([boundaryTraction2d_E{3}([x y])' 0*z]');
        boundaryTraction3d_2=@(x)boundaryTraction3d_2(x(:,1),x(:,2),x(:,3));
        
        boundaryTraction3d_3=@(x,y,z)([boundaryTraction2d_E{4}([x y])' 0*z]');
        boundaryTraction3d_3=@(x)boundaryTraction3d_3(x(:,1),x(:,2),x(:,3));
        
        boundaryTraction3d_4=@(x,y,z)([boundaryTraction2d_E{2}([x y])' 0*z]');
        boundaryTraction3d_4=@(x)boundaryTraction3d_4(x(:,1),x(:,2),x(:,3));
        
        boundaryTraction3d_5=@(x,y,z)([0*x 0*x 0*x]');
        boundaryTraction3d_5=@(x)boundaryTraction3d_5(x(:,1),x(:,2),x(:,3));
        
        boundaryTraction3d_6=@(x,y,z)([0*x 0*x 0*x]');
        boundaryTraction3d_6=@(x)boundaryTraction3d_6(x(:,1),x(:,2),x(:,3));
        
        boundaryTraction3d_E={boundaryTraction3d_1,boundaryTraction3d_2,...
            boundaryTraction3d_3,boundaryTraction3d_4,...
            boundaryTraction3d_5,boundaryTraction3d_6};
        
        isDirichlet_E=[0 0 1; 0 0 1; ... %bottom and top
            1 1 1; 0 0 1; .... %left and right
            0 0 1; 0 0 1]; %front and back

        boundaryDisplacement3d_E=@(x,y,z)([boundaryDisplacement2d_E([x y])' 0*x]'); % x,y - k-by-1, output must be numEq-by-1 or numEq-by-k
        boundaryDisplacement3d_E=@(x)(boundaryDisplacement3d_E(x(:,1),x(:,2),x(:,3)));
        %% Generation ... solution
        % ===== Mesh generation =====
        [nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect3d(domain3d, elementType,elementSize);
        boundaryIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);
        INE=IENtoINE(IEN);
        
        % misc
        numEl=size(IEN,1); %#ok<NASGU> %total number of elements in the mesh
        numNodes=size(nodeCoords,1); %total number of nodes in the mesh
        numGDoF_E=numNodes*numEq_E; %total number of global degrees of freedom (elasticity)
        
        % ===== Matrix, vector assembly =====
        K_E = formStiffnessMatrix(nodeCoords, IEN, elementType, numGP, CMatrix_E);
        Fb_E = formBodyForceVector(nodeCoords, IEN, elementType, numGP, bodyForce3d_E);
        
        % Applying Boundary Conditions
        [u_prescribed_E, Fs_E, prescribedDoF_E, freeDoF_E]=...
            formBC(nodeCoords,boundaryIEN,bndElementType,numGPbnd,...
            boundaryTraction3d_E,boundaryDisplacement3d_E,isDirichlet_E);
        
        % Preventing rigid body motifon
        [prescribedDoF_E,freeDoF_E]=preventRigidMotion(nodeCoords,prescribedDoF_E,freeDoF_E);
        F_E=Fb_E+Fs_E;
        
        % Defining the free parts of stiffness
        KK_E=K_E(freeDoF_E,freeDoF_E);
        FF_E=F_E(freeDoF_E)-K_E(freeDoF_E,prescribedDoF_E)*u_prescribed_E(prescribedDoF_E);
        
        % ===== Solution =====
        u_E=zeros(numGDoF_E,1);
        u_E(prescribedDoF_E)=u_prescribed_E(prescribedDoF_E);
        u_E(freeDoF_E)=KK_E\FF_E;
        
        u2_E=reshape(u_E,numEq_E,numNodes)';
        %% ===== Stress recovery =====
        recVals1=recoveryAvg(u2_E,nodeCoords,IEN,elementType);
        s2=CMatrix_E*recVals1;
        
        [smpEl,smpNatCoords] = meshGetNatural(nodeCoords,IEN,elementType,INE,smpPts);
        smpGrad=recoveryEvaluateGradients(u2_E,nodeCoords,IEN,elementType,smpEl,smpNatCoords);

        %===== Plots ===== %
        shift=ones(size(nodeCoords,1),1)*[0 0 1];
        drawElements(shift+nodeCoords,cell2mat(boundaryIEN),bndElementType,0,0)
        hold on;
        drawElementsInterp(shift+nodeCoords+u2_E*factor,cell2mat(boundaryIEN),bndElementType,s2(1,:),.5);
        drawNodes(shift+nodeCoords+u2_E*factor,cell2mat(boundaryIEN(3)));
    end

%% ========================= PART 2 (2D)

    function [u2_E,smpGrad,nodeCoords]=go2d(elementSize,smpPts)
            CMatrix_E=properties.CPlaneStrainFull;
%         CMatrix_E=properties.CPlaneStressFull;

        %% ===== FEM Solver Parameters ===== %
        elementType='2dQ1';
        elData=elementData(elementType);
        bndElementType=elData.bndElementType;
        numGP=elData.numGPFull;
        bndElData=elementData(bndElementType);
        numGPbnd=bndElData.numGPFull;
        
        numEq_E=2; % Number of scalar equations (elasticity)
        %% ===== Specifying Body Force and Boundary conditions =====
        
        isDirichlet_E=[0 0; 0 0; ... %bottom and right
            0 0; 1 1]; %top and left
        %%
        % ===== Mesh generation =====
        [nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect2d(domain, elementType,elementSize);
        boundaryIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);
        INE=IENtoINE(IEN);
        
        % misc
        numEl=size(IEN,1); %#ok<NASGU> %total number of elements in the mesh
        numNodes=size(nodeCoords,1); %total number of nodes in the mesh
        numGDoF_E=numNodes*numEq_E; %total number of global degrees of freedom (elasticity)
        
        % ===== Matrix, vector assembly =====
        K_E = formStiffnessMatrix(nodeCoords, IEN, elementType, numGP, CMatrix_E);
        Fb_E = formBodyForceVector(nodeCoords, IEN, elementType, numGP, bodyForce2d_E);
        
        % Applying Boundary Conditions
        [u_prescribed_E, Fs_E, prescribedDoF_E, freeDoF_E]=...
            formBC(nodeCoords,boundaryIEN,bndElementType,numGPbnd,...
            boundaryTraction2d_E,boundaryDisplacement2d_E,isDirichlet_E);
        
        %         Preventing rigid body motifon
        [prescribedDoF_E,freeDoF_E]=preventRigidMotion(nodeCoords,prescribedDoF_E,freeDoF_E);
        F_E=Fb_E+Fs_E;
        
        % Defining the free parts of stiffness
        KK_E=K_E(freeDoF_E,freeDoF_E);
        FF_E=F_E(freeDoF_E)-K_E(freeDoF_E,prescribedDoF_E)*u_prescribed_E(prescribedDoF_E);
        
        % ===== Solution =====
        u_E=zeros(numGDoF_E,1);
        u_E(prescribedDoF_E)=u_prescribed_E(prescribedDoF_E);
        u_E(freeDoF_E)=KK_E\FF_E;
        u2_E=reshape(u_E,numEq_E,numNodes)';

        %% ===== Stress recovery =====
        recVals1=recoveryAvg(u2_E,nodeCoords,IEN,elementType);
        s2=CMatrix_E*recVals1;

        [smpEl,smpNatCoords] = meshGetNatural(nodeCoords,IEN,elementType,INE,smpPts);
        smpGrad=recoveryEvaluateGradients(u2_E,nodeCoords,IEN,elementType,smpEl,smpNatCoords);
        
        %===== Plots ===== %
        drawElements(nodeCoords,IEN,elementType,0,0)
        drawElementsInterp(nodeCoords+u2_E*factor,IEN,elementType,s2(1,:));
    end

%% Main
elementSize={[20,6,3]};
factor=1e5;

for i=1:numel(elementSize)
    figure(1);clf;hold on;
    smpPts=[domain3d(1,:)+domain3d(2,:).*[.051 .5 .5]];
    [u2_3D,smpGrad3D,nodeCoords3]=go3d(elementSize{i},smpPts); %#ok<ASGLU>
    [u2_2D,smpGrad2D,nodeCoords2]=go2d(elementSize{i}(1:2),smpPts(:,1:2)); %#ok<ASGLU>
    xlabel('x');ylabel('y');
    title({'Plane stress','xx-component of stress'})
    axis equal;drawnow;
    view(3,70);
    axis equal;
    disp('Relative discrepancy in displacements:');
    norm(u2_3D(:,1:2)-repmat(u2_2D,elementSize{i}(3)+1,1)) / ...
        norm(repmat(u2_2D,elementSize{i}(3)+1,1)) %#ok<NOPRT>
%     [smpGrad3D([1 2 4 5]), smpGrad2D, smpGrad3D([1 2 4 5])-smpGrad2D]
end
end
