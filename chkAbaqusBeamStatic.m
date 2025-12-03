function chkAbaqusBeamStatic%%%% Mechanical equilibrium 3D and 2D - bending cantilever beam
% Verification with ABAQUS - Vibration of a cantilever beam
% 
% We consider a cantilever beam which is being pulled at the free end.
% Stiffness matrices as well as the resulting solution are compared.
% This file 
%   - solves the problem using our librarie
%   - generates input file for ABAQUS
%   - calls ABAQUS 
%   - extracts output data from ABAQUS and compares to our results
% ABAQUS must be installed to run. 
% Some error may occur due to unknown subtleties in the interaction between
% ABAQUS and file system. Simply running this file again often fixes the
% problem.

doPlot=1; % if we want to see plots
factor=1e5; %displacement aplification factor for plots

%% ===== Problem Parameters ===== %
domain=[0 0; 4 1]; %rectangular domain
thickness=.1;
% Material properties
% youngsModulus = 3e8%193e6;
% poissonsRatio = 0.1% 0.253;
youngsModulus = 193e6;
poissonsRatio = 0.253;
properties=elasticProperties('youngsModulus',youngsModulus,'poissonsRatio',poissonsRatio);
%% ====== BC
bodyForce2d_E=@(x,y)([0*x 0*x]'); % Elastic part
bodyForce2d_E=@(x)(bodyForce2d_E(x(:,1),x(:,2)));

traction=[0 2];
boundaryTraction2d_1=@(x,y)([0*x 0*x]');
boundaryTraction2d_1=@(x)boundaryTraction2d_1(x(:,1),x(:,2));

boundaryTraction2d_2=@(x,y)([(traction(1)+0*x) (traction(2)+0*x)]');
boundaryTraction2d_2=@(x)boundaryTraction2d_2(x(:,1),x(:,2));

boundaryTraction2d_3=@(x,y)([0*x 0*x]');
boundaryTraction2d_3=@(x)boundaryTraction2d_3(x(:,1),x(:,2));

boundaryTraction2d_4=@(x,y)([0*x 0*x]');
boundaryTraction2d_4=@(x)boundaryTraction2d_4(x(:,1),x(:,2));

boundaryTraction2d_E={boundaryTraction2d_1,boundaryTraction2d_2,...
    boundaryTraction2d_3,boundaryTraction2d_4};

boundaryDisplacement2d_E=@(x,y)([0*x 0*x]'); % x,y - k-by-1, output must be numEq-by-1 or numEq-by-k
boundaryDisplacement2d_E=@(x)(boundaryDisplacement2d_E(x(:,1),x(:,2)));

%% ============ SIMULATION SUBROUTINES ==============
    function go_C3D8R(elementSize)
        % Compute and compare stiffness matrices computed using Abaqus C3D8R (1-GP
        % 8-node brick element) and our routines.
        jobname='job-C3D8R';
        %% ===== Problem Parameters ===== %
        domain3d=[domain(1,:) 0;  domain(2,:) thickness]; %rectangular domain
        CMatrix_E=properties.CFull;
        %% ===== FEM Solver Parameters ===== %
        
        %uncomment any of the below lines
        elementType='3dQ1';
        abaqusElementType='C3D8R';
        numGP=1;
        fprintf('\n Comparing Abaqus %s vs our %s, numGP=%u \n',abaqusElementType,elementType,numGP)
        elData=elementData(elementType);
        bndElementType=elData.bndElementType;
        bndElData=elementData(bndElementType);
        numGPbnd=bndElData.numGPFull;
        
        numEq_E=3; % Number of scalar equations (elasticity)
        %% ===== Specifying Body Force and Boundary conditions =====
        bodyForce3d_E=@(x,y,z)([bodyForce2d_E([x y])' 0*z]'); % Elastic part
        bodyForce3d_E=@(x)(bodyForce3d_E(x(:,1),x(:,2),x(:,3)));
        
        boundaryTraction3d_1=@(x,y,z)([boundaryTraction2d_E{1}([x y])' 0*z]');
        boundaryTraction3d_1=@(x)boundaryTraction3d_1(x(:,1),x(:,2),x(:,3))/thickness;
        
        boundaryTraction3d_2=@(x,y,z)([boundaryTraction2d_E{3}([x y])' 0*z]');
        boundaryTraction3d_2=@(x)boundaryTraction3d_2(x(:,1),x(:,2),x(:,3))/thickness;
        
        boundaryTraction3d_3=@(x,y,z)([boundaryTraction2d_E{4}([x y])' 0*z]');
        boundaryTraction3d_3=@(x)boundaryTraction3d_3(x(:,1),x(:,2),x(:,3))/thickness;
        
        boundaryTraction3d_4=@(x,y,z)([boundaryTraction2d_E{2}([x y])' 0*z]');
        boundaryTraction3d_4=@(x)boundaryTraction3d_4(x(:,1),x(:,2),x(:,3))/thickness;
        
        boundaryTraction3d_5=@(x,y,z)([0*x 0*x 0*x]');
        boundaryTraction3d_5=@(x)boundaryTraction3d_5(x(:,1),x(:,2),x(:,3))/thickness;
        
        boundaryTraction3d_6=@(x,y,z)([0*x 0*x 0*x]');
        boundaryTraction3d_6=@(x)boundaryTraction3d_6(x(:,1),x(:,2),x(:,3))/thickness;
        
        boundaryTraction3d_E={boundaryTraction3d_1,boundaryTraction3d_2,...
            boundaryTraction3d_3,boundaryTraction3d_4,...
            boundaryTraction3d_5,boundaryTraction3d_6};
        
        isDirichlet_E=[0 0 0; 0 0 0; ... %bottom and top
            1 1 1; 0 0 0; .... %left and right
            0 0 0; 0 0 0]; %front and back
        
        boundaryDisplacement3d_E=@(x,y,z)([boundaryDisplacement2d_E([x y])' 0*x]'); % x,y - k-by-1, output must be numEq-by-1 or numEq-by-k
        boundaryDisplacement3d_E=@(x)(boundaryDisplacement3d_E(x(:,1),x(:,2),x(:,3)));
        %% ===== Our Method =====
        % ===== Mesh generation =====
        [nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect3d(domain3d, elementType,elementSize);
        BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);
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
            formBC(nodeCoords,BIEN,bndElementType,numGPbnd,...
            boundaryTraction3d_E,boundaryDisplacement3d_E,isDirichlet_E);
        
        F_E=Fb_E+Fs_E;
        
        % Defining the free parts of stiffness
        KK_E=K_E(freeDoF_E,freeDoF_E);
        FF_E=F_E(freeDoF_E)-K_E(freeDoF_E,prescribedDoF_E)*u_prescribed_E(prescribedDoF_E);
        
        %% ===== ABAQUS =====
        abaqus_inp_file(nodeCoords,IEN,elementType,boundaryElementIDs,BIEN,jobname,abaqusElementType);
        % === RUN ABAQUS ===
        cd('abaq');
        cmdstring=['abaqus job=' jobname ' ask_delete=OFF interactive'];
        tic
        disp('Running ABAQUS....')
        [~,~]=system(cmdstring);
        toc
        cd('..');
        KA=full(importMTX(['abaq' filesep jobname '_STIF2.mtx']));
        
        %% ===== Compare =====
        disp('Absolute discrepancy bw ABAQUS and our stiffness matrices is')
        disp(norm(K_E(freeDoF_E,freeDoF_E) - KA(freeDoF_E,freeDoF_E),'fro'))
        disp('Relative discrepancy bw ABAQUS and our stiffness matrices is')
        disp(norm(K_E(freeDoF_E,freeDoF_E) - KA(freeDoF_E,freeDoF_E),'fro')/norm(K_E(freeDoF_E,freeDoF_E),'fro'))
        
        %% ===== Plot =====
    end %go_C3D8R

    function go_C3D4(elementSize)
        % Compute and compare stiffness matrices computed using Abaqus C3D8R (1-GP
        % 8-node brick element) and our routines.
        jobname='job-C3D4';
        %% ===== Problem Parameters ===== %
        domain3d=[domain(1,:) 0;  domain(2,:) thickness]; %rectangular domain
        CMatrix_E=properties.CFull;
        %% ===== FEM Solver Parameters ===== %
        
        %uncomment any of the below lines
        elementType='3dP1';
        abaqusElementType='C3D4';
        numGP=1;
        fprintf('\n Comparing Abaqus %s vs our %s, numGP=%u\n',abaqusElementType,elementType,numGP)
        elData=elementData(elementType);
        bndElementType=elData.bndElementType;
        bndElData=elementData(bndElementType);
        numGPbnd=bndElData.numGPFull;
        
        numEq_E=3; % Number of scalar equations (elasticity)
        %% ===== Specifying Body Force and Boundary conditions =====
        bodyForce3d_E=@(x,y,z)([bodyForce2d_E([x y])' 0*z]'); % Elastic part
        bodyForce3d_E=@(x)(bodyForce3d_E(x(:,1),x(:,2),x(:,3)));
        
        boundaryTraction3d_1=@(x,y,z)([boundaryTraction2d_E{1}([x y])' 0*z]');
        boundaryTraction3d_1=@(x)boundaryTraction3d_1(x(:,1),x(:,2),x(:,3))/thickness;
        
        boundaryTraction3d_2=@(x,y,z)([boundaryTraction2d_E{3}([x y])' 0*z]');
        boundaryTraction3d_2=@(x)boundaryTraction3d_2(x(:,1),x(:,2),x(:,3))/thickness;
        
        boundaryTraction3d_3=@(x,y,z)([boundaryTraction2d_E{4}([x y])' 0*z]');
        boundaryTraction3d_3=@(x)boundaryTraction3d_3(x(:,1),x(:,2),x(:,3))/thickness;
        
        boundaryTraction3d_4=@(x,y,z)([boundaryTraction2d_E{2}([x y])' 0*z]');
        boundaryTraction3d_4=@(x)boundaryTraction3d_4(x(:,1),x(:,2),x(:,3))/thickness;
        
        boundaryTraction3d_5=@(x,y,z)([0*x 0*x 0*x]');
        boundaryTraction3d_5=@(x)boundaryTraction3d_5(x(:,1),x(:,2),x(:,3))/thickness;
        
        boundaryTraction3d_6=@(x,y,z)([0*x 0*x 0*x]');
        boundaryTraction3d_6=@(x)boundaryTraction3d_6(x(:,1),x(:,2),x(:,3))/thickness;
        
        boundaryTraction3d_E={boundaryTraction3d_1,boundaryTraction3d_2,...
            boundaryTraction3d_3,boundaryTraction3d_4,...
            boundaryTraction3d_5,boundaryTraction3d_6};
        
        isDirichlet_E=[0 0 0; 0 0 0; ... %bottom and top
            1 1 1; 0 0 0; .... %left and right
            0 0 0; 0 0 0]; %front and back
        
        boundaryDisplacement3d_E=@(x,y,z)([boundaryDisplacement2d_E([x y])' 0*x]'); % x,y - k-by-1, output must be numEq-by-1 or numEq-by-k
        boundaryDisplacement3d_E=@(x)(boundaryDisplacement3d_E(x(:,1),x(:,2),x(:,3)));
        %% ===== Our Method =====
        % ===== Mesh generation =====
        [nodeCoords, IEN, boundaryElementIDs, boundaryNodeLocalIDs] = meshRect3d(domain3d, elementType,elementSize);
        IEN=IEN(:,[1 3 2 4]); % this is important because Abaqus uses different node numbering and a negative volume would be obtained otherwise
        BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);
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
            formBC(nodeCoords,BIEN,bndElementType,numGPbnd,...
            boundaryTraction3d_E,boundaryDisplacement3d_E,isDirichlet_E);
        
        F_E=Fb_E+Fs_E;
        
        % Defining the free parts of stiffness
        KK_E=K_E(freeDoF_E,freeDoF_E);
        FF_E=F_E(freeDoF_E)-K_E(freeDoF_E,prescribedDoF_E)*u_prescribed_E(prescribedDoF_E);
        
        % ===== SOLUTION =====
        u=zeros(numGDoF_E,1);
        u(freeDoF_E)=KK_E\FF_E;
        u(prescribedDoF_E)=u_prescribed_E(prescribedDoF_E);
        u2=reshape(u,numEq_E,numNodes)';
        
        %% ===== ABAQUS =====
        abaqus_inp_file(nodeCoords,IEN,elementType,boundaryElementIDs,BIEN,jobname,abaqusElementType);
        % === RUN ABAQUS ===
        cd('abaq');
        cmdstring=['abaqus job=' jobname ' ask_delete=OFF interactive'];
        tic
        disp('Running ABAQUS....')
        [~,~]=system(cmdstring);
        toc
        cd('..');
        
        % === GET ABAQUS OUTPUT ===
        KA=full(importMTX(['abaq' filesep jobname '_STIF2.mtx']));
        uA2=import_displacements(['abaq\' jobname '.dat']);
        uA=uA2';
        uA=uA(:);
        
        %% ===== Compare =====
        disp('Absolute discrepancy bw ABAQUS and our stiffness matrices is')
        disp(norm(K_E(freeDoF_E,freeDoF_E) - KA(freeDoF_E,freeDoF_E),'fro'))
        disp('Relative discrepancy bw ABAQUS and our stiffness matrices is')
        disp(norm(K_E(freeDoF_E,freeDoF_E) - KA(freeDoF_E,freeDoF_E),'fro')/norm(K_E(freeDoF_E,freeDoF_E),'fro'))
        
        if norm(u)>0
            disp('Relative discrepancy bw ABAQUS and our displacements is')
            disp( norm(u-uA)/norm(u) );
        end
        
        %% ===== Plot =====
        if doPlot
            figure;
            subplot(2,1,1);
            drawElements(nodeCoords+factor*u2,{IEN,boundaryElementIDs},elementType);
            title(['Our Method, ' elementType ', numGP=' num2str(numGP)]); view(2);
            subplot(2,1,2);
            drawElements(nodeCoords+factor*uA2,{IEN,boundaryElementIDs},elementType);
            title(['Abaqus, ' abaqusElementType]); view(2);
        end
    end %go_C3D4

    function go_CPS4(elementSize)
        % Compute and compare stiffness matrices computed using Abaqus CPS4
        % (4-node plane stress quadrilateral) and our routines.
        jobname='job-CPS4';
        %% ===== Problem Parameters ===== %
        CMatrix_E=properties.CPlaneStressFull;
        %% ===== FEM Solver Parameters ===== %
        elementType='2dQ1';
        abaqusElementType='CPS4';
        numGP=4;
        fprintf('\n Comparing Abaqus %s vs our %s, numGP=%u\n',abaqusElementType,elementType,numGP)
        elData=elementData(elementType);
        bndElementType=elData.bndElementType;
        bndElData=elementData(bndElementType);
        numGPbnd=bndElData.numGPFull;
        
        numEq_E=2; % Number of scalar equations (elasticity)
        %% ===== Specifying Body Force and Boundary conditions =====
        isDirichlet_E=[0 0; 0 0; ... %bottom and right
            0 0; 1 1]; %top and left
        %% ===== Our Method =====
        % ===== Mesh generation =====
        [nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect2d(domain, elementType,elementSize);
        BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);
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
            formBC(nodeCoords,BIEN,bndElementType,numGPbnd,...
            boundaryTraction2d_E,boundaryDisplacement2d_E,isDirichlet_E);
        
        F_E=Fb_E+Fs_E;
        
        % Defining the free parts of stiffness
        KK_E=K_E(freeDoF_E,freeDoF_E);
        FF_E=F_E(freeDoF_E)-K_E(freeDoF_E,prescribedDoF_E)*u_prescribed_E(prescribedDoF_E);
        
        % ===== SOLUTION =====
        u=zeros(numGDoF_E,1);
        u(freeDoF_E)=KK_E\FF_E;
        u(prescribedDoF_E)=u_prescribed_E(prescribedDoF_E);
        u2=reshape(u,numEq_E,numNodes)';
        
        %% ===== ABAQUS =====
        abaqus_inp_file(nodeCoords,IEN,elementType,boundaryElementIDs,BIEN,jobname,abaqusElementType);
        % === RUN ABAQUS ===
        cd('abaq');
        cmdstring=['abaqus job=' jobname ' ask_delete=OFF interactive'];
        tic
        disp('Running ABAQUS....')
        [~,~]=system(cmdstring);
        toc
        cd('..');
        
        % === GET ABAQUS OUTPUT ===
        KA=full(importMTX(['abaq' filesep jobname '_STIF2.mtx']));
        uA2=import_displacements(['abaq\' jobname '.dat']);
        uA=uA2';
        uA=uA(:);
        
        %% ===== Compare =====
        disp('Absolute discrepancy bw ABAQUS and our stiffness matrices is')
        disp(norm(K_E(freeDoF_E,freeDoF_E) - KA(freeDoF_E,freeDoF_E),'fro'))
        disp('Relative discrepancy bw ABAQUS and our stiffness matrices is')
        disp(norm(K_E(freeDoF_E,freeDoF_E) - KA(freeDoF_E,freeDoF_E),'fro')/norm(K_E(freeDoF_E,freeDoF_E),'fro'))
        
        if norm(u)>0
            disp('Relative discrepancy bw ABAQUS and our displacements is')
            disp( norm(u-uA)/norm(u) );
        end
        
        %% ===== Plot =====
        if doPlot
            figure;
            subplot(2,1,1);
            drawElements(nodeCoords+factor*u2,IEN,elementType);
            title(['Our Method, ' elementType ', numGP=' num2str(numGP)]); view(2);
            subplot(2,1,2);
            drawElements(nodeCoords+factor*uA2,IEN,elementType);
            title(['Abaqus, ' abaqusElementType]); view(2);
        end
    end %go_CPS4

    function go_CPE4R(elementSize)
        % Compute and compare stiffness matrices computed using Abaqus CPS4
        % (4-node plane stress quadrilateral) and our routines.
        jobname='job-CPE4R';
        CMatrix_E=properties.CPlaneStrainFull;
        %% ===== FEM Solver Parameters ===== %
        abaqusElementType='CPE4R';
        elementType='2dQ1';
        numGP=1;
        fprintf('\n Comparing Abaqus %s vs our %s, numGP=%u\n',abaqusElementType,elementType,numGP)
        elData=elementData(elementType);
        bndElementType=elData.bndElementType;
        bndElData=elementData(bndElementType);
        numGPbnd=bndElData.numGPFull;
        
        numEq_E=2; % Number of scalar equations (elasticity)
        %% ===== Specifying Body Force and Boundary conditions =====
        isDirichlet_E=[0 0; 0 0; ... %bottom and right
            0 0; 1 1]; %top and left
        %% ===== Our Method =====
        % ===== Mesh generation =====
        [nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect2d(domain, elementType,elementSize);
        BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);
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
            formBC(nodeCoords,BIEN,bndElementType,numGPbnd,...
            boundaryTraction2d_E,boundaryDisplacement2d_E,isDirichlet_E);
        
        F_E=Fb_E+Fs_E;
        
        % Defining the free parts of stiffness
        KK_E=K_E(freeDoF_E,freeDoF_E);
        FF_E=F_E(freeDoF_E)-K_E(freeDoF_E,prescribedDoF_E)*u_prescribed_E(prescribedDoF_E);
        
        %% ===== ABAQUS =====
        abaqus_inp_file(nodeCoords,IEN,elementType,boundaryElementIDs,BIEN,jobname,abaqusElementType);
        % === RUN ABAQUS ===
        cd('abaq');
        cmdstring=['abaqus job=' jobname ' ask_delete=OFF interactive'];
        tic
        disp('Running ABAQUS....')
        [~,~]=system(cmdstring);
        toc
        cd('..');
        KA=full(importMTX(['abaq' filesep jobname '_STIF2.mtx']));
        
        %% ===== Compare =====
        disp('Absolute discrepancy bw ABAQUS and our stiffness matrices is')
        disp(norm(K_E(freeDoF_E,freeDoF_E) - KA(freeDoF_E,freeDoF_E),'fro'))
        disp('Relative discrepancy bw ABAQUS and our stiffness matrices is')
        disp(norm(K_E(freeDoF_E,freeDoF_E) - KA(freeDoF_E,freeDoF_E),'fro')/norm(K_E(freeDoF_E,freeDoF_E),'fro'))
        
        %% ===== Plot =====
    end %go_CPE4R

    function go_CPE3(elementSize)
        % Compute and compare stiffness matrices computed using Abaqus CPS4
        % (4-node plane stress quadrilateral) and our routines.
        jobname='job-CPE3';
        CMatrix_E=properties.CPlaneStrainFull;
        %% ===== FEM Solver Parameters ===== %
        elementType='2dP1';
        abaqusElementType='CPE3';
        numGP=1;
        fprintf('\n Comparing Abaqus %s vs our %s, numGP=%u\n',abaqusElementType,elementType,numGP)
        elData=elementData(elementType);
        bndElementType=elData.bndElementType;
        bndElData=elementData(bndElementType);
        numGPbnd=bndElData.numGPFull;
        
        numEq_E=2; % Number of scalar equations (elasticity)
        %% ===== Specifying Body Force and Boundary conditions =====
        isDirichlet_E=[0 0; 0 0; ... %bottom and right
            0 0; 1 1]; %top and left
        
        %% ===== Our Method =====
        % ===== Mesh generation =====
        [nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect2d(domain, elementType,elementSize);
        BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);
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
            formBC(nodeCoords,BIEN,bndElementType,numGPbnd,...
            boundaryTraction2d_E,boundaryDisplacement2d_E,isDirichlet_E);
        
        F_E=Fb_E+Fs_E;
        
        % Defining the free parts of stiffness
        KK_E=K_E(freeDoF_E,freeDoF_E);
        FF_E=F_E(freeDoF_E)-K_E(freeDoF_E,prescribedDoF_E)*u_prescribed_E(prescribedDoF_E);
        
        % ===== SOLUTION =====
        u=zeros(numGDoF_E,1);
        u(freeDoF_E)=KK_E\FF_E;
        u(prescribedDoF_E)=u_prescribed_E(prescribedDoF_E);
        u2=reshape(u,numEq_E,numNodes)';
        
        %% ===== ABAQUS =====
        abaqus_inp_file(nodeCoords,IEN,elementType,boundaryElementIDs,BIEN,jobname,abaqusElementType);
        % === RUN ABAQUS ===
        cd('abaq');
        cmdstring=['abaqus job=' jobname ' ask_delete=OFF interactive'];
        tic
        disp('Running ABAQUS....')
        [~,~]=system(cmdstring);
        toc
        cd('..');
        
        % === GET ABAQUS OUTPUT ===
        KA=full(importMTX(['abaq' filesep jobname '_STIF2.mtx']));
        uA2=import_displacements(['abaq' filesep jobname '.dat']);
        uA=uA2';
        uA=uA(:);
        
        %% ===== Compare =====
        disp('Absolute discrepancy bw ABAQUS and our stiffness matrices is')
        disp(norm(K_E(freeDoF_E,freeDoF_E) - KA(freeDoF_E,freeDoF_E),'fro'))
        disp('Relative discrepancy bw ABAQUS and our stiffness matrices is')
        disp(norm(K_E(freeDoF_E,freeDoF_E) - KA(freeDoF_E,freeDoF_E),'fro')/norm(K_E(freeDoF_E,freeDoF_E),'fro'))
        
        if norm(u)>0
            disp('Relative discrepancy bw ABAQUS and our displacements is')
            disp( norm(u-uA)/norm(u) );
        end
        
        %% ===== Plot =====
        if doPlot
            figure;
            subplot(2,1,1);
            drawElements(nodeCoords+factor*u2,IEN,elementType);
            title(['Our Method, ' elementType ', numGP=' num2str(numGP)]); view(2);
            subplot(2,1,2);
            drawElements(nodeCoords+factor*uA2,IEN,elementType);
            title(['Abaqus, ' abaqusElementType]); view(2);
        end
    end %go_CPE3

    function go_CPS3(elementSize)
        % Compute and compare stiffness matrices computed using Abaqus CPS4
        % (4-node plane stress quadrilateral) and our routines.
        jobname='job-CPS3';
        CMatrix_E=properties.CPlaneStressFull;
        %% ===== FEM Solver Parameters ===== %
        elementType='2dP1';
        abaqusElementType='CPS3';
        numGP=1;
        fprintf('\n Comparing Abaqus %s vs our %s, numGP=%u\n',abaqusElementType,elementType,numGP)
        elData=elementData(elementType);
        bndElementType=elData.bndElementType;
        bndElData=elementData(bndElementType);
        numGPbnd=bndElData.numGPFull;
        
        numEq_E=2; % Number of scalar equations (elasticity)
        %% ===== Specifying Body Force and Boundary conditions =====
        isDirichlet_E=[0 0; 0 0; ... %bottom and right
            0 0; 1 1]; %top and left
        
        %% ===== Our Method =====
        % ===== Mesh generation =====
        [nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect2d(domain, elementType,elementSize);
        BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);
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
            formBC(nodeCoords,BIEN,bndElementType,numGPbnd,...
            boundaryTraction2d_E,boundaryDisplacement2d_E,isDirichlet_E);
        
        F_E=Fb_E+Fs_E;
        
        % Defining the free parts of stiffness
        KK_E=K_E(freeDoF_E,freeDoF_E);
        FF_E=F_E(freeDoF_E)-K_E(freeDoF_E,prescribedDoF_E)*u_prescribed_E(prescribedDoF_E);
        
        % ===== SOLUTION =====
        u=zeros(numGDoF_E,1);
        u(freeDoF_E)=KK_E\FF_E;
        u(prescribedDoF_E)=u_prescribed_E(prescribedDoF_E);
        u2=reshape(u,numEq_E,numNodes)';
        
        %% ===== ABAQUS =====
        abaqus_inp_file(nodeCoords,IEN,elementType,boundaryElementIDs,BIEN,jobname,abaqusElementType);
        % === RUN ABAQUS ===
        cd('abaq');
        cmdstring=['abaqus job=' jobname ' ask_delete=OFF interactive'];
        tic
        disp('Running ABAQUS....')
        [~,~]=system(cmdstring);
        toc
        cd('..');
        
        % === GET ABAQUS OUTPUT ===
        KA=full(importMTX(['abaq' filesep jobname '_STIF2.mtx']));
        uA2=import_displacements(['abaq\' jobname '.dat']);
        uA=uA2';
        uA=uA(:);
        
        %% ===== Compare =====
        disp('Absolute discrepancy bw ABAQUS and our stiffness matrices is')
        disp(norm(K_E(freeDoF_E,freeDoF_E) - KA(freeDoF_E,freeDoF_E),'fro'))
        disp('Relative discrepancy bw ABAQUS and our stiffness matrices is')
        disp(norm(K_E(freeDoF_E,freeDoF_E) - KA(freeDoF_E,freeDoF_E),'fro')/norm(K_E(freeDoF_E,freeDoF_E),'fro'))
        
        if norm(u)>0
            disp('Relative discrepancy bw ABAQUS and our displacements is')
            disp( norm(u-uA)/norm(u) );
        end
        
        %% ===== Plot =====
        if doPlot
            figure;
            subplot(2,1,1);
            drawElements(nodeCoords+factor*u2,IEN,elementType);
            title(['Our Method, ' elementType ', numGP=' num2str(numGP)]); view(2);
            subplot(2,1,2);
            drawElements(nodeCoords+factor*uA2,IEN,elementType);
            title(['Abaqus, ' abaqusElementType]); view(2);
        end
    end %go_CPS3

%% ========= MAIN =============
elementSize={[13,6,7]};

for i=1:numel(elementSize)
    go_CPS3(elementSize{i}(1:2)); %#ok<ASGLU>
    go_CPE3(elementSize{i}(1:2)); %#ok<ASGLU>
    go_CPS4(elementSize{i}(1:2)); %#ok<ASGLU>
    go_CPE4R(elementSize{i}(1:2)); %#ok<ASGLU>
    go_C3D8R(elementSize{i}); %#ok<ASGLU>
    go_C3D4(elementSize{i}); %#ok<ASGLU>
    
end

%% ======================= MISC FUNCTIONS ==============

    function abaqus_inp_file(nodeCoords,IEN,elementType,boundaryElementIDs,BIEN,jobname,abaqusElementType)
        if size(nodeCoords,1)>1000
            warning('The mesh contains more than 1000 nodes. \nAbaqus student edition does not support more than 1000 \nFurther execution may lead to an error message.')
        end
        % generates an .inp file for abaqus
        if ~exist('abaq','dir')
            mkdir('abaq');
        end
        fid=fopen(['abaq' filesep jobname '.inp'],'w');
        disp(' ');
        disp(['Abaqus input file "abaq' filesep jobname '.inp" created.']);
        fprintf(fid,['*Heading \n' ...
            'cantileverBeam \n' ...
            '*Preprint, echo=NO, model=NO, history=NO, contact=NO\n' ...
            '*Part, name=beam-1\n']);
        
        fprintf(fid,'*Node\n');
        if size(nodeCoords,2)==2
            for ijk=1:size(nodeCoords,1)
                fprintf(fid,'%u, \t \t %20.15f, \t%20.15f \n',ijk,nodeCoords(ijk,1),nodeCoords(ijk,2));
            end
            if strcmpi(abaqusElementType,'CPS4')
                fprintf(fid,'\n*Element, type=CPS4');
            elseif strcmpi(abaqusElementType,'CPE4R')
                fprintf(fid,'\n*Element, type=CPE4R');
            elseif strcmpi(abaqusElementType,'CPE3')
                fprintf(fid,'\n*Element, type=CPE3');
            elseif strcmpi(abaqusElementType,'CPS3')
                fprintf(fid,'\n*Element, type=CPS3');
            end
        elseif size(nodeCoords,2)==3
            for ijk=1:size(nodeCoords,1)
                fprintf(fid,'%u, \t \t %20.15f, \t%20.15f, \t%20.15f \n',ijk,nodeCoords(ijk,1),nodeCoords(ijk,2),nodeCoords(ijk,3));
            end
            if strcmpi(abaqusElementType,'C3D8R')
                fprintf(fid,'\n \n *Element, type=C3D8R');
            elseif strcmpi(abaqusElementType,'C3D4')
                fprintf(fid,'\n \n *Element, type=C3D4');
            end
        end
        
        for ijk=1:size(IEN,1)
            fprintf(fid,'\n %u',ijk);
            for ii=1:size(IEN,2)
                fprintf(fid,', \t%u',IEN(ijk,ii));
            end
        end
        fprintf(fid,'\n');
        fprintf(fid,'*Nset, nset=all, generate \n');
        fprintf(fid,'%u, %u, %u\n',1,size(nodeCoords,1),1);
        fprintf(fid,'*Elset, elset=all123, generate \n');
        fprintf(fid,'%u, %u, %u\n',1,size(IEN,1),1);
        
        fprintf(fid,['*Solid Section, elset=all123, material=Steel-1, CONTROLS=EC-1 \n' ...
            ', \n' ...
            '*End Part\n' ...
            '*Assembly, name=Assembly\n' ...
            '*Instance, name=beam-1-1, part=beam-1\n' ...
            '*End Instance\n' ...
            '*Nset, nset=leftSide, instance=beam-1-1\n']);
        if size(nodeCoords,2)==2
            lsid=4;% for 2d rectangular domain, left side is in BIEN{4}
        elseif size(nodeCoords,2)==3
            lsid=3;% for 3d rectangular domain, left side is in BIEN{3}
        end
        k=0;
        for ii=1:numel(BIEN{lsid})-1
            k=k+1;
            if k<=15
                fprintf(fid,'%u, ',BIEN{lsid}(ii));
            else
                k=0;
                fprintf(fid,'%u \n ',BIEN{lsid}(ii));
            end
        end
        fprintf(fid,'%u \n',BIEN{lsid}(end));
        fprintf(fid,[...
            '*Elset, elset=_rightSide_SS, internal, instance=beam-1-1\n']);
        
        if size(nodeCoords,2)==2
            rsid=2;% for 2d quadrilateral, right side is in BIEN{2}
        elseif size(nodeCoords,2)==3
            rsid=4;% for 3d hexahedron, right side is in BIEN{4}
        end
        k=0;
        for ii=1:numel(boundaryElementIDs{rsid})-1
            k=k+1;
            if k<=15
                fprintf(fid,'%u, ',boundaryElementIDs{rsid}(ii));
            else
                k=0;
                fprintf(fid,'%u \n ',boundaryElementIDs{rsid}(ii));
            end
        end
        fprintf(fid,'%u \n',boundaryElementIDs{rsid}(end));
        fprintf(fid,['*Surface, type=ELEMENT, name=rightSide\n']);
        
        if strcmpi(elementType,'2DQ1')
            fprintf(fid,[...
                '_rightSide_SS, S2\n' ]);
        elseif strcmpi(elementType,'2DP1')
            fprintf(fid,[...
                '_rightSide_SS, S2\n' ]);
        elseif strcmpi(elementType,'3DQ1')
            fprintf(fid,[...
                '_rightSide_SS, S4\n' ]);
        elseif strcmpi(elementType,'3DP1')
            fprintf(fid,[...
                '_rightSide_SS, S3\n' ]);
        end
        fprintf(fid,[...
            '*End Assembly\n' ...
            '*Material, name=Steel-1\n' ...
            '*Elastic\n' ...
            '%f, %f\n'], youngsModulus, poissonsRatio);
        fprintf(fid,[...
            '** -------------------------\n' ...
            '*Section Controls, name=EC-1, hourglass=STIFFNESS\n' ...
            '1e-7, 0., 0.\n' ... %1st paramater is artificial hourglassing stiffness
            ... % it has no effect in fully integrated elements, e.g. CPS
            ... % the larger it is the greater the discrepancy bw stiffness matrices is
            ... % However, very small values prevent convergence
            '** -------------------------\n' ...
            '*Step, name=BeamLoad, nlgeom=NO\n' ...
            'applyTractionToBeam\n' ...
            '*Static\n' ...
            '1., 1., 1e-05, 1.\n' ...
            '*Boundary\n' ...
            'leftSide, 1, 3\n' ...
            '*Dsload, follower=NO, constant resultant=YES\n']);
        if any(strcmpi(elementType,{'2DQ1','2DP1'}))
            if norm(traction)>eps
                fprintf(fid,'rightSide, TRVEC, %f, %f, %f\n',...
                    norm(traction),traction(1)/norm(traction),traction(2)/norm(traction));
            else
                fprintf(fid,'rightSide, TRVEC, %f, %f, %f\n',...
                    0,1,0);
            end
        elseif any(strcmpi(elementType,{'3DQ1','3DP1'}))
            if norm(traction)>eps
                fprintf(fid,'rightSide, TRVEC, %f, %f, %f, %f\n',...
                    norm(traction)/thickness,traction(1)/norm(traction),traction(2)/norm(traction),0);
            else
                fprintf(fid,'rightSide, TRVEC, %f, %f, %f, %f\n',...
                    0,1,0,0);
            end
        end
        fprintf(fid,[...
            '*Restart, write, frequency=0\n' ...
            '*Output, field\n' ...
            '*Node Output\n' ...
            'CF, RT, U\n' ...
            '*Element Output, directions=YES\n' ...
            'E, S\n' ...
            '*Node Print\n'...
            'U\n'...
            '** -------------------------\n' ...
            '**ELEMENT MATRIX OUTPUT, STIFFNESS=YES, ELSET=Assembly.beam-1-1.all123\n' ...
            '** -------------------------\n' ...
            '*Output, history, frequency=0\n' ...
            '*EL PRINT, POSITION=INTEGRATION POINTS, FREQUENCY=1\n' ...
            'COORD, IVOL\n' ...
            '*End Step\n' ...
            '** -------------------------\n' ...
            '*STEP, name=exportmatrix\n' ...
            '*MATRIX GENERATE, STIFFNESS\n' ...
            '*MATRIX OUTPUT, STIFFNESS, FORMAT=MATRIX INPUT\n' ...
            '*END STEP\n' ...
            '** -------------------------            \n' ...
            '']);
        
        fprintf(fid,'\n');
        fclose(fid);
    end

end
