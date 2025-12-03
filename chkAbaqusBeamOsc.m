function chkAbaqusBeamOsc
% Verification with ABAQUS - Vibration of a cantilever beam
% 
% We consider a cantilever beam which is being pulled and then released at
% the free end. This file
%   - solves the problem using our libraries
%   - generates input file for ABAQUS
%   - calls ABAQUS 
%   - extracts output data from ABAQUS and compares to our results
% ABAQUS must be installed to run. 
% Some error may occur due to unknown subtleties of the interaction between
% ABAQUS and file system. Simply running this file again often fixes the
% problem.

%% ===== Problem Parameters ===== %
domain=[0 0; 4 1]; %rectangular domainf
thickness=.1; %thickness of a beam
tfinal=1;

% Material properties
youngsModulus = 193e6;
poissonsRatio = 0.253;
massDensity=7.999e3;
properties=elasticProperties('youngsModulus',youngsModulus,'poissonsRatio',poissonsRatio);
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

%% ============ SIMULATION SUBROUTINES ==============
    function go_CPS3_Implicit(elementSize)
        jobname='job-dyn-CPS3';
        CMatrix=properties.CPlaneStressFull;
        MMatrix=[1 0; 0 1]*massDensity;
        %% ===== FEM Solver Parameters ===== %
        elementType='2dP1';
        gamma=1/2; % Newmark parameter
        beta=1/4; % Newmark parameter
        dt = 3e-02; % Fixed timestep
        % Abaqus parameters:
        AbaqusElementType='CPS3';
        AbaqusMaxIncrements=1000;
        AbaqusAlpha=0; % Abaqus numerical dumping
        AbaqusDt=dt;  % Abaqus manual timestep
        AbaqusTimePeriod=tfinal;

        numGP=1; 
        fprintf('\n Comparing Abaqus %s vs our %s, numGP=%u\n',AbaqusElementType,elementType,numGP)
        elData=elementData(elementType);
        bndElementType=elData.bndElementType;
        bndElData=elementData(bndElementType);
        numGPbnd=bndElData.numGPFull;
        
        numEq=2; % Number of scalar equations (elasticity)
        %% ===== Our Method =====
        % ===== Mesh generation =====
        [nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect2d(domain, elementType,elementSize);
        BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);

        % id of the lower-right corner
        hlID=find(nodeCoords(:,1)==domain(2,1)&nodeCoords(:,2)==domain(1,2));
        
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
        
        outputFunction=@(i,t,u,v,a,fDoF,pDoF,M,C,K,F,u0,v0,a0,up,vp,ap){t,u,v,a};
        
        [output]=solveNewmarkLa(beta,gamma,ts,freeDoF,prescribedDoF,...
            M1,[],K,F,...
            u0,v0,a0,...
            u_prescribed,v_prescribed,a_prescribed,...
            'increment',outputFunction);
        time_my=output{1};
        uM=output{2};
        vM=output{3};
        aM=output{4};
        %% ===== Using ABAQUS =====
        % Generate Abaqus input file
        abaqus_inp_file_imp(nodeCoords,IEN,elementType,boundaryElementIDs,BIEN,jobname,AbaqusElementType,...
            hlID,AbaqusMaxIncrements, AbaqusDt, AbaqusTimePeriod, AbaqusAlpha, beta, gamma);
        % === RUN ABAQUS ===
        cd('abaq');
        cmdstring=['abaqus job=' jobname ' ask_delete=OFF interactive'];
        tic
        fprintf(1,'Starting ABAQUS: "%s" ....',cmdstring);
        [~,~]=system(cmdstring);
        fprintf(1,'done\n');
        toc
        cd('..');
        
        % === Extract ABAQUS OUTPUT ===
        recs=readAbaqFil(['abaq' filesep jobname '.fil']);
        [time_Abaqus,uA,vA,aA,rfA]=parseAbaqRecs_2(recs);
        uA=reshape(permute(uA(:,:,1:2),[1 3 2]),[size(uA,1), size(uA,2)*2]);
        vA=reshape(permute(vA(:,:,1:2),[1 3 2]),[size(vA,1), size(vA,2)*2]);
        aA=reshape(permute(aA(:,:,1:2),[1 3 2]),[size(aA,1), size(aA,2)*2]);
        MA=full(importMTX(['abaq' filesep jobname '_MASS3.mtx']));
        KA=full(importMTX(['abaq' filesep jobname '_STIF3.mtx']));
        %% ===== Compare and Plot =====
        disp('Result comparison:')
        
        fprintf(1,'U, V, A relative discrepency norms: \n %e, %e, %e \n', ...
            norm(uM-uA(1:end-1,:))/norm(uM), norm(vM-vA(1:end-1,:))/norm(vM), norm(aM-aA(1:end-1,:))/norm(aM));
        figure(1);clf;hold on;
        plot(time_my,uM(:,hlID*2),'*-');
        plot(time_Abaqus-1,uA(:,hlID*2),'.-');
        legend('Ours','Abaqus');
        title({'Deflection of the lower-right corner of the beam', ['Our dt=' num2str(dt)], ['Abaqus dt=' num2str(AbaqusDt)]});
        xlabel('time');
        ylabel('u_2');
        fprintf(1,'Timestep relative difference: %f \n' , (dt-AbaqusDt)/dt);
    end %go_CPS3_Implicit

    function go_CPS3_Explicit(elementSize)
        jobname='job-dyn-CPS3_e';
        oldjobname='job-dyn-CPS3';
        CMatrix=properties.CPlaneStressFull;
        MMatrix=[1 0; 0 1]*massDensity;
        %% ===== FEM Solver Parameters ===== %
        elementType='2dP1';
        gamma=1/2; % Newmark parameter
        beta=0; % Newmark parameter
        % timestep will be determined later based on mesh
        % Abaqus parameters:
        AbaqusElementType='CPS3';
        AbaqusDt=3e-2;  % Timestep from the implicit scheme. Used to defene 
        % number of outputs
        AbaqusTimePeriod=tfinal;

        numGP=1; 
        fprintf('\n Comparing Abaqus %s vs our %s, numGP=%u\n',AbaqusElementType,elementType,numGP)
        elData=elementData(elementType);
        bndElementType=elData.bndElementType;
        bndElData=elementData(bndElementType);
        numGPbnd=bndElData.numGPFull;
        
        numEq=2; % Number of scalar equations (elasticity)
        %% ===== Our Method =====
        % ===== Mesh generation =====
        [nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect2d(domain, elementType,elementSize);
        BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);

        % id of the lower-right corner
        hlID=find(nodeCoords(:,1)==domain(2,1)&nodeCoords(:,2)==domain(1,2));
        
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
        
        % p-wave speed
        cp=sqrt( (properties.lambda+2*properties.mu)/massDensity );
        elSizeX=(domain(2,1)-domain(1,1))/elementSize(1);
        elSizeY=(domain(2,2)-domain(1,2))/elementSize(2);
        L=min(elSizeX,elSizeY); % characteristic element length
        dt=0.9*L/cp;
        
        M1=diag(sum(M)); % Lumping mass matrix
        
        ts=0:dt:tfinal;
        outputFunction=@(i,t,u,v,a,fDoF,pDoF,M,C,K,F,u0,v0,a0,up,vp,ap){t,u,v,a};
        [output]=solveNewmarkLa(beta,gamma,ts,freeDoF,prescribedDoF,...
            M1,[],K,F,...
            u0,v0,a0,...
            u_prescribed,v_prescribed,a_prescribed,...
            'increment',outputFunction);
        time_my=output{1};
        uM=output{2};
        vM=output{3};
        aM=output{4};
        %% ===== Using ABAQUS =====
        % Generate Abaqus input file
        abaqus_inp_file_exp(nodeCoords,BIEN,jobname,...
            AbaqusDt, AbaqusTimePeriod, oldjobname);
        % === RUN ABAQUS ===
        cd('abaq');
        cmdstring=['abaqus job=' jobname ' ask_delete=OFF interactive'];
        tic
        fprintf(1,'Starting ABAQUS: "%s" ....',cmdstring);
        [~,~]=system(cmdstring);
        fprintf(1,'done\n');
        toc
        cd('..');
        % === Extract ABAQUS OUTPUT ===
        recs2=readAbaqFil(['abaq' filesep jobname '.fil']);
        [time_Abaqus,uA,vA,aA,rfA,numInc]=parseAbaqRecs_2(recs2);
        AbaqusDt=(time_Abaqus(2)-time_Abaqus(1))/(numInc(2)-numInc(1));
        uA=reshape(permute(uA(:,:,1:2),[1 3 2]),[size(uA,1), size(uA,2)*2]);
        vA=reshape(permute(vA(:,:,1:2),[1 3 2]),[size(vA,1), size(vA,2)*2]);
        aA=reshape(permute(aA(:,:,1:2),[1 3 2]),[size(aA,1), size(aA,2)*2]);
        % ===== Compare and Plot =====
        disp('Result comparison:')
        
        uMi=interp1(ts,uM,time_Abaqus(1:end-1)-1);
        vMi=interp1(ts,vM,time_Abaqus(1:end-1)-1);
        aMi=interp1(ts,aM,time_Abaqus(1:end-1)-1);
        fprintf(1,'U, V, A relative discrepency norms: \n %e, %e, %e \n', ...
            norm(uMi-uA(1:end-1,:))/norm(uMi), norm(vMi-vA(1:end-1,:))/norm(vMi), ...
            norm(aMi-aA(1:end-1,:)))/norm(aMi);
        figure(2);clf;hold on;
        plot(time_my,uM(:,hlID*2),'*-');
        plot(time_Abaqus-1,uA(:,hlID*2),'.-');
        legend('Ours','Abaqus');
        title({'Deflection of the lower-right corner of the beam',...
            ['Our dt=' num2str(dt)], ['Abaqus dt=' num2str(AbaqusDt)]});
        xlabel('time');
        ylabel('u_2');
        fprintf(1,'Timestep relative difference: %f \n' , (dt-AbaqusDt)/dt);
    end %go_CPS3_Explicit

    function go_CPE3_Implicit(elementSize)
        jobname='job-dyn-CPE3';
        CMatrix=properties.CPlaneStrainFull;
        MMatrix=[1 0; 0 1]*massDensity;
        %% ===== FEM Solver Parameters ===== %
        elementType='2dP1';
        gamma=1/2; % Newmark parameter
        beta=1/4; % Newmark parameter
        dt = 3e-02; % Fixed timestep
        % Abaqus parameters:
        AbaqusElementType='CPE3';
        AbaqusMaxIncrements=1000;
        AbaqusAlpha=0; % Abaqus numerical dumping
        AbaqusDt=dt;  % Abaqus manual timestep
        AbaqusTimePeriod=tfinal;

        numGP=1; 
        fprintf('\n Comparing Abaqus %s vs our %s, numGP=%u\n',AbaqusElementType,elementType,numGP)
        elData=elementData(elementType);
        bndElementType=elData.bndElementType;
        bndElData=elementData(bndElementType);
        numGPbnd=bndElData.numGPFull;
        
        numEq=2; % Number of scalar equations (elasticity)
        %% ===== Our Method =====
        % ===== Mesh generation =====
        [nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect2d(domain, elementType,elementSize);
        BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);

        % id of the lower-right corner
        hlID=find(nodeCoords(:,1)==domain(2,1)&nodeCoords(:,2)==domain(1,2));
        
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
        
        outputFunction=@(i,t,u,v,a,fDoF,pDoF,M,C,K,F,u0,v0,a0,up,vp,ap){t,u,v,a};
        
        [output]=solveNewmarkLa(beta,gamma,ts,freeDoF,prescribedDoF,...
            M1,[],K,F,...
            u0,v0,a0,...
            u_prescribed,v_prescribed,a_prescribed,...
            'increment',outputFunction);
        time_my=output{1};
        uM=output{2};
        vM=output{3};
        aM=output{4};
        %% ===== Using ABAQUS =====
        % Generate Abaqus input file
        abaqus_inp_file_imp(nodeCoords,IEN,elementType,boundaryElementIDs,BIEN,jobname,AbaqusElementType,...
            hlID,AbaqusMaxIncrements, AbaqusDt, AbaqusTimePeriod, AbaqusAlpha, beta, gamma);
        % === RUN ABAQUS ===
        cd('abaq');
        cmdstring=['abaqus job=' jobname ' ask_delete=OFF interactive'];
        tic
        fprintf(1,'Starting ABAQUS: "%s" ....',cmdstring);
        [~,~]=system(cmdstring);
        fprintf(1,'done\n');
        toc
        cd('..');
        
        % === Extract ABAQUS OUTPUT ===
        recs=readAbaqFil(['abaq' filesep jobname '.fil']);
        [time_Abaqus,uA,vA,aA,rfA]=parseAbaqRecs_2(recs);
        uA=reshape(permute(uA(:,:,1:2),[1 3 2]),[size(uA,1), size(uA,2)*2]);
        vA=reshape(permute(vA(:,:,1:2),[1 3 2]),[size(vA,1), size(vA,2)*2]);
        aA=reshape(permute(aA(:,:,1:2),[1 3 2]),[size(aA,1), size(aA,2)*2]);
        MA=full(importMTX(['abaq' filesep jobname '_MASS3.mtx']));
        KA=full(importMTX(['abaq' filesep jobname '_STIF3.mtx']));
        %% ===== Compare and Plot =====
        disp('Result comparison:')
        
        fprintf(1,'U, V, A relative discrepency norms: \n %e, %e, %e \n', ...
            norm(uM-uA(1:end-1,:))/norm(uM), norm(vM-vA(1:end-1,:))/norm(vM), norm(aM-aA(1:end-1,:))/norm(aM));
        figure(1);clf;hold on;
        plot(time_my,uM(:,hlID*2),'*-');
        plot(time_Abaqus-1,uA(:,hlID*2),'.-');
        legend('Ours','Abaqus');
        title({'Deflection of the lower-right corner of the beam', ['Our dt=' num2str(dt)], ['Abaqus dt=' num2str(AbaqusDt)]});
        xlabel('time');
        ylabel('u_2');
        fprintf(1,'Timestep relative difference: %f \n' , (dt-AbaqusDt)/dt);
    end %go_CPE3_Implicit

    function go_CPE3_Explicit(elementSize)
        jobname='job-dyn-CPE3_e';
        oldjobname='job-dyn-CPE3';
        CMatrix=properties.CPlaneStrainFull;
        MMatrix=[1 0; 0 1]*massDensity;
        %% ===== FEM Solver Parameters ===== %
        elementType='2dP1';
        gamma=1/2; % Newmark parameter
        beta=0; % Newmark parameter
        % timestep will be determined later based on mesh
        % Abaqus parameters:
        AbaqusElementType='CPE3';
        AbaqusDt=3e-2;  % Timestep from the implicit scheme. Used to defene 
        % number of outputs
        AbaqusTimePeriod=tfinal;

        numGP=1; 
        fprintf('\n Comparing Abaqus %s vs our %s, numGP=%u\n',AbaqusElementType,elementType,numGP)
        elData=elementData(elementType);
        bndElementType=elData.bndElementType;
        bndElData=elementData(bndElementType);
        numGPbnd=bndElData.numGPFull;
        
        numEq=2; % Number of scalar equations (elasticity)
        %% ===== Our Method =====
        % ===== Mesh generation =====
        [nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect2d(domain, elementType,elementSize);
        BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);

        % id of the lower-right corner
        hlID=find(nodeCoords(:,1)==domain(2,1)&nodeCoords(:,2)==domain(1,2));
        
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
        
        % p-wave speed
        cp=sqrt( (properties.lambda+2*properties.mu)/massDensity );
        elSizeX=(domain(2,1)-domain(1,1))/elementSize(1);
        elSizeY=(domain(2,2)-domain(1,2))/elementSize(2);
        L=min(elSizeX,elSizeY); % characteristic element length
        dt=0.8*L/cp;
        
        M1=diag(sum(M)); % Lumping mass matrix
        
        ts=0:dt:tfinal;
        outputFunction=@(i,t,u,v,a,fDoF,pDoF,M,C,K,F,u0,v0,a0,up,vp,ap){t,u,v,a};
        [output]=solveNewmarkLa(beta,gamma,ts,freeDoF,prescribedDoF,...
            M1,[],K,F,...
            u0,v0,a0,...
            u_prescribed,v_prescribed,a_prescribed,...
            'increment',outputFunction);
        time_my=output{1};
        uM=output{2};
        vM=output{3};
        aM=output{4};
        %% ===== Using ABAQUS =====
        % Generate Abaqus input file
        abaqus_inp_file_exp(nodeCoords,BIEN,jobname,...
            AbaqusDt, AbaqusTimePeriod, oldjobname);
        % === RUN ABAQUS ===
        cd('abaq');
        cmdstring=['abaqus job=' jobname ' ask_delete=OFF interactive'];
        tic
        fprintf(1,'Starting ABAQUS: "%s" ....',cmdstring);
        [~,~]=system(cmdstring);
        fprintf(1,'done\n');
        toc
        cd('..');
        % === Extract ABAQUS OUTPUT ===
        recs2=readAbaqFil(['abaq' filesep jobname '.fil']);
        [time_Abaqus,uA,vA,aA,rfA,numInc]=parseAbaqRecs_2(recs2);
        AbaqusDt=(time_Abaqus(2)-time_Abaqus(1))/(numInc(2)-numInc(1));
        uA=reshape(permute(uA(:,:,1:2),[1 3 2]),[size(uA,1), size(uA,2)*2]);
        vA=reshape(permute(vA(:,:,1:2),[1 3 2]),[size(vA,1), size(vA,2)*2]);
        aA=reshape(permute(aA(:,:,1:2),[1 3 2]),[size(aA,1), size(aA,2)*2]);
        % ===== Compare and Plot =====
        disp('Result comparison:')
        
        uMi=interp1(ts,uM,time_Abaqus(1:end-1)-1);
        vMi=interp1(ts,vM,time_Abaqus(1:end-1)-1);
        aMi=interp1(ts,aM,time_Abaqus(1:end-1)-1);
        fprintf(1,'U, V, A relative discrepency norms: \n %e, %e, %e \n', ...
            norm(uMi-uA(1:end-1,:))/norm(uMi), norm(vMi-vA(1:end-1,:))/norm(vMi), ...
            norm(aMi-aA(1:end-1,:)))/norm(aMi);
        figure(2);clf;hold on;
        plot(time_my,uM(:,hlID*2),'*-');
        plot(time_Abaqus-1,uA(:,hlID*2),'.-');
        legend('Ours','Abaqus');
        title({'Deflection of the lower-right corner of the beam',...
            ['Our dt=' num2str(dt)], ['Abaqus dt=' num2str(AbaqusDt)]});
        xlabel('time');
        ylabel('u_2');

        fprintf(1,'Timestep relative difference: %f \n' , (dt-AbaqusDt)/dt);
    end %go_CPE3_Explicit

    function go_CPS4_Implicit(elementSize)
        jobname='job-dyn-CPS4';
        CMatrix=properties.CPlaneStressFull;
        MMatrix=[1 0; 0 1]*massDensity;
        %% ===== FEM Solver Parameters ===== %
        elementType='2dQ1';
        gamma=1/2; % Newmark parameter
        beta=1/4; % Newmark parameter
        dt = 3e-02; % Fixed timestep
        % Abaqus parameters:
        AbaqusElementType='CPS4';
        AbaqusMaxIncrements=1000;
        AbaqusAlpha=0; % Abaqus numerical dumping
        AbaqusDt=dt;  % Abaqus manual timestep
        AbaqusTimePeriod=tfinal;

        numGP=4; 
        fprintf('\n Comparing Abaqus %s vs our %s, numGP=%u\n',AbaqusElementType,elementType,numGP)
        elData=elementData(elementType);
        bndElementType=elData.bndElementType;
        bndElData=elementData(bndElementType);
        numGPbnd=bndElData.numGPFull;
        
        numEq=2; % Number of scalar equations (elasticity)
        %% ===== Our Method =====
        % ===== Mesh generation =====
        [nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect2d(domain, elementType,elementSize);
        BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);

        % id of the lower-right corner
        hlID=find(nodeCoords(:,1)==domain(2,1)&nodeCoords(:,2)==domain(1,2));
        
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
        
        outputFunction=@(i,t,u,v,a,fDoF,pDoF,M,C,K,F,u0,v0,a0,up,vp,ap){t,u,v,a};
        
        [output]=solveNewmarkLa(beta,gamma,ts,freeDoF,prescribedDoF,...
            M1,[],K,F,...
            u0,v0,a0,...
            u_prescribed,v_prescribed,a_prescribed,...
            'increment',outputFunction);
        time_my=output{1};
        uM=output{2};
        vM=output{3};
        aM=output{4};
        %% ===== Using ABAQUS =====
        % Generate Abaqus input file
        abaqus_inp_file_imp(nodeCoords,IEN,elementType,boundaryElementIDs,BIEN,jobname,AbaqusElementType,...
            hlID,AbaqusMaxIncrements, AbaqusDt, AbaqusTimePeriod, AbaqusAlpha, beta, gamma);
        % === RUN ABAQUS ===
        cd('abaq');
        cmdstring=['abaqus job=' jobname ' ask_delete=OFF interactive'];
        tic
        fprintf(1,'Starting ABAQUS: "%s" ....',cmdstring);
        [~,~]=system(cmdstring);
        fprintf(1,'done\n');
        toc
        cd('..');
        
        % === Extract ABAQUS OUTPUT ===
        recs=readAbaqFil(['abaq' filesep jobname '.fil']);
        [time_Abaqus,uA,vA,aA,rfA]=parseAbaqRecs_2(recs);
        uA=reshape(permute(uA(:,:,1:2),[1 3 2]),[size(uA,1), size(uA,2)*2]);
        vA=reshape(permute(vA(:,:,1:2),[1 3 2]),[size(vA,1), size(vA,2)*2]);
        aA=reshape(permute(aA(:,:,1:2),[1 3 2]),[size(aA,1), size(aA,2)*2]);
        MA=full(importMTX(['abaq' filesep jobname '_MASS3.mtx']));
        KA=full(importMTX(['abaq' filesep jobname '_STIF3.mtx']));
        %% ===== Compare and Plot =====
        disp('Result comparison:')
        
        fprintf(1,'U, V, A relative discrepency norms: \n %e, %e, %e \n', ...
            norm(uM-uA(1:end-1,:))/norm(uM), norm(vM-vA(1:end-1,:))/norm(vM), norm(aM-aA(1:end-1,:))/norm(aM));
        figure(1);clf;hold on;
        plot(time_my,uM(:,hlID*2),'*-');
        plot(time_Abaqus-1,uA(:,hlID*2),'.-');
        legend('Ours','Abaqus');
        title({'Deflection of the lower-right corner of the beam', ['Our dt=' num2str(dt)], ['Abaqus dt=' num2str(AbaqusDt)]});
        xlabel('time');
        ylabel('u_2');
        fprintf(1,'Timestep relative difference: %f \n' , (dt-AbaqusDt)/dt);
    end %go_CPS4_Implicit

    function go_CPE4R_Implicit(elementSize)
        jobname='job-dyn-CPE4R';
        CMatrix=properties.CPlaneStrainFull;
        MMatrix=[1 0; 0 1]*massDensity;
        %% ===== FEM Solver Parameters ===== %
        elementType='2dQ1';
        gamma=1/2; % Newmark parameter
        beta=1/4; % Newmark parameter
        dt = 3e-02; % Fixed timestep
        % Abaqus parameters:
        AbaqusElementType='CPE4R';
        AbaqusMaxIncrements=1000;
        AbaqusAlpha=0; % Abaqus numerical dumping
        AbaqusDt=dt;  % Abaqus manual timestep
        AbaqusTimePeriod=tfinal;

        numGP=1; 
        fprintf('\n Comparing Abaqus %s vs our %s, numGP=%u\n',AbaqusElementType,elementType,numGP)
        elData=elementData(elementType);
        bndElementType=elData.bndElementType;
        bndElData=elementData(bndElementType);
        numGPbnd=bndElData.numGPFull;
        
        numEq=2; % Number of scalar equations (elasticity)
        %% ===== Our Method =====
        % ===== Mesh generation =====
        [nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect2d(domain, elementType,elementSize);
        BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);

        % id of the lower-right corner
        hlID=find(nodeCoords(:,1)==domain(2,1)&nodeCoords(:,2)==domain(1,2));
        
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
        
        outputFunction=@(i,t,u,v,a,fDoF,pDoF,M,C,K,F,u0,v0,a0,up,vp,ap){t,u,v,a};
        
        [output]=solveNewmarkLa(beta,gamma,ts,freeDoF,prescribedDoF,...
            M1,[],K,F,...
            u0,v0,a0,...
            u_prescribed,v_prescribed,a_prescribed,...
            'increment',outputFunction);
        time_my=output{1};
        uM=output{2};
        vM=output{3};
        aM=output{4};
        %% ===== Using ABAQUS =====
        % Generate Abaqus input file
        abaqus_inp_file_imp(nodeCoords,IEN,elementType,boundaryElementIDs,BIEN,jobname,AbaqusElementType,...
            hlID,AbaqusMaxIncrements, AbaqusDt, AbaqusTimePeriod, AbaqusAlpha, beta, gamma);
        % === RUN ABAQUS ===
        cd('abaq');
        cmdstring=['abaqus job=' jobname ' ask_delete=OFF interactive'];
        tic
        fprintf(1,'Starting ABAQUS: "%s" ....',cmdstring);
        [~,~]=system(cmdstring);
        fprintf(1,'done\n');
        toc
        cd('..');
        
        % === Extract ABAQUS OUTPUT ===
        recs=readAbaqFil(['abaq' filesep jobname '.fil']);
        [time_Abaqus,uA,vA,aA,rfA]=parseAbaqRecs_2(recs);
        uA=reshape(permute(uA(:,:,1:2),[1 3 2]),[size(uA,1), size(uA,2)*2]);
        vA=reshape(permute(vA(:,:,1:2),[1 3 2]),[size(vA,1), size(vA,2)*2]);
        aA=reshape(permute(aA(:,:,1:2),[1 3 2]),[size(aA,1), size(aA,2)*2]);
        MA=full(importMTX(['abaq' filesep jobname '_MASS3.mtx']));
        KA=full(importMTX(['abaq' filesep jobname '_STIF3.mtx']));
        %% ===== Compare and Plot =====
        disp('Result comparison:')
        
        fprintf(1,'U, V, A relative discrepency norms: \n %e, %e, %e \n', ...
            norm(uM-uA(1:end-1,:))/norm(uM), norm(vM-vA(1:end-1,:))/norm(vM), norm(aM-aA(1:end-1,:))/norm(aM));
        figure(1);clf;hold on;
        plot(time_my,uM(:,hlID*2),'*-');
        plot(time_Abaqus-1,uA(:,hlID*2),'.-');
        legend('Ours','Abaqus');
        title({'Deflection of the lower-right corner of the beam', ['Our dt=' num2str(dt)], ['Abaqus dt=' num2str(AbaqusDt)]});
        xlabel('time');
        ylabel('u_2');
        fprintf(1,'Timestep relative difference: %f \n' , (dt-AbaqusDt)/dt);
    end %go_CPE4R_Implicit

    function go_CPE4R_Explicit(elementSize)
        jobname='job-dyn-CPE4R_e';
        oldjobname='job-dyn-CPE4R';
        CMatrix=properties.CPlaneStrainFull;
        MMatrix=[1 0; 0 1]*massDensity;
        %% ===== FEM Solver Parameters ===== %
        elementType='2dQ1';
        gamma=1/2; % Newmark parameter
        beta=0; % Newmark parameter
        % timestep will be determined later based on mesh
        % Abaqus parameters:
        AbaqusElementType='CPE4R';
        AbaqusDt=3e-2;  % Timestep from the implicit scheme. Used to defene 
        % number of outputs
        AbaqusTimePeriod=tfinal;

        numGP=1; 
        fprintf('\n Comparing Abaqus %s vs our %s, numGP=%u\n',AbaqusElementType,elementType,numGP)
        elData=elementData(elementType);
        bndElementType=elData.bndElementType;
        bndElData=elementData(bndElementType);
        numGPbnd=bndElData.numGPFull;
        
        numEq=2; % Number of scalar equations (elasticity)
        %% ===== Our Method =====
        % ===== Mesh generation =====
        [nodeCoords, IEN,boundaryElementIDs, boundaryNodeLocalIDs] = meshRect2d(domain, elementType,elementSize);
        BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs);

        % id of the lower-right corner
        hlID=find(nodeCoords(:,1)==domain(2,1)&nodeCoords(:,2)==domain(1,2));
        
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
        
        % p-wave speed
        cp=sqrt( (properties.lambda+2*properties.mu)/massDensity );
        elSizeX=(domain(2,1)-domain(1,1))/elementSize(1);
        elSizeY=(domain(2,2)-domain(1,2))/elementSize(2);
        L=min(elSizeX,elSizeY); % characteristic element length
        dt=0.8*L/cp;
        
        M1=diag(sum(M)); % Lumping mass matrix
        
        ts=0:dt:tfinal;
        outputFunction=@(i,t,u,v,a,fDoF,pDoF,M,C,K,F,u0,v0,a0,up,vp,ap){t,u,v,a};
        [output]=solveNewmarkLa(beta,gamma,ts,freeDoF,prescribedDoF,...
            M1,[],K,F,...
            u0,v0,a0,...
            u_prescribed,v_prescribed,a_prescribed,...
            'increment',outputFunction);
        time_my=output{1};
        uM=output{2};
        vM=output{3};
        aM=output{4};
        %% ===== Using ABAQUS =====
        % Generate Abaqus input file
        abaqus_inp_file_exp(nodeCoords,BIEN,jobname,...
            AbaqusDt, AbaqusTimePeriod, oldjobname);
        % === RUN ABAQUS ===
        cd('abaq');
        cmdstring=['abaqus job=' jobname ' ask_delete=OFF interactive'];
        tic
        fprintf(1,'Starting ABAQUS: "%s" ....',cmdstring);
        [~,~]=system(cmdstring);
        fprintf(1,'done\n');
        toc
        cd('..');
        % === Extract ABAQUS OUTPUT ===
        recs2=readAbaqFil(['abaq' filesep jobname '.fil']);
        [time_Abaqus,uA,vA,aA,rfA,numInc]=parseAbaqRecs_2(recs2);
        AbaqusDt=(time_Abaqus(2)-time_Abaqus(1))/(numInc(2)-numInc(1));
        uA=reshape(permute(uA(:,:,1:2),[1 3 2]),[size(uA,1), size(uA,2)*2]);
        vA=reshape(permute(vA(:,:,1:2),[1 3 2]),[size(vA,1), size(vA,2)*2]);
        aA=reshape(permute(aA(:,:,1:2),[1 3 2]),[size(aA,1), size(aA,2)*2]);
        % ===== Compare and Plot =====
        disp('Result comparison:')
        
        uMi=interp1(ts,uM,time_Abaqus(1:end-1)-1);
        vMi=interp1(ts,vM,time_Abaqus(1:end-1)-1);
        aMi=interp1(ts,aM,time_Abaqus(1:end-1)-1);
        fprintf(1,'U, V, A relative discrepency norms: \n %e, %e, %e \n', ...
            norm(uMi-uA(1:end-1,:))/norm(uMi), norm(vMi-vA(1:end-1,:))/norm(vMi), ...
            norm(aMi-aA(1:end-1,:)))/norm(aMi);
        figure(2);clf;hold on;
        plot(time_my,uM(:,hlID*2),'*-');
        plot(time_Abaqus-1,uA(:,hlID*2),'.-');
        legend('Ours','Abaqus');
        title({'Deflection of the lower-right corner of the beam',...
            ['Our dt=' num2str(dt)], ['Abaqus dt=' num2str(AbaqusDt)]});
        xlabel('time');
        ylabel('u_2');

        fprintf(1,'Timestep relative difference: %f \n' , (dt-AbaqusDt)/dt);
    end %go_CPE4R_Explicit

%% ========= MAIN =============
elementSize={[13,6,7]};

for i=1:numel(elementSize)
    go_CPS3_Implicit(elementSize{i}(1:2)); %#ok<ASGLU>
    go_CPS3_Explicit(elementSize{i}(1:2))
    go_CPE3_Implicit(elementSize{i}(1:2)); %#ok<ASGLU>
    go_CPE3_Explicit(elementSize{i}(1:2))
    go_CPS4_Implicit(elementSize{i}(1:2)); %#ok<ASGLU>
%     No explicit counterpart, since CPS4 element is not available for
%     import in Abaqus/Explicit
    go_CPE4R_Implicit(elementSize{i}(1:2)); %#ok<ASGLU>
    go_CPE4R_Explicit(elementSize{i}(1:2))
end

%% ======================= MISC FUNCTIONS ==============

    function abaqus_inp_file_imp(nodeCoords,IEN,elementType,boundaryElementIDs,BIEN,jobname,abaqusElementType,...
            hlID,maxIncrementNumberAbaqus,maxTimeIncrementAbaqus, timePeriodAbaqus, alpha, beta, gamma)
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
        fprintf(fid,'*Nset, nset=tip \n');
        fprintf(fid,'%u\n',hlID);
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
            else %new line
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
            '%f, %f\n' ...
            '*Density\n' ...
            '%f\n'], youngsModulus, poissonsRatio, massDensity);
        fprintf(fid,[...
            '** -------------------------\n' ...
            '*Section Controls, name=EC-1, hourglass=STIFFNESS\n' ...
            '1e-7, 0., 0.\n' ... %1st paramater is artificial hourglassing stiffness
            ... % it has no effect in fully integrated elements, e.g. CPS
            ... % the larger it is the greater the discrepancy bw stiffness matrices is
            ... % However, very small values prevent convergence
            '** -------------------------\n']);
        
        fprintf(fid,[...
            '*Step, name=BeamLoad, nlgeom=NO\n' ...
            'applyTractionToBeam\n' ...
            '*Static\n' ...
            '1., 1., 1e-05, 1.\n' ...
            '*Boundary\n' ...
            'leftSide, 1, 3\n' ...
            '*Dsload, follower=NO, constant resultant=YES\n']);
        if any(strcmpi(elementType,{'2DQ1','2DP1'}))
            if norm(traction1)>eps
                fprintf(fid,'rightSide, TRVEC, %f, %f, %f\n',...
                    norm(traction1),traction1(1)/norm(traction1),traction1(2)/norm(traction1));
            else
                fprintf(fid,'rightSide, TRVEC, %f, %f, %f\n',...
                    0,1,0);
            end
        elseif any(strcmpi(elementType,{'3DQ1','3DP1'}))
            if norm(traction1)>eps
                fprintf(fid,'rightSide, TRVEC, %f, %f, %f, %f\n',...
                    norm(traction1)/thickness,traction1(1)/norm(traction1),traction1(2)/norm(traction1),0);
            else
                fprintf(fid,'rightSide, TRVEC, %f, %f, %f, %f\n',...
                    0,1,0,0);
            end
        end
        fprintf(fid,['*NODE FILE, FREQUENCY=1, NSET=Assembly.beam-1-1.all\n',...
            'U, V, A, RF\n']);
        fprintf(fid,[...
            '*Restart, write, FREQUENCY=1\n']);
        fprintf(fid,[...
            '*End Step\n']);
        fprintf(fid,['*STEP, name=BeamRelease, INC=%u, nlgeom=NO \n',...
            'releaseTractionAndSeeWhatHappens\n'],...
            maxIncrementNumberAbaqus);
        fprintf(fid,['*DYNAMIC, NOHAF, DIRECT, ALPHA=%f, GAMMA=%f, BETA=%f, INITIAL=NO \n',...
            '%e, %e \n'],...
            alpha,gamma,beta,maxTimeIncrementAbaqus, timePeriodAbaqus);
        fprintf(fid,['*Boundary, OP=new\n' ...
            'leftSide, 1, 3\n' ...
            '*Dsload, follower=NO, OP=new, constant resultant=YES\n']);
        if any(strcmpi(elementType,{'2DQ1','2DP1'}))
            if norm(traction2)>eps
                fprintf(fid,'rightSide, TRVEC, %f, %f, %f\n',...
                    norm(traction2),traction2(1)/norm(traction2),traction2(2)/norm(traction2));
            else
                fprintf(fid,'rightSide, TRVEC, %f, %f, %f\n',...
                    0,1,0);
            end
        elseif any(strcmpi(elementType,{'3DQ1','3DP1'}))
            if norm(traction2)>eps
                fprintf(fid,'rightSide, TRVEC, %f, %f, %f, %f\n',...
                    norm(traction2)/thickness,traction2(1)/norm(traction2),traction2(2)/norm(traction2),0);
            else
                fprintf(fid,'rightSide, TRVEC, %f, %f, %f, %f\n',...
                    0,1,0,0);
            end
        end
        fprintf(fid,['*NODE FILE, FREQUENCY=1, NSET=Assembly.beam-1-1.all\n',...
            'U, V, A, RF\n',...
            '*END STEP\n']);
        fprintf(fid,[...
            '** -------------------------\n' ...
            '*STEP, name=exportmatrix\n' ...
            '*MATRIX GENERATE, MASS, STIFFNESS\n' ...
            '*MATRIX OUTPUT, MASS, STIFFNESS, FORMAT=MATRIX INPUT\n' ...
            '*END STEP\n' ...
            '** -------------------------            \n' ...
            '']);
        fprintf(fid,'\n');
        
        fclose(fid);
    end

    function abaqus_inp_file_exp(nodeCoords,BIEN,jobname,...
            timeIncrementAbaqus, timePeriodAbaqus, oldjob)
        if ~exist('abaq','dir')
            mkdir('abaq');
        end
        fid=fopen(['abaq' filesep jobname '.inp'],'w');
        disp(' ');
        disp(['Abaqus input file "abaq' filesep jobname '.inp" created.']);
        fprintf(fid,['*Heading \n' ...
            'cantileverBeam \n']);
        fprintf(fid,[...
            '*Assembly, name=Assembly\n', ...
            '*Instance, Library=%s, instance=beam-1-1\n', ...
            '*IMPORT, UPDATE=NO, STEP=1\n',...
            '*End Instance\n'],oldjob);
        fprintf(fid,[...
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
            else %new line
                k=0;
                fprintf(fid,'%u \n ',BIEN{lsid}(ii));
            end
        end
        fprintf(fid,'%u \n',BIEN{lsid}(end));
        fprintf(fid,[...
            '*End Assembly\n']);
        
        fprintf(fid,['*STEP, name=BeamRelease, nlgeom=NO \n',...
            'releaseTractionAndSeeWhatHappens\n']);
        fprintf(fid,['*DYNAMIC, EXPLICIT, FIXED TIME INCREMENTATION\n',...
            ', %e \n'],...
            timePeriodAbaqus);
        fprintf(fid,['*Boundary, OP=new\n' ...
            'leftSide, 1, 3\n']);
        
        fprintf(fid,[...
            '*FILE OUTPUT, NUMBER INTERVAL=%u\n',...
            '*NODE FILE, NSET=Assembly.beam-1-1.all\n',...
            'U, V, A, RF\n',...
            '*END STEP\n'],ceil(timePeriodAbaqus/timeIncrementAbaqus));
        fclose(fid);
    end

end
