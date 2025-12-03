function [prescribedDoF,freeDoF]=preventRigidMotion(nodeCoords,prescribedDoF,freeDoF)
% Append the list of prescribed DoFs if necessary to remove rigid body
% motions in elasticity problems.
%
% INPUT:
% nodeCoords - #Nodes-by-#Dim array, coordinates of nodes.
% pDoF - list of IDs of prescribed (fixed) degrees of freedom.
% fDoF - list of IDs of free degrees of freedom.
%
% OUTPUT:
% pDoF - updated list of IDs of prescribed (fixed) degrees of
%   freedom. 
% fDoF - updated  list of IDs of free degrees of freedom.
%
% COMMENTS:
% If pDoF is initially empty, then no essential constraints are introduced.
% This is not guaranteed otherwise.
% 
% Example:
% [pDoF,fDoF]=preventRigidMotion(nodeCoords,pDoF,fDoF)

numDim=size(nodeCoords,2);
numGDoF=max([prescribedDoF(:);freeDoF(:)]);
numNodes=numGDoF/numDim;
if numDim==2
    [~,ll]=max(nodeCoords*[-1 -1]');
    [~,hh]=max(nodeCoords*[1 1]');
    if numel(prescribedDoF)<=2
        tofix=...%[1 2 numGDoF]'; % candidate DoFs to prescribed
            node2DoFs([ll ll hh]',numDim,[1 2 2]');
        tofix=setdiff(tofix,prescribedDoF); % DoFs to be additionally prescribed;
        prescribedDoF=[prescribedDoF; tofix];
        freeDoF=setdiff(freeDoF,prescribedDoF);
        disp(['To prevent rigid body motion we additionally fixed DoFs: ' num2str(tofix')]);
    elseif numel(prescribedDoF)>=3
        % Check if translations are possible and prevent them by fixig
        % a DoF of the first node.
        xDoFs=node2DoFs(1:numNodes,numDim,1);
        yDoFs=node2DoFs(1:numNodes,numDim,2);
        % This piece of code must be revised:
        if numel(setdiff(xDoFs,prescribedDoF))==numel(xDoFs)
            % no horizontal DoFs prescribed
            tofix=[1];
            prescribedDoF=[prescribedDoF; tofix];
            freeDoF=setdiff(freeDoF,prescribedDoF);
            disp(['To prevent rigid body motion we additionally fixed DoFs: ' num2str(tofix')]);
        elseif numel(setdiff(yDoFs,prescribedDoF))==numel(yDoFs)
            % no vertical DoFs prescribed
            tofix=[2];
            prescribedDoF=[prescribedDoF; tofix];
            freeDoF=setdiff(freeDoF,prescribedDoF);
            disp(['To prevent rigid body motion we additionally fixed DoFs: ' num2str(tofix')]);
        end
    end
elseif numDim==3
    [~,lll]=max(nodeCoords*[-1 -1 -1]');
    [~,hll]=max(nodeCoords*[1 -1 -1]');
    [~,llh]=max(nodeCoords*[-1 -1 1]');
    tofix=node2DoFs([lll lll lll hll hll llh],numDim,[1 2 3 2 3 2]); % candidate DoFs to prescribed
    if numel(prescribedDoF)<=5
        tofix=setdiff(tofix,prescribedDoF); % DoFs to be additionally prescribed;
        prescribedDoF=[prescribedDoF; tofix];
        freeDoF=setdiff(freeDoF,prescribedDoF);
        disp(['To prevent rigid body motion we additionally fixed DoFs: ' num2str(tofix')]);
    elseif numel(prescribedDoF)>=6
        % Check if translations are possible and prevent them by fixig
        % a DoF of the first node.
        xDoFs=node2DoFs(1:numNodes,numDim,1);
        yDoFs=node2DoFs(1:numNodes,numDim,2);
        zDoFs=node2DoFs(1:numNodes,numDim,3);
        % This piece of code must be revised:
        if numel(setdiff(xDoFs,prescribedDoF))==numel(xDoFs)
            % no x DoFs prescribed
            tofix=tofix(1);
            prescribedDoF=[prescribedDoF; tofix];
            freeDoF=setdiff(freeDoF,prescribedDoF);
            disp(['To prevent rigid body motion we additionally fixed DoFs: ' num2str(tofix')]);
        elseif numel(setdiff(yDoFs,prescribedDoF))==numel(yDoFs)
            % no y DoFs prescribed
            tofix=tofix(2);
            prescribedDoF=[prescribedDoF; tofix];
            freeDoF=setdiff(freeDoF,prescribedDoF);
            disp(['To prevent rigid body motion we additionally fixed DoFs: ' num2str(tofix')]);
        elseif numel(setdiff(zDoFs,prescribedDoF))==numel(zDoFs)
            % no z DoFs prescribed
            tofix=tofix(3);
            prescribedDoF=[prescribedDoF; tofix];
            freeDoF=setdiff(freeDoF,prescribedDoF);
            disp(['To prevent rigid body motion we additionally fixed DoFs: ' num2str(tofix')]);
        end
    end
    
end
