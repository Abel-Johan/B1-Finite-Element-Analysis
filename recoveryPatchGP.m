function gradu_n=recoveryPatchGP(gradu,GPCoords,nodeCoords,IEN,elementType,BIEN,varargin)
% Recover gradients at nodes from gradient values given at Gauss points
% using least square fitting (superconvergent patch recovery, SPR)
%
% INPUT:
% gradu - #Eq*#Dim-by-#El*#GP matrix containing nodal values of gradients,
%   so that gradu((i-1)*#Dim+j,(m-1)*#GP+l) corresponds to du_i/dx_j
%   evaluated at l-th Gauss point of m-th element.
% GPCoords - #El*#GP-by-#Dim matrix containing coordinates of Gauss points,
%   so that GPCoords((m-1)*#GP+l,:) are the coordinates of l-th Gauss point
%   of m-th element.
% nodeCoords - #Nodes-by-#Dim matrix, coordinates of nodes.
% IEN - #El-by-#NpE array, element-node incidence matrix. 
% elementType - single string, the type of elements.
% BIEN - a 1-by-#BR cell array containing boundary incidence matrices.
% INE (optional) - #Nodes-by-#EpN array, node-element incidence matrix, row
%   INE(i,:) corresponds to ith node and contains IDs of elements that
%   share it.
% degrees (optional) - #Nodes-by-1 matrix containing degrees (valencies) of
%   nodes, i.e. the number of elements in the mesh each node belongs to.
%
% OUTPUT:
% gradu_n - #Eq*#Dim-by-#Nodes matrix containing nodal values of gradients
%   recovered using the procedure, so that column gradu_n(:,i)
%   corresponds to the gradient evaluated at i-th node.
%
% COMMENTS:
% INE and degrees should be provided together. If not provided, they are
% computed using IENtoINE
% This function calls elementData.m in order to determine natural nodal
% coordinates, optimal number of Gauss points, the polynomial that can be
% fit exactly by the element shape functions and topology of element
% boundary.


%% parsing input
if numel(varargin)==2
    INE=varargin{1};
    degrees=varargin{2};
else
    [INE, degrees]=IENtoINE(IEN);
end

elDat=elementData(elementType);
numGP=elDat.numGPSPR;
p0=elDat.p0;
vertexLocalIDs=elDat.vertexLocalIDs;
nodeNatCoords=elDat.nodeNatCoords;
elLocBIEN=elDat.elLocBIEN;

BIEN=cell2mat(BIEN(:));

% numEq=size(u2,2);
numNodes=size(nodeCoords,1);
numDim=size(nodeCoords,2);
numEl=size(IEN,1);
numBndEl=size(BIEN,1);

if size(gradu,2)/numEl ~= numGP
    error(['Invalid numGP. Expected numGP = ' num2str(numGP)]);
end

%%
polParams=numel(p0(nodeNatCoords(1,:)));

% An important preliminary procedure is to split the domain into patches,
% such that polynomial fitting can be performed later within each patch.
% The resulting patches consist of elements adjacent to a chosen node. This
% heuristic  algorithm seems to work for many meshes, but fail for some
% (e.g. a strip of 2d quadrilaterals).

% A node is CORNER if it is not mid-side or mid-face for at least one elmt.
% Each corner node is associated with a patch, consisting of all elements,
% it is incident to (this is not always one-to-one). A node is BOUNDARY if
% it belongs to the boundary of the domain/mesh and a node is INTERIOR
% otherwise. A node is VALID if its patch have enough GPs to fit a
% polynomial of appropriate degree. A node is MASTER if it is CORNER,
% INTERIOR and VALID. We do least square fitting only in MASTER patches
% ignoring all other. Thus we rely on that the whole domain can covered by
% the master nodes' patches.

% It can happen however that there are (e.g. triangular or tetrahedral) elements
% whose every corner node is at the boundary, therefore they are not
% incident to any master node, i.e. ORPHAN. We fix this problem by
% attaching an orphan elements to an non-orphan element it shares most of
% its nodes with. We rely on on the assumption this is possible. 

% Within each patch:
% a node is OUTER if it is at the boundary of the patch, otherwise a node
% is INNER. A node is RECIPIENT if it is INNER or BOUNDARY. The value at
% recipient nodes is (partly) determined by interpolation within the
% current patch. It is expected that each node is recipient within at least
% one valid patch.

iscorner=zeros(numNodes,1);
for i=1:numEl %finding all corner nodes
    iscorner(IEN(i,vertexLocalIDs))=1;
end
isboundary=zeros(numNodes,1);
for i=1:numBndEl %finding all boundary nodes
    isboundary(BIEN(i,:))=1;
end
isvalid=(degrees*numGP)>=polParams;
ismaster=iscorner.*(1-isboundary).*isvalid;
masterNodes=find(ismaster);

%% deal with orphans
% construct master-elements incidence cell array
IME=cell(numel(masterNodes),1);
maxPatchElements=0;
for i=1:numel(masterNodes)
    IME{i}=INE(masterNodes(i),:);
    maxPatchElements=max(maxPatchElements,numel(IME{i}));
end
IMEcell=IME;
for i=1:numel(masterNodes)
    IME{i}=[IME{i} zeros(1,maxPatchElements-numel(IME{i}))];
end
IME=cell2mat(IME);
IEM=IENtoINE(IME);

orphans=setdiff(1:numEl,IME);
isorphan=zeros(numEl,1);
isorphan(orphans)=1;
for i=1:numel(orphans)
    orphanNodes=IEN(orphans(i),:);
    orphanNeighbours=INE(orphanNodes,:);
    [c,~,ic]=unique(orphanNeighbours(:));
    % remove zero
    iz=find(c==0,1);
    if ~isempty(iz)
        c(iz)=[];
        ic(ic==iz)=[];
        ic(ic>iz)=ic(ic>iz)-1;
    end
    commonCount=accumarray(ic,1); %most popular neighbour element
    [~,pref]=sort(commonCount,'descend'); %sort by popularity
    chosen=find(1-isorphan(c(pref)), 1); %choose most popular non-orphan
    if ~isempty(chosen)
        sibling=c(pref(chosen));
        adoptingPatches=nonzeros(IEM(sibling,:));
        for j=numel(adoptingPatches)
            IMEcell{adoptingPatches(j)}=[IMEcell{adoptingPatches(j)} orphans(i)];
        end
        isorphan(orphans(i))=0;
        IEM(orphans(i),:)=IEM(sibling,:); %my parents are your parents
    else
        warning('SPR failed to solve the orphan problem');
    end
end
%%
IEGP=reshape(1:numGP*numEl,numGP,numEl)'; %incidence Element --> GP
%%
gradu_n=zeros(size(gradu,1),numNodes);
cnt=zeros(1,numNodes);
for i=1:numel(masterNodes)
    patchElements=IMEcell{i};
    patchNodes=unique(IEN(patchElements,:));
    
    patchGPIDs=IEGP(patchElements,:);
    patchGPs=GPCoords(patchGPIDs,:);
    masterCoords=nodeCoords(masterNodes(i),:);

    % patch OUTER and INNER nodes
    % Outer nodes are the nodes of [element boundaries that do not repeat]
    patchElBndElIDs=repmat(patchElements,size(elLocBIEN,1),1);
    patchElBndElIDs=patchElBndElIDs(:);
    patchElBndElLocIDs=repmat(elLocBIEN,numel(patchElements),1);
    patchBIEN=IENtoBIEN(IEN,{patchElBndElIDs},{patchElBndElLocIDs});
    patchBIEN=sort(patchBIEN{1},2);
    [C,~,ic]=unique(patchBIEN,'rows');
    outerElCount=accumarray(ic,1);
    outerNodes=unique(C(ic(outerElCount(ic)==1),:)); %the tricky part
    innerNodes=setdiff(patchNodes,outerNodes);
    
    % patch BOUNDARY nodes
    patchBndNodes=patchNodes(isboundary(patchNodes)==1);
    
    % patch RECIPIENT nodes
    recipientNodes=union(innerNodes,patchBndNodes);
    recipientCoords=nodeCoords(recipientNodes,:);
    
    % least square fitting
    for j=1:size(gradu,1)
        sampleVals=gradu(j,patchGPIDs)';
        %least square fitting
        p=@(x)p0(x-repmat(masterCoords,size(x,1),1));
        pp=p(patchGPs);
        a=(pp'*pp)\(pp'*sampleVals);
        f=@(x)p(x)*a;%fitted function
        gradu_n(j,recipientNodes)=gradu_n(j,recipientNodes)+...
            f(recipientCoords)';
    end
    cnt(recipientNodes)=cnt(recipientNodes)+1;
end

if nnz(cnt==0)>0
    error('SPR failed to assign values to every nodes')
end
gradu_n=gradu_n./repmat(cnt,size(gradu,1),1);
end %function