 function [u_prescribed, Fs, prescribedDoF, freeDoF]=...
    formBC(nodeCoords,BIEN,elementType,numGP,...
    bndTractions,bndDisplacements,isDirichlets,varargin)
% Process boundary conditions and return prescribed displacements, force
% vector, and degrees of freedom lists.
%
% INPUT:
% nodeCoords - #Nodes-by-#Dim array, coordinates of nodes.
% BIEN - a single boundary incidence matrix or a #BR-by-1 cell array of such
%   matrices, where #BR is the number of boundary regions.
% elementType - a single string or #BR-by-1 cell array of strings, which
%   define element types for each boundary region. 
% numGP - a single number or a #BR-by-1 matrix or cell array, number of
%   Gauss points per element.
% boundaryTraction - a single function handle or a #BR-by-1 cell array of
%   function handles. Functions must take a k-by-#Dim matrix as input
%   (coordinates) and return a #Eq-by-k matrix as output (tractions).
% boundaryDisplacement - a single function handle or a #BR-by-1 cell
%   array of function handles. Functions must take a k-by-#Dim matrix as
%   input (coordinates) and return a #Eq-by-k matrix as output
%   (displacements). 
% isDirichlet - a #BR-by-1 or #BR-by-#Eq matrix of ones and zeros, or
%   a function handle that takes a k-by-#Dim matrix (coordinates) and
%   returns a k-by-#Eq matrix of ones and zeros. Ones and zeros indicate
%   Dirichlet and Neumann boundary conditions respectively.
% pinnedDoFs (optional) - a list of degrees of freedom IDs, which will be
%   added to the list of prescribedDoF and removed from freeDoF. Default
%   value - [].
% pinnedVals (optional) - a list of values for the pinned degrees of
%   freedom. Default value - 0.
% 
% OUTPUT:
% u_prescribed - #Nodes-by-#Eq matrix containing Dirichlet data
%   evaluated at boundary nodes and zeros at other nodes.
% Fs - #GDoF-by-1 vector containing equivalent nodal forces
%   corresponding to Neumann boundary condition.
% prescribedDoF - list of IDs of prescribed (fixed) degrees of freedom.
% freeDoF - list of IDs of free degrees of freedom.

%% Parsing input
numNodes=size(nodeCoords,1);
numDim=size(nodeCoords,2);

if numel(varargin)>=2
    pinnedDoFs=varargin{1};
    pinnedVals=varargin{2};
elseif numel(varargin)==1
    pinnedDoFs=varargin{1};
    pinnedVals=zeros(size(pinnedDoFs));
else
    pinnedDoFs=[];
    pinnedVals=0;
end

if iscell(BIEN)
    numBR=numel(BIEN);
else
    BIEN={BIEN};
    numBR=1;
end
if ~iscell(elementType)
   elementType={elementType};
end
if numel(elementType)==1
    elementType=repmat(elementType,numBR,1);
elseif numel(elementType)~=numBR
    error('Inconsistent size of elementType');
end
if numel(numGP)==1
    numGP=repmat(numGP,numBR,1);
elseif numel(numGP)~=numBR
    error('Inconsistent size of numGP');
end
if ~iscell(numGP)
   numGP=num2cell(numGP);
end
if iscell(bndTractions)
    if numel(bndTractions)~=numBR
        error('Inconsistent size of bndTractions');
    end
    numEq=size(bndTractions{1}(nodeCoords(1,:)),1);
else
    bndTractions=repmat({bndTractions},[numBR 1]);
    numEq=size(bndTractions{1}(nodeCoords(1,:)),1);
end
if ~iscell(bndDisplacements) 
    bndDisplacements={bndDisplacements};
    bndDisplacements=repmat(bndDisplacements,numBR,1);
elseif numel(bndTractions)~=numBR
    error('Inconsistent size of bndTractions');
end

numGDoF=numNodes*numEq;
Fs=zeros(numGDoF,1);

%% ===== Neumann part =====
% Integrate traction over whole boundary, no matter if Dirichlet or Neumann

for i=1:numBR
    Fs_ = formBodyForceVector(nodeCoords, BIEN{i}, elementType{i}, numGP{i},...
        bndTractions{i});
    Fs = Fs + Fs_;
end

%% ===== Dirichlet part =====
u_prescribed=zeros(numEq,numNodes);

% Find fixed DoFs
isDirichlet_=zeros(numNodes,numEq);
for i=1:numBR
    curNodeIDs=BIEN{i}(:);
    if isa(isDirichlets,'function_handle') %function_handle
        isDirichlet_(curNodeIDs,:)=isDirichlet_(curNodeIDs,:)+...
            isDirichlets{i}(nodeCoords(curNodeIDs,:));
    elseif size(isDirichlets,2)==numEq %numBR-by-numEq
        isDirichlet_(curNodeIDs,:)=isDirichlet_(curNodeIDs,:)+...
            repmat(isDirichlets(i,:),numel(curNodeIDs),1);
    else %numBR-by-1
        isDirichlet_(curNodeIDs,:)=isDirichlet_(curNodeIDs,:)+...
            isDirichlets(i);
    end
    
    % boundaryDisplacement is evaluated disregarding the type of BC
    % prescribed at nodes
    u_prescribed(:,curNodeIDs)=bndDisplacements{i}(nodeCoords(curNodeIDs,:));
end

%% get global IDs of prescribed DoFs 
prescribedDoF=[find(isDirichlet_'); pinnedDoFs(:)];
prescribedDoF=unique(prescribedDoF);
freeDoF=(1:numGDoF)';
freeDoF(prescribedDoF)=[];
u_prescribed(pinnedDoFs(:))=pinnedVals(:);

u_prescribed=u_prescribed(:);
end %function