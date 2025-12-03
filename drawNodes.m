function drawNodes(nodeCoords,varargin)
% Draw mesh nodes using MATLAB scatter3() function.
%
% INPUT:
% nodeCoords - #Nodes-by-#Dim array containing coordinates of nodes;
% IEN or {IEN,boundaryElementsIDs} or BIEN, where  
%   IEN - #El-by-#NpE array, element-node incidence matrix,
%   boundaryElementIDs - #BR-by-1 cell array, p-th cell contains l_p-by-1 
%    matrix listing boundary element IDs.
%    Providing both IEN and boundaryElementsIDs is equivalent to providing
%    IEN(unique([boundaryElementIDs{:}]),:), but element labels will differ
%    if drawn.
%   BIEN - #BR-by-1 cell array containing boundary element-node incidence
%    matrices.
% specs (optional) - string or cell array of parameters, first of which is
%   a string. The parameters are passed to scatter3() directly. Default
%   value - {'k','filled'}
%
% COMMENTS:
% Uses MATLAB's scatter3() object.
% Only nodes of elements listed in IEN (or in boundaryElementsIDs, when
% applicable) are shown. 

%
% Examples:
%   drawNodes(nodeCoords)
%   drawNodes(nodeCoords,specs)
%   drawNodes(nodeCoords,'o')
%   drawNodes(nodeCoords,IEN)
%   drawNodes(nodeCoords,IEN,specs)

%% parsing input
if isempty(nodeCoords)
    return
end;
if nargin==1 %drawNodes(nodeCoords)
    specs={'k','filled'};
    IENprovided=0;
elseif (nargin==2) && ...  % drawNodes(nodeCoords,'o') or drawNodes(nodeCoords,specs)
        ( ischar(varargin{1}) || (iscell(varargin{1})&&ischar(varargin{1}{1})) )
    specs=varargin(1);
    IENprovided=0;
elseif nargin==2
    specs={'k','filled'};
    IENprovided=1;
    IEN=varargin{1};
elseif nargin==3
    IENprovided=1;
    IEN=varargin{1};
    specs=varargin(2);
end
if iscell(specs{1})
    specs=specs{1};
end

if IENprovided==1
    if iscell(IEN)&&numel(IEN)==2&&iscell(IEN{2}) % {IEN,boundaryElementsIDs} provided in place of IEN
        [elementNumbers,~,~]=unique(vertcat(IEN{2}{:}));
        IEN=IEN{1}(elementNumbers,:);
    elseif iscell(IEN)&&size(IEN,2)==1 % BIEN
        IEN=cell2mat(IEN);
    elseif ~iscell(IEN)
        elementNumbers=1:size(IEN,1);
    else
        error('Wrong format of argument 2 (IEN)');
    end
end
numNodes=size(nodeCoords,1);

if size(nodeCoords,2)==2
    % making it 3D from 2D if necessary
    nodeCoords=[nodeCoords zeros(numNodes,1)];
elseif size(nodeCoords,2)~=3
    error('Only 2D and 3D supported');
end

%% 
% plotting
old_hold=ishold; % remember figure hold state to revert it later

if IENprovided % IEN provided
    xyz=nodeCoords(unique(nonzeros(IEN)),:);
    scatter3(xyz(:,1),xyz(:,2),xyz(:,3),[],specs{:});
    hold on;
else % IEN not provided
    scatter3(nodeCoords(:,1),nodeCoords(:,2),nodeCoords(:,3),[],specs{:});
end

if ~old_hold
    hold off
end
axis equal
end