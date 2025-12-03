function [nodeCoords2, IEN2, varargout] = meshP1toP2(nodeCoords,IEN,varargin)
% Convert '2dP1' mesh (linear, 3-node triangles) to '2dP2' mesh (quadratic,
% 6-nodes triangles) by adding a midpoint to every edge.
%
% INPUT:
% nodeCoords - #Nodes-by-#Dim matrix, coordinates of nodes.
% IEN - #El-by-#NpE array, element-node incidence matrix. #NpE=3. 
% boundaryElementIDs (optional) - #BR-by-1 cell array: p-th cell contains
%   l_p-by-1 matrix listing boundary element IDs.
% boundaryNodeLocalIDs (optional) - #BR-by-1 cell array: p-th cell contains
%   l_p-by-q1 matrix listing local IDs (within an element) of corresponding
%   boundary nodes. boundaryNodeLocalIDs(m,1) and
%   boundaryNodeLocalIDs(m,q1) define the endpoints of the boundary edge
%   of the element number boundaryElementIDs(m).
%
% OUTPUT:
% nodeCoords2 - #Nodes-by-#Dim matrix, coordinates of nodes.
% IEN2 - #El-by-#NpE array, element-node incidence matrix. 
% boundaryElementIDs2 (optional) - #BR-by-1 cell array: p-th cell contains
%   l_p-by-1 matrix listing boundary element IDs.
% boundaryNodeLocalIDs (optional) - #BR-by-1 cell array: p-th cell contains
%   l_p-by-q2 matrix listing local IDs (within an element) of corresponding
%   boundary nodes. boundaryNodeLocalIDs(m,1) and
%   boundaryNodeLocalIDs(m,q2) define the endpoints of the boundary edge
%   of the element number boundaryElementIDs(m).
%
% COMMENTS:
% boundaryElementIDs2 and boundaryNodeLocalIDs2 are computed only if 
% boundaryElementIDs and boundaryNodeLocalIDs were provided as input. In
% this implementation boundaryElementIDs2 equals boundaryElementIDs.

if size(IEN,2)~=3
    error('This only works for P1 2D mesh');
end

% we want to make a list of all edges
% Therefore we combine all edges of all elements. Most edges will be
% duplicated, as the belong to two elements. Therefore we will remove
% duplicates. In order to ease the comparison of edges (whether they are
% identical), we make sure nodes listed in the ascending order within each
% edge.

k1=1;k2=2;
pairs1=[(IEN(:,k1)<IEN(:,k2)).*IEN(:,k1) + (IEN(:,k1)>IEN(:,k2)).*IEN(:,k2), ...
        (IEN(:,k1)<IEN(:,k2)).*IEN(:,k2) + (IEN(:,k1)>IEN(:,k2)).*IEN(:,k1)];
k1=2;k2=3;
pairs2=[(IEN(:,k1)<IEN(:,k2)).*IEN(:,k1) + (IEN(:,k1)>IEN(:,k2)).*IEN(:,k2), ...
        (IEN(:,k1)<IEN(:,k2)).*IEN(:,k2) + (IEN(:,k1)>IEN(:,k2)).*IEN(:,k1)];
k1=3;k2=1;
pairs3=[(IEN(:,k1)<IEN(:,k2)).*IEN(:,k1) + (IEN(:,k1)>IEN(:,k2)).*IEN(:,k2), ...
        (IEN(:,k1)<IEN(:,k2)).*IEN(:,k2) + (IEN(:,k1)>IEN(:,k2)).*IEN(:,k1)];

pairs=[pairs1;pairs2;pairs3]; 
%now we have the list of all edges, including duplicates
%it is important that element i has edges in rows i, i+#El, i+2*#El

[unique_pairs, ~, ic]=unique(pairs,'rows'); %remove duplicates
% ic allows us to map a (unique) edge to the elements it belongs to

midpoints=.5*( nodeCoords(unique_pairs(:,1),:) + nodeCoords(unique_pairs(:,2),:) );
% list of midpoints;

nodeCoords2=[nodeCoords; midpoints]; %midpoints are added after vertices

IEN2=zeros(size(IEN,1),6);
IEN2(:,1:3)=IEN;

% now we assign a nodal ID to each midpoint. Since all existing nodes are
% vertices, midpoint IDs will be from the range
% [ #Vertices+1, ..., #Vertices + #Midpoints ] 
pairs_IDs=reshape(1:size(pairs,1),[size(pairs,1)/3 3]);
midpointNodeIDs=ic(pairs_IDs)+size(nodeCoords,1);

IEN2(:,4:6) = midpointNodeIDs;

% mesh boundary
if (nargin==4)&&(nargout==4)
    bElIDs=varargin{1};
    bNLocIDs=varargin{2};
    bNLocIDs2=cell(size(bElIDs));
    for i=1:numel(bNLocIDs)
        sorted_pairs=[(bNLocIDs{i}(:,1)>bNLocIDs{i}(:,2)).*bNLocIDs{i}(:,2) + ...
                      (bNLocIDs{i}(:,1)<bNLocIDs{i}(:,2)).*bNLocIDs{i}(:,1),...
                      (bNLocIDs{i}(:,1)>bNLocIDs{i}(:,2)).*bNLocIDs{i}(:,1) + ...
                      (bNLocIDs{i}(:,1)<bNLocIDs{i}(:,2)).*bNLocIDs{i}(:,2)];
        % Assign midpoint local IDs so that the edges are
        % 1 4 2, 2 5 3, or 3 6 1.
        midNodeLocalIDs= (sorted_pairs(:,1)==2).*5 + ...
                         (sorted_pairs(:,1)==1).*(sorted_pairs(:,2)==2).*4 + ... 
                         (sorted_pairs(:,1)==1).*(sorted_pairs(:,2)==3).*6;
        bNLocIDs2{i}=[bNLocIDs{i}(:,1) midNodeLocalIDs bNLocIDs{i}(:,2)];
    end
    varargout={bElIDs,bNLocIDs2};
end

end %function