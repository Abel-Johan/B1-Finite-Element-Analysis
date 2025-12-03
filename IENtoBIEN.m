function BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs)
% Restrict a given element-node incidence array (IEN) to given lists of
%   elements and nodes within each element. 
%
% IEN - #El-by-#NpE array, element-node incidence matrix. 
% boundaryElementIDs - #BR-by-1 cell array: p-th cell contains l_p-by-1 
%   matrix listing boundary element IDs.
% boundaryNodeLocalIDs - #BR-by-1 cell array: p-th cell contains l_p-by-q
%   matrix listing local IDs (within an element) of corresponding boundary
%   nodes.
%
% OUTPUT: 
% BIEN - a #BR-by-1 cell array. Cell BIEN{k} contains a #BEl(k)-by-#BNpE element-node
%   incidence array corresponding to boundary region k. The order of IDs
%   within row BIEN{k}(i,:) follows the order in boundaryNodeLocalIDs. 
%
% COMMENTS:
% This routine is useful for representing the boundary of a mesh as a
% sub-mesh.
%
% Examples:
% BIEN = IENtoBIEN(IEN,boundaryElementIDs,boundaryNodeLocalIDs)

numBR=numel(boundaryElementIDs);
BIEN=cell(size(boundaryElementIDs));
for i=1:numBR
BIEN{i}=IEN(...
    sub2ind(size(IEN), ...
    repmat(boundaryElementIDs{i},1,size(boundaryNodeLocalIDs{i},2)), ...
    boundaryNodeLocalIDs{i}));
end
end %function

