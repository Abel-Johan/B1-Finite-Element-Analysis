function [boundaryElementIDs, boundaryNodeLocalIDs]=edges2sublists(e,IEN,INE)
% Convert Matlab PDEToolbox (2D) edge representation to sublists of boundary
% elements and corresponding local node IDs within each boundary element. 
% 
% INPUT:
% e - 7-by-k matrix, edge representation produced by PDEToolbox initmesh()
% IEN - #El-by-#NpE array, element-node incidence matrix. 
% INE - #Nodes-by-#EpN array, node-element incidence matrix, row INE(i,:)
%   corresponds to ith node and contains IDs of elements that share it.
%
% OUTPUT:
% boundaryElementIDs - #BR-by-1 cell array: p-th cell contains l_p-by-1 
%   matrix listing boundary element IDs.
% boundaryNodeLocalIDs - #BR-by-1 cell array: p-th cell contains
%   l_p-by-#BNpE matrix  listing local IDs (within an element) of
%   corresponding boundary nodes. boundaryNodeLocalIDs(m,1) and
%   boundaryNodeLocalIDs(m,end) define the boundary edge of the element
%   number boundaryElementIDs(m). Note that l_1 + ... + l_#BR <= k.
%
% Examples:
%   [p,e,t] = initmesh(dl,'Hmax',hmax);
%   nodeCoords=p';
%   IEN=IEN(1:3,:)';
%   INE=IENtoINE(IEN);
%   [boundaryElementIDs, boundaryNodeLocalIDs]=edges2sublists(e,IEN,INE);
%
% COMMENTS:
% Only works for meshes consisting of 3-node triangular elements.
% Only external boundaries are taken into account. Interfaces between
% subdomains are disregarded.

%
b_e=e(:,((e(6,:)==0)|(e(7,:)==0))); %select boundary edges only
% thus edges at the interface in the interior are ignored
numSeg=unique(b_e(5,:)); %number of segments 

boundaryElementIDs=cell(numel(numSeg),1); 
boundaryNodeLocalIDs=cell(numel(numSeg),1);
for i=numSeg % loop over boundary segments (regions)
    curEdgeIDs=find(b_e(5,:)==i); % get the list of all edges in current segment
    curBoundElIDs=zeros(numel(curEdgeIDs),1);
    curNodeLocalIDs=zeros(numel(curEdgeIDs),2);
    for j=1:numel(curEdgeIDs) % loop through each edge in current segment
        % get the lists of elements the two ends of the edge belong to
        elementList1=nonzeros(INE(b_e(1,curEdgeIDs(j)),:));
        elementList2=nonzeros(INE(b_e(2,curEdgeIDs(j)),:));
        elID=intersect(elementList1,elementList2);
        if numel(elID)~=1
            error('Could not determine the element an edge belongs to');
        end
        curBoundElIDs(j)=elID;
        try
            curNodeLocalIDs(j,1)=find(IEN(elID,:)==b_e(1,curEdgeIDs(j)));
            curNodeLocalIDs(j,2)=find(IEN(elID,:)==b_e(2,curEdgeIDs(j)));
        catch ME
            warning('IEN and INE provided seem to be inconsistent');
            rethrow(ME);
        end %try
    end %for j=1:numel(curEdgeIDs)
    boundaryElementIDs{i}=curBoundElIDs;
    boundaryNodeLocalIDs{i}=curNodeLocalIDs;
end %for i=1:numel(segments)

end %function