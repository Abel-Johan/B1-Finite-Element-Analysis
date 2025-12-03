function DoF_IDs=node2DoFs(nodeIDs,numEq,eqIndex)
% Convert IDs of nodes to the corresponding IDs of Degrees of Freedom.
% (the two are identical in case of scalar problem, when numEq=1)
%
% INPUT:
% nodeIDs - k-by-1 vector or k-by-m array
% numEq - a single integer
% eqIndex - (optional) 1-by-m vector or k-by-m array containing non-negative
%   integers - equation numbers, i.e. IDs of local DoFs. Zeros indicate the
%   corresponding entry must not be included in the output. If k-by-m array
%   then m may equal 1 as long as condition size(eqIndex)==size(nodeIDs)
%   holds.  
% 
% OUTPUT:
% DoF_IDs - k*m-by-1 array DoF_IDs=(nodeIDs-1)*numEq + eqIndex;
%   or n-by-1 array (n<k*m) if some elements of eqIndex are zero.
%
% EXAMPLE 1: k=3, m=1;
%  Input: nodeIDs=[1 2 3]'; numEq=3; eqIndex=1;
%  Output: DoF_IDs=[1 4 7]';
%
% EXAMPLE 2: k=3, m=1;
%  Input: nodeIDs=[1 2 3]'; numEq=3; eqIndex=[1 1 2]';
%  Output: DoF_IDs=[1 4 8]';
%
% EXAMPLE 3: k=3, m=2;
%  Input: nodeIDs=[1 2 3]'; numEq=3; eqIndex=[1 1 2];
%  Output: DoF_IDs=[1 1 2 4 4 5 7 7 8]';
%
% EXAMPLE 4: k=3, m=2;
%  Input: nodeIDs=[1 1; 2 2; 1 2]; numEq=3; eqIndex=[1 2; 1 2; 2 1];
%  Output: DoF_IDs=[1 2 4 5 2 4]';
%
% EXAMPLE 5: k=3, m=2;
%  Input: nodeIDs=[1 1; 2 2; 1 2]; numEq=3; eqIndex=[1 2; 1 2; 0 1];
%  Output: DoF_IDs=[1 2 4 5 4]';
%
if nargin<=2
    eqIndex=1:numEq;
end
if size(eqIndex,1)==1
    eqIndex=repmat(eqIndex,size(nodeIDs,1),1);
end
if size(nodeIDs,2)==1
    nodeIDs=repmat(nodeIDs,1,size(eqIndex,2));
end
eqIndex=eqIndex';
nodeIDs=nodeIDs';
eqIndex=eqIndex(:);
nodeIDs=nodeIDs(:);
entries2delete=find(eqIndex==0);
eqIndex(entries2delete)=[];
nodeIDs(entries2delete)=[];

DoF_IDs=(nodeIDs-1)*numEq + eqIndex;
end