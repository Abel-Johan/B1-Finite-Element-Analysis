function [shapeFunctions, naturalDerivatives]=shapeFunctions3dP1(r,s,t)
% Shape functions for 4-node linear tetrahedral element 0<=r1,r2,r3<=1
%
% INPUT:
% r - k-by-#ElDim matrix, row r(k,:) contains natural coordinates of k-th
%   point. 
%
% OUTPUT:
% shapeFunctions - k-by-#Nodes matrix,  shapeFunctions(l,i) is the value of
%   i-th shape function at point l. 
% natDerivatives - #ElDim*k-by-#Nodes matrix, natDerivatives(k*(j-1)+l,i)
%   is the derivative of i-th shape function with respect to j-th
%   coordinate evaluated at point l.

if nargin==1
    s=r(:,2);
    t=r(:,3);
    r=r(:,1);
end

shapeFunctions=[...
    r, s, t, 1-r-s-t...
    ];
if nargout>1
    naturalDerivatives=[...
        1, 0, 0, -1;...
        0, 1, 0, -1;...
        0, 0, 1, -1;...     
       ];
end
end %function