function [shapeFunctions, naturalDerivatives]=shapeFunctions2dP2(r,s)
% Shape functions for 6-node triangular element 0<=r1,r2<=1
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
    r=r(:,1);
end
shapeFunctions=[r*(2*r-1), s*(2*s-1), (1 - r - s)*(1 - 2*r - 2*s),...
    4*r*s, 4*s*(1 - r - s), 4*r*(1 - r - s)];
if nargout>1
    naturalDerivatives=[4*r-1, 0, -3 + 4*s + 4*r,...
    4*s, -4*s, 4-8*r-4*s;
    0, 4*s-1, -3 + 4*r + 4*s,...
    4*r, -8*s+4-4*r, -4*r];
end
end %function