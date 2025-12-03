function [shapeFunctions, naturalDerivatives]=shapeFunctions2dP1(r,s)
% Shape functions for 3-node triangular element 0<=r1,r2<=1
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
shapeFunctions=[r, s, 1 - r - s];
if nargout>1
    naturalDerivatives=[1+0*r, 0+0*r, -1+0*r; ...
        0+0*r, 1+0*r, -1+0*r];
end
end %function