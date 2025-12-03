function [shapeFunctions, naturalDerivatives]=shapeFunctions2dQ1(r,s)
% Shape functions for 4-node linear quadrilateral element -1<=r1,r2<=1
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

shapeFunctions=0.25*[(1-r)*(1-s), (1+r)*(1-s), (1+r)*(1+s), (1-r)*(1+s)];
if nargout>1
    naturalDerivatives=0.25*[-(1-s), (1-s), (1+s), -(1+s);...
        -(1-r), -(1+r), (1+r), (1-r)];
end
end %function