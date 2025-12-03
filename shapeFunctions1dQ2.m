function [shapeFunctions, naturalDerivatives]=shapeFunctions1dQ2(r)
% Shape functions for 2-node quadratic 1D element -1<=r1<=1
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

shapeFunctions=[.5*(r-1)*r, 1-r^2, .5*(r+1)*r];
if nargout>1
    naturalDerivatives=[r-.5, -2*r, r+.5];
end
end %function
