function elData=elementData(elementType,varargin)
% Specific data for a requested element type.
% 
% INPUT:
% elementType - a single string specifying the type of elements used
% fieldName (optional) - a single string specifying a requested property of the
%   structure
% argument (optional) - if the requested field is a function handle, the
%   function is applied to the argument and the result is returned.
%
% OUTPUT:
% elData - elementData structure, if no fieldName specified. Otherwise,
%   the value of the requested field.
%
% Examples:
%   elData=elementData(elementType)
%   numGP=elementData(elementType,'numGPFull')
%   inElement=elementData('2dQ1','inElement') 
%   inElement([1 1]) % ans=1
%   elementData('2dQ1','inElement',[2 0]) % ans=0
%
% elementData structure FIELDS:
% 'nodeNatCoords' - #Nodes-by-#ElDim matrix, natural nodal coordinates.
% 'numGPFull' - number of Gauss points for exact integration in FE space.
% 'numGPSPR' - number of Gauss points for SPR (superconvergent patch
%   recovery).
% 'bndElementType' - string, type of boundary element.
% 'elLocBIEN' - #bndEl-by-#bndNpE matrix, elLocBIEN(k,:) lists IDs of nodes
%   incident to k-th boundary element.
% 'polTerms' - #p-to-#ElDim (#ElDim-dimensional square array with #p
%   elements in each dimension) - array of 1s and 0s representing polynomial
%   space spanned by shape functions 
% 'p0' - function handle, polynomial basis for FE space, takes natural
%   coordinates x_ (k-by-#ElDim) as input and returns functions' values p0_
%   (k-by-#Nodes) as output.
% 'vertexLocalIDs' - 1-by-l1 matrix, lists IDs of vertex nodes (as opposed
%   to midpoints)
% 'facetsLocBIEN' - #Faces-by-l2 matrix, facetsLocBIEN(k,:) lists IDs of nodes
%   incident to k-th facet. Used for element visualization.
% 'inElement' - function handle, shows if a point belongs to the element
%   natural domain. Takes natural coordinates x_ (k-by-#ElDim) as input and
%   returns k-by-1 matrix of 1s (true) and 0s (false). 
% 'inElementTol' - function handle, shows if a point belongs to the element
%   natural domain within certain tolerance. Takes natural coordinates x_
%   (k-by-#ElDim) as input and returns k-by-1 matrix of 1s (true) and 0s
%   (false).

switch elementType
    case '1dQ1'
        nodeNatCoords=[-1;1];
        polTerms=[1 1]';
        bndElementType='0d';
        elLocBIEN=[1; 2];
        numGPFull=1;
        numGPSPR=1;
        p0=@(x)[1+0*x, x];
        p0=@(x)p0(x(:,1));
        vertexLocalIDs=[1 2];
        facetsLocBIEN=[];
        inElement=@(x)...
            (-1<=x(:,1))&(x(:,1)<=1);
        inElement2=@(x,tol)...
            (-1-tol<=x(:,1))&(x(:,1)<=1+tol);
    case '1dQ2'
        nodeNatCoords=[-1; 0; 1];
        polTerms=[1 1 1]';
        bndElementType='0d';
        elLocBIEN=[1; 3];
        numGPFull=2;
        numGPSPR=2;
        p0=@(x)[1+0*x, x.^2, x];
        p0=@(x)p0(x(:,1));
        vertexLocalIDs=[1 3];        
        facetsLocBIEN=[];
        inElement=@(x)...
            (-1<=x(:,1))&(x(:,1)<=1);
        inElement2=@(x,tol)...
            (-1-tol<=x(:,1))&(x(:,1)<=1+tol);
    case '1dP1'
        nodeNatCoords=[0;1];
        polTerms=[1 1];
        bndElementType='0d';
        elLocBIEN=[1; 2];
        numGPFull=1;
        p0=@(x)[1+0*x, x];
        p0=@(x)p0(x(:,1));
        vertexLocalIDs=[1 2];     
        facetsLocBIEN=[];
        inElement=@(x)...
            (0<=x(:,1))&(x(:,1)<=1);
        inElement2=@(x,tol)...
            (-tol<=x(:,1))&(x(:,1)<=1+tol);
    case '2dQ1'
        nodeNatCoords=[-1 -1; 1 -1; 1 1; -1 1];
        polTerms=[1 1; 1 1];
        bndElementType='1dQ1';
        elLocBIEN=[1 2; 2 3; 3 4; 4 1];
        numGPFull=4;
        numGPSPR=1;
        p0=@(x,y)[1+0*x, x, y, x.*y];
        p0=@(x)p0(x(:,1),x(:,2));
        vertexLocalIDs=[1 2 3 4];
        facetsLocBIEN=[1 2 3 4];
        inElement=@(x)...
            (-1<=x(:,1))&(x(:,1)<=1)&...
            (-1<=x(:,2))&(x(:,2)<=1);
        inElement2=@(x,tol)...
            (-1-tol<=x(:,1))&(x(:,1)<=1+tol)&...
            (-1-tol<=x(:,2))&(x(:,2)<=1+tol);
    case '2dQ2r'
        nodeNatCoords=[-1 -1; 1 -1; 1 1; -1 1; 0 -1; 1 0; 0 1; -1 0];
        polTerms=[1 1 1; 1 1 1; 1 1 0];
        bndElementType='1dQ2';
        elLocBIEN=[1 5 2; 2 6 3; 3 7 4; 4 8 1];
        numGPFull=4;
        numGPSPR=4;
        p0=@(x,y)[1+0*x, x, y, x.^2, x.*y, y.^2, x.^2.*y, x.*y.^2, x.^2.*y.^2];
        p0=@(x)p0(x(:,1),x(:,2));
        vertexLocalIDs=[1 2 3 4];
        facetsLocBIEN=[1 5 2 6 3 7 4 8];
        inElement=@(x)...
            (-1<=x(:,1))&(x(:,1)<=1)&...
            (-1<=x(:,2))&(x(:,2)<=1);
        inElement2=@(x,tol)...
            (-1-tol<=x(:,1))&(x(:,1)<=1+tol)&...
            (-1-tol<=x(:,2))&(x(:,2)<=1+tol);
    case '2dQ2'
        nodeNatCoords=[-1 -1; 1 -1; 1 1; -1 1; 0 -1; 1 0; 0 1; -1 0; 0 0];
        polTerms=[1 1 1; 1 1 1; 1 1 1];
        bndElementType='1dQ2';
        elLocBIEN=[1 5 2; 2 6 3; 3 7 4; 4 8 1];
        numGPFull=4;
        numGPSPR=4;
        p0=@(x,y)[1+0*x, x, y, x.^2, x.*y, y.^2, x.^2.*y, x.*y.^2, x.^2.*y.^2];
        p0=@(x)p0(x(:,1),x(:,2));
        vertexLocalIDs=[1 2 3 4];
        facetsLocBIEN=[1 5 2 6 3 7 4 8];
        inElement=@(x)...
            (-1<=x(:,1))&(x(:,1)<=1)&...
            (-1<=x(:,2))&(x(:,2)<=1);
        inElement2=@(x,tol)...
            (-1-tol<=x(:,1))&(x(:,1)<=1+tol)&...
            (-1-tol<=x(:,2))&(x(:,2)<=1+tol);
    case '2dP1'
        nodeNatCoords=[1 0; 0 1; 0 0];
        polTerms=[1 1; 1 0];
        bndElementType='1dQ1';
        elLocBIEN=[1 2; 2 3; 3 1];
        numGPFull=1;
        numGPSPR=1;
        p0=@(x,y)[1+0*x, x, y];
        p0=@(x)p0(x(:,1),x(:,2));
        vertexLocalIDs=[1 2 3];
        facetsLocBIEN=[1 2 3];
        inElement=@(x)...
            (0<=x(:,1))&(x(:,1)<=1)&...
            (0<=x(:,2))&(x(:,2)<=1)&...
            (0<=x(:,2)+x(:,1))&(x(:,1)+x(:,2)<=1);
        inElement2=@(x,tol)...
            (-tol<=x(:,1))&(x(:,1)<=1+tol)&...
            (-tol<=x(:,2))&(x(:,2)<=1+tol)&...
            (-tol<=x(:,2)+x(:,1))&(x(:,1)+x(:,2)<=1+tol);
    case '2dP2'
        nodeNatCoords=[1 0; 0 1; 0 0; .5 .5; 0 .5; .5 0];
        polTerms=[1 1 1; 1 1 0; 1 0 0];
        elLocBIEN=[1 4 2; 2 5 3; 3 6 1];
        bndElementType='1dQ2';
        numGPFull=3;
        numGPSPR=3;
        p0=@(x,y)[1+0*x, x, y, x.^2, x.*y, y.^2];
        p0=@(x)p0(x(:,1),x(:,2));        
        vertexLocalIDs=[1 2 3];
        facetsLocBIEN=[1 4 2 5 3 6];
        inElement=@(x)...
            (0<=x(:,1))&(x(:,1)<=1)&...
            (0<=x(:,2))&(x(:,2)<=1)&...
            (0<=x(:,2)+x(:,1))&(x(:,1)+x(:,2)<=1);
        inElement2=@(x,tol)...
            (-tol<=x(:,1))&(x(:,1)<=1+tol)&...
            (-tol<=x(:,2))&(x(:,2)<=1+tol)&...
            (-tol<=x(:,2)+x(:,1))&(x(:,1)+x(:,2)<=1+tol);
    case '3dQ1' % 8-node brick
        nodeNatCoords=[-1 -1 -1; 1 -1 -1; 1 1 -1; -1 1 -1;...
                       -1 -1 +1; 1 -1 +1; 1 1 +1; -1 1 +1];
        polTerms=[1 1; 1 1];
        polTerms(:,:,2)=[1 1; 1 1];
        bndElementType='2dQ1';
        elLocBIEN=[1 2 3 4; 8 7 6 5; ...
            1 4 8 5; 3 2 6 7; ...
            2 1 5 6; 4 3 7 8];
        numGPFull=8;
        numGPSPR=1;
        p00=@(x,y,z)[1+0*x, x, y, z, x.*y, y.*z, z.*x, x.*y.*z];
        p0=@(x)p00(x(:,1),x(:,2),x(:,3));        
        vertexLocalIDs=[1 2 3 4 5 6 7 8];
        facetsLocBIEN=[1 2 3 4; 8 7 6 5; ...
            1 4 8 5; 3 2 6 7; ...
            2 1 5 6; 4 3 7 8];
        % used for drawing
        inElement=@(x)...
            (-1<=x(:,1))&(x(:,1)<=1)&...
            (-1<=x(:,2))&(x(:,2)<=1)&...
            (-1<=x(:,3))&(x(:,3)<=1);
        inElement2=@(x,tol)...
            (-1-tol<=x(:,1))&(x(:,1)<=1+tol)&...
            (-1-tol<=x(:,2))&(x(:,2)<=1+tol)&...
            (-1-tol<=x(:,3))&(x(:,3)<=1+tol);
    case '3dP1' % 4-node tetrahedral
        nodeNatCoords=[1 0 0; 0 1 0; 0 0 1; 0 0 0];
        polTerms=[1 1; 1 0];
        polTerms(:,:,2)=[1 0; 0 0];
        bndElementType='2dP1';
        elLocBIEN=[2 3 4; 3 1 4; 1 2 4; 1 3 2];
        numGPFull=1;
        numGPSPR=1;
        p00=@(x,y,z)[1+0*x, x, y, z];
        p0=@(x)p00(x(:,1),x(:,2),x(:,3));        
        vertexLocalIDs=[1 2 3 4];
        facetsLocBIEN=[2 3 4; 3 1 4; 1 2 4; 1 3 2;];
        % used for drawing
        inElement=@(x)...
            (0<=x(:,1))&(x(:,1)<=1)&...
            (0<=x(:,2))&(x(:,2)<=1)&...
            (0<=x(:,3))&(x(:,3)<=1)&...
            (0<=x(:,3)+x(:,2)+x(:,1))&(x(:,3)+x(:,2)+x(:,1)<=1);
        inElement2=@(x,tol)...
            (0-tol<=x(:,1))&(x(:,1)<=1+tol)&...
            (0-tol<=x(:,2))&(x(:,2)<=1+tol)&...
            (0-tol<=x(:,3))&(x(:,3)<=1+tol)&...
            (0-tol<=x(:,3)+x(:,2)+x(:,1))&(x(:,3)+x(:,2)+x(:,1)<=1+tol);
    otherwise
        elData=[];
        warning(['Element type ' elementType ' not supported']);
        return;
end

elData=struct(...
    'nodeNatCoords',nodeNatCoords,...
    'bndElementType',bndElementType,...
    'elLocBIEN',elLocBIEN,...
    'numGPFull',numGPFull,...
    'numGPSPR',numGPSPR,...
    'polTerms',polTerms,...
    'p0',p0,...
    'vertexLocalIDs',vertexLocalIDs,...
    'facetsLocBIEN',facetsLocBIEN,...
    'inElement',inElement,...
    'inElementTol',inElement2...
    );
if numel(varargin)>=1
    elData=elData.(varargin{1});
end
if numel(varargin)==2
    elData=elData(varargin{2});
end

    
end