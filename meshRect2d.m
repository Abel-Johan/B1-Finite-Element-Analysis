function [nodeCoords, IEN, boundaryElementIDs, boundaryNodeLocalIDs, e] = meshRect2d(domain, elementType, elementSize)
% Create a mesh for a rectangular 2D domain.
%
% INPUT:
% domain (optional) - 1-by-4 matrix [x1 y1 x2 y2] or 2-by-2 [x1 y1; x2 y2]
%   with coordinates defining the domain. Default value - [-2 -2 2 2].
% elementType (optional) - single string, the type of elements. Default
%   value - '2dQ1'.
% elementSize (optional) - a single real or 2-by-1 matrix. A single real
%   defines the maximum diameter of elements. A 2-by-1 matrix defines the
%   number of "cells" in each row and column of the rectangular grid. Each
%   "cell" is a quadrilateral element or a pair of triangular elements.
%   Default value - sqrt(2). 
%
% OUTPUT:
% nodeCoords - #Nodes-by-#Dim matrix, coordinates of nodes. #Dim=2.
% IEN - #El-by-#NpE array, element-node incidence matrix. 
% boundaryElementIDs - #BR-by-1 cell array: p-th cell contains l_p-by-1
%   matrix listing boundary element IDs. #BR=4.
% boundaryNodeLocalIDs - #BR-by-1 cell array: p-th cell contains l_p-by-q
%   matrix listing local IDs (within an element) of corresponding boundary
%   nodes. boundaryNodeLocalIDs(m,1) and boundaryNodeLocalIDs(m,end) define
%   the endpoints of the boundary edge of the element number
%   boundaryElementIDs(m).
% e - edges, the data structure supported by PDEToolbox, defines the
%   boundaries of the domain.
%
% COMMENTS:
% The supported element types are '2dQ1', '2dQ2r', '2dQ2', '2dP1', '2dP2'.
% Within each row IEN(i,:), vertices are listed before the midpoints. All
% nodes are listed anti-clockwise.
% The four cells in boundaryElementIDs and boundaryNodeLocalIDs correspond
% respectively to LOWER, RIGHT, UPPER and LEFT boundaries of the
% rectangular domain.
%
% Examples:
% [nodeCoords, IEN, bndElIDs, bndNodeLocIDs] = meshRect2d(domain,...
%   elementType, elementSize) 
% [nodeCoords, IEN, bndElIDs, bndNodeLocIDs] = meshRect2d([0 0; 1 1],...
%   '2dQ1', [10 10]) 
% [nodeCoords, IEN, bndElIDs, bndNodeLocIDs] = meshRect2d([0 0; 1 1],...
%   '2dQ1', .1) 
%%
% Assigning default values if some input is missing
if nargin < 1
    domain = [-2 -2 2 2];
elseif size(domain,1)==2 && size(domain,2)==2
    domain=domain';
end
if nargin < 2
    elementType = '2dQ1';
end
if nargin < 3
    elementSize = sqrt(2);
end

% Checking validity of the input
if ~( (numel(elementSize)==1 && elementSize>0) || ...
        (numel(elementSize)==2 && elementSize(1) >= 1 && elementSize(2) >= 1) )
    error('Element size must be a positive real number or a pair of positive integer real numbers');
end

if (domain(3)-domain(1))*(domain(4)-domain(2))==0
    error('The specified domain is degenerate');
end

% NumX, NumY - number of cells in each direction, chosen in a way that
% guarantees a bound on the element diameter. Each cell corresponds to one
% rectangular or to two triangular elements.
x1=domain(1);y1=domain(2);
x2=domain(3);y2=domain(4);
if numel(elementSize)==1
    NumX = ceil(2^.5*abs(x2-x1)/elementSize);
    NumY = ceil(2^.5*abs(y2-y1)/elementSize);
elseif numel(elementSize)==2
    NumX=elementSize(1);
    NumY=elementSize(2);
end

%% nodeCoords and IEN
switch elementType
    case '2dQ1' %4-node quadrilateral element
        k=0;
        nodeCoords=zeros((NumX+1)*(NumY+1),2);
        nodeCoordsX=linspace(x1,x2,NumX+1);
        nodeCoordsY=linspace(y1,y2,NumY+1);
        for j=1:NumY+1
            for i=1:NumX+1
                k=k+1;
                nodeCoords(k,:)=[nodeCoordsX(i),nodeCoordsY(j)];
            end
        end
        
        IEN=zeros(NumX*NumY,4);
        k=0;
        for j=1:NumY
            for i=1:NumX
                k=k+1;
                IEN(k,:)=[(j-1)*(NumX+1)+i, (j-1)*(NumX+1)+i+1,...
                    j*(NumX+1)+i+1, j*(NumX+1)+i];
            end
        end
    case '2dQ2r' %8-node quadrilateral element
        k=0;
        nodeCoords=zeros((2*NumX+1)*(2*NumY+1)-NumX*NumY,2);
        nodeCoordsX=linspace(x1,x2,2*NumX+1);
        nodeCoordsY=linspace(y1,y2,2*NumY+1);
        for j=1:2:2*NumY-1 %iterate over all but last rows grouping by two
            for i=1:2*NumX+1 %vertex and midpoints
                k=k+1;
                nodeCoords(k,:)=[nodeCoordsX(i),nodeCoordsY(j)];
            end
            for i=1:2:2*NumX+1 %mipoints only
                k=k+1;
                nodeCoords(k,:)=[nodeCoordsX(i),nodeCoordsY(j+1)];
            end
        end
        for i=1:2*NumX+1 %vertex and midpoints, last row
            k=k+1;
            nodeCoords(k,:)=[nodeCoordsX(i),nodeCoordsY(end)];
        end
        
        IEN=zeros(NumX*NumY,8);
        k=0;
        for j=1:NumY
            for i=1:NumX
                k=k+1;
                IEN(k,:)=[(j-1)*(3*NumX+2)+2*i-1,   (j-1)*(3*NumX+2)+2*i+1,...
                    j*(3*NumX+2)+2*i+1,             j*(3*NumX+2)+2*i-1,...
                    (j-1)*(3*NumX+2)+2*i,       (j-1)*(3*NumX+2)+2*NumX+2+i,...
                    j*(3*NumX+2)+2*i,           (j-1)*(3*NumX+2)+2*NumX+1+i];
            end
        end
    case '2dQ2' %9-node quadrilateral element
        k=0;
        nodeCoords=zeros((2*NumX+1)*(2*NumY+1),2);
        nodeCoordsX=linspace(x1,x2,2*NumX+1);
        nodeCoordsY=linspace(y1,y2,2*NumY+1);
        for j=1:2*NumY+1 %iterate over all rows
            for i=1:2*NumX+1 %iterate over all columns
                k=k+1;
                nodeCoords(k,:)=[nodeCoordsX(i),nodeCoordsY(j)];
            end
        end
        
        IEN=zeros(NumX*NumY,9);
        k=0;
        for j=1:NumY
            for i=1:NumX
                k=k+1;
                IEN(k,:)=[(2*j-2)*(2*NumX+1)+2*i-1,   (2*j-2)*(2*NumX+1)+2*i+1,...
                    (2*j)*(2*NumX+1)+2*i+1,     (2*j)*(2*NumX+1)+2*i-1,...
                    (2*j-2)*(2*NumX+1)+2*i,     (2*j-1)*(2*NumX+1)+2*i+1,...
                    (2*j)*(2*NumX+1)+2*i,       (2*j-1)*(2*NumX+1)+2*i-1,...
                    (2*j-1)*(2*NumX+1)+2*i];
            end
        end
    case '2dP1' %3-node triangular element
        k=0;
        nodeCoords=zeros((NumX+1)*(NumY+1),2);
        nodeCoordsX=linspace(x1,x2,NumX+1);
        nodeCoordsY=linspace(y1,y2,NumY+1);
        for j=1:NumY+1
            for i=1:NumX+1
                k=k+1;
                nodeCoords(k,:)=[nodeCoordsX(i),nodeCoordsY(j)];
            end
        end
        
        IEN=zeros(2*NumX*NumY,3);
        k=0;
        for j=1:NumY
            for i=1:NumX
                k=k+1;
                IEN(k,:)=[(j-1)*(NumX+1)+i, ...
                    j*(NumX+1)+i+1, j*(NumX+1)+i];
                k=k+1;
                IEN(k,:)=[(j-1)*(NumX+1)+i, (j-1)*(NumX+1)+i+1,...
                    j*(NumX+1)+i+1];
            end
        end
    case '2dP2' %6-node triangular element
        k=0;
        nodeCoords=zeros((2*NumX+1)*(2*NumY+1),2);
        nodeCoordsX=linspace(x1,x2,2*NumX+1);
        nodeCoordsY=linspace(y1,y2,2*NumY+1);
        for j=1:2*NumY+1 %iterate over all rows
            for i=1:2*NumX+1 %iterate over all columns
                k=k+1;
                nodeCoords(k,:)=[nodeCoordsX(i),nodeCoordsY(j)];
            end
        end
        
        IEN=zeros(2*NumX*NumY,6);
        k=0;
        for j=1:NumY
            for i=1:NumX
                k=k+1;
                IEN(k,:)=[(2*j-2)*(2*NumX+1)+2*i-1, ...
                    (2*j)*(2*NumX+1)+2*i+1,     (2*j)*(2*NumX+1)+2*i-1,...
                    (2*j-1)*(2*NumX+1)+2*i,     (2*j)*(2*NumX+1)+2*i,...
                    (2*j-1)*(2*NumX+1)+2*i-1];
                k=k+1;
                IEN(k,:)=[(2*j-2)*(2*NumX+1)+2*i-1,     (2*j-2)*(2*NumX+1)+2*i+1,...
                    (2*j)*(2*NumX+1)+2*i+1,         ...
                    (2*j-2)*(2*NumX+1)+2*i,         (2*j-1)*(2*NumX+1)+2*i+1,...
                    (2*j-1)*(2*NumX+1)+2*i];
            end
        end
    otherwise
        error(['Wrong element type ' elementType ]);
end %switch

%% boundaryElementIDs and boundaryNodeLocalIDs
boundaryElementIDs=cell(4,1);
boundaryNodeLocalIDs=cell(4,1);
switch elementType
    case '2dQ1'
        boundaryElementIDs{1}=1:NumX;
        boundaryNodeLocalIDs{1}=repmat([1 2],[NumX,1]);
        boundaryElementIDs{2}=NumX:NumX:NumY*NumX;
        boundaryNodeLocalIDs{2}=repmat([2 3],[NumY,1]);
        boundaryElementIDs{3}=(NumY-1)*NumX+1:1:NumY*NumX;
        boundaryNodeLocalIDs{3}=repmat([3 4],[NumX,1]);
        boundaryElementIDs{4}=1:NumX:(NumY-1)*NumX+1;
        boundaryNodeLocalIDs{4}=repmat([4 1],[NumY,1]);
    case {'2dQ2r', '2dQ2'}
        boundaryElementIDs{1}=1:NumX;
        boundaryNodeLocalIDs{1}=repmat([1 5 2],[NumX,1]);
        boundaryElementIDs{2}=NumX:NumX:NumY*NumX;
        boundaryNodeLocalIDs{2}=repmat([2 6 3],[NumY,1]);
        boundaryElementIDs{3}=(NumY-1)*NumX+1:1:NumY*NumX;
        boundaryNodeLocalIDs{3}=repmat([3 7 4],[NumX,1]);
        boundaryElementIDs{4}=1:NumX:(NumY-1)*NumX+1;
        boundaryNodeLocalIDs{4}=repmat([4 8 1],[NumY,1]);
    case '2dP1'
        boundaryElementIDs{1}=2:2:2*NumX;
        boundaryNodeLocalIDs{1}=repmat([1 2],[NumX,1]);
        boundaryElementIDs{2}=2*NumX:2*NumX:2*NumX*NumY;
        boundaryNodeLocalIDs{2}=repmat([2 3],[NumY,1]);
        boundaryElementIDs{3}=(NumY-1)*2*NumX+1:2:2*NumX*NumY-1;
        boundaryNodeLocalIDs{3}=repmat([3 2],[NumX,1]);
        boundaryElementIDs{4}=1:2*NumX:(NumY-1)*2*NumX+1;
        boundaryNodeLocalIDs{4}=repmat([3 1],[NumY,1]);
    case '2dP2'
        boundaryElementIDs{1}=2:2:2*NumX;
        boundaryNodeLocalIDs{1}=repmat([1 4 2],[NumX,1]);
        boundaryElementIDs{2}=2*NumX:2*NumX:2*NumX*NumY;
        boundaryNodeLocalIDs{2}=repmat([2 5 3],[NumY,1]);
        boundaryElementIDs{3}=(NumY-1)*2*NumX+1:2:2*NumX*NumY-1;
        boundaryNodeLocalIDs{3}=repmat([3 5 2],[NumX,1]);
        boundaryElementIDs{4}=1:2*NumX:(NumY-1)*2*NumX+1;
        boundaryNodeLocalIDs{4}=repmat([3 6 1],[NumY,1]);
end %switch

for i=1:4
    boundaryElementIDs{i}=boundaryElementIDs{i}';
end

%% edges e
e=zeros(7,4);
e(1,1)=find((nodeCoords(:,1)==domain(1)).*((nodeCoords(:,2)==domain(2))),1,'first');
e(1,2)=find((nodeCoords(:,1)==domain(3)).*((nodeCoords(:,2)==domain(2))),1,'first');
e(2,1)=e(1,2);
e(2,2)=find((nodeCoords(:,1)==domain(3)).*((nodeCoords(:,2)==domain(4))),1,'first');
e(3,1)=e(2,2);
e(3,2)=find((nodeCoords(:,1)==domain(1)).*((nodeCoords(:,2)==domain(4))),1,'first');
e(4,1)=e(3,2);
e(4,2)=e(1,1);
e(5,:)=[1 2 3 4];
e(6,:)=e(6,:)+1;

end %function