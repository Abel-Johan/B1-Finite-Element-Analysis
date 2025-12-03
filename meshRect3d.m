function [nodeCoords, IEN, boundaryElementIDs, boundaryNodeLocalIDs] = meshRect3d(domain, elementType, elementSize)
% Create a mesh for a cuboid 3D domain.
%
% INPUT:
% domain (optional) - 1-by-4 matrix [x1 y1 z1 x2 y2 z2] or 2-by-2 
%   [x1 y1 z1; x2 y2 z2] with coordinates defining the domain. Default
%   value - [-2 -2 -2 2 2 2]. 
% elementType (optional) - single string, the type of elements. Default
%   value - '3dQ1'.
% elementSize (optional) - a single real or 3-by-1 matrix. A single real
%   defines the maximum diameter of elements. A 3-by-1 matrix defines the
%   number of "cells" in each dimension of the grid. Each "cell" is a
%   hexahedron or 5 tetrahedra depending on the element type.
%   Default value - sqrt(2). 
%
% OUTPUT:
% nodeCoords - #Nodes-by-#Dim matrix, coordinates of nodes. #Dim=3.
% IEN - #El-by-#NpE array, element-node incidence matrix. 
% boundaryElementIDs - 6-by-1 cell array: p-th cell contains l_p-by-1
%   matrix listing boundary element IDs. #BR=6.
% boundaryNodeLocalIDs - 6-by-1 cell array: p-th cell contains l_p-by-q
%   matrix listing local IDs (within an element) of corresponding boundary
%   nodes.
%
% COMMENTS:
% The supported element types are '3dQ1', '3dP1'.
% Within each row IEN(i,:), vertices are listed before the midpoints.
% The six cells in boundaryElementIDs and boundaryNodeLocalIDs correspond
% respectively to BOTTOM, TOP, LEFT, RIGHT, FRONT and BACK boundaries of
% the domain, i.e. corresponding to -z, +z, -x, +x, -y, +y directions.
%
% Examples:
% [nodeCoords, IEN, bndElIDs, bndNodeLocIDs] = ...
%   meshRect3d(domain, elementType, elementSize)
% [nodeCoords, IEN, bndElIDs, bndNodeLocIDs] = ... 
%   meshRect3d([0 0 0; 1 1 1], '3dQ1', [10 10 10]) 
% [nodeCoords, IEN, bndElIDs, bndNodeLocIDs] = ...
%   meshRect3d([0 0 0 1 1 1], '3dQ1', .1) 

%% Parsing Input
% Assigning default values if some input is missing
if nargin < 1
    domain = [-2 -2 -2 2 2 2];
elseif (size(domain,1)==2)&&(size(domain,2)==3)
    domain=domain'; % for linear indexing.
end
if nargin < 2
    elementType = '3dQ1';
end
if nargin < 3
    elementSize = sqrt(2);
end

% Checking validity of the input
if ~( (numel(elementSize)==1 && elementSize>0) || ...
        (numel(elementSize)==3 && elementSize(1) >= 1 ...
        && elementSize(2) >= 1 && elementSize(3) >= 1) )
    error('Element size must be a positive real number or three positive integer real numbers');
end

switch elementType
    case {'3dQ1'}
    case {'3dP1'}
    otherwise
        error('Element type is not supported');
end

if (domain(4)-domain(1))*(domain(5)-domain(2))*(domain(6)-domain(3))==0
    error('The specified domain is degenerate');
end

%%
elData=elementData(elementType);
elLocBIEN=elData.elLocBIEN;
nodeNatCoords=elData.nodeNatCoords;

%% Coordinates
% The following code is only applicable to rectangular domain
% NumX, NumY - number of cells in each direction, chosen in a way that
% guarantees a bound on the element diameter. Each cell corresponds to one
% rectangular or to two triangular elements.
x1=domain(1);y1=domain(2);z1=domain(3);
x2=domain(4);y2=domain(5);z2=domain(6);
if numel(elementSize)==1
    NumX = ceil(2^.5*abs(x2-x1)/elementSize);
    NumY = ceil(2^.5*abs(y2-y1)/elementSize);
    NumZ = ceil(2^.5*abs(z2-z1)/elementSize);
elseif numel(elementSize)==3
    NumX=elementSize(1);
    NumY=elementSize(2);
    NumZ=elementSize(3);
end

%% Nodes coordinates and IEN
switch elementType
    case '3dQ1' %8-node hexahedral element
        s=0;
        nodeCoords=zeros((NumX+1)*(NumY+1)*(NumZ+1),3);
        nodeCoordsX=linspace(x1,x2,NumX+1);
        nodeCoordsY=linspace(y1,y2,NumY+1);
        nodeCoordsZ=linspace(z1,z2,NumZ+1);
        for k=1:NumZ+1
            for j=1:NumY+1
                for i=1:NumX+1
                    s=s+1;
                    nodeCoords(s,:)=...
                        [nodeCoordsX(i),nodeCoordsY(j),nodeCoordsZ(k)];
                end
            end
        end
        
        IEN=zeros(NumX*NumY*NumZ,size(nodeNatCoords,1));
        s=0;
for k=1:NumZ
for j=1:NumY
for i=1:NumX
    s=s+1;
    IEN(s,:)=[(k-1)*(NumY+1)*(NumX+1)+(j-1)*(NumX+1)+i,...
        (k-1)*(NumY+1)*(NumX+1)+(j-1)*(NumX+1)+i+1,...
        (k-1)*(NumY+1)*(NumX+1)+(j)*(NumX+1)+i+1,...
        (k-1)*(NumY+1)*(NumX+1)+(j)*(NumX+1)+i,...
        (k)*(NumY+1)*(NumX+1)+(j-1)*(NumX+1)+i,...
        (k)*(NumY+1)*(NumX+1)+(j-1)*(NumX+1)+i+1,...
        (k)*(NumY+1)*(NumX+1)+(j)*(NumX+1)+i+1,...
        (k)*(NumY+1)*(NumX+1)+(j)*(NumX+1)+i];
end
end
end
% end case '3dQ1'
    case '3dP1' %4-node tetrahedral element
        s=0;
        nodeCoords=zeros((NumX+1)*(NumY+1)*(NumZ+1),3);
        nodeCoordsX=linspace(x1,x2,NumX+1);
        nodeCoordsY=linspace(y1,y2,NumY+1);
        nodeCoordsZ=linspace(z1,z2,NumZ+1);
        for k=1:NumZ+1
            for j=1:NumY+1
                for i=1:NumX+1
                    s=s+1;
                    nodeCoords(s,:)=...
                        [nodeCoordsX(i),nodeCoordsY(j),nodeCoordsZ(k)];
                end
            end
        end
        
        IEN=zeros(5*NumX*NumY*NumZ,size(nodeNatCoords,1));
        s=0;
for k=1:NumZ
for j=1:NumY
for i=1:NumX
    IEN(s+1,:)=[...
        (k-1)*(NumY+1)*(NumX+1)+(j-1)*(NumX+1)+i+1,...%2
        (k-1)*(NumY+1)*(NumX+1)+(j)*(NumX+1)+i,...%4
        (k)*(NumY+1)*(NumX+1)+(j-1)*(NumX+1)+i,...%5
        (k-1)*(NumY+1)*(NumX+1)+(j-1)*(NumX+1)+i,...%1
        ];
    IEN(s+2,:)=[...
        (k-1)*(NumY+1)*(NumX+1)+(j)*(NumX+1)+i,...%4
        (k-1)*(NumY+1)*(NumX+1)+(j-1)*(NumX+1)+i+1,...%2
        (k)*(NumY+1)*(NumX+1)+(j)*(NumX+1)+i+1,...%7
        (k-1)*(NumY+1)*(NumX+1)+(j)*(NumX+1)+i+1,...%3
        ];
    IEN(s+3,:)=[...
        (k)*(NumY+1)*(NumX+1)+(j)*(NumX+1)+i+1,...%7
        (k)*(NumY+1)*(NumX+1)+(j-1)*(NumX+1)+i,...%5
        (k-1)*(NumY+1)*(NumX+1)+(j)*(NumX+1)+i,...%4
        (k)*(NumY+1)*(NumX+1)+(j)*(NumX+1)+i...%8
        ];
    IEN(s+4,:)=[...
        (k)*(NumY+1)*(NumX+1)+(j-1)*(NumX+1)+i,...%5
        (k)*(NumY+1)*(NumX+1)+(j)*(NumX+1)+i+1,...%7
        (k-1)*(NumY+1)*(NumX+1)+(j-1)*(NumX+1)+i+1,...%2
        (k)*(NumY+1)*(NumX+1)+(j-1)*(NumX+1)+i+1,...%6
        ];
    IEN(s+5,:)=[...
        (k)*(NumY+1)*(NumX+1)+(j-1)*(NumX+1)+i,...%5
        (k-1)*(NumY+1)*(NumX+1)+(j)*(NumX+1)+i,...%4
        (k-1)*(NumY+1)*(NumX+1)+(j-1)*(NumX+1)+i+1,...%2
        (k)*(NumY+1)*(NumX+1)+(j)*(NumX+1)+i+1,...%7
        ];
    s=s+5;
end
end
end
end %switch

%% Boundary
boundaryElementIDs=cell(6,1);
boundaryNodeLocalIDs=cell(6,1);
switch elementType
    case '3dQ1'
        boundaryElementIDs{1}=1:NumX*NumY;
        boundaryNodeLocalIDs{1}=repmat(elLocBIEN(1,:),[NumX*NumY,1]);
        boundaryElementIDs{2}=(1:1:NumX*NumY)+NumX*NumY*(NumZ-1);
        boundaryNodeLocalIDs{2}=repmat(elLocBIEN(2,:),[NumX*NumY,1]);

        boundaryElementIDs{3}=1:NumX:NumX*NumY*NumZ;
        boundaryNodeLocalIDs{3}=repmat(elLocBIEN(3,:),[NumY*NumZ,1]);
        boundaryElementIDs{4}=(1:NumX:NumX*NumY*NumZ)+(NumX-1);
        boundaryNodeLocalIDs{4}=repmat(elLocBIEN(4,:),[NumY*NumZ,1]);

        temp1=repmat((1:NumX),NumZ,1);
        temp2=repmat((0:NumX*NumY:NumX*NumY*(NumZ-1))',1,NumX);
        temp3=temp1+temp2;
        boundaryElementIDs{5}=temp3(:)';
        boundaryNodeLocalIDs{5}=repmat(elLocBIEN(5,:),[NumX*NumZ,1]);

        temp3=temp1+temp2+NumX*(NumY-1);
        boundaryElementIDs{6}=temp3(:)';
        boundaryNodeLocalIDs{6}=repmat(elLocBIEN(6,:),[NumX*NumZ,1]);
    case '3dP1'
        shft=1:NumX*NumY; %bottom
        boundaryElementIDs{1}=[1+5*(shft-1);2+5*(shft-1)];
        boundaryElementIDs{1}=boundaryElementIDs{1}(:)';
        boundaryNodeLocalIDs{1}=repmat(elLocBIEN(3,:),[2*NumX*NumY,1]);
        
        shft=(1:1:NumX*NumY)+NumX*NumY*(NumZ-1); %top
        boundaryElementIDs{2}=[3+5*(shft-1);4+5*(shft-1)];
        boundaryElementIDs{2}=boundaryElementIDs{2}(:)';
        boundaryNodeLocalIDs{2}=repmat(elLocBIEN(3,:),[2*NumX*NumY,1]);

        shft=1:NumX:NumX*NumY*NumZ; %left
        boundaryElementIDs{3}=[3+5*(shft-1);1+5*(shft-1)];
        boundaryElementIDs{3}=boundaryElementIDs{3}(:)';
        boundaryNodeLocalIDs{3}=repmat(elLocBIEN(1,:),[2*NumY*NumZ,1]);

        shft=(1:NumX:NumX*NumY*NumZ)+(NumX-1); %right
        boundaryElementIDs{4}=[2+5*(shft-1);4+5*(shft-1)];
        boundaryElementIDs{4}=boundaryElementIDs{4}(:)';
        boundaryNodeLocalIDs{4}=repmat(elLocBIEN(1,:),[2*NumY*NumZ,1]);

        temp1=repmat((1:NumX),NumZ,1); %front
        temp2=repmat((0:NumX*NumY:NumX*NumY*(NumZ-1))',1,NumX);
        temp3=temp1+temp2;
        shft=temp3(:)';
        boundaryElementIDs{5}=[1+5*(shft-1);4+5*(shft-1)];
        boundaryElementIDs{5}=boundaryElementIDs{5}(:)';
        boundaryNodeLocalIDs{5}=repmat(elLocBIEN(2,:),[2*NumX*NumZ,1]);

        temp3=temp1+temp2+NumX*(NumY-1); %back
        shft=temp3(:)';
        boundaryElementIDs{6}=[3+5*(shft-1);2+5*(shft-1)];
        boundaryElementIDs{6}=boundaryElementIDs{6}(:)';
        boundaryNodeLocalIDs{6}=repmat(elLocBIEN(2,:),[2*NumX*NumZ,1]);
end %switch

for i=1:6
    boundaryElementIDs{i}=boundaryElementIDs{i}';
end

end %function