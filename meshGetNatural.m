function [el, ptNatCoords] = meshGetNatural(nodeCoords,IEN,elementType,INE,ptCoords)
% Map global coordinates to natural coordinates.
%
% INPUT:
% nodeCoords - #Nodes-by-#Dim matrix, coordinates of nodes.
% IEN - #El-by-#NpE array, element-node incidence matrix. 
% elementType - single string, the type of elements.
% ptCoords - k-by-#Dim matrix, row ptCoords(i,:) contains global
%   coordinates of the query point i. 
%
% OUTPUT:
% el - k-by-1 matrix, element IDs of the query points.
% ptNatCoords - k-by-#ElDim matrix, natural coordinates of the query points.
%
% COMMENTS:
% This routine uses Newton's iterative method to find the pre-image in
% higher-order elements. If a preimage of some particular point was not
% found, a warning message is generated and the corresponding entry in
% the output array el is zeros.
% This routine uses the 'inElementTol' property of an element, which allows
% to tell whether a requested point belongs to the element's natural domain
% up to some tolerance. See elementData.m for details

% The following parameters may be tweaked to improve performance
maxiter=16; %maximum number of iterations in Newton method
atol=1e-4; %absolute tolerance
%% 
numDim=size(nodeCoords,2);
numNodes=size(nodeCoords,1);
numPts=size(ptCoords,1);

if exist(['shapeFunctions' elementType],'file')==2
    shapeFunctions=str2func(['shapeFunctions' elementType]);
else
    error(['Can''t find shape functions for specified element type: "' ...
        elementType '"' ]);
end
inElement=elementData(elementType,'inElementTol');
%% subroutine
    function [curNatCoords,aerr,i]=getNatural(elementNodeCoords,curPt)
        % Trying to solve for natural coordinates within a given element
        if elementType(1)=='2'
            xi=[0 0];
        elseif elementType(1)=='3'
            xi=[0 0 0];
        elseif elementType(1)=='1'
            xi=[0];
        end
        aerr=Inf;
        i=0;
        % Newton's iteration
        while (i<=maxiter)&&(atol<=aerr)
            [shapeFunctionsVals,naturalDerivatives]=shapeFunctions(xi);
            J=elementNodeCoords'*naturalDerivatives';
            F=shapeFunctionsVals*elementNodeCoords-curPt;
            dx=(J\F')';
            xi=xi-dx;
            i=i+1;
            aerr=norm(dx);
        end
        curNatCoords=xi;
    end %subfunction getNatural()

%% main 
el=zeros(numPts,1);
ptNatCoords=zeros(numPts,numDim);
% We iterate over the given points, whose preimages are to be found. For
% each point we find the node closest to it. Most likely the point belongs
% to an element containing the closest node. Less likely the point belongs
% to an element adjacent to an element containing the closest node. In rare
% cases neither of these hold (e.g. in distorted meshes or when the point
% is near a boundary where a self contact is about to happen), but we
% ignore this and only consider the first two possibilities. Within each of
% candidate element we attempt to find the preimage using Newton's
% iteration.

extended=0; % only look among elements containing the node closest to the given point 
i=1;
fail=0;
while i<=numPts
    if extended==0
        curPt=ptCoords(i,:);
        % determine candidate elements - those incident to the closest
        % point
        dist=sum((nodeCoords-ones(numNodes,1)*curPt).^2,2);
        [~,clNodeID]=min(dist);
        cand=nonzeros(INE(clNodeID,:));
    else
        % consider neighbouring elements of the elements that contain the
        % closest node
        cand2=unique(nonzeros(INE(unique(IEN(cand,:)),:)));
        cand=setdiff(cand2,cand); %excluding candidates that did not work
    end
    for j=1:numel(cand)
        curCand=cand(j);
        elementNodeCoords=nodeCoords(IEN(curCand,:),:);
        %Newton's method to compute natural coordinates
        [curNatCoords,aerr,~]=getNatural(elementNodeCoords,curPt);
        
        if (aerr<atol)&&(inElement(curNatCoords,atol)) %found it!
            break; 
        end
    end
    if (aerr<atol)&&(inElement(curNatCoords,atol)) %found it
        el(i)=cand(j);
        ptNatCoords(i,:)=curNatCoords;
        i=i+1;
        extended=0; % reset to search withing closest elements
    elseif extended==0 
        extended=1;
    else %extended=1;
        el(i)=0; %found none
        fail=1;
        i=i+1;
        extended=0;
    end

    if fail==1
        warning('Couldn''t find preimages for some points');
    end
end %loop over pts

end%function