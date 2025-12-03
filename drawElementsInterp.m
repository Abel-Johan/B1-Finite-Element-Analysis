function drawElementsInterp(nodeCoords,IEN,elementType,varargin)
% Draw the mesh using MATLAB patch() function. Values may be specified for
% each node to define the color. See also drawElements function.
%
% INPUT:
% nodeCoords - #Nodes-by-#Dim array containing coordinates of nodes;
% IEN or {IEN,boundaryElementsIDs} or BIEN, where  
%   IEN - #El-by-#NpE array, element-node incidence matrix,
%   boundaryElementIDs - #BR-by-1 cell array, p-th cell contains l_p-by-1 
%    matrix listing boundary element IDs.
%    Providing both IEN and boundaryElementsIDs is equivalent to providing
%    IEN(unique([boundaryElementIDs{:}]),:), but element labels will differ
%    if drawn.
%   BIEN - #BR-by-1 cell array containing boundary element-node incidence
%    matrices.
% elementType - a single string specifying the type of elements used;
% vals (optional) - a scalar or #Nodes-by-1 vector containing scalar values
%   to display as colors;
% alphas (optional) - a scalar or #Nodes-by-1 vector containing scalar
%   values to display as opacity;
% options (optional) - a string 'Labels', 'elementLabels',
%   'nodeLabels' or 'nodeLocalLabels'. If specified, corresponding IDs are
%   displayed.
%
% COMMENTS:
% Uses MATLAB's patch() object.
% The edges of elements are piece-wise linear and do not represent true
% geometry of higher-order elements.
% Only elements listed in IEN (or in boundaryElementsIDs, when applicable)
% are shown. 
% When BIEN is passed as the second argument, elementType should
% correspond to the boundary element type. {IEN,boundaryElementsIDs} and
% BIEN are told apart based on the fact that BIEN is a coloumn of cells.
% 
% Examples:
% drawElementsInterp(nodeCoords,IEN,elementType) 
% drawElementsInterp(nodeCoords,IEN,elementType, vals) 
% drawElementsInterp(nodeCoords,IEN,elementType, vals, alphas) 
% drawElementsInterp(nodeCoords,IEN,elementType, 0, alphas) 
% drawElementsInterp(nodeCoords,IEN,elementType,0,0) %plots transparent elements
% drawElementsInterp(nodeCoords,{IEN,boundaryElementIDs},...) 
% drawElementsInterp(nodeCoords,{IEN,{[1;2];[5;6]}},...) 


%% Presentation parameters
gcf;gca;
linestyle='-';
edgealpha=.3;
maxElLabels=99;
maxNodeLabels=99;
shiftfactor=.6; %how much to shift local node labels towards the centroid

%% Parsing input
numNodes=size(nodeCoords,1);

if iscell(IEN)&&numel(IEN)==2&&iscell(IEN{2}) % {IEN,boundaryElementsIDs} provided in place of IEN
    [elementNumbers,~,~]=unique(vertcat(IEN{2}{:}));
    IEN=IEN{1}(elementNumbers,:);
elseif iscell(IEN)&&size(IEN,2)==1 % BIEN
    IEN=cell2mat(IEN);
elseif ~iscell(IEN)
    elementNumbers=1:size(IEN,1);
else
    error('Wrong format of argument 2 (IEN)');
end
numElements=size(IEN,1);

if size(nodeCoords,2)==2
    nodeCoords=[nodeCoords zeros(numNodes,1)];
elseif size(nodeCoords,2)==3
else
    error('Only 2D and 3D supported');
end

if numel(varargin)>=1&&ischar(varargin{1})
    str=varargin{1};
    varargin(1)=[];
elseif numel(varargin)>=2&&ischar(varargin{2})
    str=varargin{2};
    varargin(2)=[];
elseif numel(varargin)>=3&&ischar(varargin{3})
    str=varargin{3};
    varargin(3)=[];
else
    str=[];
end

if ~isempty(str)&&strcmpi(str,'elementlabels')
    drawElementLabels=1;
    drawNodeLabels=0;
    drawNodeLocalLabels=0;
elseif ~isempty(str)&&strcmpi(str,'nodelabels')
    drawElementLabels=0;
    drawNodeLabels=1;
    drawNodeLocalLabels=0;
elseif ~isempty(str)&&strcmpi(str,'nodelocallabels')
    drawElementLabels=0;
    drawNodeLabels=0;
    drawNodeLocalLabels=1;
elseif ~isempty(str)&&strcmpi(str,'labels')
    drawElementLabels=1;
    drawNodeLabels=1;
    drawNodeLocalLabels=1;
else
    drawElementLabels=0;
    drawNodeLabels=0;
    drawNodeLocalLabels=0;
end


if ~isempty(varargin)
    vals=varargin{1};
    vals=vals(:);
    if isscalar(vals)
        vals=vals*ones(size(nodeCoords,1),1);
    elseif numel(vals)~=size(nodeCoords,1)
        error('Number of values does not agree with the number of nodes');
    end
else
    vals=zeros(numNodes,1);
end

if numel(varargin)>1
    alphas=varargin{2};
    if size(vals)~=size(alphas)
        error('vals and alphas must be of the same size');
    end
else
    alphas=ones(size(vals));
    alphas(isnan(vals))=0;
end

%% 
% perm=elementData(elementType,'facetsLocBIEN');
perm=elementData(elementType,'facetsLocBIEN');
if size(perm,1)>1 %strcmpi(elementType(1:2),'3d')
    faces=arrayfun(@(i){IEN(:,perm(i,:))},1:size(perm,1))';
    faces=cell2mat(faces);
else %strcmpi(elementType(1:2),'2d')
    faces=IEN(:,perm);
end

norm_vals=(vals-min(vals))/(max(vals)-min(vals)); %normalize to [0 1]
cmap=colormap(jet);
norm_vals(isnan(norm_vals))=0;
color_vals=interp1(linspace(0,1,length(cmap)),cmap,norm_vals);

if ~ishold %if hold off then erase current axis
    cla;
end 

patch('Faces',faces,...
    'Vertices',nodeCoords,...
    'FaceVertexCData',color_vals,...
    'FaceColor','interp',...
    'FaceVertexAlphaData',alphas,...
    'FaceAlpha','flat',...
    'LineStyle',linestyle,...
    'EdgeAlpha',edgealpha);

colorbar;
if max(vals)>min(vals)
    if max(vals)-min(vals)>1e-12
        colorbar('Ticks',linspace(min(vals),max(vals),10))
    end
    caxis([min(vals) max(vals)])
end

%% labels
cCoords=zeros(numElements,3);
for i=1:numElements
    cCoords(i,:)=sum(nodeCoords(IEN(i,:),:),1)/numel(IEN(i,:));
end

if drawElementLabels==1
    for i=1:numElements
        if i>maxElLabels
            disp('Maximum number of element labels reached - no more labels plotted');
            break;
        end
        text(cCoords(i,1),cCoords(i,2),cCoords(i,3),num2str(elementNumbers(i)),...
            'Color','b','FontSize',7);
    end
end

if drawNodeLabels==1
    for i=1:numNodes
        if i>maxNodeLabels
            disp('Maximum number of node labels reached - no more labels plotted');
            break;
        end
        text(nodeCoords(i,1),nodeCoords(i,2),nodeCoords(i,3),num2str(i),...
            'Color','k','FontSize',7);
    end
end

s=0;
if drawNodeLocalLabels==1
    for i=1:numElements
        if s>maxNodeLabels
            break;
        end
        for k=1:numel(IEN(i,:))
            if s>maxNodeLabels
                disp('Maximum number of node labels reached - no more labels plotted');
                break;
            end
            s=s+1;
            text(nodeCoords(IEN(i,k),1)*shiftfactor+cCoords(i,1)*(1-shiftfactor),...
                nodeCoords(IEN(i,k),2)*shiftfactor+cCoords(i,2)*(1-shiftfactor),...
                nodeCoords(IEN(i,k),3)*shiftfactor+cCoords(i,3)*(1-shiftfactor),...
                num2str(k),...
                'Color','r','FontSize',5,'HorizontalAlignment','center');
        end
    end
end
%%
axis equal
end