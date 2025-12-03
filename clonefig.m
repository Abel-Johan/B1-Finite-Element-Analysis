function clonefig(varargin)
% Clone a figure by creating an exact copy of it
%
% INPUT:
% i (optional) - graphic object identifier. See copyobj reference for more
%   details.
%
% Examples:
%   clonefig % clones current figure;
%   clonefig(i) % clones figure i;
%   clonefig(0) % clones all figures;

if ~isempty(varargin)
    if isempty(varargin{1})||varargin{1}==0 %clone all figures
        gr=groot;
        ch=gr.Children;
        copyobj(ch,groot)
    else %clone
        copyobj(varargin{1},groot)
    end
else %clone current figure
    copyobj(gcf,groot)
end
end%function