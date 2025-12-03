function snapfig
% Take a snapshot of the current figure and save it in a .png file
% 
% Examples:
%   snapfig;

%determine the first available filename
i=1;
x=dir('*.png');
if ~isempty(x)
    x=struct2table(x);
    x=cellstr(x.name);
    while nnz(  strcmp(  x , ['snapshot' num2str(i,'%3.0f') '.png']  )  )
        i=i+1;
    end
end
%make a snapshot
set(gcf, 'PaperPositionMode', 'auto');
print('-dpng', ['snapshot' num2str(i,'%3.0f') '.png']);
end