function u=import_displacements(filename)
% This ad hoc procedure parses an Abaqus data file (.dat) and extracts
% nodal displacements. 
% This procedure does not work when all displacements are zeros, because
% the dimensionality of the displacements (2d or 3d) is determined from
% non-zero rows.  
% This file worked successfully with Abaqus 6.14 Student Edition.
%
% Disclaimer: The author is not associated with owners of Abaqus
% trademark in any way.

fid=fopen(filename,'r');
    k=[];
    s='';
    while isempty(k)&&~isnumeric(s); %look for the header
        s=fgetl(fid);
        k=strfind(s,'N O D E   O U T P U T');
        l=strfind(s,'NUMBER OF NODES IS');
        if ~isempty(l)
            numN=sscanf(s(l+19:end),'%u'); % read the number of nodes.
        end
    end
    for i=1:9 %skip 9 lines
        s=fgetl(fid);
    end
    s=fgetl(fid);
    A1=sscanf(s,'%u %f %f %f');
    if numel(A1)==3
        A = fscanf(fid,'%u %f %f');
        A=reshape(A(:),3,numel(A)/3)';
        A=[A1(1), A1(2), A1(3); A];
    elseif numel(A1)==4
        A = fscanf(fid,'%u %f %f %f');
        A=reshape(A(:),4,numel(A)/4)';
        A=[A1(1), A1(2), A1(3), A1(4); A];
    end
    u=zeros(numN,size(A,2)-1);
    u(A(:,1),:)=A(:,2:end);
fclose(fid);
end
