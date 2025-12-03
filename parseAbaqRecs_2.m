function [t,u,v,a,rf,incNums]=parseAbaqRecs_2(recs)
% Extracts time and displacements from Abaqus output records. It is
%   anticipated that at each increment output contains displacements,
%   velocities and accelerations at all nodes.
% INPUT:
% recs - 1-by-k structure array:
%   recs.Pos - uint64, position of record in a processed file '.fi_';
%   recs.NW  - uint64, lenght of record;
%   recs.Key - uint64, record key;
%   recs.NW  - n(k)-by-1 uint64, record content.
% See Abaqus Analysis User's Guide Chapter 5.1.2 for details
%
% Disclaimer: The author is not associated with owners of Abaqus
% trademark in any way.

%% functions
    function res=todouble(x) 
    % used to reintepret uint64 as double without changing its binary data
        res=typecast(x,'double');
    end
    function res=tochar(x)
        res=typecast(x,'uint8');
        res=char(res)';
    end
%%
fprintf(1,'Parsing Abaqus output....');

Poss=cell2mat({recs.Pos}'); %k-by-1 integer array of records' start positions
NWs=cell2mat({recs.NW}'); %k-by-1 integer array of records' length
Keys=cell2mat({recs.Key}'); %k-by-1 integer array of records' keys

numNodes=sum(Keys==uint64(1901));
incStartIDs=find(Keys==uint64(2000));
numInc=numel(incStartIDs);

t=zeros(numInc,1);
incNums=zeros(numInc,1);
u=zeros(numInc,numNodes,3);
v=zeros(numInc,numNodes,3);
a=zeros(numInc,numNodes,3);
rf=zeros(numInc,numNodes,3);

for i=1:numInc % loop over increment records
    % Extract "Total time" from Increment Start (2000-type) record
    cont=recs(incStartIDs(i)).Cont;
    t(i)=todouble(cont(1,:));
    incNums(i)=cont(7,:);
    % Find corresponding Increment End (2001-type) record
    incEndID=find( (Keys==uint64(2001) & ... % Increment End record
        Poss >= recs(incStartIDs(i)).Pos), 1, 'first'); %that follows this
    % Increment Start record
    
    % Find all nodal output records bw 2000 and 2001 type records
    nodeOutIDs=find( (Keys==uint64(101)) & ... %nodal output DISPLACEMENTS
        (Poss>=recs(incStartIDs(i)).Pos) & ... %after this increment starts
        (Poss<=recs(incEndID).Pos) );
    for j=1:numel(nodeOutIDs)
        nodeOutID=nodeOutIDs(j);
        % Parse Node Output (101-type) record
        if NWs(nodeOutID)==5 %2D case
            cont=recs(nodeOutID).Cont;
            u(i,j,1)=todouble(cont(2,:));
            u(i,j,2)=todouble(cont(3,:));
        elseif NWs(nodeOutID)==6 %3D case
            cont=recs(nodeOutID).Cont;
            u(i,j,1)=todouble(cont(2,:));
            u(i,j,2)=todouble(cont(3,:));
            u(i,j,3)=todouble(cont(4,:));
        else
            error('Unexpected length of 101-type record');
        end
    end
    
    % Find all nodal output records bw 2000 and 2001 type records
    nodeOutIDs=find( (Keys==uint64(102)) & ... %nodal output VELOCITIES
        (Poss>=recs(incStartIDs(i)).Pos) & ... %after this increment starts
        (Poss<=recs(incEndID).Pos) );
    for j=1:numel(nodeOutIDs)
        nodeOutID=nodeOutIDs(j);
        % Parse Node Output (101-type) record
        if NWs(nodeOutID)==5 %2D case
            cont=recs(nodeOutID).Cont;
            v(i,j,1)=todouble(cont(2,:));
            v(i,j,2)=todouble(cont(3,:));
        elseif NWs(nodeOutID)==6 %3D case
            cont=recs(nodeOutID).Cont;
            v(i,j,1)=todouble(cont(2,:));
            v(i,j,2)=todouble(cont(3,:));
            v(i,j,3)=todouble(cont(4,:));
        else
            error('Unexpected length of 101-type record');
        end
    end
    
    % Find all nodal output records bw 2000 and 2001 type records
    nodeOutIDs=find( (Keys==uint64(103)) & ... %nodal output ACCELERATIONS
        (Poss>=recs(incStartIDs(i)).Pos) & ... %after this increment starts
        (Poss<=recs(incEndID).Pos) );
    for j=1:numel(nodeOutIDs)
        nodeOutID=nodeOutIDs(j);
        % Parse Node Output (101-type) record
        if NWs(nodeOutID)==5 %2D case
            cont=recs(nodeOutID).Cont;
            a(i,j,1)=todouble(cont(2,:));
            a(i,j,2)=todouble(cont(3,:));
        elseif NWs(nodeOutID)==6 %3D case
            cont=recs(nodeOutID).Cont;
            a(i,j,1)=todouble(cont(2,:));
            a(i,j,2)=todouble(cont(3,:));
            a(i,j,3)=todouble(cont(4,:));
        else
            error('Unexpected length of 101-type record');
        end
    end

    % Find all nodal output records bw 2000 and 2001 type records
    nodeOutIDs=find( (Keys==uint64(104)) & ... %nodal output REACTION FORCES
        (Poss>=recs(incStartIDs(i)).Pos) & ... %after this increment starts
        (Poss<=recs(incEndID).Pos) );
    for j=1:numel(nodeOutIDs)
        nodeOutID=nodeOutIDs(j);
        % Parse Node Output (101-type) record
        if NWs(nodeOutID)==5 %2D case
            cont=recs(nodeOutID).Cont;
            rf(i,j,1)=todouble(cont(2,:));
            rf(i,j,2)=todouble(cont(3,:));
        elseif NWs(nodeOutID)==6 %3D case
            cont=recs(nodeOutID).Cont;
            rf(i,j,1)=todouble(cont(2,:));
            rf(i,j,2)=todouble(cont(3,:));
            rf(i,j,3)=todouble(cont(4,:));
        else
            error('Unexpected length of 101-type record');
        end
    end

end % loop over increment records
fprintf(1,'done\n');

end