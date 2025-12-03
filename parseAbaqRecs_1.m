function [t,u]=parseAbaqRecs_1(recs)
% Extracts time and displacements from Abaqus output records. It is
%   anticipated that at each increment only one displacements output is
%   produced, i.e. output was requested at a single node.
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

%%
if nargin==0
    recs=readAbaqFil();
end
Poss=cell2mat({recs.Pos}'); %k-by-1 integer array of records' start positions
NWs=cell2mat({recs.NW}'); %k-by-1 integer array of records' length
Keys=cell2mat({recs.Key}'); %k-by-1 integer array of records' keys

incStartIDs=find(Keys==uint64(2000));
numInc=numel(incStartIDs);

t=zeros(numInc,1);
u=zeros(numInc,3);

for i=1:numInc % loop over increment records
    % Extract "Total time" from Increment Start (2000-type) record
    cont=recs(incStartIDs(i)).Cont;
    t(i)=todouble(cont(1,:));
    % Find corresponding Increment End (2001-type) record
    incEndID=find( (Keys==uint64(2001) & ... % Increment End record
        Poss >= recs(incStartIDs(i)).Pos), 1, 'first'); %that follows this
    % Increment Start record
    
    % Find all nodal output records bw 2000 and 2001 type records
    nodeOutID=find( (Keys==uint64(101)) & ... %nodal output
        (Poss>=recs(incStartIDs(i)).Pos) & ... %after this increment starts
        (Poss<=recs(incEndID).Pos) );
    
    % We expect a single node output record
    if numel(nodeOutID)>1
        error('Unexpected multiple outputs for current nodal output');
    end
    
    % Parse Node Output (101-type) record
    if NWs(nodeOutID)==5 %2D case
        cont=recs(nodeOutID).Cont;
        u(i,1)=todouble(cont(2,:));
        u(i,2)=todouble(cont(3,:));
    elseif NWs(nodeOutID)==6 %3D case
        cont=recs(nodeOutID).Cont;
        u(i,1)=todouble(cont(2,:));
        u(i,2)=todouble(cont(3,:));
        u(i,3)=todouble(cont(4,:));
    else
        error('Unexpected length of 101-type record');
    end
end % loop over increment records
end