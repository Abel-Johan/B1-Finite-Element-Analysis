function recs=readAbaqFil(filename)
% Reads records from binary Abaqus Results File (.fil). Produces auxiliary
%   .fi_ file.
% INPUT:
% filename - string, file name.
%
% OUTPUT:
% recs - 1-by-k structure array:
%   recs.Pos - uint64, position of record in a processed file '.fi_';
%   recs.NW  - uint64, lenght of record;
%   recs.Key - uint64, record key;
%   recs.NW  - n(k)-by-1 uint64, record content.
% See Abaqus Analysis User's Guide Chapter 5.1.2 for details
%
% Disclaimer: The author is not associated with owners of Abaqus
% trademark in any way.

%% create a temporary file '.fi_' without binary frame marks
fprintf(1,'Reading Abaqus output from %s ....', filename);

fid=fopen(filename,'r');
fname2=[filename(1:end-1) '_'];
fid2=fopen(fname2,'w');

x=1;
while ~isempty(x)
    framesize=fread(fid,1,'*uint32');
    if isempty(framesize)
        break
    end
    x=fread(fid,framesize,'*uint8');
    fwrite(fid2,x,'*uint8');
    framesize2=fread(fid,1,'*uint32');
    if framesize2~=framesize
        error('Unexpected binary file format');
    end
end
fclose(fid);
fclose(fid2);

%% Reading records from the temporary file '_.fil'
fid=fopen(fname2,'r');

% initialise struct array and index
recs=struct('Pos',[],'NW',[],'Key',[],'Cont',[]);
i=0;

status=0;
while status>=0 %while not EoF or something like that
    pos=ftell(fid)+1;
    recNW=fread(fid,1,'*uint64');
    % if unable to read further (i.e. EoF already reached), fread returns empty
    % array. If can read, but not all requested length (i.e. EoF reached in
    % the process), fread returns a truncated array (such check is
    % performed with recCont).
    if isempty(recNW)
        status=-1;
        continue;
    end
    recKey=fread(fid,1,'*uint64');
    if isempty(recNW)
        status=-1;
        continue;
    end
    % read record
    recCont=fread(fid,recNW-2,'*uint64');
    if size(recCont)<recNW-2
        status=-1;
        continue;
    end
    i=i+1;
    recs(i).Pos=pos;
    recs(i).NW=recNW;
    recs(i).Key=recKey;
    recs(i).Cont=recCont;
end
fclose(fid);
fprintf(1,'done\n');
fprintf(1, ['File %s created\n'], fname2);

end