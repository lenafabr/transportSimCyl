function bindposdist = readbindpos(fname)
%%% reads the cargo binding positions from the file located at fname
% fname must be a *.bindpos.out file

fileID = fopen(fname,'r');
siminfo = textscan(fileID,'%f %f',1);
npart = siminfo{1};
nread = siminfo{2}*3;

fmtstr = strtrim(repmat('%f ',1,npart));
data = textscan(fileID,fmtstr,nread);

bindposdist = cat(2,data{:})';
bindposdist = reshape(bindposdist,npart,3,siminfo{2});

fclose(fileID);

end