function fptdist = readfpt(fname)
%%% reads the cargo capture time distribution from the file located at fname
% fname must be a *.fpt.out file

fileID = fopen(fname,'r');
siminfo = textscan(fileID,'%f %f',1);
npart = siminfo{1};
ntrials = siminfo{2};

fmtstr = strtrim(repmat('%f ',1,npart));
data = textscan(fileID,fmtstr,ntrials);
fptdist = cat(2,data{:});
fclose(fileID);


end