clear all
close all
clc

% open file
fid = fopen('metadata_raw.txt');

% read line-by-line
l = 1;
tline{l} = fgetl(fid);
while ischar(tline{l})
    if strcmpi(tline{l},''); tline{l} = '-9999'; end
    l = l+1;
    tline{l} = fgetl(fid);
end
tline(end) = [];

% close file
fclose(fid);

% number of lines
Nl = length(tline);
Ns = Nl/8;

% reshape array
tline = reshape(tline,[8,Nl/8])';
tline(:,2) = [];

% write new file
fname = 'metadata_organized.csv';
format = strcat(repmat('%s,',1,size(tline,2)),'\n');
fid = fopen(fname,'w');
for s = 1:Ns
    fprintf(fid,format,tline{s,:});
end
fclose(fid);