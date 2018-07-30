function [latlon,igbp,net] = read_metadata(sname)

% init storage
latlon = zeros(2,1)./0;
igbp = {'missing'};
net = {'missing'};

fname = 'metadata_fluxnet_ameriflux.csv';
fid = fopen(fname);

while ~feof(fid)
   tline = fgetl(fid);
   tline = strsplit(tline,'.');
   if strcmpi(sname,tline{1})
       latlon(1) = str2num(tline{2});
       latlon(2) = str2num(tline{3});
       igbp      = str2num(tline(4));
       net       = str2num(tline(5));
   end
end

fclose(fid);