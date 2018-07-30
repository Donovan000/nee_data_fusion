function [latlon,igbp] = read_metadata(sname,network)

% init storage
latlon = zeros(2,1)./0;
igbp = {'missing'};
net = {'missing'};

fname = 'metadata_fluxnet_ameriflux.csv';
fid = fopen(fname);

while ~feof(fid)
   tline = fgetl(fid);
   tline = strsplit(tline,',');
   if strcmpi(sname,tline{1}) && strcmpi(network,tline{5})
       latlon(1) = str2num(tline{2});
       latlon(2) = str2num(tline{3});
       igbp      = tline{4};
       break
   end
end

fclose(fid);