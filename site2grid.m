function [lat, lon] = site2grid(s)

latlon = csvread('./data/latlon_fn_af.csv');

lat = latlon(s, 1);
lon = latlon(s, 2);

%fprintf('\nLatitude: %g   Longitude: %g \n', lat, lon);

end