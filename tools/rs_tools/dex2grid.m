function [lat, lon] = dex2grid(r, c)
assert(r <= 180);
assert(r > 0);
assert(c <= 360);
assert(c > 0);

if r <= 90
    lat = 90 - r;
else
    lat = (r - 90)*-1;
end

if c <= 180
    lon = c - 180;
else
    lon = (180 - c)*-1;
end

lat;
lon;

