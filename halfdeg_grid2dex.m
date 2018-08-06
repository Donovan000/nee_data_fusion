function [r, c] = halfdeg_grid2dex(lat, lon)
assert(lat <= 90);
assert(lat >= -90);
assert(lon <= 180);
assert(lon >= -180)

if lat <= 90
    r = 180 - (lat*2);
    r = round(r);
elseif lat <= 0
    r = abs(lat*2) + 180;
    r = round(r);
end

if lon <= 0
    c = 360 - abs(lon*2);
    c = round(c);
elseif lon > 0
    c = 360 + (lon*2); 
    c = round(c);    
end

r;
c;