function [r, c] = grid2dex(lat, lon)
assert(lat <= 90);
assert(lat >= -90);
assert(lon <= 180);
assert(lon >= -180)

if lat <= 90
    r = 90 - lat;
    r = round(r);
elseif lat <= 0
    r = abs(lat) + 90;
    r = round(r);
end

if lon <= 0
    c = 180 - abs(lon);
    c = round(c);
elseif lon > 0
    c = 180 + lon; 
    c = round(c);    
end

r;
c;