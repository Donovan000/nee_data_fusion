function [r, c] = site2dex(s)

[lat, lon] = site2grid(s);
[r, c] = grid2dex(lat, lon);

r;
c;
