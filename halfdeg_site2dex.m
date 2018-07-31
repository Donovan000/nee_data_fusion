function [r, c] = halfdeg_site2dex(s)

[lat, lon] = site2grid(s);
[r, c] = halfdeg_grid2dex(lat, lon);

r;
c;