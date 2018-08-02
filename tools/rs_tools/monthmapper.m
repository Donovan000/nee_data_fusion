function [month] = monthmapper(year, day)

if mod(year, 4) == 0
    day = day - 1;
end

if day <= 31
    month = 1;
elseif day <= 59
    month = 2;
elseif day <= 90
    month = 3;
elseif day <= 120
    month = 4;
elseif day <= 151
    month = 5;
elseif day <= 181
    month = 6;
elseif day <= 212
    month = 7;
elseif day <= 243
    month = 8;
elseif day <= 273
    month = 9;
elseif day <= 304
    month = 10;
elseif day <= 334
    month = 11;
elseif day <= 365
    month = 12;
end