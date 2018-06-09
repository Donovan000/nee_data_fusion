function [mon,doy] = doy2month(y,d)

% Define year in months and days
Ndays = [0,31,28,31,30,31,30,31,31,30,31,30,31];
if rem(y,4) == 0 
    Ndays(3) = 29; 
end

% Assert a real date
assert(y >= 1901 && y <= 2099);
assert(y == floor(y));
assert(d == floor(d));
assert(d <= sum(Ndays));


% Define unique months
mon = find(cumsum(Ndays)<d,1,'last');

% Identify day of month
doy = d-sum(Ndays(1:mon));
