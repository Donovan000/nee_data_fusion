function [di,ef] = budyko(X)

% grandmas
cols = 3:9;
I = all(~isnan(X(:,cols)'));
assert(length(I) == size(X,1));
I = find(I);

% extract variables
Pp = X(I,3);            % precipitation                 [mm/d]
Ta = X(I,4);            % air temperature               [C]
Pa = X(I,5);            % air pressure                  [kPa]
Ws = X(I,6);            % wind speed                    [m/s]
Vp = X(I,7);            % vapor pressure deficit        [hPa]
Rn = X(I,8);            % net radiation                 [W/m2]
Qe = X(I,9);            % latent heat                   [W/m2]

% slope of saturation vapor pressure curve
expterm = (17.27*Ta) ./ (Ta+237.3);
numerat = 4098 * (0.6108 * exp(expterm));
denomin = (Ta+237.3).^2;
m  = numerat./denomin;  %                               [~]

% potential evap - penman
Rn = Rn * 0.0864;       % net radiation                 [MJ/m2/d]
gv = 2264.705 /1e3;     % latent heat of vaporization   [MJ/kg]
g  = 0.0016286*Pa/gv;   % psychormetric constant        [Pa/K]
de = Vp/1e1;            % vapor pressure deficit        [kPa]

Ep = m.*Rn + g.*6.43.*(1+0.536.*Ws).*de;
Ep = Ep ./ gv.*(m+g);    %                              [mm/d]

% actual evap
Ea = Qe * 0.0864;       %                               [MJ/m2/d]
Ea = Ea / 2.45;         %                               [mm/d]

% dryness index
pp = sum(Pp);
ep = sum(Ep);
di = ep/pp;

% evaporative fraction
ea = sum(Ea);
ef = ea/pp;