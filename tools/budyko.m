function [di,ef] = budyko(X)

% grandmas
cols = 1:8;
I = all(~isnan(X(:,cols)'));
assert(length(I) == size(X,1));
I = find(I);
if length(I) < 2*365
    ef = 0/0;
    di = 0/0;
    return
end

% extract variables
Pp = X(I,1);            % precipitation                 [mm/d]
Ta = X(I,2);            % air temperature               [C]
Pa = X(I,3);            % air pressure                  [kPa]
Ws = X(I,4);            % wind speed                    [m/s]
Vp = X(I,5);            % vapor pressure deficit        [hPa]
Rn = X(I,6);            % net radiation                 [W/m2]
Qe = X(I,7);            % latent heat                   [W/m2]
Qh = X(I,8);            % sensible heat                 [W/m2]

% % slope of saturation vapor pressure curve
% expterm = (17.27*Ta) ./ (Ta+237.3);
% numerat = 4098 * (0.6108 * exp(expterm));
% denomin = (Ta+237.3).^2;
% m  = numerat./denomin;  %                               [~]
% 
% % potential evap - penman
% Rn = Rn * 0.0864;       % net radiation                 [MJ/m2/d]
% gv = 2264.705 /1e3;     % latent heat of vaporization   [MJ/kg]
% g  = 0.0016286*Pa/gv;   % psychormetric constant        [Pa/K]
% de = Vp/1e1;            % vapor pressure deficit        [kPa]
% 
% Ep = m.*Rn + g.*6.43.*(1+0.536.*Ws).*de;
% Ep = Ep ./ gv.*(m+g);    %                              [mm/d]

% potential evap from Milly & Dunne (2016)
Ep = 0.8*(Qe+Qh);
Ep = Ep * 0.0864;       %                               [MJ/m2/d]
Ep = Ep / 2.45;         %                               [mm/d]

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