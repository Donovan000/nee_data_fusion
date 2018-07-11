function [Fsi,Tsi] = sobol_combinatorial(I,Y)








% % dimensions
% [Ns,Np] = size(Y12);
% Np = Np/2;
% 
% % homogeneous stats - for later rcalculations
% mu1 = mean(Y1); mu2 = mean(Y2);
% sg1 = var(Y1);  sg2 = var(Y1);
% 
% % init storage
% D1 = zeros(Np,2);
% D2 = zeros(Np,2);
% 
% % variance integrals
% for p = 1:Np
%     for n = 1:Ns
%         D1(p,1) = D1(p,1) + Y1(n)*Y12(n,p);
%         D1(p,2) = D1(p,2) + Y1(n)*Y21(n,p);
%         D2(p,1) = D2(p,1) + Y2(n)*Y21(n,p);
%         D2(p,2) = D2(p,2) + Y2(n)*Y12(n,p);
%     end
% end
% 
% % normalize by sample size
% D1 = D1./Ns - fo^2;
% D2 = D2./Ns - fo^2;
% 
% % remove mean
% D1 = D1 - mu1^2;
% D2 = D2 - mu2^2;
% 
% % correction factor
% % Dcorrection = mean(Y(:,1).*Y(:,2)) - fo^2;
% 
% % main effects
% Fsi(:,1) = D1(:,1)./sg1; %+ Dcorrection
% Fsi(:,2) = D2(:,1)./sg2; %+ Dcorrection
% 
% % total effects
% Tsi(:,1) = 1 - D1(:,2)./sg1; %+ Dcorrection
% Tsi(:,2) = 1 - D2(:,2)./sg2; %+ Dcorrection
