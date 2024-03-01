function [alp_, bet_] = calcCoeff(phiBound, thetaBound, thetaN, n, lam, mu)
% ----------------------------------------------------------------------
%   calculate 1st-order approx. coefficients
%    20210124  y.yoshimura
%    Inputs:
%   Outputs:
%   related function files:
%   note:
%   cf:
%   revisions;
%
%   (c) 2021 yasuhiro yoshimura
%----------------------------------------------------------------------

%% pre-allocation
nFlag = zeros(size(phiBound)-1); 
a = zeros(size(phiBound)-1); 
b = zeros(size(phiBound)-1); 
alp_ = zeros(size(phiBound)-1); 
bet_ = zeros(size(phiBound)-1); 

%% calc
for i = 1:size(phiBound,1)-1
    for k = 1:size(phiBound,2)-1
        nFlag(i,k) = (phiBound(i,k) <= pi/2) * (pi/2 <= phiBound(i,k+1));
        nFlag(i,k) = nFlag(i,k) * (thetaBound(i,k) <= thetaN) * (thetaN <= thetaBound(i+1,k));
        [a(i,k), b(i,k)] = intBound([phiBound(i,k),phiBound(i,k+1)], [thetaBound(i,k),thetaBound(i+1,k)], n, nFlag(i,k));
        p = (a(i,k) + b(i,k)) / 2;
        q = (b(i,k) - a(i,k)) / 2;
        r0 = mu * exp(lam*(p - 1)) * sinh(lam * q) /(lam * q);
        r1 = 3 * mu * exp(lam*(p - 1)) * (lam * q * cosh(lam*q) - sinh(lam*q)) / (lam^2 * q^2);
        alp_(i,k) = r0 - r1 * (a(i,k) + b(i,k)) / (b(i,k) - a(i,k)); % alpha
        bet_(i,k) = 2 * r1 / (b(i,k) - a(i,k)); % beta
    end
end

end