function [a, b] = intBound(phi, theta, n, lobeIn)
% ----------------------------------------------------------------------
%   calculte the integration bound
%    20201201  y.yoshimura
%    Inputs:
%   Outputs:
%   related function files:
%   note:
%   cf:
%   revisions;
%
%   (c) 2020 yasuhiro yoshimura
%----------------------------------------------------------------------

phi0 = phi(1);
phi1 = phi(2);
theta0 = theta(1);
theta1 = theta(2);

c0 = [sin(theta0)*cos(phi0), sin(theta0)*sin(phi0), cos(theta0)];
c1 = [sin(theta0)*cos(phi1), sin(theta0)*sin(phi1), cos(theta0)];
c2 = [sin(theta1)*cos(phi1), sin(theta1)*sin(phi1), cos(theta1)];
c3 = [sin(theta1)*cos(phi0), sin(theta1)*sin(phi0), cos(theta1)];
a = min([c0*n', c1*n', c2*n', c3*n']);

if lobeIn == 1
    b = 1;
elseif lobeIn == -1
    a = -1;
    b = max([c0*n', c1*n', c2*n', c3*n']);
else
    b = max([c0*n', c1*n', c2*n', c3*n']);
end

end