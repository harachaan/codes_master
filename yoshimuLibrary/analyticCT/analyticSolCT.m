function [Axyz, Bxyz] = analyticSolCT(thetaN, alp, bet, phi, theta)
% ----------------------------------------------------------------------
%    analytic solution for 1st order approximation
%    20210131  y.yoshimura
%    Inputs:
%   Outputs: Axyz: 1x3
%            Bxyz: 1x3
%   related function files:
%   note:
%   cf:
%   revisions;
%   
%   (c) 2021 yasuhiro yoshimura
%----------------------------------------------------------------------

% integration range
phi0 = phi(1);
phi1 = phi(2);

theta0 = theta(1);
theta1 = theta(2);
 

Axyz = [((sin(phi0) - sin(phi1))*(4*theta0 - 4*theta1 - sin(4*theta0) + sin(4*theta1)))/4, -((cos(phi0) - cos(phi1))*(4*theta0 - 4*theta1 - sin(4*theta0) + sin(4*theta1)))/4, -((cos(4*theta0) - cos(4*theta1))*(phi0 - phi1))/4];
Axyz = alp .* Axyz;

Bxyz = [(64*cos(theta0/2)^3*sin(theta0/2)*sin(phi1))/3 - (64*cos(theta0/2)^3*sin(theta0/2)*sin(phi0))/3 + (64*cos(theta1/2)^3*sin(theta1/2)*sin(phi0))/3 + (1088*cos(theta0/2)^5*sin(theta0/2)*sin(phi0))/15 - (64*cos(theta1/2)^3*sin(theta1/2)*sin(phi1))/3 - (1088*cos(theta0/2)^5*sin(theta0/2)*sin(phi1))/15 - (1088*cos(theta1/2)^5*sin(theta1/2)*sin(phi0))/15 - (512*cos(theta0/2)^7*sin(theta0/2)*sin(phi0))/5 + (1088*cos(theta1/2)^5*sin(theta1/2)*sin(phi1))/15 + (512*cos(theta0/2)^7*sin(theta0/2)*sin(phi1))/5 + (512*cos(theta1/2)^7*sin(theta1/2)*sin(phi0))/5 + (256*cos(theta0/2)^9*sin(theta0/2)*sin(phi0))/5 - (512*cos(theta1/2)^7*sin(theta1/2)*sin(phi1))/5 - (256*cos(theta0/2)^9*sin(theta0/2)*sin(phi1))/5 - (256*cos(theta1/2)^9*sin(theta1/2)*sin(phi0))/5 + (256*cos(theta1/2)^9*sin(theta1/2)*sin(phi1))/5 + 32*cos(theta0/2)^4*cos(thetaN/2)*sin(thetaN/2)*cos(phi0)^2 - 32*cos(theta0/2)^4*cos(thetaN/2)*sin(thetaN/2)*cos(phi1)^2 - 32*cos(theta1/2)^4*cos(thetaN/2)*sin(thetaN/2)*cos(phi0)^2 + 32*cos(theta1/2)^4*cos(thetaN/2)*sin(thetaN/2)*cos(phi1)^2 - (320*cos(theta0/2)^6*cos(thetaN/2)*sin(thetaN/2)*cos(phi0)^2)/3 + (320*cos(theta0/2)^6*cos(thetaN/2)*sin(thetaN/2)*cos(phi1)^2)/3 + (320*cos(theta1/2)^6*cos(thetaN/2)*sin(thetaN/2)*cos(phi0)^2)/3 - (320*cos(theta1/2)^6*cos(thetaN/2)*sin(thetaN/2)*cos(phi1)^2)/3 + 128*cos(theta0/2)^8*cos(thetaN/2)*sin(thetaN/2)*cos(phi0)^2 - 128*cos(theta0/2)^8*cos(thetaN/2)*sin(thetaN/2)*cos(phi1)^2 - 128*cos(theta1/2)^8*cos(thetaN/2)*sin(thetaN/2)*cos(phi0)^2 + 128*cos(theta1/2)^8*cos(thetaN/2)*sin(thetaN/2)*cos(phi1)^2 - (256*cos(theta0/2)^10*cos(thetaN/2)*sin(thetaN/2)*cos(phi0)^2)/5 + (256*cos(theta0/2)^10*cos(thetaN/2)*sin(thetaN/2)*cos(phi1)^2)/5 + (256*cos(theta1/2)^10*cos(thetaN/2)*sin(thetaN/2)*cos(phi0)^2)/5 - (256*cos(theta1/2)^10*cos(thetaN/2)*sin(thetaN/2)*cos(phi1)^2)/5 + (128*cos(theta0/2)^3*cos(thetaN/2)^2*sin(theta0/2)*sin(phi0))/3 - (128*cos(theta0/2)^3*cos(thetaN/2)^2*sin(theta0/2)*sin(phi1))/3 - (128*cos(theta1/2)^3*cos(thetaN/2)^2*sin(theta1/2)*sin(phi0))/3 - (2176*cos(theta0/2)^5*cos(thetaN/2)^2*sin(theta0/2)*sin(phi0))/15 + (128*cos(theta1/2)^3*cos(thetaN/2)^2*sin(theta1/2)*sin(phi1))/3 + (2176*cos(theta0/2)^5*cos(thetaN/2)^2*sin(theta0/2)*sin(phi1))/15 + (2176*cos(theta1/2)^5*cos(thetaN/2)^2*sin(theta1/2)*sin(phi0))/15 + (1024*cos(theta0/2)^7*cos(thetaN/2)^2*sin(theta0/2)*sin(phi0))/5 - (2176*cos(theta1/2)^5*cos(thetaN/2)^2*sin(theta1/2)*sin(phi1))/15 - (1024*cos(theta0/2)^7*cos(thetaN/2)^2*sin(theta0/2)*sin(phi1))/5 - (1024*cos(theta1/2)^7*cos(thetaN/2)^2*sin(theta1/2)*sin(phi0))/5 - (512*cos(theta0/2)^9*cos(thetaN/2)^2*sin(theta0/2)*sin(phi0))/5 + (1024*cos(theta1/2)^7*cos(thetaN/2)^2*sin(theta1/2)*sin(phi1))/5 + (512*cos(theta0/2)^9*cos(thetaN/2)^2*sin(theta0/2)*sin(phi1))/5 + (512*cos(theta1/2)^9*cos(thetaN/2)^2*sin(theta1/2)*sin(phi0))/5 - (512*cos(theta1/2)^9*cos(thetaN/2)^2*sin(theta1/2)*sin(phi1))/5, (sin(2*phi1)*sin(thetaN)*((4*cos(theta0)^3*(3*cos(theta0)^2 - 5))/15 - (4*cos(theta1)^3*(3*cos(theta1)^2 - 5))/15))/2 - (sin(2*phi0)*sin(thetaN)*((4*cos(theta0)^3*(3*cos(theta0)^2 - 5))/15 - (4*cos(theta1)^3*(3*cos(theta1)^2 - 5))/15))/2 + phi0*sin(thetaN)*((4*cos(theta0)^3*(3*cos(theta0)^2 - 5))/15 - (4*cos(theta1)^3*(3*cos(theta1)^2 - 5))/15) - phi1*sin(thetaN)*((4*cos(theta0)^3*(3*cos(theta0)^2 - 5))/15 - (4*cos(theta1)^3*(3*cos(theta1)^2 - 5))/15) + 2*cos(phi0)*cos(thetaN)*((4*sin(theta0)^3*(3*sin(theta0)^2 - 5))/15 - (4*sin(theta1)^3*(3*sin(theta1)^2 - 5))/15) - 2*cos(phi1)*cos(thetaN)*((4*sin(theta0)^3*(3*sin(theta0)^2 - 5))/15 - (4*sin(theta1)^3*(3*sin(theta1)^2 - 5))/15), 4*phi0*cos(thetaN)*(cos(theta0)^3/3 - cos(theta1)^3/3) - 4*phi1*cos(thetaN)*(cos(theta0)^3/3 - cos(theta1)^3/3) - 8*phi0*cos(thetaN)*(cos(theta0)^5/5 - cos(theta1)^5/5) + 8*phi1*cos(thetaN)*(cos(theta0)^5/5 - cos(theta1)^5/5) + 4*cos(phi0)*sin(thetaN)*(sin(theta0)^3/3 - sin(theta1)^3/3) - 4*cos(phi1)*sin(thetaN)*(sin(theta0)^3/3 - sin(theta1)^3/3) + 8*cos(phi0)*sin(thetaN)*((sin(theta0)^3*(3*sin(theta0)^2 - 5))/15 - (sin(theta1)^3*(3*sin(theta1)^2 - 5))/15) - 8*cos(phi1)*sin(thetaN)*((sin(theta0)^3*(3*sin(theta0)^2 - 5))/15 - (sin(theta1)^3*(3*sin(theta1)^2 - 5))/15)];
Bxyz = bet .* Bxyz; 

end