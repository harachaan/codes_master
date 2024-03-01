function appMag = mag(fobs, d)
% ----------------------------------------------------------------------
%   apparent magnitude of object with distance d [m]
%    20201112  y.yoshimura
%    Inputs: fobs: c_total * (s'*n) * A *  (v'*n), Nx1 vector
%            d: distance between object and observer, m, Nx1 vector
%   Outputs: appMag: apparent magnitude of object
%   related function files:
%   note:
%   cf: 'Light Curve Approximation Using an Attitude Model 
%        of Solar Sail Spacecraft' for the definitoin of fobs      
%   revisions;
%   
%   (c) 2020 yasuhiro yoshimura
%----------------------------------------------------------------------


msun = -26.7; % apparent magnitude of Sun
appMag = msun - 2.5 .* log10(fobs ./ d.^2);


end