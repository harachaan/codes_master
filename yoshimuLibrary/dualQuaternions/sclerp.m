function dqt = sclerp(t, scalar, dq1, dq2)
% ----------------------------------------------------------------------
%   SCrew Linear Interpolation of quaternions (slerp)
%    20210223  y.yoshimura
%    Inputs: t: normalzied time, i.e., 0<=t<=1
%            scalar: specifies the definition of the quaternion
%            scalar == 0,  q0:= cos(theta/2), q = [q0, q1, q2, q3]
%            scalar == 4,  q4:= cos(theta/2), q = [q1, q2, q3, q4]         
%            dq1: starting dual quaternion, 1x8
%            dq2: ending dual quaternion, 1x8
%   Outputs: qt: interpolated quaternions, Nx8
%   related function files:
%   note:
%   cf:
%   revisions;
%   function qt = sclerp(t, scalar, dq1, dq2)
%   (c) 2021 yasuhiro yoshimura
%----------------------------------------------------------------------

t = t(:); % column vector

% rotation
dqtmp = dqMult(scalar, 0, dqInv(scalar,dq1), dq2);

eTheta = (scalar == 0) .* acos(dqtmp(1)) ...
    + (scalar == 4) .* acos(dqtmp(4));  % Euler angle of rotation
eTheta = eTheta * 2;
eAxis = (scalar == 0) .* dqtmp(2:4) ...  % Euler axis of rotation
    + (scalar == 4) .* dqtmp(1:3);
eAxis = eAxis ./ sin(eTheta/2);

% translation
d = (scalar == 0) .* dqtmp(5) ...
    + (scalar == 4) .* dqtmp(8);
d = -2 .* d ./ sin(eTheta/2);

dqVec = (scalar == 0) .* dqtmp(6:8) ...
    + (scalar == 4) .* dqtmp(5:7);
p = (dqVec - d ./ 2 .* cos(eTheta/2) .* eAxis) ./ sin(eTheta/2);

dq1tmp = kron(ones(size(t,1),1), dq1);

dqt = (scalar == 0) .* dqMult(scalar, 0, dq1tmp, [cos(t.*eTheta/2), eAxis.*sin(t.*eTheta./2), -t.*d./2.*sin(t.*eTheta./2), p.*sin(t.*eTheta./2)+t.*d./2.*cos(t.*eTheta./2).*eAxis]) ...
    + (scalar == 4) .* dqMult(scalar, 0, dq1tmp, [eAxis.*sin(t.*eTheta./2), cos(t.*eTheta/2), p.*sin(t.*eTheta./2)+t.*d./2.*cos(t.*eTheta./2).*eAxis, -t.*d./2.*sin(t.*eTheta./2)]);

end