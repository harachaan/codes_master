function [dqa, dqb] = ctrlDq(dtp, dq1, dq2, w1, v1, w2, v2)
% ----------------------------------------------------------------------
%   calculate control dual quaternions for C1 interpolation
%    20210310  y.yoshimura
%    Inputs: dtp: time step, s
%            dq1, dq2: dual quaternions intepolated from dq1 to dq2 (for
%            each)
%            w1, w2: corresponding angular rates
%            v1, v2: corresponding translational velocities
%   Outputs: dqa, dqb: control dual quaternions, 1x8 (for each)
%   related function files:
%   note:
%   cf:
%   revisions;
%
%   (c) 2021 yasuhiro yoshimura
%----------------------------------------------------------------------

% for dqa
tmp = [dtp/6.*w1, dtp/6.*v1];
eTheta = 2 .* norm(tmp(1:3));
eAxis = tmp(1:3) ./ (eTheta/2);
dqtmp_r = [eAxis.*sin(eTheta/2) cos(eTheta/2)];
dqtmp_d = 0.5 .* qMult(4, 0, [tmp(4:6) 0], dqtmp_r);
dq_tmp = [dqtmp_r dqtmp_d];
dqa = dqMult(4, 0, dq1, dq_tmp);

% for dqb
tmp = [-dtp/6.*w2, -dtp/6.*v2];
eTheta = 2 .* norm(tmp(1:3));
eAxis = tmp(1:3) ./ (eTheta/2);
dqtmp_r = [eAxis.*sin(eTheta/2) cos(eTheta/2)];
dqtmp_d = 0.5 .* qMult(4, 0, [tmp(4:6) 0], dqtmp_r);
dq_tmp = [dqtmp_r dqtmp_d];
dqb = dqMult(4, 0, dq2, dq_tmp);


end