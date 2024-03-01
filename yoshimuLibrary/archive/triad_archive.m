function R = triad(v1_i, v2_i, w1_b, w2_b)
% ---------------------------------------------------------------------
%  calculate Diretional cosine matrix using TRIAD method
%    Input: v1_i, v2_i: reference vector w.r.t inertial frame
%           w1_b, w2_b: referecne vector w.r.t b-frame(e.g., body-fixed frame) 
%   Output: R: directional cosine matrix from inertial frame to b-frame
%   R  = triad(v1_i, v2_i, w1_b, w2_b)
%   (c) 2016 yasuhiro yoshimura
%----------------------------------------------------------------------

% normalization
v1_i = v1_i ./ norm(v1_i);
v2_i = v2_i ./ norm(v2_i);
w1_b = w1_b ./ norm(w1_b);
w2_b = w2_b ./ norm(w2_b);

r1 = v1_i;
r2 = cross(v1_i, v2_i);
r3 = cross(r1, r2);

s1 = w1_b;
s2 = cross(w1_b, w2_b);
s3 = cross(s1, s2);

R = [s1 s2 s3] * [r1 r2 r3]';




end