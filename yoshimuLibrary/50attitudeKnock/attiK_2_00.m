% ----------------------------------------------------------------------
%   2-00 1-00で設定した慣性テンソルJの主慣性軸への回転行列に対応する1-2-3オイラー角を求めよ
%    20210510  y.yoshimura
%    Inputs:
%   Outputs:
%   related function files:
%   note:
%   cf:
%   revisions;
%   
%   (c) 2021 yasuhiro yoshimura
%----------------------------------------------------------------------

clc
clear
close all

%% 慣性テンソル
J = [10 2 -3
    2 14 1
    -3 1 12];

%% 主慣性モーメント
[v, d] = eig(J); % 固有ベクトルv, 固有値dを求めるeig関数．

% 固有ベクトルを並べたものは慣性テンソルを次式で対角化する．
% テンソルの座標変換でもあるので，vそのものが回転行列
v^(-1) * J * v

% phi -> theta -> psi の順とする．回転行列は
% R = [cos(theta)*cos(psi) sin(theta)*cos(psi)*sin(phi)+sin(psi)*cos(phi) sin(psi)*sin(phi)-sin(theta)*cos(psi)*cos(phi)
%     -cos(theta)*sin(psi) cos(psi)*cos(phi)-sin(theta)*sin(psi)*sin(phi) sin(theta)*sin(psi)*cos(phi)+cos(psi)*sin(phi)
%     sin(theta) -cos(theta)*sin(phi) cos(theta)*cos(phi)];
% なので，R(3,1), R(3,2), R(3,3),R(1,1),R(2,1)から逆算

R = v;

phi = atan2(-R(3,2), R(3,3)) % rad
theta = atan2(R(3,1), sqrt(R(1,1).^2 + R(2,1).^2))
psi = atan2(-R(2,1), R(1,1))