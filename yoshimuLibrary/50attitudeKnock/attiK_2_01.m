% ----------------------------------------------------------------------
%   2-01 2-00で求めた1-2-3オイラー角を3-2-1オイラー角へ変換せよ
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

%% phi -> theta -> psi の順とする．ここでは適当に値を設定
phi = 0.1;
theta = 0.2;
psi = 0.3;
R = [cos(theta)*cos(psi) sin(theta)*cos(psi)*sin(phi)+sin(psi)*cos(phi) sin(psi)*sin(phi)-sin(theta)*cos(psi)*cos(phi)
    -cos(theta)*sin(psi) cos(psi)*cos(phi)-sin(theta)*sin(psi)*sin(phi) sin(theta)*sin(psi)*cos(phi)+cos(psi)*sin(phi)
    sin(theta) -cos(theta)*sin(phi) cos(theta)*cos(phi)];

%% 3-2-1オイラー角の回転行列
% phi(z軸） -> theta(y軸） -> psi(x軸） の順とする．
% 資料によってはphiがx軸周りだったり，定義が異なる場合があるので注意．
% R321 = [cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta)
%     sin(theta)*cos(phi)*sin(psi)-sin(phi)*cos(psi) sin(theta)*sin(phi)*sin(psi)+cos(phi)*cos(psi) cos(theta)*sin(psi)
%     sin(theta)*cos(phi)*cos(psi)+sin(phi)*sin(psi) sin(theta)*sin(phi)*cos(psi)-cos(phi)*sin(psi) cos(theta)*cos(psi)];
% 回転行列はどのような姿勢表現パラメータを使っても一意なので，1-2-3オイラー角から
% 回転行列を計算し，そのR(3,2), R(3,3),R(1,1),R(1,2)から3-2-1オイラー角を逆算

psi = atan2(R(2,3), R(3,3))
phi = atan2(R(1,2), R(1,1))
theta = atan2(-R(1,3), sqrt(R(2,3).^2 + R(3,3).^2))
