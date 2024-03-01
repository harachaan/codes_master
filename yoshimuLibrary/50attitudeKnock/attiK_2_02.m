% ----------------------------------------------------------------------
%   2-02 2-01で求めた3-2-1オイラー角をクォータニオンへ変換せよ
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

%% 3-2-1オイラー角の回転行列
% phi(z軸） -> theta(y軸） -> psi(x軸） の順とする．
% ここでは値を適当に設定．
% 資料によってはphiがx軸周りだったり，定義が異なる場合があるので注意．
phi = 0.1;
theta = 0.2;
psi = 0.3;

R = [cos(theta)*cos(phi) cos(theta)*sin(phi) -sin(theta)
     sin(theta)*cos(phi)*sin(psi)-sin(phi)*cos(psi) sin(theta)*sin(phi)*sin(psi)+cos(phi)*cos(psi) cos(theta)*sin(psi)
     sin(theta)*cos(phi)*cos(psi)+sin(phi)*sin(psi) sin(theta)*sin(phi)*cos(psi)-cos(phi)*sin(psi) cos(theta)*cos(psi)];

%% 回転行列から逆算する方法（もあるが面倒）
% 同様にクォータニオンの回転行列が
% R = [q(1)^2 + q(2)^2 - q(3)^2 - q(4)^2, 2*(q(2)*q(3) + q(1)*q(4)), 2*(q(2)*q(4) - q(1)*q(3))
%     2*(q(2) * q(3) - q(1) * q(4)), q(1)^2 - q(2)^2 + q(3)^2 - q(4)^2, 2*(q(3)*q(4) + q(1)*q(2))
%     2*(q(2)*q(4) + q(1)*q(3)), 2*(q(3)*q(4) - q(1)*q(2)), q(1)^2 - q(2)^2 - q(3)^2 + q(4)^2];
% なので，回転行列からクォータニオンを逆算．

%% オイラー角の回転をそれぞれクォータニオンに対応させる方法
% この方法だとオイラー角の時間履歴をfor文を使わずに一気に変換できる．

% q4がスカラー部の場合を想定
n = length(phi); % 3-2-1オイラー角の時間履歴にも対応させる．
q_phi = [zeros(n,2), sin(phi/2), cos(phi/2)]; % z軸まわり，nx4 matrix
q_theta = [zeros(n,1), sin(theta/2), zeros(n,1), cos(theta/2)];
q_psi = [sin(psi/2), zeros(n,2), cos(psi/2)];
q = qMult(4,1, q_psi, qMult(4,1, q_theta, q_phi))