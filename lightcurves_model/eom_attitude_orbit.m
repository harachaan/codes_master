% 姿勢，角速度，軌道位置・速度を求める常微分方程式のfunction
function dydt = eom_attitude_orbit(t, y, J, mu)

% トルク [N・m]
tau = [0
       0
       0];
w1 = y(1,1); w2 = y(2,1); w3 = y(3,1);
w = [w1
     w2
     w3];
q1 = y(4,1); q2 = y(5,1); q3 = y(6,1); q4 = y(7,1);
q = [q1
     q2
     q3
     q4];
r1 = y(8,1); r2 = y(9,1); r3 = y(10,1);
r = [r1
     r2
     r3];
v1 = y(11,1); v2 = y(12,1); v3 = y(13,1);
v = [v1
     v2
     v3];

% 角運動量(機体固定座標系)
h_b = J * w;

% euler's equation 
% dydt(1:3, 1) = inv(J) * (tau - cross(w, hb)); % 下の方が計算の効率がいいらしい？
dydt(1:3, 1) = J \ (tau - cross(w, h_b)); 
% quaternions' kinematics
w_q = [w
       0]; % 便宜上，quaternionsに合わせて4行ベクトルにした
dydt(4:7, 1) = 1/2 * q_pro(w_q, q);

% Cowell's Formulation
FI = 0; FJ = 0; FK = 0;
F = [FI
     FJ
     FK];
dydt(8:10) = v;
dydt(11:13) = -mu / norm(r)^3 * r + F;

end