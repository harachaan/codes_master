% light_curvesの時間履歴を出力したいけど，
% とりあえずパラメータ全部適当に決めて点のスカラー量を出した


% 必要なパラメータ--------------------------
% BRDF f_r
% 表面反射率 s
% 鏡面反射成分 R_s
% 全体 R_spec
% 
%     i番目の面の法線ベクトル u_n_i
%     二等分単位ベクトル u_h
%         太陽光ベクトル u_sun
%         観測者ベクトル u_obs
%     
%     面の反射特性を表す定数 z
%         n_u
%         n_v
%     定数 k_1
%     
%     フレネル係数 F_reflect
%         フレネル反射率 F_o
%     
% 表面拡散率 d
% 拡散反射成分 R_d
% 全体 R_diff
%     拡散反射率 rho
%     定数 k_2
% 
% 
% i番目の面に入射して吸収されなかった太陽光の割合 F_sun_i
%     i番目の面のBRDF f_r_i
%     物体表面に当たる可視光域太陽光の強さ C_sun W/m^2
% i番目の面で反射した太陽光の割合 F_obs_i
%     i番目の面の面積 A_i
%     宇宙機の高度 h_t
% 
% N面からなる宇宙機の等級 m_app
%     太陽光の見かけの等級 -26.7
% 
% 
% 平板モデル　(z_b方向の片面のみ反射するものとする？)
%  形状パラメータ
%     寸法(2辺) a, b
%     質量中心　c_g
%  表面特性パラメータ
%     s, d, F_o, rho, n_u, n_v
%---------------------------------------------
clc
clear
close all


s = 0.7;
d = 1-s;
F_o = 0.3;
rho = 0.3;
surface_reflection = [1000 1000];
    n_u = surface_reflection(1, 1); n_v = surface_reflection(1, 2);

r_earth = 6378.14;
h_t = 600 ; % km

% 物体固定座標系単位ベクトル
u_u = [1
       0
       0];
u_v = [0
       1
       0];
u_n = [0
       0
       1];

u_sun = [0
         0
         1];
u_obs = [0
         1 / sqrt(2)
         1 / sqrt(2)]; % 太陽と観測者同じでいいの？ 
    
shape = [1.0 1.0];
    a = shape(1, 1); b = shape(1, 2);
    
C_sun = 455; % W/m^2
%----------------------------------------------
A = a * b;
u_h = (u_sun + u_obs); u_h = u_h ./ norm(u_h);
% phi_h = atan(u_h(2, 1) / u_h(1, 1));
k_1 = sqrt((n_u+1) * (n_v+1)) / (8*pi);
k_2 = (28*rho / 23*pi) * (1 - s*F_o);
% zは２つの表し方があって，こっちだと分子と分母が０になってNaNになっちゃうからもう一つの方に変更
% しようと思ったけど，そもそも幾何学的に特異点になりそうだから，どっちでも良さそう．
z = (n_u*dot(u_h, u_u)^2 + n_v*dot(u_h, u_v)^2) / (1 - dot(u_h, u_n)^2);
% z = n_u * cos(phi_h)^2 + n_v * sin(phi_h)^2;
F_reflect = s*F_o + (1 - s*F_o)*(1 - dot(u_sun, u_h))^5;

% Ashikhmin-shirley Model
R_s = k_1 * dot(u_h, u_n)^z / (dot(u_sun, u_h) ...
    * max([dot(u_obs, u_n) dot(u_sun, u_n)])) * F_reflect;
R_d = k_2 * (1 - (1 - dot(u_obs, u_n)/2)^5) ...
    * (1 - (1 - dot(u_sun, u_n)/2)^5);
f_r = s*R_s + d*R_d; % f_r(theta, phi)
% F_r = integral2(f_r, )

F_sun = C_sun * f_r * dot(u_sun, u_n);

F_obs = (F_sun * A * dot(u_obs, u_n)) / (h_t^2);

m_app = -26.7 - 2.5 * log10(abs(F_obs / C_sun))


