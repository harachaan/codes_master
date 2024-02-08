% Error quaternionと等価なrotation vectorを回帰
% ~/ku/master_research/codes_master/regression

clc
clear
close all
% -------------------------------------------------------------------------

curdir = pwd;
addpath(strcat(curdir, '/../hara_functions/'));
addpath(strcat(curdir, '/../gpryui/'));

% input values directly ------------------------------------------------------

validScalar = 2; % どのテストデータを使うか

tau = log(1);
sigma = log(0.3);
eta = log(0.1);
params = [tau, sigma, eta]; % hyper parameters

% read train data ---------------------------------------------------------
% [t_, mApp, q1, q2, q3, q4, w1, w2, w3], shape=(tspan, 9)

cd '../data/traindata/';
trainStruct = dir('*.csv'); % 構造体にtraindata fileを格納
fName = trainStruct(1).name; df = readmatrix(fName); % for acquiring the size of a file
trainDF = zeros(size(df, 1), size(df, 2), size(trainStruct, 1)); % 事前割り当て, shape=(tspan, 9, valid fileの数)
for i = 1:1:size(trainStruct, 1)
    fName = trainStruct(i).name;
    df = readmatrix(fName);
    trainDF(:, :, i) = df;
end

% read valid data ----------------------------------------------------------

cd '../validdata'
validStruct = dir('*.csv'); % Structureにvaliddata fileを格納
fName = validStruct(1).name; df = readmatrix(fName); % for acquiring the size of a file
validDF = zeros(size(df, 1), size(df, 2), size(validStruct, 1)); % 事前割り当て, shape=(tspan, 9, valid fileの数)
for i = 1:1:size(validStruct, 1)
    fName = validStruct(i).name;
    df = readmatrix(fName);
    validDF(:, :, i) = df;
end

cd '../../regression';

% data handling -----------------------------------------------------------
% 姿勢の差分はquaternion errorで計算して，それと等価な回転軸ベクトルを出力のデータセットにする

% train data (x_y_train: shape=(:, inputDim+outputDim), input:quaternion and w, output:mApp, delta_rv and delta_w 学習データ全部繋げてる)
inputDim = 7; outputDim = 7;
x_train = zeros(size(trainDF, 1)-1, inputDim, size(trainDF, 3)); 
y_train = zeros(size(trainDF, 1)-1, outputDim, size(trainDF, 3)); % preallocation (出力を差分にすること考慮)
x_y_train = [];
for i = 1:1:size(trainDF, 3) % fileごとに読み込んで操作
    mApp = trainDF(:, 2, i);
    q = trainDF(:, 3:6, i);
    w = trainDF(:, 7:9, i);
    delta_q = zeros(size(q, 1)-1, size(q, 2)); 
    delta_rv = zeros(size(delta_q, 1), size(delta_q, 2)-1);
    delta_w = zeros(size(w, 1)-1, size(w, 2)); % preallocation
    for j = 1:1:(size(q, 1)-1)
        delta_q(j, :) = qMult_h(q(j+1, :), qInv_h(q(j, :))); % error quaternion
        delta_rv(j, :) = q2rv(delta_q(j, :)); % convert delta_q to rotation vector
        delta_w(j, :) = w(j+1, :) - w(j, :);
    end
    x_train(:, :, i) = [trainDF(1:end-1, 3:6, i), w(1:end-1, :)]; % 各ファイルごとの入力と出力を作成
    y_train(:, :, i) = [mApp(1:end-1, :), delta_rv, delta_w];
    x_y_train = [x_y_train; x_train(:, :, i), y_train(:, :, i)]; % 1つの行列にまとめる
    x_y_train(any(isinf(x_y_train)'), :) = []; % infのある行を削除
end
x_train = x_y_train(:, 1:7); y_train = x_y_train(:, 8:14);
% train dataの各スカラー出力の平均を0にする
y_train_mean = mean(y_train); % shape=(1, outputDim)
y_train = y_train - y_train_mean; % broadcasting (y_trainの各行ベクトルについてy_train_meanを引く)

% valid data (t_x_y_valid: (t_, mapp, quaternion, w), shape=(:,9), テストデータは繋げずに変数にそれぞれ格納
t_x_y_valid = validDF(:, :, validScalar);
t_valid = t_x_y_valid(:, 1);
t_x_y_valid(isinf(t_x_y_valid)) = NaN;
t_x_y_valid = fillmissing(t_x_y_valid, 'previous'); % mAppのinfを直前の値に置換
% t_x_y_valid(any(isinf(t_x_y_valid)'), :) = []; % infのある行を削除
t_valid = t_valid(1:size(t_x_y_valid, 1), 1); x_valid = t_x_y_valid(:, 3:9); y_valid = t_x_y_valid(:, 2:end);
% for i = 1:1:size(validDF, 3)
%     assignin('base', sprintf('t_x_y_valid%d', i), validDF(:, :, i)); % 変数't_x_y_valid%d'にvaliddataを格納．(infの数が違うと3次元配列にできなかった)
% end


% 逐次的にregression -------------------------------------------------------

tStart_gpr = tic;
y_reg_q = zeros(size(y_valid)); % shape=(:, outputDim+1), rvをqにするから列を+1
y_reg_mu = zeros(size(y_reg_q, 1), size(y_reg_q, 2)-1);
y_reg_var = zeros(size(y_reg_mu));

y_reg_q(1,2:8) = y_valid(1,2:8); % initial condittion of q and w

[K, Kinv] = kernelMatrix(x_train, params);

for i = 1:1:size(y_reg_q, 1)
    if i < size(y_reg_q, 1)
        for j = 1:1:size(y_reg_mu, 2)
            [y_reg_mu(i, j), y_reg_var(i, j)] = gprNormal(y_reg_q(i, 2:8), x_train, y_train(:, j), params, Kinv); % 入力(:,7): q,w, 出力(:,6): delta_rv, w
        end
        y_reg_mu(i, :) = y_reg_mu(i,:) + y_train_mean; % 出力の平均を0にしたのを戻す
        y_reg_q(i, 1) = y_reg_mu(i, 1); % mAppはそのまま
        y_reg_q(i+1, 2:5) = qMult_h(rv2q(y_reg_mu(i, 2:4)), y_reg_q(i, 2:5)); % 回帰されたdelta_rvをdelta_qに変換してquaternion積で次時刻の姿勢を計算
        y_reg_q(i+1, 6:8) = y_reg_q(i, 6:8) + y_reg_mu(i, 5:7); % 回帰されたdelta_wを時刻iのwに足す．
    elseif i == size(y_reg_q, 1) % 最後の時刻のmAppを回帰する
        y_reg_mu(i, 1) = gprNormal(y_reg_q(i, 2:8), x_train, y_train(:, 1), params, Kinv);
        y_reg_mu(i, 1) = y_reg_mu(i, 1) + y_train_mean(1, 1);
        y_reg_q(i,1) = y_reg_mu(i, 1); % mAppはそのまま, 最後の時刻の回帰mAppを代入
    end
end
% 角速度をdegに変換
y_reg_q(:, 6:8) = rad2deg(y_reg_q(:, 6:8));
t_x_y_valid(:, 7:9) = rad2deg(t_x_y_valid(:, 7:9));
y_valid(:, 6:8) = rad2deg(y_valid(:, 6:8));

tEnd_gpr = toc(tStart_gpr) / 60;

% export to csv ----------------------------------------------------------

savedir = strcat(curdir, '/results/');
savename = strcat(savedir, 't_mApp_q_w_reg.csv'); writematrix([t_valid, y_reg_q], savename);
savename = strcat(savedir, 't_mApp_q_w_valid.csv'); writematrix(t_x_y_valid, savename);


% plot --------------------------------------------------------------------

savedir = strcat(curdir, '/results/figs/');
fignames = ["mApp", "q1", "q2", "q3", "q4", "w1", "w2", "w3"];
x_label = "time [s]";
y_labels = ["magnitude", "q1", "q2", "q3", "q4", "\omega_1 [deg/s]", "\omega_2 [deg/s]", "\omega_3 [deg/s]"];

for i = 1:1:size(fignames, 2)
    figure('Name', fignames(i));
    plot_gpr(t_valid, y_reg_q(:, i), y_valid(:, i), fignames(i), x_label, y_labels(i), savedir)
end

% -------------------------------------------------------------------------








