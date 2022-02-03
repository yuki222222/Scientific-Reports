%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% 読み込んだ呼吸のピーク抽出を行い，20Hzでリサンプリングしたピーク間隔波形を出力する
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clear all
close all


%% データ読込，フィルタ処理
% 0mg/L:a, c, e
% 1mg/L:a, b, c
% 50mg/L:a, b, d
% 100mg/L:a, d, f, i, j
fishnum = 'caf100g';
% loadname1 = append('D:\harada\研究\matlab\呼吸波心電位解析\短時間計測実験データ\呼吸波\luminance_ROI_', fishnum, '.txt.mat');
loadname1 = append('D:\harada\研究\matlab\呼吸波心電位解析\短時間計測実験データ\心電位\sig_', fishnum, '.txt.mat');
loadname2 = append('D:\harada\研究\matlab\呼吸波心電位解析\短時間計測実験データ\心電位\フィルタ処理後\wave_ecg_', fishnum, '.mat');
load(loadname1)
load(loadname2)


% wave_ecg = data_sig(:,1);
locs_ecg = locs;

% 出力ファイル名
name_fish = append('D:\harada\研究\matlab\呼吸波心電位解析\短時間計測実験データ\間隔波形データ\因果解析用間隔データ\pks_', fishnum, '.mat');

% サンプリング周波数
% Fs_resp = 100;
Fs_heart = 1000;
Fs_resp = 1000;

% カットオフ周波数
resp_fl = 1;
resp_fh = 6;

% wave_resp = BPF_but(int_list(1,1:120*Fs_resp),Fs_resp, resp_fl, resp_fh);
wave_resp = BPF_but(data_sig(1:120*Fs_resp, 2),Fs_resp, resp_fl, resp_fh);

%% ピーク抽出(呼吸)
% ピーク抽出条件
mindistance = 100; % ピーク間の最低距離
pheight = 0;   % ピークの最低値
minpro = 0; 
click_resp = 0;

% ピーク抽出
[~,locs_resp,~,~] = findpeaks(wave_resp(1500:end), 'MinPeakDistance', mindistance,  'MinPeakProminence', minpro, 'MinPeakHeight', pheight);
findpeaks(wave_resp(1500:end), 'MinPeakDistance', mindistance,  'MinPeakProminence', minpro, 'MinPeakHeight', pheight)

% locs_resp = locs_resp';

%%  PP間隔計算
RawPP = [locs_resp(2:end)/Fs_resp, diff(locs_resp)]; % 1列目：時刻(「秒」に変換)，2列目：RR間隔
figure();
plot(RawPP(:,1), RawPP(:,2)); % RRデータの描画
xlabel('Time[s]'); ylabel('PP interval[ms]');
title('PP interval');
% ylim([0 200])

%%  RR間隔計算
RawRR = [locs_ecg(2:end)/Fs_heart, diff(locs_ecg)]; % 1列目：時刻(「秒」に変換)，2列目：RR間隔
figure();
plot(RawRR(:,1), RawRR(:,2)); % RRデータの描画
xlabel('Time[s]'); ylabel('RR interval[ms]');
title('RR interval');
ylim([0 200])
%% PP間隔データをリサンプリング
newFs = 20; % リサンプリングの時間間隔
newtind = RawPP(1, 1):1/newFs:RawPP(end,1); % 等間隔の時刻データを定義 1/newFs刻みにしている
PP = [newtind; interp1(RawPP(:, 1), RawPP(:, 2), newtind)]';   % リサンプリング後データ　列に転置している
figure();
plot(RawPP(:,1), RawPP(:, 2));  % リサンプリング前のRR間隔
hold on;
plot(PP(:,1), PP(:,2), 'r');    % リサンプリング後のRR間隔
xlabel('Time [s]');ylabel('PP interval [ms]');
title('PP interval');

% 外れ値を検出し，平均値に置き換える
% [B,TF] = rmoutliers(PP(:,2),'mean');
% PP(TF,2) = mean(PP(:,2)); 

figure();
plot(PP(:,2))

%% RR間隔データをリサンプリング
newFs = 20; % リサンプリングの時間間隔
newtind = RawRR(1, 1):1/newFs:RawRR(end,1); % 等間隔の時刻データを定義 1/newFs刻みにしている
RR = [newtind; interp1(RawRR(:, 1), RawRR(:, 2), newtind)]';   % リサンプリング後データ　列に転置している
figure();
plot(RawRR(:,1), RawRR(:, 2));  % リサンプリング前のRR間隔
hold on;
plot(RR(:,1), RR(:,2), 'r');    % リサンプリング後のRR間隔
xlabel('Time [s]');ylabel('RR interval [ms]');
title('RR interval');

% 外れ値を検出し，平均値に置き換える
% [B,TF] = rmoutliers(RR(:,2),'mean');
% RR(TF,2) = mean(RR(:,2)); 

% 120sに揃えるため，それ以降を削除
idx = knnsearch(RR(:,1), 120);
RR(idx:end, :) = [];

figure();
plot(RR(:,2))


%% matファイルに出力
save(name_fish)