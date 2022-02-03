%% q秒ごとの変動係数を計算する
clear all
close all

%% データ読込
fishnum = 'caf100g';
loadname1 = append('D:\harada\研究\matlab\呼吸波心電位解析\短時間計測実験データ\間隔波形データ\因果解析用間隔データ\pks_', fishnum, '.mat');
load(loadname1)
% 0mg/L:a, c, e
% 1mg/L:a, b, c
% 50mg/L:a, b, d
% 100mg/L:a, d, f, i, j

% サンプリング周波数
Fs = 20;

% 20Hzでリサンプリングした間隔波形
PP = PP(:,2)/1000;
RR = RR(:,2)/1000;

% figure();
% plot(PP, '-b');
% hold on;
% plot(RR, '-r');

flag_save = 1;

%% 変動係数計算
% 窓幅(心拍)
qh = 10*Fs;
% 窓幅（呼吸）
qr = 5*Fs;
% 分割数
Mh= fix(length(PP)/qh);
Mr= fix(length(PP)/qr);
% CV計算結果
resp_cv = zeros(Mr,1);
ecg_cv = zeros(Mh,1);
% 平均値
resp_av = mean(PP);
ecg_av = mean(RR);

for i = 1:Mh
    ecg_cv(i) = var(RR((i-1)*qh+1:i*qh))/ecg_av;
end

for i = 1:Mr
    resp_cv(i) = var(PP((i-1)*qr+1:i*qr))/resp_av;
end

%% 変動係数の時系列変化
time_h = (1:Mh)*qh/60;
time_r = (1:Mr)*qr/60;
% figure();
% plot(time_r, resp_cv, '-b');
% hold on;
% plot(time_h, ecg_cv, '-r');
% ylim([0 8]);
% xline(15);

% 移動平均で滑らかに
maf_resp = MAF(resp_cv, 10);
maf_ecg = MAF(ecg_cv, 10);

time_maf_h = (1:Mh-8)*qh/60;
time_maf_r = (1:Mr-8)*qr/60;
% figure();
% plot(time_maf_r, maf_resp, '-b');
% hold on 
% plot(time_maf_h, maf_ecg, '-r');
% xline(15);
% ylim([0 8]);

%% 結果を保存
if flag_save == 1
     savefile = append('D:\harada\研究\matlab\呼吸波心電位解析\短時間計測実験データ\変動係数結果\再実験後データ\cv_', fishnum, '.mat');
    savevar1 = 'resp_cv';
    savevar2 = 'ecg_cv';
    save(savefile, savevar1, savevar2); 
end
