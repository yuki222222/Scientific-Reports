

% ----------------------------------------------------------------
close all
clear all

%% 2次元の平均輝度値データint_list(line:部位, row:フレーム数)と呼吸波を読み込み，バンドパスフィルタで濾波

fishnum = 'caf100a';
loadname1 = append('D:\harada\研究\matlab\呼吸波心電位解析\短時間計測実験データ\呼吸波\luminance_ROI_', fishnum, '.mat');
% loadname1 = append('D:\harada\研究\matlab\呼吸波心電位解析\短時間計測実験データ\心電位\sig_', fishnum, '.txt.mat');
loadname2 = append('D:\harada\研究\matlab\呼吸波心電位解析\短時間計測実験データ\心電位\フィルタ処理後\wave_ecg_', fishnum, '.mat');
savename = append('D:\harada\研究\matlab\呼吸波心電位解析\短時間計測実験データ\ピーク周波数結果\AIC10_50\freq_', fishnum, '.mat');
load(loadname1)
load(loadname2)

% 保存するかどうか
flag_save = 0;

% サンプリング周波数
Fs = 100;

% カットオフ周波数
brth_l = 0.5;
brth_h = 6;

% 描画フラグ
flag_t = 0;
flag_hl = 0;
flag_ha = 0;
flag_freq = 1;
flag_filt = 0;

%バンドパスフィルタで濾波
wave_resp = BPF_but(int_list(1,1:fr_num-1),Fs, brth_l, brth_h);
wave_resp = wave_resp(1:120*Fs);
% wave_resp = BPF_but(data_sig(1:120*Fs, 2),Fs, brth_l, brth_h);

% データ長
len = length(wave_resp);
% 
% % 心拍と呼吸をプロット
% % time_all = (1:length(wave_resp))/Fs/60;
% figure();
% plot(ave_wave_ecg_high);
% title('心拍');
% figure();
% plot(wave_resp);
% xline(15*60*Fs);
% title('呼吸');
% pbaspect([8 1 1])
% xlim([1 2400000]);
% ylim([-500 500]);
% 
% パワポ用
figure();
plot([wave_resp])
title('拡大');
xlim([50*Fs 60*Fs]);
pbaspect([5.63 1 1]);
ylim([-0.1 0.1]);

figure();
plot(wave_resp);
%% 時間区間win[s]でwave_ecgを分割

% 時間区間
 win = 5;
win_fr = win*Fs;

% 分割数
Mb = fix(len/win_fr);

% wave_brを分割
div_wave_resp = zeros(Mb, win_fr);

cnt = 1;
    for i=1:Mb
        for j=1:win_fr
            div_wave_resp(i,j) = wave_resp(cnt);
        
            if cnt >= len
                break;
            end
            cnt = cnt+1;
        end
    end

%% 各時間区間で振幅を正規分布に標準正規化

 for i=1:Mb
     div_wave_resp(i,:) = zscore(div_wave_resp(i,:));
  
 end
 
 
% 結果を描画
if flag_t == 1
%     Mb = 30;
    for i=1:Mb
    graphname = [num2str((i-1)*win_fr), '～',num2str(i*win_fr),'フレームにおける振幅の標準偏差'];
    figure('Name', graphname) 
    % 呼吸波描画
     plot(div_wave_resp(i,:),'-b');
    
%      pbaspect([9 2 1]) 
%     pbaspect([8 1 1])

  
    
    
     xticklabels({})
     yticklabels({})
%     xlabel('フレーム数')
%     ylabel('振幅の標準偏差')
%     legend('1番目','2番目','3番目','4番目')
    end
end


%% 周波数解析
% ARモデルの係数決定にはユール・ウォーカー法，次数決定にはAICを使用

pAR = zeros(Mb,2); % 各区間のAR次数を格納　1列目：呼吸　2列目：心電図
dft= 512;   % dft点数を設定(偶数にすること！)　

% PSD値と離散周波数を格納
resp_psd = zeros(Mb, dft/2+1);
resp_F = zeros(Mb, dft/2+1);

% ecg_psd = zeros(Mb, dft/2+1);
% ecg_F = zeros(Mb, dft/2+1);

for i=1:Mb
    % 適切なARモデル次数の設定(呼吸)
    for mo = 10:50
        ARmodel = ar(div_wave_resp(i,:), mo, 'yw');
        AIC(mo) = aic(ARmodel); % AICの計算
    end
    [~, morder] = min(AIC);           % AICが最小となる次数を抽出
    pAR(i, 1) = morder;
%     
%     % 適切なARモデル次数の設定(心電図)
%     for mo = 1:20
%         ARmodel = ar(div_wave_ecg(i,:), mo, 'yw');
%         AIC(mo) = aic(ARmodel); % AICの計算
%     end
%     [tmp, morder] = min(AIC);           % AICが最小となる次数を抽出
%     pAR(i, 2) = morder;
    
    % PSD推定
%     [Pr, Fr] = pyulear(div_wave_resp(i,:)-mean(div_wave_resp(i,:)), 10, dft, Fs);
    [Pr, Fr] = pburg(div_wave_resp(i,:)-mean(div_wave_resp(i,:)), pAR(i,1), dft, Fs);
%     [Pe, Fe] = pburg(div_ave_wave_ecg_high(i,:)-mean(div_ave_wave_ecg_high(i,:)), 10, dft, Fs);
    
    % 結果を格納
    resp_psd(i,:) = Pr;
    resp_F(i,:) = Fr;
%     ecg_psd(i,:) = Pe;
%     ecg_F(i,:) = Fe;
     
end

resp_psd_dm = resp_psd;
% ecg_psd_dm = ecg_psd;

%psdの正規化
    for i=1:Mb
        for j=1:length(Pr)
            resp_psd(i,j) = resp_psd_dm(i,j)/sum(resp_psd_dm(i,:));
        end
    end

%      for i=1:Mb
%         for j=1:length(Pr)
%             ecg_psd(i,j) = ecg_psd_dm(i,j)/sum(ecg_psd_dm(i,:));
%         end
%     end

% 結果を描画
if flag_freq == 1
    for i=1:Mb
    graphname = [num2str((i-1)*win_fr), '～',num2str(i*win_fr),'フレームにおけるPSD'];
    figure('Name', graphname)
        plot(resp_F(i, :), resp_psd(i, :), '-b')
%         ylim([0 0.8])
%         ylim([0 0.15])
%          hold on
       
        
%         plot(ecg_F(i, :), ecg_psd(i, :), '-r')
        xlim([0 10])
        ylim([0 0.6])
        % スライド用
        pbaspect([10 1 1])
%          xticklabels({})
          yticklabels({})
        
    end

end

% ピーク周波数計算
pks_fr_resp = zeros(Mb, 1);
% pks_fr_ecg = zeros(Mb, 1);

for i = 1:Mb
    [~, I1] = max(resp_psd(i, :));
%     [~, I2] = max(ecg_psd(i, :));
    pks_fr_resp(i, 1) = resp_F(1, I1);
%     pks_fr_ecg(i, 1) = resp_F(1, I2);
end 

pks_fr_resp = pks_fr_resp.';
% pks_fr_ecg = pks_fr_ecg.';

% zeroindex = find(pks_fr_resp == 0);
% pks_fr_resp(zeroindex) = pks_fr_resp(zeroindex-1);
%% ピーク周波数の時間変化を描画
figure();
plot(pks_fr_resp, '-b');
ylim([0 9]);
% pks_fr_resp(19*Fs:40*Fs) = NaN;
% pks_fr_resp(51*Fs:71*Fs) = NaN;

%% 保存
pks_fr_ecg = pks_fr_ecg(1:23);
if flag_save == 1
%     savefile = append('D:\harada\研究\matlab\呼吸波心電位解析\短時間計測実験データ\ピーク周波数結果\freq_', fishnum, '.mat');
    savevar1 = 'pks_fr_resp';
%     savevar2 = 'pks_fr_ecg';
    save(savename, savevar1);
 
end