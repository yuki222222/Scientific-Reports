

% ----------------------------------------------------------------
close all
clear all


%% 2次元の平均輝度値データint_list(line:部位, row:フレーム数)と呼吸波を読み込み，バンドパスフィルタで濾波

load('D:\harada\研究\matlab\心電図計測\matファイル\luminance_ROI_caf0za.mat')
load('D:\harada\研究\matlab\心電図計測\心電_呼吸波データ\sig_caf0zz.txt.mat')
Fs_brth = 100;
Fs_ecg = 1000;


% 呼吸：青　，心電図：オレンジ

brth_l = 0.5;
brth_h = 5;
ecg_l = 0.5;
ecg_h = 5;

% brth_l = 0.5;
% brth_h = 20;
% ecg_l = 0.5;
% ecg_h = 20;

% brth_l = 2;
% brth_h = 7;
% ecg_l = 0.5;
% ecg_h = 40;

%--------波形描画について-------------%
% 時系列データを描画
flag_t = 0;
 % 周波数データを描画
 flag_f = 1;

save_flag = 0;
savename = 'm_br.mat';

%part =1:呼吸
part = 1;

%オレンジ：顎，青：心臓，緑：呼吸波
    orange = [1 102/255 0];
    blue = [0 102/255, 204/255];
    color = orange;
    green = [51/255, 153/255, 102/255];

%バンドパスフィルタで濾波
lum_brth = BPF_but(int_list(part,1:fr_num-1),Fs_brth, brth_l, brth_h);
wave_ecg = BPF_but(data_sig(:,2),Fs_ecg, ecg_l, ecg_h);

%心電図を1000Hzから100Hzにサンプリング
[P,Q] = rat(Fs_brth/Fs_ecg);
wave_ecg_all = resample(wave_ecg, P, Q);

%輝度データを心電図の長さに合わせる
if length(lum_brth)>length(wave_ecg_all)
   lum_brth(length(wave_ecg_all)+1:end) = []; 
else
    wave_ecg_all(length(lum_brth)+1:end) = [];
end

 fr_num = length(wave_ecg_all); % データ数
wave_ecg = wave_ecg_all;
 
% データ長を120秒間にカット
wave_ecg(12001:end) = [];
lum_brth(12001:end) = [];
% fr_num = 12000;

figure();
plot(wave_ecg)

figure();
plot(wave_ecg)
xlim([15*Fs_brth 25*Fs_brth])
ylim([-200 200])
pbaspect([24 3 1])

figure();
plot(lum_brth)
xlim([15*Fs_brth 25*Fs_brth])
ylim([-0.1 0.1])
pbaspect([24 3 1])

%% 位相

 xcorr_on = 0;
 % 位相のずれを解消
 if xcorr_on == 1
     %相互相関を計算
     [R, lag] = xcorr(lum_brth(:),wave_ecg(:));
     %相関が一番大きいところがずれを示す
     [M1,I1] = max(R);
     t_lag = lag(I1);
     %t_lag<0であればwave_br(:)の方が進んでいる
     if(t_lag < 0)
         wave_ecg(:) = circshift(wave_ecg(:), t_lag);
         for i=0:abs(t_lag)-1
            wave_ecg(length(wave_ecg)-i) = 0;
         end
     else
         lum_brth(:) = circshift(lum_brth(:), -t_lag);
         for i=0:t_lag-1
            lum_brth(fr_num-i) = 0; 
         end
     end
     
     
    
 end




%% 時間区間win[s]でlum_brthを分割

% 時間区間
 win = 5;
% win = 2.5;
win_fr = win*Fs_brth;

% 分割数
M = fix(fr_num/win_fr);

% lum_bpfを分割
div_lum_brth = zeros(M, win_fr);

 cnt = 1;
    for i=1:M
        for j=1:win_fr
            div_lum_brth(i,j) = lum_brth(cnt);
            if cnt >= fr_num
                break;
            end
            cnt = cnt+1;
        end
    end
    
    
    
    
    %% 時間区間win[s]でwave_ecgを分割

% 時間区間
 win = 5;
% win = 2.5;
win_fr = win*Fs_brth;

% 分割数
Mb = fix(length(wave_ecg)/win_fr);

% wave_brを分割
div_wave_ecg = zeros(Mb, win_fr);

 cnt = 1;
    for i=1:Mb
        for j=1:win_fr
            div_wave_ecg(i,j) = wave_ecg(cnt);
            if cnt >= length(wave_ecg)
                break;
            end
            cnt = cnt+1;
        end
    end

%% 各時間区間で振幅を正規分布に標準正規化

 for i=1:M
    div_lum_brth(i,:) = zscore(div_lum_brth(i,:));
 end

 for i=1:Mb
    div_wave_ecg(i,:) = zscore(div_wave_ecg(i,:));
 end

% 結果を描画
 
    if(part == 1)
        color = blue;
    else
        color = orange;
    end

if(M < Mb)
    min_M = M;
else
    min_M = Mb;
end


if flag_t == 1
for i=1:min_M
    graphname = [num2str((i-1)*win_fr), '～',num2str(i*win_fr),'フレームにおける振幅の標準偏差'];
    figure('Name', graphname)  
    %顎or心臓描画
     plot(div_lum_brth(i,:), 'Color', color)
% %     
     hold on
    
    
    %心電図描画
    plot(div_wave_ecg(i,:))
    pbaspect([3 1 1])
%      pbaspect([9 2 1]) 
%     pbaspect([8 1 1])

%     hold on
    ylim([-3 3])
    
    
     xticklabels({})
     yticklabels({})
%     xlabel('フレーム数')
%     ylabel('振幅の標準偏差')
%     legend('1番目','2番目','3番目','4番目')
end

end

%% -------------------------------------------周波数解析------------------------------------------------------------------
%% 各分割区間のPSDを，ARモデルを用いて推定
% ARモデルの係数決定にはユール・ウォーカー法，次数決定にはAICを使用

pAR = zeros(M,2); % 各区間のAR次数を格納　1列目：呼吸　2列目：心電図
dft= 512*7;   % dft点数を設定(偶数にすること！)　

% PSD値と離散周波数を格納
lum_psd = zeros(M, dft/2+1);
lum_F = zeros(M, dft/2+1);

ecg_psd = zeros(Mb, dft/2+1);
ecg_F = zeros(Mb, dft/2+1);

for i=1:M
    % 適切なARモデル次数の設定(輝度)
    for mo = 10:50
        ARmodel = ar(div_lum_brth(i,:), mo, 'yw');
        AIC(mo) = aic(ARmodel); % AICの計算
    end
    [~, morder] = min(AIC);           % AICが最小となる次数を抽出
    pAR(i, 1) = morder;
    
    % 適切なARモデル次数の設定(心電図)
    for mo = 10:50
        ARmodel = ar(div_wave_ecg(i,:), mo, 'yw');
        AIC(mo) = aic(ARmodel); % AICの計算
    end
    [tmp, morder] = min(AIC);           % AICが最小となる次数を抽出
    pAR(i, 2) = morder;
    
    % PSD推定
    [Pb, Fb] = pburg(div_lum_brth(i,:)-mean(div_lum_brth(i,:)), 30, dft, Fs_brth);
    [Pe, Fe] = pburg(div_wave_ecg(i,:)-mean(div_wave_ecg(i,:)), 30, dft, Fs_brth);
    
    % 結果を格納
    lum_psd(i,:) = Pb;
    lum_F(i,:) = Fb;
    ecg_psd(i,:) = Pe;
    ecg_F(i,:) = Fe;
     
end

lum_psd_dm = lum_psd;
ecg_psd_dm = ecg_psd;

% [Pxx, F] = pburg(div_lum_brth(1,:)-mean(div_lum_brth(1,:)), 33, 512, Fs_brth);
% [Pxxb, Fb] = pburg(div_wave_ecg(1,:)-mean(div_wave_ecg(1,:)), 33, 512, Fs_brth);
% 
% % [Pxx, F] = pburg(div_lum_bpf(1,:,1)-mean(div_lum_bpf(1,:,1)), 50, 64, Fs);
% 
% % PSD値と離散周波数を格納
% lum_psd = zeros(M, length(Pxx));
% lum_F = zeros(M, length(F));
% 
% br_psd = zeros(Mb, length(Pxxb));
% br_F = zeros(Mb, length(Fb));
% 
%     for i=1:M
%         [Pxx, F] = pburg(div_lum_brth(i,:)-mean(div_lum_brth(i,:)), 33, 512, Fs_brth);
%         lum_psd(i,1:length(Pxx)) = Pxx;
%         lum_F(i,1:length(Pxx)) = F;
%     end
% 
%      for i=1:Mb
%         [Pxxb, Fb] = pburg(div_wave_ecg(i,:)-mean(div_wave_ecg(i,:)), 33, 512, Fs_brth);
%         br_psd(i,1:length(Pxxb)) = Pxxb;
%         br_F(i,1:length(Pxxb)) = Fb;
%      end



%psdの正規化
    for i=1:M
        for j=1:length(Pb)
            lum_psd(i,j) = lum_psd_dm(i,j)/sum(lum_psd_dm(i,:));
        end
    end

     for i=1:Mb
        for j=1:length(Pb)
            ecg_psd(i,j) = ecg_psd_dm(i,j)/sum(ecg_psd_dm(i,:));
        end
    end

%     for i=1:M
%         for j=1:length(Pxx)
%             lum_psd(i,j) = lum_psd_dm(i,j)-min(lum_psd_dm(i,:))/(max(lum_psd_dm(i,:))-min(lum_psd_dm(i,:)));
%         end
%     end
% 
%      for i=1:Mb
%         for j=1:length(Pxxb)
%             br_psd(i,j) = br_psd_dm(i,j)-min(br_psd_dm(i,:))/(max(br_psd_dm(i,:))-min(br_psd_dm(i,:)));
%         end
%     end



% 結果を描画
if flag_f == 1
for i=1:M
    graphname = [num2str((i-1)*win_fr), '～',num2str(i*win_fr),'フレームにおけるPSD'];
    figure('Name', graphname)
        plot(lum_F(i, :), lum_psd(i, :), 'Color',color)
%         ylim([0 0.8])
%         ylim([0 0.15])
         hold on
       
        
        plot(ecg_F(i, :), ecg_psd(i, :), 'Color', green)
        hold on
        xlim([0 10])
        ylim([0 0.6])
        pbaspect([9 1 1])
         xticklabels({})
          yticklabels({})
        
%     xlabel('離散周波数')
%     ylabel('パワースペクトル密度')
%     legend('1番目','2番目','3番目','4番目')
    
%     pbaspect([2.5 1 1])
%     pbaspect([4.5 1 1])
end

end
% 平均PSD計算
% mean_Pe = zeros(1,dft/2+1);
% mean_Pb = zeros(1,dft/2+1);
% for i=1:dft/2+1
%     mean_Pe(1,i) = mean(ecg_psd(:,i));
%     mean_Pb(1,i) = mean(lum_psd(:,i));
% end
% 
% % 平均PSD描画
% figure()
% plot(lum_F(1,:), mean_Pb(1,:))
% hold on
% plot(lum_F(1,:), mean_Pe(1,:), 'Color', green)
% title('平均PSD')
% ylim([0 0.15])

%% 2つのROIのピーク周波数比を求める(呼吸/心拍)

fr_rat = zeros(M, 1);
pks_fr_br = zeros(M, 1);
pks_fr_ecg = zeros(M, 1);

for i = 1:M
    [~, I1] = max(lum_psd(i, :));
    [~, I2] = max(ecg_psd(i, :));
    pks_fr_br(i, 1) = lum_F(1, I1);
    pks_fr_ecg(i, 1) = lum_F(1, I2);
    fr_rat(i, 1) = pks_fr_br(i, 1)/pks_fr_ecg(i, 1);
end 

pks_fr_br = pks_fr_br.';
pks_fr_ecg = pks_fr_ecg.';
fr_rat = fr_rat.';
% caf0f_fr_rat = fr_rat;
% save('caf0f_ratio.mat')
% 
%% 各時間帯におけるピーク周波数をプロットした散布図を描画
% m_pk_br_psd = zeros(1,M);
% m_pk_lum_psd = zeros(1,M);
% 
% for i = 1:M
%     [~, I1] = max(ecg_psd(i, :));
%     [~, I2] = max(lum_psd(i, :));
%     m_pk_br_psd(1,i) = ecg_F(i, I1);
%     m_pk_lum_psd(1,i) = lum_F(i, I2);
% end 
% 
% graphname = ['各時間帯におけるピーク周波数'];
% figure('Name', graphname)
% scatter(m_pk_br_psd, m_pk_lum_psd);
% % ylim([2 5])
% % xlim([2 5])
% 
% %% matファイルに出力
% if save_flag == 1
%     save(savename)
% end
% 
% %% 心電と呼吸を重ねて描画
% graphname = ['全時間領域における心電図と呼吸'];
% figure('Name', graphname)
% plot(zscore(lum_brth))
% hold on 
% plot(zscore(wave_ecg))
% 
%  [Pxx_br, F_br] = pburg(lum_brth-mean(lum_brth), 33, 512, Fs_brth);
%  figure()
% plot(F_br, Pxx_br)
%  
% 
% %% 心電のみ描画（スライド用）
% % figure()
% % plot(wave_ecg)
% % xlim([3000 5000])
% % pbaspect([9 2 1])
% % 
% % figure()
% % plot(wave_ecg)
% % xlim([4200 4400])
% % pbaspect([6 2 1])
% 
% figure()
% plot(zscore(wave_ecg))
% hold on
% plot(zscore(lum_brth))
%  ylim([-3 7])
% xlim([1000 2000])

%% 対応有の両側t検定
% % 各フレームにおけるピーク周波数を算出し，全時間領域で比較．
% % 同等性検定を行う．
% 
% cut_flag = 1;
% rd_fr = 0;
% 
% % おかしいフレームを取り除く
% if cut_flag == 1
%     % カットしたい秒
%     %小さい順に記述すること！
%     % (22秒がおかしい)
%     cut_fr = [22];
%     
%     %削除によるずれを埋める変数
%     cut_count = 0;
%     
%     %指定したフレームをカット
%     for i = 1:length(cut_fr)
%         cut_id = fix(cut_fr(i)/win)+1;
%         lum_psd(cut_id-cut_count, :) = [];
%         br_psd(cut_id-cut_count, :) = [];
%         cut_count = cut_count + 1;
%     end 
%     rd_fr = length(cut_fr);
% end
% 
% 
% 
% %1列目に鰓の各フレームにおけるピーク周波数を，2列目に呼吸波の各フレームにおけるピーク周波数を格納
% f_comp = zeros(M-rd_fr, 2);
% 
% for i = 1:M-rd_fr
%     [~, I] = max(lum_psd(i, :));
%     f_comp(i, 1) = lum_F(i, I);
%     [~, I] = max(br_psd(i, :));
%     f_comp(i, 2) = lum_F(i, I);   
% end
% 
% [p1, p2, CI] = TOST(f_comp(:,1), f_comp(:,2), -0.5, 0.5, 0.05);
% 
% % % 各離散周波数におけるPSDの平均と分散を各離散周波数ごとに求める
% % F_max = 10;
% % cut_flag = 1;
% % rd_fr = 0;
% % 
% % if cut_flag == 1
% %     % カットしたい秒
% %     %小さい順に記述すること！
% %     % (22秒，52秒がおかしい)
% %     cut_fr = [22 52];
% %     
% %     %削除によるずれを埋める変数
% %     cut_count = 0;
% %     
% %     %指定したフレームをカット
% %     for i = 1:length(cut_fr)
% %         cut_id = fix(cut_fr(i)/win)+1;
% %         lum_psd(cut_id-cut_count, :) = [];
% %         br_psd(cut_id-cut_count, :) = [];
% %         cut_count = cut_count + 1;
% %     end 
% %     rd_fr = length(cut_fr);
% % end
% 
% % 各離散周波数のPSDを格納
% f_comp2 = zeros(M-rd_fr, F_max);
% 
% %1[Hz]～F_max[Hz]まで
% for F_unit=1:F_max
%     %整数F_unitの近傍にある離散周波数Fを取り出す
%     unit = knnsearch(F, F_unit);
%     %
%     f_comp2(:, F_unit) = lum_psd(:, unit);
%     f_comp2(:, F_unit, 2) = br_psd(:, unit);
% end
% 
% % エクセルに出力
% filename = '各離散周波数における鰓ROIと呼吸波のpsd.xlsx';
% writetable(array2table(f_comp2(:,:,1)),filename,'Sheet',1,'Range','A1:Z100')
% writetable(array2table(f_comp2(:,:,2)),filename,'Sheet',2,'Range','A1:Z100')
% 
% %各離散周波数でのpsdの両側t検定
% h = zeros(1, F_max);
% pt = zeros(1, F_max);
% 
% for i=1:F_max
%    [h(i),pt(i)] = ttest(f_comp2(:,i,1), f_comp2(:,i,2)); 
% end
% 
% f_num = length(lum_psd);    % 離散周波数の数
% 
% % 顎ROIと呼吸波における各離散周波数のpsd平均と分散
% psd_mean = zeros(f_num, 1);
% psd_var = zeros(f_num, 1);
% psd_mean_br = zeros(f_num, 1);
% psd_var_br = zeros(f_num, 1);
% 
% %鰓ROIの各離散周波数におけるpsdの平均と分散を計算
% for i=1:f_num
%         psd_mean(i,1) = mean(lum_psd(:,i));
%         psd_var(i,1) = var(lum_psd(:,i));
% end
% 
% %呼吸波の各離散周波数におけるpsdの平均と分散を計算
% for i=1:f_num
%         psd_mean_br(i,1) = mean(br_psd(:,i));
%         psd_var_br(i,1) = var(br_psd(:,i));
% end
% 
% 
% % [h_mean,p_mean] = ttest(psd_mean(:,1), psd_mean_br(:,1));
% % [h_var,p_var] = ttest(psd_var(:,1), psd_var_br(:,1));
% 
% % 平均と分散の集団をエクセルに出力
% filename = '平均と分散の計算結果(鰓ROIと呼吸波).xlsx';
% writetable(array2table(psd_mean(:,1)),filename,'Sheet',1,'Range','A1')
% writetable(array2table(psd_mean_br(:,1)),filename,'Sheet',1,'Range','B1')
% writetable(array2table(psd_var(:,1)),filename,'Sheet',2,'Range','A1')
% writetable(array2table(psd_var_br(:,1)),filename,'Sheet',2,'Range','B1')




