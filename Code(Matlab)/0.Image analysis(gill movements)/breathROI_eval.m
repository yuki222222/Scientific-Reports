

% ----------------------------------------------------------------
close all
clear all


%% 2�����̕��ϋP�x�l�f�[�^int_list(line:����, row:�t���[����)�ƌċz�g��ǂݍ��݁C�o���h�p�X�t�B���^���h�g

load('D:\harada\����\matlab\�S�d�}�v��\mat�t�@�C��\luminance_ROI_caf0za.mat')
load('D:\harada\����\matlab\�S�d�}�v��\�S�d_�ċz�g�f�[�^\sig_caf0zz.txt.mat')
Fs_brth = 100;
Fs_ecg = 1000;


% �ċz�F�@�C�S�d�}�F�I�����W

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

%--------�g�`�`��ɂ���-------------%
% ���n��f�[�^��`��
flag_t = 0;
 % ���g���f�[�^��`��
 flag_f = 1;

save_flag = 0;
savename = 'm_br.mat';

%part =1:�ċz
part = 1;

%�I�����W�F�{�C�F�S���C�΁F�ċz�g
    orange = [1 102/255 0];
    blue = [0 102/255, 204/255];
    color = orange;
    green = [51/255, 153/255, 102/255];

%�o���h�p�X�t�B���^���h�g
lum_brth = BPF_but(int_list(part,1:fr_num-1),Fs_brth, brth_l, brth_h);
wave_ecg = BPF_but(data_sig(:,2),Fs_ecg, ecg_l, ecg_h);

%�S�d�}��1000Hz����100Hz�ɃT���v�����O
[P,Q] = rat(Fs_brth/Fs_ecg);
wave_ecg_all = resample(wave_ecg, P, Q);

%�P�x�f�[�^��S�d�}�̒����ɍ��킹��
if length(lum_brth)>length(wave_ecg_all)
   lum_brth(length(wave_ecg_all)+1:end) = []; 
else
    wave_ecg_all(length(lum_brth)+1:end) = [];
end

 fr_num = length(wave_ecg_all); % �f�[�^��
wave_ecg = wave_ecg_all;
 
% �f�[�^����120�b�ԂɃJ�b�g
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

%% �ʑ�

 xcorr_on = 0;
 % �ʑ��̂��������
 if xcorr_on == 1
     %���ݑ��ւ��v�Z
     [R, lag] = xcorr(lum_brth(:),wave_ecg(:));
     %���ւ���ԑ傫���Ƃ��낪���������
     [M1,I1] = max(R);
     t_lag = lag(I1);
     %t_lag<0�ł����wave_br(:)�̕����i��ł���
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




%% ���ԋ��win[s]��lum_brth�𕪊�

% ���ԋ��
 win = 5;
% win = 2.5;
win_fr = win*Fs_brth;

% ������
M = fix(fr_num/win_fr);

% lum_bpf�𕪊�
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
    
    
    
    
    %% ���ԋ��win[s]��wave_ecg�𕪊�

% ���ԋ��
 win = 5;
% win = 2.5;
win_fr = win*Fs_brth;

% ������
Mb = fix(length(wave_ecg)/win_fr);

% wave_br�𕪊�
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

%% �e���ԋ�ԂŐU���𐳋K���z�ɕW�����K��

 for i=1:M
    div_lum_brth(i,:) = zscore(div_lum_brth(i,:));
 end

 for i=1:Mb
    div_wave_ecg(i,:) = zscore(div_wave_ecg(i,:));
 end

% ���ʂ�`��
 
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
    graphname = [num2str((i-1)*win_fr), '�`',num2str(i*win_fr),'�t���[���ɂ�����U���̕W���΍�'];
    figure('Name', graphname)  
    %�{or�S���`��
     plot(div_lum_brth(i,:), 'Color', color)
% %     
     hold on
    
    
    %�S�d�}�`��
    plot(div_wave_ecg(i,:))
    pbaspect([3 1 1])
%      pbaspect([9 2 1]) 
%     pbaspect([8 1 1])

%     hold on
    ylim([-3 3])
    
    
     xticklabels({})
     yticklabels({})
%     xlabel('�t���[����')
%     ylabel('�U���̕W���΍�')
%     legend('1�Ԗ�','2�Ԗ�','3�Ԗ�','4�Ԗ�')
end

end

%% -------------------------------------------���g�����------------------------------------------------------------------
%% �e������Ԃ�PSD���CAR���f����p���Đ���
% AR���f���̌W������ɂ̓��[���E�E�H�[�J�[�@�C��������ɂ�AIC���g�p

pAR = zeros(M,2); % �e��Ԃ�AR�������i�[�@1��ځF�ċz�@2��ځF�S�d�}
dft= 512*7;   % dft�_����ݒ�(�����ɂ��邱�ƁI)�@

% PSD�l�Ɨ��U���g�����i�[
lum_psd = zeros(M, dft/2+1);
lum_F = zeros(M, dft/2+1);

ecg_psd = zeros(Mb, dft/2+1);
ecg_F = zeros(Mb, dft/2+1);

for i=1:M
    % �K�؂�AR���f�������̐ݒ�(�P�x)
    for mo = 10:50
        ARmodel = ar(div_lum_brth(i,:), mo, 'yw');
        AIC(mo) = aic(ARmodel); % AIC�̌v�Z
    end
    [~, morder] = min(AIC);           % AIC���ŏ��ƂȂ鎟���𒊏o
    pAR(i, 1) = morder;
    
    % �K�؂�AR���f�������̐ݒ�(�S�d�})
    for mo = 10:50
        ARmodel = ar(div_wave_ecg(i,:), mo, 'yw');
        AIC(mo) = aic(ARmodel); % AIC�̌v�Z
    end
    [tmp, morder] = min(AIC);           % AIC���ŏ��ƂȂ鎟���𒊏o
    pAR(i, 2) = morder;
    
    % PSD����
    [Pb, Fb] = pburg(div_lum_brth(i,:)-mean(div_lum_brth(i,:)), 30, dft, Fs_brth);
    [Pe, Fe] = pburg(div_wave_ecg(i,:)-mean(div_wave_ecg(i,:)), 30, dft, Fs_brth);
    
    % ���ʂ��i�[
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
% % PSD�l�Ɨ��U���g�����i�[
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



%psd�̐��K��
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



% ���ʂ�`��
if flag_f == 1
for i=1:M
    graphname = [num2str((i-1)*win_fr), '�`',num2str(i*win_fr),'�t���[���ɂ�����PSD'];
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
        
%     xlabel('���U���g��')
%     ylabel('�p���[�X�y�N�g�����x')
%     legend('1�Ԗ�','2�Ԗ�','3�Ԗ�','4�Ԗ�')
    
%     pbaspect([2.5 1 1])
%     pbaspect([4.5 1 1])
end

end
% ����PSD�v�Z
% mean_Pe = zeros(1,dft/2+1);
% mean_Pb = zeros(1,dft/2+1);
% for i=1:dft/2+1
%     mean_Pe(1,i) = mean(ecg_psd(:,i));
%     mean_Pb(1,i) = mean(lum_psd(:,i));
% end
% 
% % ����PSD�`��
% figure()
% plot(lum_F(1,:), mean_Pb(1,:))
% hold on
% plot(lum_F(1,:), mean_Pe(1,:), 'Color', green)
% title('����PSD')
% ylim([0 0.15])

%% 2��ROI�̃s�[�N���g��������߂�(�ċz/�S��)

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
%% �e���ԑтɂ�����s�[�N���g�����v���b�g�����U�z�}��`��
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
% graphname = ['�e���ԑтɂ�����s�[�N���g��'];
% figure('Name', graphname)
% scatter(m_pk_br_psd, m_pk_lum_psd);
% % ylim([2 5])
% % xlim([2 5])
% 
% %% mat�t�@�C���ɏo��
% if save_flag == 1
%     save(savename)
% end
% 
% %% �S�d�ƌċz���d�˂ĕ`��
% graphname = ['�S���ԗ̈�ɂ�����S�d�}�ƌċz'];
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
% %% �S�d�̂ݕ`��i�X���C�h�p�j
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

%% �Ή��L�̗���t����
% % �e�t���[���ɂ�����s�[�N���g�����Z�o���C�S���ԗ̈�Ŕ�r�D
% % ������������s���D
% 
% cut_flag = 1;
% rd_fr = 0;
% 
% % ���������t���[������菜��
% if cut_flag == 1
%     % �J�b�g�������b
%     %���������ɋL�q���邱�ƁI
%     % (22�b����������)
%     cut_fr = [22];
%     
%     %�폜�ɂ�邸��𖄂߂�ϐ�
%     cut_count = 0;
%     
%     %�w�肵���t���[�����J�b�g
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
% %1��ڂ��҂̊e�t���[���ɂ�����s�[�N���g�����C2��ڂɌċz�g�̊e�t���[���ɂ�����s�[�N���g�����i�[
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
% % % �e���U���g���ɂ�����PSD�̕��ςƕ��U���e���U���g�����Ƃɋ��߂�
% % F_max = 10;
% % cut_flag = 1;
% % rd_fr = 0;
% % 
% % if cut_flag == 1
% %     % �J�b�g�������b
% %     %���������ɋL�q���邱�ƁI
% %     % (22�b�C52�b����������)
% %     cut_fr = [22 52];
% %     
% %     %�폜�ɂ�邸��𖄂߂�ϐ�
% %     cut_count = 0;
% %     
% %     %�w�肵���t���[�����J�b�g
% %     for i = 1:length(cut_fr)
% %         cut_id = fix(cut_fr(i)/win)+1;
% %         lum_psd(cut_id-cut_count, :) = [];
% %         br_psd(cut_id-cut_count, :) = [];
% %         cut_count = cut_count + 1;
% %     end 
% %     rd_fr = length(cut_fr);
% % end
% 
% % �e���U���g����PSD���i�[
% f_comp2 = zeros(M-rd_fr, F_max);
% 
% %1[Hz]�`F_max[Hz]�܂�
% for F_unit=1:F_max
%     %����F_unit�̋ߖT�ɂ��闣�U���g��F�����o��
%     unit = knnsearch(F, F_unit);
%     %
%     f_comp2(:, F_unit) = lum_psd(:, unit);
%     f_comp2(:, F_unit, 2) = br_psd(:, unit);
% end
% 
% % �G�N�Z���ɏo��
% filename = '�e���U���g���ɂ�������ROI�ƌċz�g��psd.xlsx';
% writetable(array2table(f_comp2(:,:,1)),filename,'Sheet',1,'Range','A1:Z100')
% writetable(array2table(f_comp2(:,:,2)),filename,'Sheet',2,'Range','A1:Z100')
% 
% %�e���U���g���ł�psd�̗���t����
% h = zeros(1, F_max);
% pt = zeros(1, F_max);
% 
% for i=1:F_max
%    [h(i),pt(i)] = ttest(f_comp2(:,i,1), f_comp2(:,i,2)); 
% end
% 
% f_num = length(lum_psd);    % ���U���g���̐�
% 
% % �{ROI�ƌċz�g�ɂ�����e���U���g����psd���ςƕ��U
% psd_mean = zeros(f_num, 1);
% psd_var = zeros(f_num, 1);
% psd_mean_br = zeros(f_num, 1);
% psd_var_br = zeros(f_num, 1);
% 
% %��ROI�̊e���U���g���ɂ�����psd�̕��ςƕ��U���v�Z
% for i=1:f_num
%         psd_mean(i,1) = mean(lum_psd(:,i));
%         psd_var(i,1) = var(lum_psd(:,i));
% end
% 
% %�ċz�g�̊e���U���g���ɂ�����psd�̕��ςƕ��U���v�Z
% for i=1:f_num
%         psd_mean_br(i,1) = mean(br_psd(:,i));
%         psd_var_br(i,1) = var(br_psd(:,i));
% end
% 
% 
% % [h_mean,p_mean] = ttest(psd_mean(:,1), psd_mean_br(:,1));
% % [h_var,p_var] = ttest(psd_var(:,1), psd_var_br(:,1));
% 
% % ���ςƕ��U�̏W�c���G�N�Z���ɏo��
% filename = '���ςƕ��U�̌v�Z����(��ROI�ƌċz�g).xlsx';
% writetable(array2table(psd_mean(:,1)),filename,'Sheet',1,'Range','A1')
% writetable(array2table(psd_mean_br(:,1)),filename,'Sheet',1,'Range','B1')
% writetable(array2table(psd_var(:,1)),filename,'Sheet',2,'Range','A1')
% writetable(array2table(psd_var_br(:,1)),filename,'Sheet',2,'Range','B1')




