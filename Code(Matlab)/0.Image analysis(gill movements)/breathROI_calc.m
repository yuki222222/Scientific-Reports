%% 読み込んだ動画の関心領域を複数個指定する．それらの平均輝度値を計算し，比較する．

% ----------------------------------------------------------------

%% 第一フレームの読み取り
close all;clear;clc

Videoname = '計測動画\caf0za.avi';
% Videoname = '無麻酔解析動画\opfisha.avi';

ROIpngname = 'ROI_caf0za.png';
luminance_ROIname = 'matファイル\luminance_ROI_caf0za.mat';

stTime = 0; % 動画の読み取りを開始する時間[s]

Video = VideoReader(Videoname);
videoplay = 0;
img = readFrame(Video); % 動画1フレーム目の画像を格納
figure;imshow(img)

% 動画の再生y
if videoplay == 1
    implay(Videoname)
end
%% 関心領域を複数個指定（左クリックで対象の多角形を作成＝＞右クリック＝＞マスクの作成）

ROI = roipoly;  %関心領域をマウスで指定
index_ROI = zeros(10000, 10000);
index_ROI(1:length(find(ROI>=1)),1) =find(ROI>=1); %指定した関心領域を取り出す(座標情報)
page = 1;

while 1
    prompt = '"y"+"enter"押下でROIの指定終了, "enter"押下で継続\n';
    str = input(prompt,'s');
    if str =='y'
       break 
    end
    page = page + 1;
    ROI(:, :, page)=roipoly;
    index_ROI(1:length(find(ROI(:,:,page)>=1)), page) = find(ROI(:,:,page)>=1);
end

%% ROIを枠線にする処理

flag = 0;
sz = size(ROI(:,:,1));

% ROIの塗りつぶしを無しにする処理
for i = 1:page
    for ln = 1:sz(1,1)
        for rw = 1:sz(1,2)
            if flag == 1 && ROI(ln, rw+1,i) ~= 0
                ROI(ln, rw, i) = 0;
    
            elseif flag == 1 && ROI(ln, rw+1,i) == 0
                flag = 0;
                
            elseif ROI(ln, rw, i) >= 1
                flag = 1;
            end
        end
    end
end

% ROIの値を反転
for i = 1:page
    for ln = 1:sz(1,1)
        for rw = 1:sz(1,2)
            if(ROI(ln,rw,i) > 0)
                ROI(ln,rw,i) = 0;
            else
                ROI(ln,rw,i) = 1;
            end
        end
    end
end

% 全ROI(ROI_all)を重ね描く
ROI_all = ROI(:,:,1);

for i = 1:page-1
    for In = 1:sz(1,1)
        for rw = 1:sz(1,2)
            if ROI(In,rw,i+1) == 0
                ROI_all(In,rw) = 0;
            end
        end
    end
end


imshow(imfuse(ROI_all, img))
imwrite(ROI_all,ROIpngname)
%% 各フレームを読み取り=>そのフレームの輝度を計算

fr_num=1;
Video.CurrentTime=stTime; % 1フレーム目から読み取り
fr_cnt = round(Video.Duration)*round(Video.FrameRate);
int_list=zeros(2,fr_cnt);    % 計算した平均輝度値を格納
int_in_ROI = zeros(10000,10000);    % あるフレームにおけるROI内の輝度値を格納


while hasFrame(Video)
   img = readFrame(Video);
   int=mean(img,3);     % RGB平均
   
   for i=1:page
       %i番目のROIの座標情報を格納
       index_ROI_i=index_ROI(:,i);
       
        % index_ROI(:,i)のゼロ成分を削除
       for j = 1:length(index_ROI_i)
           if index_ROI_i(j)==0
               break
           end
       end
       
       index_ROI_i(j:end,:)=[];
       
       int_in_ROI(1:length(int(index_ROI_i)), i)=int(index_ROI_i);
       
       % フレームfr_numにおける輝度値の平均を計算
       int_list(i, fr_num)=mean(int_in_ROI(:,i)); 
   end
   
   fr_num=fr_num+1;
   if fr_num > fr_cnt
      break 
   end
   
end

% 平均輝度値を正規化
% for i=1:page
%     for j=1:fr_num-1
%         int_list(i,j) = (int_list(i,j) - min(int_list(i,:)))/(max(int_list(i,:)) - min(int_list(i,:)));
%     end
% end

%% グラフの表示

% frame_st = 2;
% 
% for i = 1:page
%     graphname = [num2str(i),'番目のROIにおける各フレームの平均輝度'];
%     figure('Name', graphname)
%     plot(1:fr_num-1,int_list(i,1:fr_num-1)')
%     xlim([frame_st fr_num-1])
%     xlabel('フレーム数')
%     ylabel('平均輝度')
% end
% 
% figure;
% for i = 1:page
%     plot(1:fr_num-1,int_list(i,1:fr_num-1)')
%     xlim([frame_st fr_num-1])
%     hold on
% end
% 
% xlabel('フレーム数')
% ylabel('平均輝度')
% title('各フレームにおけるROIの平均輝度')
% legend('1番目','2番目','3番目','4番目')

%% matファイルに出力
 save(luminance_ROIname)
 
%% 周波数解析（同期性の確認）

% Fs=100;   
% 
% figure('Name', '呼吸波のパワースペクトル');
% for i = 1:page
%     list_bpf=BPF_but(int_list(i,1:fr_num-1), Fs, 1, 6); 
%     [Pxx,F] = pburg(list_bpf-mean(list_bpf),50,1024,Fs);
%     plot(F,10*log10(Pxx))
%     hold on
% end
% 
% 
