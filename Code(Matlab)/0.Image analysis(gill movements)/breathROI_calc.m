%% �ǂݍ��񂾓���̊֐S�̈�𕡐��w�肷��D�����̕��ϋP�x�l���v�Z���C��r����D

% ----------------------------------------------------------------

%% ���t���[���̓ǂݎ��
close all;clear;clc

Videoname = '�v������\caf0za.avi';
% Videoname = '��������͓���\opfisha.avi';

ROIpngname = 'ROI_caf0za.png';
luminance_ROIname = 'mat�t�@�C��\luminance_ROI_caf0za.mat';

stTime = 0; % ����̓ǂݎ����J�n���鎞��[s]

Video = VideoReader(Videoname);
videoplay = 0;
img = readFrame(Video); % ����1�t���[���ڂ̉摜���i�[
figure;imshow(img)

% ����̍Đ�y
if videoplay == 1
    implay(Videoname)
end
%% �֐S�̈�𕡐��w��i���N���b�N�őΏۂ̑��p�`���쐬�����E�N���b�N�����}�X�N�̍쐬�j

ROI = roipoly;  %�֐S�̈���}�E�X�Ŏw��
index_ROI = zeros(10000, 10000);
index_ROI(1:length(find(ROI>=1)),1) =find(ROI>=1); %�w�肵���֐S�̈�����o��(���W���)
page = 1;

while 1
    prompt = '"y"+"enter"������ROI�̎w��I��, "enter"�����Ōp��\n';
    str = input(prompt,'s');
    if str =='y'
       break 
    end
    page = page + 1;
    ROI(:, :, page)=roipoly;
    index_ROI(1:length(find(ROI(:,:,page)>=1)), page) = find(ROI(:,:,page)>=1);
end

%% ROI��g���ɂ��鏈��

flag = 0;
sz = size(ROI(:,:,1));

% ROI�̓h��Ԃ��𖳂��ɂ��鏈��
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

% ROI�̒l�𔽓]
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

% �SROI(ROI_all)���d�˕`��
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
%% �e�t���[����ǂݎ��=>���̃t���[���̋P�x���v�Z

fr_num=1;
Video.CurrentTime=stTime; % 1�t���[���ڂ���ǂݎ��
fr_cnt = round(Video.Duration)*round(Video.FrameRate);
int_list=zeros(2,fr_cnt);    % �v�Z�������ϋP�x�l���i�[
int_in_ROI = zeros(10000,10000);    % ����t���[���ɂ�����ROI���̋P�x�l���i�[


while hasFrame(Video)
   img = readFrame(Video);
   int=mean(img,3);     % RGB����
   
   for i=1:page
       %i�Ԗڂ�ROI�̍��W�����i�[
       index_ROI_i=index_ROI(:,i);
       
        % index_ROI(:,i)�̃[���������폜
       for j = 1:length(index_ROI_i)
           if index_ROI_i(j)==0
               break
           end
       end
       
       index_ROI_i(j:end,:)=[];
       
       int_in_ROI(1:length(int(index_ROI_i)), i)=int(index_ROI_i);
       
       % �t���[��fr_num�ɂ�����P�x�l�̕��ς��v�Z
       int_list(i, fr_num)=mean(int_in_ROI(:,i)); 
   end
   
   fr_num=fr_num+1;
   if fr_num > fr_cnt
      break 
   end
   
end

% ���ϋP�x�l�𐳋K��
% for i=1:page
%     for j=1:fr_num-1
%         int_list(i,j) = (int_list(i,j) - min(int_list(i,:)))/(max(int_list(i,:)) - min(int_list(i,:)));
%     end
% end

%% �O���t�̕\��

% frame_st = 2;
% 
% for i = 1:page
%     graphname = [num2str(i),'�Ԗڂ�ROI�ɂ�����e�t���[���̕��ϋP�x'];
%     figure('Name', graphname)
%     plot(1:fr_num-1,int_list(i,1:fr_num-1)')
%     xlim([frame_st fr_num-1])
%     xlabel('�t���[����')
%     ylabel('���ϋP�x')
% end
% 
% figure;
% for i = 1:page
%     plot(1:fr_num-1,int_list(i,1:fr_num-1)')
%     xlim([frame_st fr_num-1])
%     hold on
% end
% 
% xlabel('�t���[����')
% ylabel('���ϋP�x')
% title('�e�t���[���ɂ�����ROI�̕��ϋP�x')
% legend('1�Ԗ�','2�Ԗ�','3�Ԗ�','4�Ԗ�')

%% mat�t�@�C���ɏo��
 save(luminance_ROIname)
 
%% ���g����́i�������̊m�F�j

% Fs=100;   
% 
% figure('Name', '�ċz�g�̃p���[�X�y�N�g��');
% for i = 1:page
%     list_bpf=BPF_but(int_list(i,1:fr_num-1), Fs, 1, 6); 
%     [Pxx,F] = pburg(list_bpf-mean(list_bpf),50,1024,Fs);
%     plot(F,10*log10(Pxx))
%     hold on
% end
% 
% 
