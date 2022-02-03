function BPF_v=BPF_but(data,fs, low, high)
% バタワース2*n次フィルタ
% 配列dataの値をpassbandの領域でバンドパス

n = 3;                     % n:BPF(BandPassFilter)の次数
passband=[low high]/(fs/2);    %[low high] pass filter


[ButB,ButA] = butter(n,passband,'bandpass');
assignin('base', 'ButB', ButB);
assignin('base', 'ButA', ButA);
BPF_v=filtfilt(ButB,ButA,data); 