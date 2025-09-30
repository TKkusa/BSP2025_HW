% --- 載入數據 ---
% EEG data (fs = 100 Hz)
f3 = load('eeg1-f3.dat');
o1 = load('eeg1-o1.dat');
o2 = load('eeg1-o2.dat');
fs_eeg = 100; % 取樣率 100 Hz 

% ECG data (fs = 1000 Hz)
ecg = load('ecg_hfn.dat');
fs_ecg = 1000; % 取樣率 1000 Hz

% --- Task 1.1: Autocorrelation of f3 and o1 ---
% 定義時間與索引
time_start_acf = 4.2;
time_end_acf = 4.96;
index_start_acf = time_start_acf * fs_eeg;
index_end_acf = time_end_acf * fs_eeg;

%截取訊號片段
f3_segment = f3(index_start_acf:index_end_acf);
o1_segment = o1(index_start_acf:index_end_acf);

% 計算自相關 (使用 'coeff' 選項進行正規化)
[acf_f3, lags_f3] = xcorr(f3_segment, 'coeff');
[acf_o1, lags_o1] = xcorr(o1_segment, 'coeff');

% 轉換 lag 索引為秒
delay_f3 = lags_f3 / fs_eeg;
delay_o1 = lags_o1 / fs_eeg;

% 繪圖
figure;
subplot(2,1,1);
plot(delay_f3, acf_f3);
title('ACF of f3 channel (4.2-4.96s)');
xlabel('Delay (seconds)');
ylabel('Autocorrelation');

subplot(2,1,2);
plot(delay_o1, acf_o1);
title('ACF of o1 channel (4.2-4.96s)');
xlabel('Delay (seconds)');
ylabel('Autocorrelation');