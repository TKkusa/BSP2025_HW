%% 1. 初始化與定義參數
clear;          
close all;      
clc;            

fs = 1000;      
filter_order = 8; 
passband_freq = [0.5, 50]; 

%% 2. 設計 Butterworth Bandpass Filter (使用 SOS 格式)
nyquist_freq = fs / 2;
Wn = passband_freq / nyquist_freq; 

N = filter_order / 2; % N 現在是 4

[sos, g] = butter(N, Wn, 'bandpass'); 

fprintf('Butterworth Filter 設計完成 (%d 階, SOS 格式)。\n', filter_order);

%% 3. 載入並應用濾波器到 ECG 訊號
ecg_data = load('ecg_lfn.dat');
ecg_original = ecg_data(:, 1); 

ecg_filtered = filtfilt(sos, g, ecg_original); 

fprintf('濾波器已應用到 ECG 訊號。\n');

% 檢查濾波結果 
min_val = min(ecg_filtered);
max_val = max(ecg_filtered);
fprintf('濾波後訊號範圍: [%f, %f]\n', min_val, max_val);
if any(isnan(ecg_filtered)) || any(isinf(ecg_filtered)) || max_val > 1e6 || min_val < -1e6
    warning('濾波結果可能仍有數值問題 (NaN/Inf 或 幅度過大)。');
else
    fprintf('濾波結果數值看起來在合理範圍內。\n');
end

%% 4. 繪製預覽比較圖
figure;
t = (0:length(ecg_original)-1) / fs; 
plot(t, ecg_original, 'b-', 'DisplayName', '原始 ECG');
hold on;
plot(t, ecg_filtered, 'r-', 'LineWidth', 1.5, 'DisplayName', '濾波後 ECG');
hold off;
title(sprintf('[預覽] 原始 ECG vs. %d階濾波後 ECG (使用 SOS)', filter_order));
xlabel('時間 (seconds)');
ylabel('Amplitude');
legend;
grid on;
xlim([0, 5]); 
ylim([-3, 3]); % 保持合理的 Y 軸範圍