% =========================================================================
% 生醫訊號處理 HW3 - Task 1: Pan-Tompkins 演算法
% 完整版 - 包含所有中間步驟的獨立圖表
% =========================================================================

%% 1. 初始化與定義參數
clear;          
close all;      
clc;            

fs = 200;      % 取樣率 (Hz)
filter_order = 8; 
passband_freq = [5, 15]; 

%% 2. 設計 Butterworth Bandpass Filter (使用 SOS 格式)
nyquist_freq = fs / 2;
Wn = passband_freq / nyquist_freq; 
N = filter_order / 2; 

[sos_bp, g_bp] = butter(N, Wn, 'bandpass'); 

fprintf('Pan-Tompkins 帶通濾波器設計完成 (%d 階, [5-15] Hz, SOS 格式)。\n', filter_order);

%% 3. 載入並應用帶通濾波器 (以 ECG3.dat 為例)
ecg_data = load('ECG3.dat');
ecg_original = ecg_data(:, 1); 
t = (0:length(ecg_original)-1) / fs; % 建立時間軸

ecg_bandpassed = filtfilt(sos_bp, g_bp, ecg_original); 

fprintf('帶通濾波器已應用到 ECG3.dat。\n');


%% 4. 步驟 2: 設計並應用微分濾波器
b_deriv = [1, -1];
a_deriv = 1;
ecg_differentiated = filtfilt(b_deriv, a_deriv, ecg_bandpassed);

fprintf('步驟 2: 微分濾波器已應用。\n');


%% 5. 步驟 3: 逐點平方
ecg_squared = ecg_differentiated .^ 2;

fprintf('步驟 3: 逐點平方已完成。\n');


%% 6. 步驟 4: 設計並應用移動積分窗口
window_duration_ms = 150;
window_length_samples = round((window_duration_ms / 1000) * fs); 
b_integral = ones(1, window_length_samples);
a_integral = 1;
ecg_integrated = filtfilt(b_integral, a_integral, ecg_squared);

fprintf('步驟 4: 移動窗口積分已完成 (窗口 %d 點)。\n', window_length_samples);

% =========================================================================
% 步驟 5: 閾值偵測與參數計算
% =========================================================================

%% 7. 偵測 QRS 峰值
simple_threshold = max(ecg_integrated) / 2;
min_peak_distance_sec = 0.25;
min_peak_distance_samples = round(min_peak_distance_sec * fs);

[pks, locs, w, p] = findpeaks(ecg_integrated, ...
    'MinPeakHeight', simple_threshold, ...
    'MinPeakDistance', min_peak_distance_samples);

fprintf('QRS 峰值偵測完成，共找到 %d 個心跳。\n', length(locs));

%% 8. 計算平均心率與 QRS 寬度
% --- 計算平均心率 ---
rr_intervals_sec = diff(locs) / fs;
bpm = 60 ./ rr_intervals_sec;
avg_bpm = mean(bpm);

% --- 計算平均 QRS 寬度 ---
avg_qrs_width_samples = mean(w);
avg_qrs_width_ms = (avg_qrs_width_samples / fs) * 1000;

% --- 顯示結果 ---
fprintf('--- 分析結果 (ECG3.dat) ---\n');
fprintf('平均心率 (Averaged Heart Rate): %.1f BPM\n', avg_bpm);
fprintf('平均 QRS 寬度 (Averaged QRS Width): %.1f ms\n', avg_qrs_width_ms);

% =========================================================================
% 步驟 6: 視覺化所有圖表
% =========================================================================

%% 9. 繪製所有中間步驟的比較圖 (每個步驟一個獨立視窗)
% 圖 1: 原始訊號
figure; 
plot(t, ecg_original, 'b');
title('原始 ECG 訊號 (ECG3.dat)');
xlabel('時間 (seconds)');
ylabel('Amplitude');
grid on;
axis tight;

% 圖 2: 帶通濾波後
figure; 
plot(t, ecg_bandpassed, 'r');
title('步驟 1: 帶通濾波後 (5-15 Hz)');
xlabel('時間 (seconds)');
ylabel('Amplitude');
grid on;
axis tight;

% 圖 3: 微分後
figure; 
plot(t, ecg_differentiated, 'g');
title('步驟 2: 微分後');
xlabel('時間 (seconds)');
ylabel('Amplitude');
grid on;
axis tight;

% 圖 4: 平方後
figure; 
plot(t, ecg_squared, 'm'); 
title('步驟 3: 平方後');
xlabel('時間 (seconds)');
ylabel('Amplitude');
grid on;
axis tight;

% 圖 5: 移動積分後
figure; 
plot(t, ecg_integrated, 'c'); 
title('步驟 4: 移動窗口積分後');
xlabel('時間 (seconds)');
ylabel('Amplitude');
grid on;
axis tight;

%% 10. 視覺化最終驗證
% 將偵測到的峰值標記回原始 ECG 訊號上
figure;
plot(t, ecg_original, 'b', 'DisplayName', '原始 ECG 訊號');
hold on;
% 標記 R 峰值位置
plot(t(locs), ecg_original(locs), 'rv', 'MarkerFaceColor', 'r', 'DisplayName', '偵測到的 R 波');
title('HW3-1: Pan-Tompkins QRS 偵測結果 (ECG3.dat)');
xlabel('時間 (seconds)');
ylabel('Amplitude');
legend;
grid on;
axis tight;

fprintf('所有 6 張圖表已繪製完成。\n');