% =========================================================================
% 生醫訊號處理 HW3 - Problem 2: 睡眠 ECG 分析 (v4 - 增加積分圖)
% =========================================================================
%% 1. 初始化與定義參數
clear;          
close all;      
clc;            

fs_original = 512; % 原始訊號取樣率 (Hz)
fs = 200;      % 我們的 Pan-Tompkins 演算法的目標取樣率 (Hz)

file_to_analyze_mat = 'ECG_REM.mat'; % 要分析的檔案 (ECG_Wake.mat, ECG_N2.mat, etc.)

% Pan-Tompkins 演算法參數 (for fs=200Hz)
filter_order = 8; 
passband_freq = [5, 15]; 
window_duration_ms = 150;
min_peak_distance_sec = 0.25;

fprintf('開始分析 %s...\n', file_to_analyze_mat);

%% 2. 載入並降採樣 (Resample)
try
    data_struct = load(file_to_analyze_mat);
    fnames = fieldnames(data_struct);
    
    % 找出 .mat 檔案中長度最長的變數 (這一定是訊號)
    max_len = 0;
    signal_var_name = '';
    for i = 1:length(fnames)
        current_len = numel(data_struct.(fnames{i}));
        if current_len > max_len
            max_len = current_len;
            signal_var_name = fnames{i};
        end
    end
    
    if isempty(signal_var_name)
        error('在 .mat 檔案中找不到任何數據變數。');
    end
    
    ecg_original = data_struct.(signal_var_name);
    fprintf('成功從 .mat 檔案中載入變數: %s\n', signal_var_name);
    
catch ME
    fprintf('錯誤：無法讀取 %s。\n', file_to_analyze_mat);
    fprintf('請確保 .mat 檔案 (例如 ECG_Wake.mat) 與 .m 腳本在同一個資料夾中。\n');
    rethrow(ME);
end

% 確保訊號是 N x 1 的欄位向量
if size(ecg_original, 1) == 1 && size(ecg_original, 2) > 1
    ecg_original = ecg_original';
elseif size(ecg_original, 1) > 1 && size(ecg_original, 2) > 1
    ecg_original = ecg_original(:, 1);
end

fprintf('成功載入原始訊號，長度為: %d 點。\n', length(ecg_original));

if length(ecg_original) < fs_original * 10
    error('載入的訊號長度過短 (%d 點)，請檢查 .mat 檔案內容。', length(ecg_original));
end

% --- !!! 關鍵步驟：降採樣 !!! ---
ecg_resampled = resample(ecg_original, fs, fs_original);
fprintf('訊號已從 %d Hz 降採樣至 %d Hz。新訊號長度: %d 點。\n', fs_original, fs, length(ecg_resampled));

% 建立新的時間軸 (for 200Hz)
t = (0:length(ecg_resampled)-1) / fs; 

% =========================================================================
% 步驟 3: Pan-Tompkins 演算法 (fs=200Hz)
% =========================================================================
% 步驟 3.1: 帶通濾波
nyquist_freq = fs / 2;
Wn = passband_freq / nyquist_freq; 
N = filter_order / 2; 
[sos_bp, g_bp] = butter(N, Wn, 'bandpass'); 
ecg_bandpassed = filtfilt(sos_bp, g_bp, ecg_resampled); 
fprintf('步驟 1: 帶通濾波完成。\n');

% 步驟 3.2: 微分
b_deriv = [1, -1];
a_deriv = 1;
ecg_differentiated = filter(b_deriv, a_deriv, ecg_bandpassed);
fprintf('步驟 2: 微分完成。\n');

% 步驟 3.3: 平方
ecg_squared = ecg_differentiated .^ 2;
fprintf('步驟 3: 平方完成。\n');

% 步驟 3.4: 移動積分窗口
window_length_samples = round((window_duration_ms / 1000) * fs); 
b_integral = ones(1, window_length_samples);
a_integral = 1;
ecg_integrated = filter(b_integral, a_integral, ecg_squared);
fprintf('步驟 4: 移動積分完成。\n');

% =========================================================================
% 步驟 4: 閾值偵測與參數計算
% =========================================================================
% --- 閾值策略：使用百分位數 (Percentile) ---
simple_threshold = prctile(ecg_integrated, 80); % 使用第 80 百分位數作為閾值

min_peak_distance_samples = round(min_peak_distance_sec * fs);
[pks, locs, w, p] = findpeaks(ecg_integrated, ...
    'MinPeakHeight', simple_threshold, ...
    'MinPeakDistance', min_peak_distance_samples);

% 補償 filter 函數引入的延遲 (Group Delay)
delay_deriv = (length(b_deriv) - 1) / 2; % 0.5
delay_integ = (length(b_integral) - 1) / 2; % 14.5 (假設 150ms 窗口 @ 200Hz = 30 點)
total_delay = round(delay_deriv + delay_integ); % 15
locs_corrected = locs - total_delay; 

% --- 過濾無效索引 ---
valid_indices_mask = (locs_corrected > 1) & (locs_corrected <= length(ecg_resampled));
locs_corrected = locs_corrected(valid_indices_mask);
w = w(valid_indices_mask); 
fprintf('QRS 峰值偵測完成，共找到 %d 個心跳 (已應用 %d 點延遲校正)。\n', length(locs_corrected), total_delay);

%% 5. 計算平均心率與 QRS 寬度
rr_intervals_sec = diff(locs_corrected) / fs; 
bpm = 60 ./ rr_intervals_sec;
avg_bpm = mean(bpm);
avg_qrs_width_samples = mean(w);
avg_qrs_width_ms = (avg_qrs_width_samples / fs) * 1000;

% --- 顯示結果 ---
fprintf('--- 分析結果 (%s) ---\n', file_to_analyze_mat);
fprintf('平均心率 (Averaged Heart Rate): %.1f BPM\n', avg_bpm);
fprintf('平均 QRS 寬度 (Averaged QRS Width): %.1f ms\n', avg_qrs_width_ms);

% =========================================================================
% 步驟 6: 視覺化最終驗證
% =========================================================================

% --- !!! 新增的積分圖 (用於偵錯) !!! ---
figure; 
plot(t, ecg_integrated, 'c'); 
hold on;
yline(simple_threshold, 'r--', 'LineWidth', 2, 'DisplayName', '偵測閾值 (80%)');
title('步驟 4: 移動窗口積分後 (用於閾值調整)');
xlabel('時間 (seconds)');
ylabel('Amplitude');
grid on;
axis tight;
legend;
% --- !!! 偵錯圖結束 !!! ---


figure;
plot(t, ecg_resampled, 'b', 'DisplayName', '降採樣後 ECG 訊號 (200Hz)');
hold on;
plot(t(locs_corrected), ecg_resampled(locs_corrected), 'rv', 'MarkerFaceColor', 'r', 'DisplayName', '偵測到的 R 波');
title(sprintf('HW3-2: Pan-Tompkins QRS 偵測結果 (%s)', file_to_analyze_mat));
xlabel('時間 (seconds)');
ylabel('Amplitude');
legend;
grid on;
axis tight;

fprintf('最終驗證圖表已繪製完成。\n');