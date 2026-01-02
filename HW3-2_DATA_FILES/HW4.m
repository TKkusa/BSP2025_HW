% --- MATLAB HW #4 - 心率變異度 (HRV) 分析完整腳本 (標準版 v3.0) ---
%
% 前提：需安裝 Signal Processing Toolbox
% 包含：
% 1. R 峰偵測：HW3 v2.1 演算法 (使用官方 findpeaks)
% 2. HRV 分析：SDNN, LF, HF (使用官方 pwelch)
%
% -----------------------------------------------------------------

clear; clc; close all;

%% --- 步驟 0: 檢查資料檔 ---
files_list = {'ECG_Wake.mat', 'ECG_N2.mat', 'ECG_N3.mat', 'ECG_REM.mat'}; 

fs_orig = 512;    % 原始取樣率
fs_target = 200;  % 目標取樣率 (HW4 要求)

% 儲存結果的表格
results_table = table();

for i = 1:length(files_list)
    filename = files_list{i};
    disp(['處理中: ', filename]);
    
    % 1. 載入資料 (智慧搜尋變數)
    try
        mat_data = load(filename);
        vars = fieldnames(mat_data);
        
        % 找出長度最長的數值變數
        max_len = 0;
        data_var_name = '';
        for k = 1:length(vars)
            curr_var = mat_data.(vars{k});
            if isnumeric(curr_var) && numel(curr_var) > max_len
                max_len = numel(curr_var);
                data_var_name = vars{k};
            end
        end
        
        if isempty(data_var_name)
            error('在檔案中找不到數值訊號變數');
        end
        
        data_orig = mat_data.(data_var_name);
        
    catch ME
        disp(['錯誤：無法載入或解析 ', filename]);
        disp(['訊息: ' ME.message]);
        continue;
    end
    
    data_orig = data_orig(:); % 確保為行向量

    % 2. 降採樣 (512 -> 200 Hz)
    %    (使用官方 resample，需 Signal Processing Toolbox)
    if exist('resample', 'file')
        data_resampled = resample(data_orig, fs_target, fs_orig);
    else
        error('錯誤: 找不到 resample 函數。請安裝 Signal Processing Toolbox。');
    end
    
    % 3. 執行 R 峰偵測 (使用 HW3 v2.1 演算法)
    try
        [r_peaks, ~] = hw4_pan_tompkins(data_resampled, fs_target);
    catch ME
        error(['呼叫 hw4_pan_tompkins 失敗。\n錯誤訊息: ' ME.message]);
    end
    
    % 4. 執行 HRV 分析 (HW4 Task 2)
    [sdnn, lf, hf, lf_hf_ratio, rr_sec] = analyze_hrv_simple(r_peaks, fs_target, filename);
    
    % 5. 顯示單檔結果
    fprintf('  R 峰數: %d\n', length(r_peaks));
    if ~isempty(rr_sec)
        fprintf('  平均心率: %.2f BPM\n', 60 / mean(rr_sec));
    else
        fprintf('  平均心率: NaN\n');
    end
    fprintf('  SDNN:   %.2f ms\n', sdnn);
    fprintf('  LF:     %.2f ms^2\n', lf);
    fprintf('  HF:     %.2f ms^2\n', hf);
    fprintf('  LF/HF:  %.2f\n', lf_hf_ratio);
    fprintf('--------------------------------------\n');
    
    % 儲存到總表
    if isempty(rr_sec)
        hr_val=NaN; 
    else
        hr_val=60/mean(rr_sec); 
    end
    new_row = {filename, length(r_peaks), hr_val, sdnn, lf, hf, lf_hf_ratio};
    results_table = [results_table; new_row];
end

% 設定表格欄位名稱
results_table.Properties.VariableNames = {'File', 'R_Peaks', 'HeartRate', 'SDNN', 'LF', 'HF', 'LF_HF_Ratio'};

disp(' ');
disp('=== HW4 最終結果總表 ===');
disp(results_table);
disp('--------------------------');


%% ==========================================================
%% ⚠️ 請務必複製以下所有函數，直到最後一個 end！
%% ==========================================================

%% --- 附屬函數 1: HRV 分析函數 (HW4 Task 2) ---
function [sdnn, lf, hf, ratio, rr_sec] = analyze_hrv_simple(r_peaks, fs, label_name)
    % 計算 Time Domain (SDNN) 和 Frequency Domain (LF, HF) 指標
    
    if length(r_peaks) < 2
        sdnn=NaN; lf=NaN; hf=NaN; ratio=NaN; rr_sec=[];
        return;
    end

    % 1. 計算 RR 間期 (秒)
    rr_samples = diff(r_peaks);
    rr_sec = rr_samples / fs;
    
    % 簡單過濾異常的 RR 間期 (Artifact Removal)
    % 排除 < 300ms (200 BPM) 或 > 2000ms (30 BPM) 的極端值
    valid_idx = find(rr_sec > 0.3 & rr_sec < 2.0);
    rr_sec = rr_sec(valid_idx);
    
    if length(rr_sec) < 5
        sdnn=NaN; lf=NaN; hf=NaN; ratio=NaN;
        return;
    end
    
    % --- Time Domain: SDNN ---
    % SDNN: RR 間期的標準差 (單位: ms)
    sdnn = std(rr_sec) * 1000;
    
    % --- Frequency Domain: LF, HF ---
    
    % 建立時間軸
    beat_times = cumsum(rr_sec);
    beat_times = beat_times - beat_times(1);
    
    % 重新取樣設定 (4 Hz)
    fs_hrv = 4; 
    t_interp = 0 : 1/fs_hrv : max(beat_times);
    
    if length(beat_times) < 2
        lf = NaN; hf = NaN; ratio = NaN;
        return;
    end
    
    % 插值
    rr_interp = interp1(beat_times, rr_sec, t_interp, 'spline');
    
    % 去除直流分量 (Detrend)
    rr_interp_detrend = rr_interp - mean(rr_interp);
    
    % 計算功率頻譜密度 (PSD) - 使用 Welch 法 (需 Signal Processing Toolbox)
    nfft = 1024;
    window_len = min(256, length(rr_interp_detrend));
    
    if window_len == 0
        lf=NaN; hf=NaN; ratio=NaN; 
        return; 
    end
    
    if exist('pwelch', 'file')
        [pxx, f] = pwelch(rr_interp_detrend, hamming(window_len), [], nfft, fs_hrv);
    else
        error('找不到 pwelch 函數。請安裝 Signal Processing Toolbox。');
    end
    
    % 定義頻帶 (標準 HRV 定義)
    lf_idx = find(f >= 0.04 & f <= 0.15);
    hf_idx = find(f >= 0.15 & f <= 0.40);
    
    % 積分計算功率 (單位換算：s^2 -> ms^2，需乘 10^6)
    lf_power = trapz(f(lf_idx), pxx(lf_idx)) * 1e6; 
    hf_power = trapz(f(hf_idx), pxx(hf_idx)) * 1e6; 
    
    lf = lf_power;
    hf = hf_power;
    ratio = lf / hf;
    
    % --- 繪圖 (PSD) ---
    figure('Name', ['HRV PSD Analysis: ' label_name]);
    plot(f, 10*log10(pxx)); grid on;
    xlabel('Frequency (Hz)'); ylabel('PSD (dB/Hz)');
    title(['PSD Spectrum: ' strrep(label_name, '_', '\_')]);
    xlim([0 0.5]); 
    
    hold on;
    xline(0.04, 'r--'); xline(0.15, 'r--'); xline(0.40, 'r--');
    
    msg = sprintf('LF/HF = %.2f', ratio);
    text(0.25, max(10*log10(pxx))*0.9, msg, 'Color', 'blue', 'FontSize', 12, 'FontWeight', 'bold');

end % analyze_hrv_simple 結束


%% --- 附屬函數 2: Pan-Tompkins 演算法 (HW3 v2.1 最終版) ---
function [R_locs_final, avg_QRS_width] = hw4_pan_tompkins(ecg_signal, fs)
    % HW3 v2.1 演算法：靜態閾值 (0.25) + 鎖定期刪除 (0.2s)
    
    ecg_signal = ecg_signal(:);
    
    % 1. 帶通濾波
    b1 = [1, zeros(1,6), -1]; a1 = [1, -1];
    filter1_hp = filter(b1, a1, ecg_signal);
    b2 = [1, zeros(1,12), -2, zeros(1,12), 1]; a2 = [1, -2, 1];
    filter1_bp = filter(b2, a2, filter1_hp);
    delay1 = 18; filter1 = [filter1_bp(delay1+1:end); zeros(delay1,1)];

    % 2. 微分
    b_der = [1, 2, 0, -2, -1] / 8; a_der = 1;
    filter2_der = filter(b_der, a_der, filter1);
    delay2 = 2; filter2 = [filter2_der(delay2+1:end); zeros(delay2,1)];

    % 3. 平方
    filter3 = filter2 .^ 2; 

    % 4. 移動窗口積分
    N = 31; b_mwi = ones(1, N) / N; a_mwi = 1;
    filter4_mwi = filter(b_mwi, a_mwi, filter3);
    delay4 = 15; filter4 = [filter4_mwi(delay4+1:end); zeros(delay4,1)];

    % 5. 偵測
    learn_samples = 2 * fs;
    if length(filter4) < learn_samples
        R_locs_final = []; avg_QRS_width = 0; 
        return;
    end
    
    % v2.1 參數
    threshold_coefficient = 0.25; 
    peak_estimate = max(filter4(1:learn_samples));
    noise_estimate = mean(filter4(1:learn_samples));
    thres = noise_estimate + threshold_coefficient * (peak_estimate - noise_estimate);
    
    if thres <= 0
        thres = 0.1 * peak_estimate; 
    end
    
    % 使用官方 findpeaks (需 Signal Processing Toolbox)
    if exist('findpeaks', 'file')
        [~, R_locs] = findpeaks(filter4, 'MinPeakHeight', thres);
    else
        error('找不到 findpeaks 函數。請安裝 Signal Processing Toolbox。');
    end
    
    if isempty(R_locs)
        R_locs_final = []; avg_QRS_width = 0; 
        return; 
    end

    refractory_period = 0.2 * fs;
    R_locs_final = R_locs(1);
    for i = 2:length(R_locs)
        if (R_locs(i) - R_locs_final(end)) > refractory_period
            R_locs_final = [R_locs_final; R_locs(i)];
        end
    end
    R_locs_final = R_locs_final(:);
    
    total_delay = delay1 + delay2;
    R_locs_final = R_locs_final - total_delay;
    
    buffer_zone = 0.5 * fs;
    R_locs_final = R_locs_final(R_locs_final > buffer_zone & R_locs_final < (length(ecg_signal) - buffer_zone));
    avg_QRS_width = 0; 

end % hw4_pan_tompkins 結束