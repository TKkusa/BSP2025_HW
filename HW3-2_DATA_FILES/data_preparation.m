% =========================================================================
% 生醫訊號處理 HW4 - 資料準備 (基於 HW3 的偵測結果)
% =========================================================================
clear; close all; clc;

% 檔案列表
files = {'ECG_Wake.mat', 'ECG_N2.mat', 'ECG_N3.mat', 'ECG_REM.mat'};
fs = 200; % HW3 的降採樣後頻率
fs_original = 512;

% 儲存結果的結構
rr_data = struct();

fprintf('正在重新執行 HW3 的偵測流程以獲取 R 波位置...\n');

for i = 1:length(files)
    filename = files{i};
    fprintf('處理 %s ... ', filename);
    
    % --- 1. 載入與降採樣 (HW3 邏輯) ---
    try
        data_struct = load(filename);
        fnames = fieldnames(data_struct);
        % 尋找最長的變數
        max_len = 0; sig_name = '';
        for k=1:length(fnames), len=numel(data_struct.(fnames{k})); if len>max_len, max_len=len; sig_name=fnames{k}; end; end
        ecg_orig = data_struct.(sig_name);
        if size(ecg_orig,1)==1, ecg_orig=ecg_orig'; elseif size(ecg_orig,2)>1, ecg_orig=ecg_orig(:,1); end
        
        % 降採樣
        ecg_resampled = resample(ecg_orig, fs, fs_original);
        
        % --- 2. Pan-Tompkins 演算法 (簡化版) ---
        % 帶通濾波 [5-25 Hz] (針對高雜訊調整)
        [sos, g] = butter(4, [5, 25]/(fs/2), 'bandpass');
        ecg_bp = filtfilt(sos, g, ecg_resampled);
        
        % 微分 & 平方
        ecg_diff = filter([1, -1], 1, ecg_bp);
        ecg_sq = ecg_diff .^ 2;
        
        % 移動積分 (150ms)
        win_len = round(0.15 * fs);
        ecg_int = filter(ones(1, win_len), 1, ecg_sq);
        
        % 閾值偵測 (Mean + 2*STD)
        threshold = mean(ecg_int) + 2 * std(ecg_int);
        min_dist = round(0.25 * fs);
        [~, locs] = findpeaks(ecg_int, 'MinPeakHeight', threshold, 'MinPeakDistance', min_dist);
        
        % 延遲校正 (Total delay = 15 samples)
        delay = round((length([1, -1])-1)/2 + (win_len-1)/2);
        locs = locs - delay;
        locs = locs(locs > 0 & locs <= length(ecg_resampled));
        
        % --- 3. 計算 RR 間期 (Task 1) ---
        % 將位置轉換為時間 (秒)
        r_times = locs / fs;
        % 計算 RR 間期 (秒)
        rr_intervals = diff(r_times);
        
        % 儲存到結構中
        short_name = filename(1:end-4); % 去掉 .mat
        rr_data.(short_name).rr = rr_intervals;
        rr_data.(short_name).r_times = r_times(2:end); % 對應 RR 的時間軸
        
        fprintf('完成。偵測到 %d 個心跳。\n', length(locs));
        
    catch ME
        fprintf('失敗: %s\n', ME.message);
    end
end

% 儲存資料供後續分析
save('HW4_RR_Data.mat', 'rr_data');
fprintf('\n所有 RR 間期資料已儲存至 HW4_RR_Data.mat\n');
fprintf('現在可以開始進行 Task 2 (HRV 分析) 了。\n');