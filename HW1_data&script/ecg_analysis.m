
%% 1. 初始化與載入資料
clear;
close all;
clc;

ecg_data = load('ecg_hfn.dat');
ecg = ecg_data(:, 1);

fs_ecg = 1000;
t_ecg = (0:length(ecg)-1) / fs_ecg;

%% 2. 選取 QRS 波樣板 (Template)
template_start_index = 250;
template_length = 86;
template_end_index = template_start_index + template_length - 1;
qrs_template = ecg(template_start_index:template_end_index);

figure;
plot(qrs_template);
title('您選取的 QRS 樣板波形');
xlabel('取樣點 (Samples)');
ylabel('Amplitude');

%% 3. 使用互相關進行心跳偵測
[correlation_output, ~] = xcorr(ecg, qrs_template);
correlation_output = correlation_output / max(correlation_output);

% --- 可調整參數 ---
detection_threshold = 0.8; % 您可以根據需要調整此閾值
min_beat_distance = 0.5 * fs_ecg;
% -------------------------

% --- !!! 新增的圖表 !!! ---
% 視覺化互相關結果與您設定的閾值
figure;
plot(correlation_output);
hold on;
yline(detection_threshold, 'r--', 'LineWidth', 2);
title('正規化互相關結果與偵測閾值');
xlabel('取樣點 (Samples)');
ylabel('相似度分數 (Normalized Correlation)');
legend('互相關結果', sprintf('偵測閾值 = %.2f', detection_threshold));
grid on;
axis tight; % 自動調整座標軸範圍
hold off;
% --- !!! 新增結束 !!! ---


[pks, locs] = findpeaks(correlation_output, 'MinPeakHeight', detection_threshold, 'MinPeakDistance', min_beat_distance);

% 將 findpeaks 在長向量中找到的位置(locs)，映射回原始 ecg 訊號的正確索引
locs_corrected = locs - (length(ecg) - 1);

% 手動微調標記點位置
manual_shift = -45;  
locs_corrected = locs_corrected - manual_shift;

%% 4. 計算 R-R 間期與 BPM
rr_intervals_samples = diff(locs_corrected);
rr_intervals_ms = (rr_intervals_samples / fs_ecg) * 1000;
bpm = 60 ./ (rr_intervals_samples / fs_ecg);

%% 5. 視覺化最終結果
figure;
plot(t_ecg, ecg);
hold on;
if ~isempty(locs_corrected)
    plot(t_ecg(locs_corrected), ecg(locs_corrected), 'rv', 'MarkerFaceColor', 'r', 'MarkerSize', 8);
    for i = 1:length(bpm)
        label_pos_x = t_ecg(locs_corrected(i)) + (rr_intervals_samples(i) / fs_ecg / 2);
        label_text = sprintf('%.0f ms\n%d BPM', rr_intervals_ms(i), round(bpm(i)));
        text(label_pos_x, max(ecg)*0.9, label_text, 'HorizontalAlignment', 'center', 'FontSize', 10, 'FontWeight', 'bold');
    end
end
hold off;

title('ECG 心跳偵測結果 (Beat Detection Results)');
xlabel('時間 (seconds)');
ylabel('ECG Amplitude');
legend('原始 ECG 訊號', '偵測到的 R 波');
grid on;
axis tight;

fprintf('分析完成。\n');
fprintf('總共偵測到 %d 個心跳。\n', length(locs_corrected));
if ~isempty(bpm)
    fprintf('平均心率為: %.1f BPM。\n', mean(bpm));
else
    fprintf('平均心率為: NaN BPM。\n');
end