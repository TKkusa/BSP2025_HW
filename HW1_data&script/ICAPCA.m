%% HW5: Biomedical Signal Processing - ICA & PCA Blind Source Separation
% Date: 2025/12/08
% Modified for Dark Mode Visualization
clear; close all; clc;

%% 1. 數據載入與前處理 (Data Loading & Pre-processing)
if exist('eeg1-o1.dat', 'file') && exist('eeg2-f3.dat', 'file') && exist('ecg_hfn.dat', 'file')
    disp('Loading real data files...');
    data_eeg1 = load('eeg1-o1.dat'); 
    data_eeg2 = load('eeg2-f3.dat');
    data_ecg_high = load('ecg_hfn.dat');
    
    % [Source: 141] Downsample ECG
    data_ecg_low = resample(data_ecg_high, 1, 10); 
    
    % [Source: 142, 143] 定義訊號源
    s1 = data_eeg1(461:660)';      
    s2 = data_eeg2(1:200)';        
    s3 = data_ecg_low(1:200)';     
else
    warning('找不到 .dat 檔案，正在生成模擬數據...');
    N = 200; t = (1:N)/100;
    s1 = 0.5*sin(2*pi*10*t) + 0.1*randn(1,N); 
    s2 = 0.3*sawtooth(2*pi*2*t) + 0.1*randn(1,N); 
    s3 = zeros(1,N); s3(20:50:end) = 1; 
    s3 = conv(s3, hamming(5), 'same');
end

% 確保格式正確
if size(s1,1) > size(s1,2), s1 = s1'; end
if size(s2,1) > size(s2,2), s2 = s2'; end
if size(s3,1) > size(s3,2), s3 = s3'; end

S = [s1; s2; s3];
N = length(s1);
time_axis = (0:N-1) * (1000/100); % msec

%% 2. 繪製原始訊號 (Original Components) - 深色模式 / 淺藍線條
% 背景設為黑色 ('k')，線條設為淺藍 ([0.4 0.7 1])
figure('Name', 'Original Components', 'Color', 'k'); 
titles = {'Original s1 (EEG1)', 'Original s2 (EEG2)', 'Original s3 (ECG)'};

for i = 1:3
    subplot(3,1,i); 
    % 使用 RGB 顏色 [0.4 0.7 1] 代表淺藍色
    plot(time_axis, S(i,:), 'Color', [0.4 0.7 1], 'LineWidth', 1.2); 
    title(titles{i}, 'Color', 'w'); % 標題改為白色
    ylabel('Amp', 'Color', 'w'); 
    grid on;
    % 設定座標軸為黑色背景、白色文字
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w', 'GridAlpha', 0.3);
end
xlabel('Time (msec)', 'Color', 'w');

%% 3. 訊號混合 (Signal Mixing) - 深色模式 / 白色線條
A = [0.5, 0.5, 0.5;
     0.2, 0.7, 0.7;
     0.7, 0.4, 0.2;
    -0.5, 0.2, -0.6;
     0.7, -0.5, -0.4];
X = A * S;

% 背景設為黑色 ('k')，線條設為白色 ('w')
figure('Name', 'Mixed Signals', 'Color', 'k');
offset = max(max(abs(X))) * 1.5; 
hold on;
for i = 1:5
    plot(time_axis, X(i,:) + (i-1)*offset, 'Color', 'w', 'LineWidth', 1.2);
end
title('Mixed Signals (Observations X)', 'Color', 'w');
xlabel('Time (msec)', 'Color', 'w'); 
ylabel('Channel Offset', 'Color', 'w');
% 設定座標軸樣式
set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w', 'GridAlpha', 0.3);
grid on; hold off;

%% 4. 任務 1: 執行 ICA (Independent Component Analysis)
disp('Performing ICA...');
num_IC = 3; 

if exist('jadeR', 'file')
    W = jadeR(X, num_IC);
    IC_signals = W * X;
elseif exist('rica', 'file')
    Mdl = rica(X', num_IC); 
    IC_signals = transform(Mdl, X')'; 
elseif exist('fastica', 'file')
    [IC_signals, ~] = fastica(X, 'numOfIC', num_IC, 'verbose', 'off'); 
else
    error('錯誤: 找不到 ICA 函數');
end

% 繪製 ICA 結果 (保持紅色，但背景改深色以統一風格)
figure('Name', 'ICA Results', 'Color', 'k');
for i = 1:num_IC
    subplot(num_IC, 1, i);
    plot(time_axis, IC_signals(i,:), 'r', 'LineWidth', 1.2); % 保持紅色 'r'
    title(['ICA Component ', num2str(i)], 'Color', 'w');
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w', 'GridAlpha', 0.3);
    grid on;
end
xlabel('Time (msec)', 'Color', 'w');

%% 5. 任務 2: 執行 PCA (Principal Component Analysis)
disp('Performing PCA...');
X_centered = X - mean(X, 2);
[U, Sigma, V] = svd(X_centered, 'econ');
PCA_signals = U' * X_centered;
PCA_signals = PCA_signals(1:3, :);

% 繪製 PCA 結果 (保持洋紅色，但背景改深色以統一風格)
figure('Name', 'PCA Results', 'Color', 'k');
for i = 1:3
    subplot(3, 1, i);
    plot(time_axis, PCA_signals(i,:), 'm', 'LineWidth', 1.2); % 保持洋紅色 'm'
    title(['PCA Component ', num2str(i)], 'Color', 'w');
    set(gca, 'Color', 'k', 'XColor', 'w', 'YColor', 'w', 'GridColor', 'w', 'GridAlpha', 0.3);
    grid on;
end
xlabel('Time (msec)', 'Color', 'w');

%% 6. 結果討論
disp('--- Analysis Complete ---');