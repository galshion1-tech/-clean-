clear; clc; close all;

%% 1. 配置与信号生成
cfg = sonar_config();

[s, mf, ~] = generate_lfm(cfg);
[paths, cfg] = generate_virtual_sources(cfg);
[phi_true, ~] = generate_phase_error(cfg, 1);

[X, ~, ~] = simulate_received_signal(cfg, s, paths, phi_true, 1);
[Y, mf_data, ~] = matched_filter_bank(cfg, X, mf, paths);

%% 2. Baseline CLEAN
[~, clean_base, ~] = clean_baseline(cfg, Y, mf_data, s, mf, paths);

%% 3. 相位误差估计
[phi_hat, ~, ~] = estimate_phase_error(cfg, Y, mf_data, s, mf, paths);

%% 4. Calibration-aware CLEAN
[~, clean_cal, ~] = clean_calibration_aware(cfg, Y, mf_data, s, mf, paths, phi_hat);

%% 5. Oracle CLEAN
[~, clean_oracle, ~] = clean_calibration_aware(cfg, Y, mf_data, s, mf, paths, phi_true);

%% 6. 统一指标计算
metrics = compute_metrics(cfg, Y, mf_data, paths, ...
    clean_base, clean_cal, clean_oracle, phi_true, phi_hat);

%% 7. 检查 metrics 是否真的生成
fprintf('\nmetrics variable generated successfully.\n');

fprintf('\nFields in metrics.table:\n');
disp(metrics.table.Properties.VariableNames.');

fprintf('\nMetrics table:\n');
disp(metrics.table);

%% 8. 检查 target_pdr_dB 字段
if ismember('target_pdr_dB', metrics.table.Properties.VariableNames)
    fprintf('\ntarget_pdr_dB exists.\n');
else
    error('metrics.table 中没有 target_pdr_dB，请检查 compute_metrics.m 是否为最新版。');
end

%% 9. 打印关键结果
fprintf('\nKey metrics:\n');
fprintf('Baseline clean mask rel   = %.6g\n', metrics.summary.baseline_clean_mask_rel);
fprintf('Calibrated clean mask rel = %.6g\n', metrics.summary.calibrated_clean_mask_rel);
fprintf('Oracle clean mask rel     = %.6g\n', metrics.summary.oracle_clean_mask_rel);

fprintf('\nPDR:\n');
fprintf('Baseline PDR   = %.3f dB\n', metrics.summary.baseline_pdr_dB);
fprintf('Calibrated PDR = %.3f dB\n', metrics.summary.calibrated_pdr_dB);
fprintf('Oracle PDR     = %.3f dB\n', metrics.summary.oracle_pdr_dB);

fprintf('\nPhase RMSE active = %.6f deg\n', metrics.phase.rmse_active_deg);

%% 10. 画 PDR 图
figure;
bar(metrics.table{:, 'target_pdr_dB'});
set(gca, 'XTickLabel', metrics.table.method);
ylabel('Target PDR / dB');
title('Target Peak Dominance Ratio');
grid on;