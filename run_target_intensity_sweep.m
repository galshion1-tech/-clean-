function sweep_result = run_target_intensity_sweep()
%RUN_TARGET_INTENSITY_SWEEP 扫描目标回波相对直达波强度 Ae。
%
% 用法:
%   sweep_result = run_target_intensity_sweep();
%
% 物理设置:
%   固定直达波功率 P_direct；
%   固定噪声功率 P_noise；
%   改变目标回波功率 P_target；
%
%   Ae = 10log10(P_target / P_direct)
%
% 这样更接近“目标回波强弱变化”的实验，而不是单纯改变直达波强度。
%
% 输出:
%   sweep_result.detail_table
%   sweep_result.summary_table

clc;
close all;

%% 1. 基础配置
cfg0 = sonar_config();

% 目标相对直达波强度，单位 dB。
% 对应目标功率比直达波弱 33 dB 到 15 dB。
Ae_dB_list = -33:3:-15;

% Monte Carlo 次数。调试建议 3~10；正式论文图建议 50 或 100。
num_mc = 3;

% 固定阵列相位误差标准差。默认使用 sonar_config() 中的设置。
fixed_phase_error_std_deg = cfg0.phase_error_std * 180/pi;

% 是否保存结果和图片。
save_results = true;
save_figures = true;

% 运行选项。
opts = struct();
opts.verbose = false;
opts.store_signals = false;
opts.store_clean_data = false;
opts.run_baseline = true;
opts.run_calibrated = true;
opts.run_oracle = true;

%% 2. 固定直达波功率和噪声功率
% 默认配置下:
%   target_power = 1
%   INR = 30 dB
%   direct_power_ref = 1000
%   noise_power_ref = 0.1
%
% 扫描 Ae 时，保持 direct_power_ref 和 noise_power_ref 不变，
% 只改变 target_power。
direct_power_ref = cfg0.target_power * 10^(cfg0.INR/10);
noise_power_ref = cfg0.noise_power;

fprintf('\n========== Target Intensity Sweep ==========\n');
fprintf('Ae list / dB:\n');
disp(Ae_dB_list);
fprintf('num_mc = %d\n', num_mc);
fprintf('fixed phase error std = %.3f deg\n', fixed_phase_error_std_deg);
fprintf('fixed direct power = %.6g\n', direct_power_ref);
fprintf('fixed noise power  = %.6g\n', noise_power_ref);

%% 3. 创建输出目录
if ~exist(cfg0.result_dir, 'dir')
    mkdir(cfg0.result_dir);
end

if ~exist(cfg0.figure_dir, 'dir')
    mkdir(cfg0.figure_dir);
end

%% 4. 预分配详细结果
num_Ae = numel(Ae_dB_list);
num_runs = num_Ae * num_mc;

Ae_dB = zeros(num_runs, 1);
target_power = zeros(num_runs, 1);
INR_dB = zeros(num_runs, 1);
SNR_dB = zeros(num_runs, 1);
mc_index_col = zeros(num_runs, 1);

original_pdr_dB = NaN(num_runs, 1);
baseline_pdr_dB = NaN(num_runs, 1);
calibrated_pdr_dB = NaN(num_runs, 1);
oracle_pdr_dB = NaN(num_runs, 1);

baseline_mask_rel = NaN(num_runs, 1);
calibrated_mask_rel = NaN(num_runs, 1);
oracle_mask_rel = NaN(num_runs, 1);

baseline_interference_suppression_dB = NaN(num_runs, 1);
calibrated_interference_suppression_dB = NaN(num_runs, 1);
oracle_interference_suppression_dB = NaN(num_runs, 1);

clean_mask_gain_dB = NaN(num_runs, 1);
pdr_gain_dB = NaN(num_runs, 1);
interference_gain_dB = NaN(num_runs, 1);

cal_oracle_mask_gap_dB = NaN(num_runs, 1);
cal_oracle_pdr_gap_dB = NaN(num_runs, 1);

phase_rmse_active_deg = NaN(num_runs, 1);
phase_rmse_all_deg = NaN(num_runs, 1);
phase_max_abs_deg = NaN(num_runs, 1);

baseline_target_peak_ratio = NaN(num_runs, 1);
calibrated_target_peak_ratio = NaN(num_runs, 1);
oracle_target_peak_ratio = NaN(num_runs, 1);

baseline_false_amp = NaN(num_runs, 1);
calibrated_false_amp = NaN(num_runs, 1);
oracle_false_amp = NaN(num_runs, 1);

baseline_target_amp = NaN(num_runs, 1);
calibrated_target_amp = NaN(num_runs, 1);
oracle_target_amp = NaN(num_runs, 1);

baseline_iter = NaN(num_runs, 1);
calibrated_iter = NaN(num_runs, 1);
oracle_iter = NaN(num_runs, 1);

runtime_sec = NaN(num_runs, 1);
ok = false(num_runs, 1);

baseline_stop_reason = cell(num_runs, 1);
calibrated_stop_reason = cell(num_runs, 1);
oracle_stop_reason = cell(num_runs, 1);

%% 5. 主扫描循环
row = 0;

for ia = 1:num_Ae
    cur_Ae_dB = Ae_dB_list(ia);

    % 根据 Ae 设置目标功率:
    %   P_target = P_direct * 10^(Ae/10)
    cur_target_power = direct_power_ref * 10^(cur_Ae_dB/10);

    % 为了让 main_single_run 的功率派生逻辑保持一致:
    %   INR = 10log10(P_direct/P_target) = -Ae
    cur_INR_dB = -cur_Ae_dB;

    % 固定噪声功率:
    %   noise_power = target_power / 10^(SNR/10)
    % 所以:
    %   SNR = 10log10(target_power/noise_power_ref)
    cur_SNR_dB = 10 * log10(cur_target_power / noise_power_ref);

    fprintf('\n--- Ae = %.2f dB, target_power = %.4g, INR = %.2f dB, SNR = %.2f dB ---\n', ...
        cur_Ae_dB, cur_target_power, cur_INR_dB, cur_SNR_dB);

    for mc = 1:num_mc
        row = row + 1;

        overrides = struct();
        overrides.power_control_mode = 'hybrid_physical_loss';
        overrides.max_order = 2;
        overrides.target_power = cur_target_power;
        overrides.INR = cur_INR_dB;
        overrides.SNR = cur_SNR_dB;
        overrides.phase_error_std = fixed_phase_error_std_deg * pi/180;

        try
            result = main_single_run(cfg0, mc, overrides, opts);

            s = result.metrics.summary;
            tbl = result.metrics.table;

            Ae_dB(row) = cur_Ae_dB;
            target_power(row) = cur_target_power;
            INR_dB(row) = cur_INR_dB;
            SNR_dB(row) = cur_SNR_dB;
            mc_index_col(row) = mc;

            original_pdr_dB(row) = get_method_value(tbl, 'Original', 'target_pdr_dB');
            baseline_pdr_dB(row) = s.baseline_pdr_dB;
            calibrated_pdr_dB(row) = s.calibrated_pdr_dB;
            oracle_pdr_dB(row) = s.oracle_pdr_dB;

            baseline_mask_rel(row) = s.baseline_clean_mask_rel;
            calibrated_mask_rel(row) = s.calibrated_clean_mask_rel;
            oracle_mask_rel(row) = s.oracle_clean_mask_rel;

            baseline_interference_suppression_dB(row) = s.baseline_interference_suppression_dB;
            calibrated_interference_suppression_dB(row) = s.calibrated_interference_suppression_dB;
            oracle_interference_suppression_dB(row) = s.oracle_interference_suppression_dB;

            clean_mask_gain_dB(row) = s.calibrated_vs_baseline_clean_mask_gain_dB;
            pdr_gain_dB(row) = s.calibrated_vs_baseline_pdr_gain_dB;
            interference_gain_dB(row) = s.calibrated_vs_baseline_interference_gain_dB;

            cal_oracle_mask_gap_dB(row) = s.calibrated_oracle_gap_clean_mask_dB;
            cal_oracle_pdr_gap_dB(row) = s.calibrated_oracle_pdr_gap_dB;

            phase_rmse_active_deg(row) = s.phase_rmse_active_deg;
            phase_rmse_all_deg(row) = s.phase_rmse_all_deg;
            phase_max_abs_deg(row) = s.phase_max_abs_deg;

            baseline_target_peak_ratio(row) = get_method_value(tbl, 'Baseline', 'target_peak_ratio');
            calibrated_target_peak_ratio(row) = get_method_value(tbl, 'Calibrated', 'target_peak_ratio');
            oracle_target_peak_ratio(row) = get_method_value(tbl, 'Oracle', 'target_peak_ratio');

            baseline_false_amp(row) = get_method_value(tbl, 'Baseline', 'false_amp_for_pdr');
            calibrated_false_amp(row) = get_method_value(tbl, 'Calibrated', 'false_amp_for_pdr');
            oracle_false_amp(row) = get_method_value(tbl, 'Oracle', 'false_amp_for_pdr');

            baseline_target_amp(row) = get_method_value(tbl, 'Baseline', 'target_amp_for_pdr');
            calibrated_target_amp(row) = get_method_value(tbl, 'Calibrated', 'target_amp_for_pdr');
            oracle_target_amp(row) = get_method_value(tbl, 'Oracle', 'target_amp_for_pdr');

            baseline_iter(row) = get_method_value(tbl, 'Baseline', 'num_iter');
            calibrated_iter(row) = get_method_value(tbl, 'Calibrated', 'num_iter');
            oracle_iter(row) = get_method_value(tbl, 'Oracle', 'num_iter');

            baseline_stop_reason{row} = get_method_text(tbl, 'Baseline', 'stop_reason');
            calibrated_stop_reason{row} = get_method_text(tbl, 'Calibrated', 'stop_reason');
            oracle_stop_reason{row} = get_method_text(tbl, 'Oracle', 'stop_reason');

            runtime_sec(row) = result.runtime_sec;
            ok(row) = result.ok;

            fprintf('MC %3d/%3d: Base PDR = %7.3f dB, Cal PDR = %7.3f dB, PDR gain = %7.3f dB\n', ...
                mc, num_mc, baseline_pdr_dB(row), calibrated_pdr_dB(row), pdr_gain_dB(row));

        catch ME
            Ae_dB(row) = cur_Ae_dB;
            target_power(row) = cur_target_power;
            INR_dB(row) = cur_INR_dB;
            SNR_dB(row) = cur_SNR_dB;
            mc_index_col(row) = mc;
            ok(row) = false;

            baseline_stop_reason{row} = ['error: ', ME.message];
            calibrated_stop_reason{row} = ['error: ', ME.message];
            oracle_stop_reason{row} = ['error: ', ME.message];

            warning('Run failed at Ae %.2f dB, MC %d: %s', ...
                cur_Ae_dB, mc, ME.message);
        end
    end
end

%% 6. 详细结果表
detail_table = table();

detail_table.Ae_dB = Ae_dB;
detail_table.target_power = target_power;
detail_table.INR_dB = INR_dB;
detail_table.SNR_dB = SNR_dB;
detail_table.mc_index = mc_index_col;
detail_table.ok = ok;
detail_table.runtime_sec = runtime_sec;

detail_table.original_pdr_dB = original_pdr_dB;
detail_table.baseline_pdr_dB = baseline_pdr_dB;
detail_table.calibrated_pdr_dB = calibrated_pdr_dB;
detail_table.oracle_pdr_dB = oracle_pdr_dB;

detail_table.baseline_mask_rel = baseline_mask_rel;
detail_table.calibrated_mask_rel = calibrated_mask_rel;
detail_table.oracle_mask_rel = oracle_mask_rel;

detail_table.baseline_interference_suppression_dB = baseline_interference_suppression_dB;
detail_table.calibrated_interference_suppression_dB = calibrated_interference_suppression_dB;
detail_table.oracle_interference_suppression_dB = oracle_interference_suppression_dB;

detail_table.clean_mask_gain_dB = clean_mask_gain_dB;
detail_table.pdr_gain_dB = pdr_gain_dB;
detail_table.interference_gain_dB = interference_gain_dB;

detail_table.cal_oracle_mask_gap_dB = cal_oracle_mask_gap_dB;
detail_table.cal_oracle_pdr_gap_dB = cal_oracle_pdr_gap_dB;

detail_table.phase_rmse_active_deg = phase_rmse_active_deg;
detail_table.phase_rmse_all_deg = phase_rmse_all_deg;
detail_table.phase_max_abs_deg = phase_max_abs_deg;

detail_table.baseline_target_peak_ratio = baseline_target_peak_ratio;
detail_table.calibrated_target_peak_ratio = calibrated_target_peak_ratio;
detail_table.oracle_target_peak_ratio = oracle_target_peak_ratio;

detail_table.baseline_false_amp = baseline_false_amp;
detail_table.calibrated_false_amp = calibrated_false_amp;
detail_table.oracle_false_amp = oracle_false_amp;

detail_table.baseline_target_amp = baseline_target_amp;
detail_table.calibrated_target_amp = calibrated_target_amp;
detail_table.oracle_target_amp = oracle_target_amp;

detail_table.baseline_iter = baseline_iter;
detail_table.calibrated_iter = calibrated_iter;
detail_table.oracle_iter = oracle_iter;

detail_table.baseline_stop_reason = baseline_stop_reason;
detail_table.calibrated_stop_reason = calibrated_stop_reason;
detail_table.oracle_stop_reason = oracle_stop_reason;

%% 7. 汇总统计
summary_table = make_Ae_summary(detail_table, Ae_dB_list);

%% 8. 输出结构体
sweep_result = struct();

sweep_result.model = 'target_intensity_Ae_sweep';
sweep_result.Ae_dB_list = Ae_dB_list;
sweep_result.num_mc = num_mc;
sweep_result.fixed_phase_error_std_deg = fixed_phase_error_std_deg;
sweep_result.direct_power_ref = direct_power_ref;
sweep_result.noise_power_ref = noise_power_ref;
sweep_result.detail_table = detail_table;
sweep_result.summary_table = summary_table;
sweep_result.cfg0 = cfg0;
sweep_result.opts = opts;

fprintf('\n========== Target Intensity Sweep Summary ==========\n');
disp(summary_table);

%% 9. 保存结果
if save_results
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');

    mat_file = fullfile(cfg0.result_dir, ...
        ['target_intensity_sweep_', timestamp, '.mat']);
    detail_csv = fullfile(cfg0.result_dir, ...
        ['target_intensity_sweep_detail_', timestamp, '.csv']);
    summary_csv = fullfile(cfg0.result_dir, ...
        ['target_intensity_sweep_summary_', timestamp, '.csv']);

    save(mat_file, 'sweep_result');
    writetable(detail_table, detail_csv);
    writetable(summary_table, summary_csv);

    fprintf('\nSaved results:\n');
    fprintf('MAT     : %s\n', mat_file);
    fprintf('Detail  : %s\n', detail_csv);
    fprintf('Summary : %s\n', summary_csv);
end

%% 10. 绘图
plot_target_intensity_sweep(summary_table);

if save_figures
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    fig_handles = findall(0, 'Type', 'figure');

    for k = 1:numel(fig_handles)
        fig = fig_handles(k);
        fig_name = fullfile(cfg0.figure_dir, ...
            sprintf('target_intensity_sweep_fig%d_%s.png', fig.Number, timestamp));
        saveas(fig, fig_name);
    end
end

fprintf('\nTarget intensity sweep finished.\n');
end

%% ========================================================================
% 局部函数
% ========================================================================

function val = get_method_value(tbl, method_name, var_name)
%GET_METHOD_VALUE 从 metrics.table 中读取指定方法的数值字段。

idx = find(strcmp(tbl.method, method_name), 1);

if isempty(idx) || ~ismember(var_name, tbl.Properties.VariableNames)
    val = NaN;
else
    val = tbl.(var_name)(idx);
end
end

function txt = get_method_text(tbl, method_name, var_name)
%GET_METHOD_TEXT 从 metrics.table 中读取指定方法的文本字段。

idx = find(strcmp(tbl.method, method_name), 1);

if isempty(idx) || ~ismember(var_name, tbl.Properties.VariableNames)
    txt = 'missing';
    return;
end

x = tbl.(var_name)(idx);

if iscell(x)
    txt = x{1};
elseif ischar(x)
    txt = x;
else
    txt = sprintf('%s', x);
end
end

function summary_table = make_Ae_summary(detail_table, Ae_dB_list)
%MAKE_AE_SUMMARY 对每个 Ae 进行 Monte Carlo 统计。

num_Ae = numel(Ae_dB_list);

summary_table = table();
summary_table.Ae_dB = Ae_dB_list(:);
summary_table.num_mc = zeros(num_Ae, 1);
summary_table.num_ok = zeros(num_Ae, 1);

summary_table.target_power_mean = NaN(num_Ae, 1);
summary_table.INR_dB_mean = NaN(num_Ae, 1);
summary_table.SNR_dB_mean = NaN(num_Ae, 1);

summary_table.original_pdr_dB_mean = NaN(num_Ae, 1);
summary_table.baseline_pdr_dB_mean = NaN(num_Ae, 1);
summary_table.calibrated_pdr_dB_mean = NaN(num_Ae, 1);
summary_table.oracle_pdr_dB_mean = NaN(num_Ae, 1);

summary_table.original_pdr_dB_std = NaN(num_Ae, 1);
summary_table.baseline_pdr_dB_std = NaN(num_Ae, 1);
summary_table.calibrated_pdr_dB_std = NaN(num_Ae, 1);
summary_table.oracle_pdr_dB_std = NaN(num_Ae, 1);

summary_table.baseline_mask_rel_mean = NaN(num_Ae, 1);
summary_table.calibrated_mask_rel_mean = NaN(num_Ae, 1);
summary_table.oracle_mask_rel_mean = NaN(num_Ae, 1);

summary_table.baseline_mask_rel_dB_mean = NaN(num_Ae, 1);
summary_table.calibrated_mask_rel_dB_mean = NaN(num_Ae, 1);
summary_table.oracle_mask_rel_dB_mean = NaN(num_Ae, 1);

summary_table.clean_mask_gain_dB_mean = NaN(num_Ae, 1);
summary_table.clean_mask_gain_dB_std = NaN(num_Ae, 1);
summary_table.pdr_gain_dB_mean = NaN(num_Ae, 1);
summary_table.pdr_gain_dB_std = NaN(num_Ae, 1);
summary_table.interference_gain_dB_mean = NaN(num_Ae, 1);

summary_table.cal_oracle_mask_gap_dB_mean = NaN(num_Ae, 1);
summary_table.cal_oracle_pdr_gap_dB_mean = NaN(num_Ae, 1);

summary_table.phase_rmse_active_deg_mean = NaN(num_Ae, 1);
summary_table.phase_rmse_active_deg_std = NaN(num_Ae, 1);

summary_table.baseline_false_amp_mean = NaN(num_Ae, 1);
summary_table.calibrated_false_amp_mean = NaN(num_Ae, 1);
summary_table.oracle_false_amp_mean = NaN(num_Ae, 1);

summary_table.baseline_target_amp_mean = NaN(num_Ae, 1);
summary_table.calibrated_target_amp_mean = NaN(num_Ae, 1);
summary_table.oracle_target_amp_mean = NaN(num_Ae, 1);

summary_table.baseline_iter_mean = NaN(num_Ae, 1);
summary_table.calibrated_iter_mean = NaN(num_Ae, 1);
summary_table.oracle_iter_mean = NaN(num_Ae, 1);

for ia = 1:num_Ae
    Ae = Ae_dB_list(ia);

    idx = detail_table.Ae_dB == Ae;
    idx_ok = idx & detail_table.ok;

    summary_table.num_mc(ia) = sum(idx);
    summary_table.num_ok(ia) = sum(idx_ok);

    summary_table.target_power_mean(ia) = nanmean_local(detail_table.target_power(idx_ok));
    summary_table.INR_dB_mean(ia) = nanmean_local(detail_table.INR_dB(idx_ok));
    summary_table.SNR_dB_mean(ia) = nanmean_local(detail_table.SNR_dB(idx_ok));

    summary_table.original_pdr_dB_mean(ia) = nanmean_local(detail_table.original_pdr_dB(idx_ok));
    summary_table.baseline_pdr_dB_mean(ia) = nanmean_local(detail_table.baseline_pdr_dB(idx_ok));
    summary_table.calibrated_pdr_dB_mean(ia) = nanmean_local(detail_table.calibrated_pdr_dB(idx_ok));
    summary_table.oracle_pdr_dB_mean(ia) = nanmean_local(detail_table.oracle_pdr_dB(idx_ok));

    summary_table.original_pdr_dB_std(ia) = nanstd_local(detail_table.original_pdr_dB(idx_ok));
    summary_table.baseline_pdr_dB_std(ia) = nanstd_local(detail_table.baseline_pdr_dB(idx_ok));
    summary_table.calibrated_pdr_dB_std(ia) = nanstd_local(detail_table.calibrated_pdr_dB(idx_ok));
    summary_table.oracle_pdr_dB_std(ia) = nanstd_local(detail_table.oracle_pdr_dB(idx_ok));

    summary_table.baseline_mask_rel_mean(ia) = nanmean_local(detail_table.baseline_mask_rel(idx_ok));
    summary_table.calibrated_mask_rel_mean(ia) = nanmean_local(detail_table.calibrated_mask_rel(idx_ok));
    summary_table.oracle_mask_rel_mean(ia) = nanmean_local(detail_table.oracle_mask_rel(idx_ok));

    summary_table.baseline_mask_rel_dB_mean(ia) = ...
        10 * log10(summary_table.baseline_mask_rel_mean(ia));
    summary_table.calibrated_mask_rel_dB_mean(ia) = ...
        10 * log10(summary_table.calibrated_mask_rel_mean(ia));
    summary_table.oracle_mask_rel_dB_mean(ia) = ...
        10 * log10(summary_table.oracle_mask_rel_mean(ia));

    summary_table.clean_mask_gain_dB_mean(ia) = nanmean_local(detail_table.clean_mask_gain_dB(idx_ok));
    summary_table.clean_mask_gain_dB_std(ia) = nanstd_local(detail_table.clean_mask_gain_dB(idx_ok));
    summary_table.pdr_gain_dB_mean(ia) = nanmean_local(detail_table.pdr_gain_dB(idx_ok));
    summary_table.pdr_gain_dB_std(ia) = nanstd_local(detail_table.pdr_gain_dB(idx_ok));
    summary_table.interference_gain_dB_mean(ia) = nanmean_local(detail_table.interference_gain_dB(idx_ok));

    summary_table.cal_oracle_mask_gap_dB_mean(ia) = nanmean_local(detail_table.cal_oracle_mask_gap_dB(idx_ok));
    summary_table.cal_oracle_pdr_gap_dB_mean(ia) = nanmean_local(detail_table.cal_oracle_pdr_gap_dB(idx_ok));

    summary_table.phase_rmse_active_deg_mean(ia) = nanmean_local(detail_table.phase_rmse_active_deg(idx_ok));
    summary_table.phase_rmse_active_deg_std(ia) = nanstd_local(detail_table.phase_rmse_active_deg(idx_ok));

    summary_table.baseline_false_amp_mean(ia) = nanmean_local(detail_table.baseline_false_amp(idx_ok));
    summary_table.calibrated_false_amp_mean(ia) = nanmean_local(detail_table.calibrated_false_amp(idx_ok));
    summary_table.oracle_false_amp_mean(ia) = nanmean_local(detail_table.oracle_false_amp(idx_ok));

    summary_table.baseline_target_amp_mean(ia) = nanmean_local(detail_table.baseline_target_amp(idx_ok));
    summary_table.calibrated_target_amp_mean(ia) = nanmean_local(detail_table.calibrated_target_amp(idx_ok));
    summary_table.oracle_target_amp_mean(ia) = nanmean_local(detail_table.oracle_target_amp(idx_ok));

    summary_table.baseline_iter_mean(ia) = nanmean_local(detail_table.baseline_iter(idx_ok));
    summary_table.calibrated_iter_mean(ia) = nanmean_local(detail_table.calibrated_iter(idx_ok));
    summary_table.oracle_iter_mean(ia) = nanmean_local(detail_table.oracle_iter(idx_ok));
end
end

function plot_target_intensity_sweep(summary_table)
%PLOT_TARGET_INTENSITY_SWEEP 绘制目标强度扫描结果。

x = summary_table.Ae_dB;

figure;
plot(x, summary_table.original_pdr_dB_mean, '-k', 'LineWidth', 1.2); hold on;
plot(x, summary_table.baseline_pdr_dB_mean, '-o', 'LineWidth', 1.5);
plot(x, summary_table.calibrated_pdr_dB_mean, '-s', 'LineWidth', 1.5);
plot(x, summary_table.oracle_pdr_dB_mean, '-^', 'LineWidth', 1.5);
grid on;
xlabel('Target-to-direct ratio A_e / dB');
ylabel('Target PDR / dB');
title('Target Peak Dominance vs Target Intensity');
legend('Original', 'Baseline', 'Calibrated', 'Oracle', 'Location', 'best');

figure;
plot(x, summary_table.baseline_mask_rel_dB_mean, '-o', 'LineWidth', 1.5); hold on;
plot(x, summary_table.calibrated_mask_rel_dB_mean, '-s', 'LineWidth', 1.5);
plot(x, summary_table.oracle_mask_rel_dB_mean, '-^', 'LineWidth', 1.5);
grid on;
xlabel('Target-to-direct ratio A_e / dB');
ylabel('Residual mask energy ratio / dB');
title('CLEAN Residual Energy vs Target Intensity');
legend('Baseline', 'Calibrated', 'Oracle', 'Location', 'best');

figure;
plot(x, summary_table.clean_mask_gain_dB_mean, '-o', 'LineWidth', 1.5); hold on;
plot(x, summary_table.pdr_gain_dB_mean, '-s', 'LineWidth', 1.5);
grid on;
xlabel('Target-to-direct ratio A_e / dB');
ylabel('Gain over baseline / dB');
title('Calibration Gain vs Target Intensity');
legend('Clean-mask gain', 'PDR gain', 'Location', 'best');

figure;
plot(x, summary_table.baseline_false_amp_mean, '-o', 'LineWidth', 1.5); hold on;
plot(x, summary_table.calibrated_false_amp_mean, '-s', 'LineWidth', 1.5);
plot(x, summary_table.oracle_false_amp_mean, '-^', 'LineWidth', 1.5);
grid on;
xlabel('Target-to-direct ratio A_e / dB');
ylabel('False peak amplitude');
title('Maximum False Peak vs Target Intensity');
legend('Baseline', 'Calibrated', 'Oracle', 'Location', 'best');

figure;
plot(x, summary_table.baseline_target_amp_mean, '-o', 'LineWidth', 1.5); hold on;
plot(x, summary_table.calibrated_target_amp_mean, '-s', 'LineWidth', 1.5);
plot(x, summary_table.oracle_target_amp_mean, '-^', 'LineWidth', 1.5);
grid on;
xlabel('Target-to-direct ratio A_e / dB');
ylabel('Target peak amplitude');
title('Target Peak Amplitude vs Target Intensity');
legend('Baseline', 'Calibrated', 'Oracle', 'Location', 'best');

figure;
plot(x, summary_table.phase_rmse_active_deg_mean, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Target-to-direct ratio A_e / dB');
ylabel('Phase estimation RMSE / deg');
title('Phase Estimation Error vs Target Intensity');
end

function m = nanmean_local(x)
%NANMEAN_LOCAL 不依赖统计工具箱的 NaN 均值。

x = x(isfinite(x));

if isempty(x)
    m = NaN;
else
    m = mean(x);
end
end

function s = nanstd_local(x)
%NANSTD_LOCAL 不依赖统计工具箱的 NaN 标准差。

x = x(isfinite(x));

if numel(x) <= 1
    s = NaN;
else
    s = std(x);
end
end
