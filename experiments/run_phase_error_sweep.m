function sweep_result = run_phase_error_sweep()
%RUN_PHASE_ERROR_SWEEP 扫描阵列相位误差强度，评估校准感知 CLEAN 性能。
%
% 用法:
%   sweep_result = run_phase_error_sweep();
%
% 输出:
%   sweep_result.detail_table   每个相位误差强度、每次 MC 的详细结果。
%   sweep_result.summary_table  每个相位误差强度下的均值/标准差统计。
%
% 物理含义:
%   phase_error_std 越大，真实阵列流形越偏离名义导向矢量。
%   预期现象:
%       1. Baseline CLEAN 残余能量随相位误差增大而上升；
%       2. Calibrated CLEAN 由于估计并补偿相位误差，应接近 Oracle；
%       3. PDR 应在 Calibrated CLEAN 后明显改善。

clc;
close all;

%% 0. 添加项目目录到 MATLAB 路径
root_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(root_dir);
addpath(fullfile(root_dir, 'lib'));
addpath(fullfile(root_dir, 'adaptive'));

%% 1. 基础配置
cfg0 = sonar_config();

% 扫描相位误差标准差，单位 deg。
phase_std_deg_list = [0, 2, 5, 10, 15, 20, 25, 30];

% Monte Carlo 次数。调试建议 3~10；正式论文图建议 50 或 100。
num_mc = 3;

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

%% 2. 创建输出目录
if ~exist(cfg0.result_dir, 'dir')
    mkdir(cfg0.result_dir);
end

if ~exist(cfg0.figure_dir, 'dir')
    mkdir(cfg0.figure_dir);
end

%% 3. 预分配详细结果
num_phase = numel(phase_std_deg_list);
num_runs = num_phase * num_mc;

phase_std_deg = zeros(num_runs, 1);
mc_index_col = zeros(num_runs, 1);

baseline_mask_rel = NaN(num_runs, 1);
calibrated_mask_rel = NaN(num_runs, 1);
oracle_mask_rel = NaN(num_runs, 1);

baseline_pdr_dB = NaN(num_runs, 1);
calibrated_pdr_dB = NaN(num_runs, 1);
oracle_pdr_dB = NaN(num_runs, 1);

baseline_interference_suppression_dB = NaN(num_runs, 1);
calibrated_interference_suppression_dB = NaN(num_runs, 1);
oracle_interference_suppression_dB = NaN(num_runs, 1);

baseline_direct_suppression_dB = NaN(num_runs, 1);
calibrated_direct_suppression_dB = NaN(num_runs, 1);
oracle_direct_suppression_dB = NaN(num_runs, 1);

clean_mask_gain_dB = NaN(num_runs, 1);
pdr_gain_dB = NaN(num_runs, 1);
interference_gain_dB = NaN(num_runs, 1);
direct_gain_dB = NaN(num_runs, 1);

cal_oracle_mask_gap_dB = NaN(num_runs, 1);
cal_oracle_pdr_gap_dB = NaN(num_runs, 1);

phase_rmse_active_deg = NaN(num_runs, 1);
phase_rmse_all_deg = NaN(num_runs, 1);
phase_max_abs_deg = NaN(num_runs, 1);

baseline_iter = NaN(num_runs, 1);
calibrated_iter = NaN(num_runs, 1);
oracle_iter = NaN(num_runs, 1);

runtime_sec = NaN(num_runs, 1);

baseline_stop_reason = cell(num_runs, 1);
calibrated_stop_reason = cell(num_runs, 1);
oracle_stop_reason = cell(num_runs, 1);

ok = false(num_runs, 1);

%% 4. 主扫描循环
fprintf('\n========== Phase Error Sweep ==========\n');
fprintf('num_phase = %d, num_mc = %d, total runs = %d\n', ...
    num_phase, num_mc, num_runs);
fprintf('phase std list / deg:\n');
disp(phase_std_deg_list);

row = 0;

for ip = 1:num_phase
    cur_phase_deg = phase_std_deg_list(ip);
    cur_phase_rad = cur_phase_deg * pi / 180;

    fprintf('\n--- Phase std = %.2f deg ---\n', cur_phase_deg);

    for mc = 1:num_mc
        row = row + 1;

        overrides = struct();
        overrides.phase_error_std = cur_phase_rad;

        try
            result = main_single_run(cfg0, mc, overrides, opts);

            s = result.metrics.summary;
            tbl = result.metrics.table;

            phase_std_deg(row) = cur_phase_deg;
            mc_index_col(row) = mc;

            baseline_mask_rel(row) = s.baseline_clean_mask_rel;
            calibrated_mask_rel(row) = s.calibrated_clean_mask_rel;
            oracle_mask_rel(row) = s.oracle_clean_mask_rel;

            baseline_pdr_dB(row) = s.baseline_pdr_dB;
            calibrated_pdr_dB(row) = s.calibrated_pdr_dB;
            oracle_pdr_dB(row) = s.oracle_pdr_dB;

            baseline_interference_suppression_dB(row) = s.baseline_interference_suppression_dB;
            calibrated_interference_suppression_dB(row) = s.calibrated_interference_suppression_dB;
            oracle_interference_suppression_dB(row) = s.oracle_interference_suppression_dB;

            baseline_direct_suppression_dB(row) = s.baseline_direct_suppression_dB;
            calibrated_direct_suppression_dB(row) = s.calibrated_direct_suppression_dB;
            oracle_direct_suppression_dB(row) = s.oracle_direct_suppression_dB;

            clean_mask_gain_dB(row) = s.calibrated_vs_baseline_clean_mask_gain_dB;
            pdr_gain_dB(row) = s.calibrated_vs_baseline_pdr_gain_dB;
            interference_gain_dB(row) = s.calibrated_vs_baseline_interference_gain_dB;
            direct_gain_dB(row) = s.calibrated_vs_baseline_direct_gain_dB;

            cal_oracle_mask_gap_dB(row) = s.calibrated_oracle_gap_clean_mask_dB;
            cal_oracle_pdr_gap_dB(row) = s.calibrated_oracle_pdr_gap_dB;

            phase_rmse_active_deg(row) = s.phase_rmse_active_deg;
            phase_rmse_all_deg(row) = s.phase_rmse_all_deg;
            phase_max_abs_deg(row) = s.phase_max_abs_deg;

            baseline_iter(row) = get_method_value(tbl, 'Baseline', 'num_iter');
            calibrated_iter(row) = get_method_value(tbl, 'Calibrated', 'num_iter');
            oracle_iter(row) = get_method_value(tbl, 'Oracle', 'num_iter');

            baseline_stop_reason{row} = get_method_text(tbl, 'Baseline', 'stop_reason');
            calibrated_stop_reason{row} = get_method_text(tbl, 'Calibrated', 'stop_reason');
            oracle_stop_reason{row} = get_method_text(tbl, 'Oracle', 'stop_reason');

            runtime_sec(row) = result.runtime_sec;
            ok(row) = result.ok;

            fprintf('MC %3d/%3d: gain = %7.3f dB, PDR gain = %7.3f dB, RMSE = %.4f deg\n', ...
                mc, num_mc, clean_mask_gain_dB(row), pdr_gain_dB(row), ...
                phase_rmse_active_deg(row));

        catch ME
            phase_std_deg(row) = cur_phase_deg;
            mc_index_col(row) = mc;
            ok(row) = false;

            baseline_stop_reason{row} = ['error: ', ME.message];
            calibrated_stop_reason{row} = ['error: ', ME.message];
            oracle_stop_reason{row} = ['error: ', ME.message];

            warning('Run failed at phase %.2f deg, MC %d: %s', ...
                cur_phase_deg, mc, ME.message);
        end
    end
end

%% 5. 详细结果表
detail_table = table();

detail_table.phase_std_deg = phase_std_deg;
detail_table.mc_index = mc_index_col;
detail_table.ok = ok;
detail_table.runtime_sec = runtime_sec;

detail_table.baseline_mask_rel = baseline_mask_rel;
detail_table.calibrated_mask_rel = calibrated_mask_rel;
detail_table.oracle_mask_rel = oracle_mask_rel;

detail_table.baseline_pdr_dB = baseline_pdr_dB;
detail_table.calibrated_pdr_dB = calibrated_pdr_dB;
detail_table.oracle_pdr_dB = oracle_pdr_dB;

detail_table.baseline_interference_suppression_dB = baseline_interference_suppression_dB;
detail_table.calibrated_interference_suppression_dB = calibrated_interference_suppression_dB;
detail_table.oracle_interference_suppression_dB = oracle_interference_suppression_dB;

detail_table.baseline_direct_suppression_dB = baseline_direct_suppression_dB;
detail_table.calibrated_direct_suppression_dB = calibrated_direct_suppression_dB;
detail_table.oracle_direct_suppression_dB = oracle_direct_suppression_dB;

detail_table.clean_mask_gain_dB = clean_mask_gain_dB;
detail_table.pdr_gain_dB = pdr_gain_dB;
detail_table.interference_gain_dB = interference_gain_dB;
detail_table.direct_gain_dB = direct_gain_dB;

detail_table.cal_oracle_mask_gap_dB = cal_oracle_mask_gap_dB;
detail_table.cal_oracle_pdr_gap_dB = cal_oracle_pdr_gap_dB;

detail_table.phase_rmse_active_deg = phase_rmse_active_deg;
detail_table.phase_rmse_all_deg = phase_rmse_all_deg;
detail_table.phase_max_abs_deg = phase_max_abs_deg;

detail_table.baseline_iter = baseline_iter;
detail_table.calibrated_iter = calibrated_iter;
detail_table.oracle_iter = oracle_iter;

detail_table.baseline_stop_reason = baseline_stop_reason;
detail_table.calibrated_stop_reason = calibrated_stop_reason;
detail_table.oracle_stop_reason = oracle_stop_reason;

%% 6. 按相位误差强度汇总统计
summary_table = make_phase_summary(detail_table, phase_std_deg_list);

%% 7. 输出结构体
sweep_result = struct();

sweep_result.model = 'phase_error_std_sweep';
sweep_result.phase_std_deg_list = phase_std_deg_list;
sweep_result.num_mc = num_mc;
sweep_result.detail_table = detail_table;
sweep_result.summary_table = summary_table;
sweep_result.cfg0 = cfg0;
sweep_result.opts = opts;

fprintf('\n========== Sweep Summary ==========\n');
disp(summary_table);

%% 8. 保存结果
if save_results
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');

    mat_file = fullfile(cfg0.result_dir, ...
        ['phase_error_sweep_', timestamp, '.mat']);
    detail_csv = fullfile(cfg0.result_dir, ...
        ['phase_error_sweep_detail_', timestamp, '.csv']);
    summary_csv = fullfile(cfg0.result_dir, ...
        ['phase_error_sweep_summary_', timestamp, '.csv']);

    save(mat_file, 'sweep_result');
    writetable(detail_table, detail_csv);
    writetable(summary_table, summary_csv);

    fprintf('\nSaved results:\n');
    fprintf('MAT     : %s\n', mat_file);
    fprintf('Detail  : %s\n', detail_csv);
    fprintf('Summary : %s\n', summary_csv);
end

%% 9. 绘图
plot_phase_error_sweep(summary_table);

if save_figures
    timestamp = datestr(now, 'yyyymmdd_HHMMSS');
    fig_handles = findall(0, 'Type', 'figure');

    for k = 1:numel(fig_handles)
        fig = fig_handles(k);
        fig_name = fullfile(cfg0.figure_dir, ...
            sprintf('phase_error_sweep_fig%d_%s.png', fig.Number, timestamp));
        saveas(fig, fig_name);
    end
end

fprintf('\nPhase error sweep finished.\n');
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

function summary_table = make_phase_summary(detail_table, phase_std_deg_list)
%MAKE_PHASE_SUMMARY 对每个相位误差强度进行 Monte Carlo 统计。

num_phase = numel(phase_std_deg_list);

summary_table = table();
summary_table.phase_std_deg = phase_std_deg_list(:);
summary_table.num_mc = zeros(num_phase, 1);
summary_table.num_ok = zeros(num_phase, 1);

summary_table.baseline_mask_rel_mean = NaN(num_phase, 1);
summary_table.calibrated_mask_rel_mean = NaN(num_phase, 1);
summary_table.oracle_mask_rel_mean = NaN(num_phase, 1);

summary_table.baseline_mask_rel_std = NaN(num_phase, 1);
summary_table.calibrated_mask_rel_std = NaN(num_phase, 1);
summary_table.oracle_mask_rel_std = NaN(num_phase, 1);

summary_table.baseline_mask_rel_dB_mean = NaN(num_phase, 1);
summary_table.calibrated_mask_rel_dB_mean = NaN(num_phase, 1);
summary_table.oracle_mask_rel_dB_mean = NaN(num_phase, 1);

summary_table.baseline_pdr_dB_mean = NaN(num_phase, 1);
summary_table.calibrated_pdr_dB_mean = NaN(num_phase, 1);
summary_table.oracle_pdr_dB_mean = NaN(num_phase, 1);

summary_table.baseline_pdr_dB_std = NaN(num_phase, 1);
summary_table.calibrated_pdr_dB_std = NaN(num_phase, 1);
summary_table.oracle_pdr_dB_std = NaN(num_phase, 1);

summary_table.clean_mask_gain_dB_mean = NaN(num_phase, 1);
summary_table.clean_mask_gain_dB_std = NaN(num_phase, 1);
summary_table.pdr_gain_dB_mean = NaN(num_phase, 1);
summary_table.pdr_gain_dB_std = NaN(num_phase, 1);

summary_table.interference_gain_dB_mean = NaN(num_phase, 1);
summary_table.direct_gain_dB_mean = NaN(num_phase, 1);

summary_table.cal_oracle_mask_gap_dB_mean = NaN(num_phase, 1);
summary_table.cal_oracle_pdr_gap_dB_mean = NaN(num_phase, 1);

summary_table.phase_rmse_active_deg_mean = NaN(num_phase, 1);
summary_table.phase_rmse_active_deg_std = NaN(num_phase, 1);
summary_table.phase_max_abs_deg_mean = NaN(num_phase, 1);

summary_table.baseline_iter_mean = NaN(num_phase, 1);
summary_table.calibrated_iter_mean = NaN(num_phase, 1);
summary_table.oracle_iter_mean = NaN(num_phase, 1);

for ip = 1:num_phase
    ph = phase_std_deg_list(ip);

    idx = detail_table.phase_std_deg == ph;
    idx_ok = idx & detail_table.ok;

    summary_table.num_mc(ip) = sum(idx);
    summary_table.num_ok(ip) = sum(idx_ok);

    summary_table.baseline_mask_rel_mean(ip) = nanmean_local(detail_table.baseline_mask_rel(idx_ok));
    summary_table.calibrated_mask_rel_mean(ip) = nanmean_local(detail_table.calibrated_mask_rel(idx_ok));
    summary_table.oracle_mask_rel_mean(ip) = nanmean_local(detail_table.oracle_mask_rel(idx_ok));

    summary_table.baseline_mask_rel_std(ip) = nanstd_local(detail_table.baseline_mask_rel(idx_ok));
    summary_table.calibrated_mask_rel_std(ip) = nanstd_local(detail_table.calibrated_mask_rel(idx_ok));
    summary_table.oracle_mask_rel_std(ip) = nanstd_local(detail_table.oracle_mask_rel(idx_ok));

    summary_table.baseline_mask_rel_dB_mean(ip) = ...
        10 * log10(summary_table.baseline_mask_rel_mean(ip));
    summary_table.calibrated_mask_rel_dB_mean(ip) = ...
        10 * log10(summary_table.calibrated_mask_rel_mean(ip));
    summary_table.oracle_mask_rel_dB_mean(ip) = ...
        10 * log10(summary_table.oracle_mask_rel_mean(ip));

    summary_table.baseline_pdr_dB_mean(ip) = nanmean_local(detail_table.baseline_pdr_dB(idx_ok));
    summary_table.calibrated_pdr_dB_mean(ip) = nanmean_local(detail_table.calibrated_pdr_dB(idx_ok));
    summary_table.oracle_pdr_dB_mean(ip) = nanmean_local(detail_table.oracle_pdr_dB(idx_ok));

    summary_table.baseline_pdr_dB_std(ip) = nanstd_local(detail_table.baseline_pdr_dB(idx_ok));
    summary_table.calibrated_pdr_dB_std(ip) = nanstd_local(detail_table.calibrated_pdr_dB(idx_ok));
    summary_table.oracle_pdr_dB_std(ip) = nanstd_local(detail_table.oracle_pdr_dB(idx_ok));

    summary_table.clean_mask_gain_dB_mean(ip) = nanmean_local(detail_table.clean_mask_gain_dB(idx_ok));
    summary_table.clean_mask_gain_dB_std(ip) = nanstd_local(detail_table.clean_mask_gain_dB(idx_ok));

    summary_table.pdr_gain_dB_mean(ip) = nanmean_local(detail_table.pdr_gain_dB(idx_ok));
    summary_table.pdr_gain_dB_std(ip) = nanstd_local(detail_table.pdr_gain_dB(idx_ok));

    summary_table.interference_gain_dB_mean(ip) = nanmean_local(detail_table.interference_gain_dB(idx_ok));
    summary_table.direct_gain_dB_mean(ip) = nanmean_local(detail_table.direct_gain_dB(idx_ok));

    summary_table.cal_oracle_mask_gap_dB_mean(ip) = nanmean_local(detail_table.cal_oracle_mask_gap_dB(idx_ok));
    summary_table.cal_oracle_pdr_gap_dB_mean(ip) = nanmean_local(detail_table.cal_oracle_pdr_gap_dB(idx_ok));

    summary_table.phase_rmse_active_deg_mean(ip) = nanmean_local(detail_table.phase_rmse_active_deg(idx_ok));
    summary_table.phase_rmse_active_deg_std(ip) = nanstd_local(detail_table.phase_rmse_active_deg(idx_ok));
    summary_table.phase_max_abs_deg_mean(ip) = nanmean_local(detail_table.phase_max_abs_deg(idx_ok));

    summary_table.baseline_iter_mean(ip) = nanmean_local(detail_table.baseline_iter(idx_ok));
    summary_table.calibrated_iter_mean(ip) = nanmean_local(detail_table.calibrated_iter(idx_ok));
    summary_table.oracle_iter_mean(ip) = nanmean_local(detail_table.oracle_iter(idx_ok));
end
end

function plot_phase_error_sweep(summary_table)
%PLOT_PHASE_ERROR_SWEEP 绘制相位误差扫描结果。

x = summary_table.phase_std_deg;

figure;
plot(x, summary_table.baseline_mask_rel_dB_mean, '-o', 'LineWidth', 1.5); hold on;
plot(x, summary_table.calibrated_mask_rel_dB_mean, '-s', 'LineWidth', 1.5);
plot(x, summary_table.oracle_mask_rel_dB_mean, '-^', 'LineWidth', 1.5);
grid on;
xlabel('Phase error std / deg');
ylabel('Residual mask energy ratio / dB');
title('CLEAN Residual Energy vs Phase Error');
legend('Baseline', 'Calibrated', 'Oracle', 'Location', 'best');

figure;
plot(x, summary_table.clean_mask_gain_dB_mean, '-o', 'LineWidth', 1.5); hold on;
plot(x, summary_table.pdr_gain_dB_mean, '-s', 'LineWidth', 1.5);
grid on;
xlabel('Phase error std / deg');
ylabel('Gain over baseline / dB');
title('Calibration Gain vs Phase Error');
legend('Clean-mask gain', 'PDR gain', 'Location', 'best');

figure;
plot(x, summary_table.baseline_pdr_dB_mean, '-o', 'LineWidth', 1.5); hold on;
plot(x, summary_table.calibrated_pdr_dB_mean, '-s', 'LineWidth', 1.5);
plot(x, summary_table.oracle_pdr_dB_mean, '-^', 'LineWidth', 1.5);
grid on;
xlabel('Phase error std / deg');
ylabel('Target PDR / dB');
title('Target Peak Dominance vs Phase Error');
legend('Baseline', 'Calibrated', 'Oracle', 'Location', 'best');

figure;
plot(x, summary_table.phase_rmse_active_deg_mean, '-o', 'LineWidth', 1.5);
grid on;
xlabel('Phase error std / deg');
ylabel('Phase estimation RMSE / deg');
title('Phase Estimation Error vs Phase Error');

figure;
plot(x, summary_table.baseline_iter_mean, '-o', 'LineWidth', 1.5); hold on;
plot(x, summary_table.calibrated_iter_mean, '-s', 'LineWidth', 1.5);
plot(x, summary_table.oracle_iter_mean, '-^', 'LineWidth', 1.5);
grid on;
xlabel('Phase error std / deg');
ylabel('Mean iteration count');
title('CLEAN Iterations vs Phase Error');
legend('Baseline', 'Calibrated', 'Oracle', 'Location', 'best');
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
