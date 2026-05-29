function metrics = compute_metrics(cfg, Y, mf_data, paths, clean_base, clean_cal, clean_oracle, clean_mvdr, clean_lms, clean_rls, phi_true, phi_hat)
%COMPUTE_METRICS 统一计算所有方法的直达波抑制效果和相位估计指标。
%
% 用法:
%   metrics = compute_metrics(cfg, Y, mf_data, paths, clean_base, clean_cal);
%   metrics = compute_metrics(cfg, Y, mf_data, paths, clean_base, clean_cal, clean_oracle);
%   metrics = compute_metrics(cfg, Y, mf_data, paths, ...
%       clean_base, clean_cal, clean_oracle, clean_mvdr, clean_lms, clean_rls);
%   metrics = compute_metrics(cfg, Y, mf_data, paths, ...
%       clean_base, clean_cal, clean_oracle, clean_mvdr, clean_lms, clean_rls, ...
%       phi_true, phi_hat);
%
% 输入:
%   cfg          sonar_config() 配置结构体。
%   Y            原始匹配滤波输出，M x N_mf。
%   mf_data      matched_filter_bank() 输出结构体，至少包含 t_mf。
%   paths        generate_virtual_sources() 输出路径结构体。
%   clean_base   clean_baseline() 输出 clean_data。
%   clean_cal    clean_calibration_aware() 输出 clean_data。
%   clean_oracle 可选，Oracle CLEAN 输出 clean_data。没有则传 []。
%   clean_mvdr   可选，MVDR 波束形成器输出 clean_data。没有则传 []。
%   clean_lms    可选，LMS 对消器输出 clean_data。没有则传 []。
%   clean_rls    可选，RLS 对消器输出 clean_data。没有则传 []。
%   phi_true     可选，真实阵列相位误差。
%   phi_hat      可选，估计阵列相位误差。
%
% 输出:
%   metrics      指标结构体。
%
% 指标说明:
%   direct window:
%       直达波时延附近的局部残余能量。
%
%   interference window:
%       direct / surface / bottom / multipath 等强干扰路径附近的联合窗口残余能量。
%
%   target window:
%       目标路径附近的区域能量和峰值变化。该指标表示目标区域变化，
%       不严格等于纯目标能量保持率，因为窗口中可能包含噪声和旁瓣。
%
%   phase RMSE:
%       phi_hat 与 phi_true 的相对相位误差 RMSE，active RMSE 不计参考阵元。
%
%   PDR:
%       target peak dominance ratio，定义为
%           PDR = 20log10(A_t / A_f)
%       其中 A_t 是目标窗口内阵列非相干幅度峰值，A_f 是目标保护窗外
%       的最大虚警峰值。

%% 1. 输入检查
if nargin < 7
    clean_oracle = [];
end

if nargin < 8
    clean_mvdr = [];
end

if nargin < 9
    clean_lms = [];
end

if nargin < 10
    clean_rls = [];
end

if nargin < 11
    phi_true = [];
end

if nargin < 12
    phi_hat = [];
end

required_fields = { ...
    'M', ...
    'direct_window', ...
    'multipath_window', ...
    'target_window'};

for k = 1:numel(required_fields)
    field_name = required_fields{k};
    assert(isfield(cfg, field_name), 'cfg.%s 是必需字段。', field_name);
end

assert(isfield(mf_data, 't_mf'), 'mf_data.t_mf 是必需字段。');
assert(~isempty(paths), 'paths 不能为空。');

M = cfg.M;
t_mf = mf_data.t_mf(:);
N_mf = numel(t_mf);

assert(size(Y, 1) == M, 'Y 的行数必须等于 cfg.M。');
assert(size(Y, 2) == N_mf, 'Y 的列数必须等于 numel(mf_data.t_mf)。');
assert(all(isfinite(Y(:))), 'Y 中存在非有限数。');

%% 2. 构造物理评价窗口
[masks, window_info] = build_metric_masks(cfg, t_mf, paths);

%% 3. 整理待评价方法
method_names = {};
method_signals = {};
method_clean_data = {};

method_names{end+1, 1} = 'Original';
method_signals{end+1, 1} = Y;
method_clean_data{end+1, 1} = [];

if ~isempty(clean_base)
    assert(isfield(clean_base, 'Y_residual'), ...
        'clean_base 必须包含 Y_residual 字段。');
    method_names{end+1, 1} = 'Baseline';
    method_signals{end+1, 1} = clean_base.Y_residual;
    method_clean_data{end+1, 1} = clean_base;
end

if ~isempty(clean_cal)
    assert(isfield(clean_cal, 'Y_residual'), ...
        'clean_cal 必须包含 Y_residual 字段。');
    method_names{end+1, 1} = 'Calibrated';
    method_signals{end+1, 1} = clean_cal.Y_residual;
    method_clean_data{end+1, 1} = clean_cal;
end

if ~isempty(clean_oracle)
    assert(isfield(clean_oracle, 'Y_residual'), ...
        'clean_oracle 必须包含 Y_residual 字段。');
    method_names{end+1, 1} = 'Oracle';
    method_signals{end+1, 1} = clean_oracle.Y_residual;
    method_clean_data{end+1, 1} = clean_oracle;
end

if ~isempty(clean_mvdr)
    assert(isfield(clean_mvdr, 'Y_residual'), ...
        'clean_mvdr 必须包含 Y_residual 字段。');
    method_names{end+1, 1} = 'MVDR';
    method_signals{end+1, 1} = clean_mvdr.Y_residual;
    method_clean_data{end+1, 1} = clean_mvdr;
end

if ~isempty(clean_lms)
    assert(isfield(clean_lms, 'Y_residual'), ...
        'clean_lms 必须包含 Y_residual 字段。');
    method_names{end+1, 1} = 'LMS';
    method_signals{end+1, 1} = clean_lms.Y_residual;
    method_clean_data{end+1, 1} = clean_lms;
end

if ~isempty(clean_rls)
    assert(isfield(clean_rls, 'Y_residual'), ...
        'clean_rls 必须包含 Y_residual 字段。');
    method_names{end+1, 1} = 'RLS';
    method_signals{end+1, 1} = clean_rls.Y_residual;
    method_clean_data{end+1, 1} = clean_rls;
end

num_methods = numel(method_names);

%% 4. 计算每个方法的窗口能量和峰值
direct_energy = zeros(num_methods, 1);
interference_energy = zeros(num_methods, 1);
target_energy = zeros(num_methods, 1);
total_energy = zeros(num_methods, 1);

direct_peak = zeros(num_methods, 1);
interference_peak = zeros(num_methods, 1);
target_peak = zeros(num_methods, 1);
target_pdr_dB = NaN(num_methods, 1);
target_amp_for_pdr = NaN(num_methods, 1);
false_amp_for_pdr = NaN(num_methods, 1);
target_peak_delay_for_pdr = NaN(num_methods, 1);
false_peak_delay_for_pdr = NaN(num_methods, 1);

clean_mask_rel = NaN(num_methods, 1);
clean_total_rel = NaN(num_methods, 1);
num_iter = NaN(num_methods, 1);
stop_reason = cell(num_methods, 1);
uses_phase_calibration = false(num_methods, 1);

for i = 1:num_methods
    Yi = method_signals{i};

    assert(size(Yi, 1) == M && size(Yi, 2) == N_mf, ...
        '方法 %s 的 Y_residual 尺寸与原始 Y 不一致。', method_names{i});
    assert(all(isfinite(Yi(:))), ...
        '方法 %s 的结果中存在非有限数。', method_names{i});

    direct_energy(i) = energy_in_mask(Yi, masks.direct);
    interference_energy(i) = energy_in_mask(Yi, masks.interference);
    target_energy(i) = energy_in_mask(Yi, masks.target);
    total_energy(i) = energy_sum(Yi);

    direct_peak(i) = incoherent_peak_in_mask(Yi, masks.direct);
    interference_peak(i) = incoherent_peak_in_mask(Yi, masks.interference);
    target_peak(i) = incoherent_peak_in_mask(Yi, masks.target);

    [target_pdr_dB(i), target_amp_for_pdr(i), false_amp_for_pdr(i), ...
        target_peak_delay_for_pdr(i), false_peak_delay_for_pdr(i)] = ...
        compute_pdr(Yi, t_mf, masks.target, masks.false_peak);

    cdata = method_clean_data{i};

    if isempty(cdata)
        clean_mask_rel(i) = NaN;
        clean_total_rel(i) = NaN;
        num_iter(i) = 0;
        stop_reason{i} = 'not_cleaned';
        uses_phase_calibration(i) = false;
    else
        clean_mask_rel(i) = get_field_or_nan(cdata, 'final_mask_energy_rel');
        clean_total_rel(i) = get_field_or_nan(cdata, 'final_total_energy_rel');
        num_iter(i) = get_field_or_nan(cdata, 'num_iter');

        if isfield(cdata, 'stop_reason')
            stop_reason{i} = cdata.stop_reason;
        else
            stop_reason{i} = 'unknown';
        end

        uses_phase_calibration(i) = isfield(cdata, 'phi_hat');
    end
end

%% 5. 以 Original 为基准计算抑制比和目标区域变化
original_direct_energy = direct_energy(1);
original_interference_energy = interference_energy(1);
original_target_energy = target_energy(1);

original_direct_peak = direct_peak(1);
original_interference_peak = interference_peak(1);
original_target_peak = target_peak(1);

direct_suppression_dB = zeros(num_methods, 1);
interference_suppression_dB = zeros(num_methods, 1);
direct_peak_suppression_dB = zeros(num_methods, 1);
interference_peak_suppression_dB = zeros(num_methods, 1);

target_energy_ratio = zeros(num_methods, 1);
target_peak_ratio = zeros(num_methods, 1);
target_energy_change_dB = zeros(num_methods, 1);
target_peak_change_dB = zeros(num_methods, 1);

for i = 1:num_methods
    direct_suppression_dB(i) = safe_db_ratio(original_direct_energy, direct_energy(i));
    interference_suppression_dB(i) = safe_db_ratio(original_interference_energy, interference_energy(i));

    direct_peak_suppression_dB(i) = safe_db_ratio(original_direct_peak, direct_peak(i));
    interference_peak_suppression_dB(i) = safe_db_ratio(original_interference_peak, interference_peak(i));

    target_energy_ratio(i) = safe_ratio(target_energy(i), original_target_energy);
    target_peak_ratio(i) = safe_ratio(target_peak(i), original_target_peak);
    target_energy_change_dB(i) = safe_db_ratio(target_energy(i), original_target_energy);
    target_peak_change_dB(i) = safe_db_ratio(target_peak(i), original_target_peak);
end

%% 6. 相位估计误差
phase_metrics = struct();

if ~isempty(phi_true) && ~isempty(phi_hat)
    phi_true = phi_true(:);
    phi_hat = phi_hat(:);

    assert(numel(phi_true) == M, 'phi_true 长度必须等于 cfg.M。');
    assert(numel(phi_hat) == M, 'phi_hat 长度必须等于 cfg.M。');

    phase_err = angle(exp(1j * (phi_hat - phi_true)));

    if isfield(cfg, 'ref_sensor')
        ref_sensor = cfg.ref_sensor;
    else
        ref_sensor = 1;
    end

    active_idx = setdiff((1:M).', ref_sensor);

    phase_metrics.available = true;
    phase_metrics.ref_sensor = ref_sensor;
    phase_metrics.active_idx = active_idx;
    phase_metrics.phase_err = phase_err;
    phase_metrics.phase_err_deg = phase_err * 180/pi;

    phase_metrics.rmse_all_rad = sqrt(mean(phase_err.^2));
    phase_metrics.rmse_all_deg = phase_metrics.rmse_all_rad * 180/pi;

    phase_metrics.rmse_active_rad = sqrt(mean(phase_err(active_idx).^2));
    phase_metrics.rmse_active_deg = phase_metrics.rmse_active_rad * 180/pi;

    phase_metrics.max_abs_rad = max(abs(phase_err));
    phase_metrics.max_abs_deg = phase_metrics.max_abs_rad * 180/pi;

    phase_metrics.phi_true_deg = phi_true * 180/pi;
    phase_metrics.phi_hat_deg = phi_hat * 180/pi;
else
    phase_metrics.available = false;
    phase_metrics.ref_sensor = NaN;
    phase_metrics.active_idx = [];
    phase_metrics.rmse_all_rad = NaN;
    phase_metrics.rmse_all_deg = NaN;
    phase_metrics.rmse_active_rad = NaN;
    phase_metrics.rmse_active_deg = NaN;
    phase_metrics.max_abs_rad = NaN;
    phase_metrics.max_abs_deg = NaN;
end

%% 7. 汇总表
metrics_table = table();

metrics_table.method = method_names;
metrics_table.uses_phase_calibration = uses_phase_calibration;
metrics_table.num_iter = num_iter;
metrics_table.stop_reason = stop_reason;

metrics_table.direct_energy = direct_energy;
metrics_table.interference_energy = interference_energy;
metrics_table.target_energy = target_energy;
metrics_table.total_energy = total_energy;

metrics_table.direct_suppression_dB = direct_suppression_dB;
metrics_table.interference_suppression_dB = interference_suppression_dB;

metrics_table.direct_peak = direct_peak;
metrics_table.interference_peak = interference_peak;
metrics_table.target_peak = target_peak;
metrics_table.target_pdr_dB = target_pdr_dB;
metrics_table.target_amp_for_pdr = target_amp_for_pdr;
metrics_table.false_amp_for_pdr = false_amp_for_pdr;
metrics_table.target_peak_delay_for_pdr = target_peak_delay_for_pdr;
metrics_table.false_peak_delay_for_pdr = false_peak_delay_for_pdr;

metrics_table.direct_peak_suppression_dB = direct_peak_suppression_dB;
metrics_table.interference_peak_suppression_dB = interference_peak_suppression_dB;

metrics_table.target_energy_ratio = target_energy_ratio;
metrics_table.target_peak_ratio = target_peak_ratio;
metrics_table.target_energy_change_dB = target_energy_change_dB;
metrics_table.target_peak_change_dB = target_peak_change_dB;

metrics_table.clean_mask_rel = clean_mask_rel;
metrics_table.clean_total_rel = clean_total_rel;

%% 8. 输出结构体
metrics = struct();

metrics.model = 'direct_path_suppression_metrics';
metrics.note = 'target window metrics describe target-region change, not pure target-only preservation';

metrics.table = metrics_table;
metrics.method_names = method_names;

metrics.masks = masks;
metrics.window_info = window_info;
metrics.phase = phase_metrics;

metrics.reference.original_direct_energy = original_direct_energy;
metrics.reference.original_interference_energy = original_interference_energy;
metrics.reference.original_target_energy = original_target_energy;
metrics.reference.original_direct_peak = original_direct_peak;
metrics.reference.original_interference_peak = original_interference_peak;
metrics.reference.original_target_peak = original_target_peak;

metrics.summary = make_summary(metrics_table, phase_metrics);
end

%% ========================================================================
% 局部辅助函数
% ========================================================================

function [masks, window_info] = build_metric_masks(cfg, t_mf, paths)
%BUILD_METRIC_MASKS 构造 direct/interference/target 评价窗口。

N_mf = numel(t_mf);

masks = struct();
masks.direct = false(N_mf, 1);
masks.interference = false(N_mf, 1);
masks.target = false(N_mf, 1);
masks.target_guard = false(N_mf, 1);
masks.false_peak = true(N_mf, 1);

window_info = struct();
window_info.entries = struct('type', {}, 'delay', {}, 'half_window', {}, ...
    'role', {}, 'idx_start', {}, 'idx_end', {});

if isfield(cfg, 'target_guard_window')
    target_guard_window = cfg.target_guard_window;
else
    target_guard_window = 2 * cfg.target_window;
end

for p = 1:numel(paths)
    assert(isfield(paths(p), 'delay'), 'paths(%d).delay 缺失。', p);
    assert(isfield(paths(p), 'type'), 'paths(%d).type 缺失。', p);

    tau0 = paths(p).delay;
    is_direct = strcmp(paths(p).type, 'direct');

    if isfield(paths(p), 'is_interference')
        is_interference = paths(p).is_interference;
    else
        is_interference = ~strcmp(paths(p).type, 'target');
    end

    if isfield(paths(p), 'is_target')
        is_target = paths(p).is_target;
    else
        is_target = strcmp(paths(p).type, 'target');
    end

    if is_direct
        local_mask = abs(t_mf - tau0) <= cfg.direct_window;
        masks.direct = masks.direct | local_mask;
        window_info.entries(end+1) = make_window_entry(paths(p).type, tau0, ...
            cfg.direct_window, 'direct', local_mask); %#ok<AGROW>
    end

    if is_interference
        if isfield(cfg, 'interference_window')
            half_window = cfg.interference_window;
        elseif is_direct
            half_window = cfg.direct_window;
        else
            half_window = cfg.multipath_window;
        end

        local_mask = abs(t_mf - tau0) <= half_window;
        masks.interference = masks.interference | local_mask;

        window_info.entries(end+1) = make_window_entry(paths(p).type, tau0, ...
            half_window, 'interference', local_mask); %#ok<AGROW>
    end

    if is_target
        local_mask = abs(t_mf - tau0) <= cfg.target_window;
        masks.target = masks.target | local_mask;

        window_info.entries(end+1) = make_window_entry(paths(p).type, tau0, ...
            cfg.target_window, 'target', local_mask); %#ok<AGROW>

        guard_mask = abs(t_mf - tau0) <= target_guard_window;
        masks.target_guard = masks.target_guard | guard_mask;

        window_info.entries(end+1) = make_window_entry(paths(p).type, tau0, ...
            target_guard_window, 'target_guard', guard_mask); %#ok<AGROW>
    end
end

masks.false_peak = ~masks.target_guard;

window_info.num_direct_bins = sum(masks.direct);
window_info.num_interference_bins = sum(masks.interference);
window_info.num_target_bins = sum(masks.target);
window_info.num_target_guard_bins = sum(masks.target_guard);
window_info.num_false_peak_bins = sum(masks.false_peak);
window_info.target_guard_window = target_guard_window;

if window_info.num_direct_bins == 0
    warning('direct 评价窗口为空。');
end

if window_info.num_interference_bins == 0
    warning('interference 评价窗口为空。');
end

if window_info.num_target_bins == 0
    warning('target 评价窗口为空。');
end

if window_info.num_false_peak_bins == 0
    warning('false peak 评价窗口为空，PDR 将无法计算。');
end
end

function entry = make_window_entry(type, delay, half_window, role, local_mask)
%MAKE_WINDOW_ENTRY 构造窗口信息条目。

idx = find(local_mask);

entry = struct();
entry.type = type;
entry.delay = delay;
entry.half_window = half_window;
entry.role = role;

if isempty(idx)
    entry.idx_start = NaN;
    entry.idx_end = NaN;
else
    entry.idx_start = idx(1);
    entry.idx_end = idx(end);
end
end

function e = energy_sum(X)
%ENERGY_SUM 计算复矩阵总能量。

e = sum(abs(X(:)).^2);
end

function e = energy_in_mask(Y, mask)
%ENERGY_IN_MASK 计算指定时延窗口内的阵列总能量。

if ~any(mask)
    e = 0;
else
    e = energy_sum(Y(:, mask));
end
end

function p = incoherent_peak_in_mask(Y, mask)
%INCOHERENT_PEAK_IN_MASK 计算指定窗口内的非相干峰值。

if ~any(mask)
    p = NaN;
    return;
end

power_incoh = sum(abs(Y(:, mask)).^2, 1);
p = max(power_incoh);
end

function [pdr_dB, At, Af, tau_t, tau_f] = compute_pdr(Y, t_mf, target_mask, false_mask)
%COMPUTE_PDR 计算目标峰值占优度 PDR。
%
% 阵列非相干幅度谱:
%   A(tau) = sqrt(sum_m |Y_m(tau)|^2)
%
% 目标峰值占优度:
%   PDR = 20log10(A_t/A_f)
%
% 其中 A_t 是目标窗口内最大峰值，A_f 是目标保护窗外最大虚警峰值。

amp_profile = sqrt(sum(abs(Y).^2, 1)).';

if ~any(target_mask) || ~any(false_mask)
    pdr_dB = NaN;
    At = NaN;
    Af = NaN;
    tau_t = NaN;
    tau_f = NaN;
    return;
end

target_amp = amp_profile;
target_amp(~target_mask) = -inf;
[At, idx_t] = max(target_amp);

false_amp = amp_profile;
false_amp(~false_mask) = -inf;
[Af, idx_f] = max(false_amp);

if At <= 0 || Af <= 0 || ~isfinite(At) || ~isfinite(Af)
    pdr_dB = NaN;
else
    pdr_dB = 20 * log10(At / Af);
end

tau_t = t_mf(idx_t);
tau_f = t_mf(idx_f);
end

function r = safe_ratio(num, den)
%SAFE_RATIO 安全比值。

if ~isfinite(num) || ~isfinite(den) || abs(den) <= eps
    r = NaN;
else
    r = num / den;
end
end

function v = safe_db_ratio(num, den)
%SAFE_DB_RATIO 计算 10log10(num/den)。

if ~isfinite(num) || ~isfinite(den) || num <= 0 || den <= 0
    v = NaN;
else
    v = 10 * log10(num / den);
end
end

function v = get_field_or_nan(s, field_name)
%GET_FIELD_OR_NAN 从结构体读取字段，不存在则返回 NaN。

if isstruct(s) && isfield(s, field_name)
    v = s.(field_name);
else
    v = NaN;
end
end

function summary = make_summary(metrics_table, phase_metrics)
%MAKE_SUMMARY 生成常用摘要指标。

summary = struct();

method_list = metrics_table.method;

idx_base = find(strcmp(method_list, 'Baseline'), 1);
idx_cal = find(strcmp(method_list, 'Calibrated'), 1);
idx_oracle = find(strcmp(method_list, 'Oracle'), 1);

summary.has_baseline = ~isempty(idx_base);
summary.has_calibrated = ~isempty(idx_cal);
summary.has_oracle = ~isempty(idx_oracle);

if summary.has_baseline
    summary.baseline_interference_suppression_dB = ...
        metrics_table.interference_suppression_dB(idx_base);
    summary.baseline_direct_suppression_dB = ...
        metrics_table.direct_suppression_dB(idx_base);
    summary.baseline_clean_mask_rel = ...
        metrics_table.clean_mask_rel(idx_base);
    summary.baseline_pdr_dB = ...
        metrics_table.target_pdr_dB(idx_base);
else
    summary.baseline_interference_suppression_dB = NaN;
    summary.baseline_direct_suppression_dB = NaN;
    summary.baseline_clean_mask_rel = NaN;
    summary.baseline_pdr_dB = NaN;
end

if summary.has_calibrated
    summary.calibrated_interference_suppression_dB = ...
        metrics_table.interference_suppression_dB(idx_cal);
    summary.calibrated_direct_suppression_dB = ...
        metrics_table.direct_suppression_dB(idx_cal);
    summary.calibrated_clean_mask_rel = ...
        metrics_table.clean_mask_rel(idx_cal);
    summary.calibrated_pdr_dB = ...
        metrics_table.target_pdr_dB(idx_cal);
else
    summary.calibrated_interference_suppression_dB = NaN;
    summary.calibrated_direct_suppression_dB = NaN;
    summary.calibrated_clean_mask_rel = NaN;
    summary.calibrated_pdr_dB = NaN;
end

if summary.has_oracle
    summary.oracle_interference_suppression_dB = ...
        metrics_table.interference_suppression_dB(idx_oracle);
    summary.oracle_direct_suppression_dB = ...
        metrics_table.direct_suppression_dB(idx_oracle);
    summary.oracle_clean_mask_rel = ...
        metrics_table.clean_mask_rel(idx_oracle);
    summary.oracle_pdr_dB = ...
        metrics_table.target_pdr_dB(idx_oracle);
else
    summary.oracle_interference_suppression_dB = NaN;
    summary.oracle_direct_suppression_dB = NaN;
    summary.oracle_clean_mask_rel = NaN;
    summary.oracle_pdr_dB = NaN;
end

if summary.has_baseline && summary.has_calibrated
    summary.calibrated_vs_baseline_interference_gain_dB = ...
        metrics_table.interference_suppression_dB(idx_cal) - ...
        metrics_table.interference_suppression_dB(idx_base);

    summary.calibrated_vs_baseline_direct_gain_dB = ...
        metrics_table.direct_suppression_dB(idx_cal) - ...
        metrics_table.direct_suppression_dB(idx_base);

    summary.calibrated_vs_baseline_clean_mask_gain_dB = ...
        10*log10(metrics_table.clean_mask_rel(idx_base) / ...
        max(metrics_table.clean_mask_rel(idx_cal), eps));

    summary.calibrated_vs_baseline_pdr_gain_dB = ...
        metrics_table.target_pdr_dB(idx_cal) - ...
        metrics_table.target_pdr_dB(idx_base);
else
    summary.calibrated_vs_baseline_interference_gain_dB = NaN;
    summary.calibrated_vs_baseline_direct_gain_dB = NaN;
    summary.calibrated_vs_baseline_clean_mask_gain_dB = NaN;
    summary.calibrated_vs_baseline_pdr_gain_dB = NaN;
end

if summary.has_oracle && summary.has_calibrated
    summary.calibrated_oracle_gap_interference_dB = ...
        metrics_table.interference_suppression_dB(idx_oracle) - ...
        metrics_table.interference_suppression_dB(idx_cal);

    summary.calibrated_oracle_gap_clean_mask_dB = ...
        10*log10(metrics_table.clean_mask_rel(idx_cal) / ...
        max(metrics_table.clean_mask_rel(idx_oracle), eps));

    summary.calibrated_oracle_pdr_gap_dB = ...
        metrics_table.target_pdr_dB(idx_oracle) - ...
        metrics_table.target_pdr_dB(idx_cal);
else
    summary.calibrated_oracle_gap_interference_dB = NaN;
    summary.calibrated_oracle_gap_clean_mask_dB = NaN;
    summary.calibrated_oracle_pdr_gap_dB = NaN;
end

% --- 自适应对比算法摘要指标 ---

idx_mvdr = find(strcmp(method_list, 'MVDR'), 1);
idx_lms = find(strcmp(method_list, 'LMS'), 1);
idx_rls = find(strcmp(method_list, 'RLS'), 1);

summary.has_mvdr = ~isempty(idx_mvdr);
summary.has_lms = ~isempty(idx_lms);
summary.has_rls = ~isempty(idx_rls);

adaptive_fields = { ...
    'interference_suppression_dB', ...
    'direct_suppression_dB', ...
    'clean_mask_rel', ...
    'pdr_dB'};

adaptive_table_fields = { ...
    'interference_suppression_dB', ...
    'direct_suppression_dB', ...
    'clean_mask_rel', ...
    'target_pdr_dB'};

adaptive_methods = {'mvdr', 'lms', 'rls'};
adaptive_indices = {idx_mvdr, idx_lms, idx_rls};

for am = 1:numel(adaptive_methods)
    method_prefix = adaptive_methods{am};
    idx_val = adaptive_indices{am};

    if ~isempty(idx_val)
        summary.([method_prefix, '_interference_suppression_dB']) = ...
            metrics_table.interference_suppression_dB(idx_val);
        summary.([method_prefix, '_direct_suppression_dB']) = ...
            metrics_table.direct_suppression_dB(idx_val);
        summary.([method_prefix, '_clean_mask_rel']) = ...
            metrics_table.clean_mask_rel(idx_val);
        summary.([method_prefix, '_pdr_dB']) = ...
            metrics_table.target_pdr_dB(idx_val);
        summary.([method_prefix, '_num_iter']) = ...
            metrics_table.num_iter(idx_val);
    else
        for af = 1:numel(adaptive_table_fields)
            summary.([method_prefix, '_', adaptive_fields{af}]) = NaN;
        end
        summary.([method_prefix, '_num_iter']) = NaN;
    end

    % 自适应方法 vs Baseline gain
    if summary.has_baseline && ~isempty(idx_val)
        summary.([method_prefix, '_vs_baseline_interference_gain_dB']) = ...
            metrics_table.interference_suppression_dB(idx_val) - ...
            metrics_table.interference_suppression_dB(idx_base);
        summary.([method_prefix, '_vs_baseline_pdr_gain_dB']) = ...
            metrics_table.target_pdr_dB(idx_val) - ...
            metrics_table.target_pdr_dB(idx_base);
        summary.([method_prefix, '_vs_baseline_clean_mask_gain_dB']) = ...
            10*log10(metrics_table.clean_mask_rel(idx_base) / ...
            max(metrics_table.clean_mask_rel(idx_val), eps));
    else
        summary.([method_prefix, '_vs_baseline_interference_gain_dB']) = NaN;
        summary.([method_prefix, '_vs_baseline_pdr_gain_dB']) = NaN;
        summary.([method_prefix, '_vs_baseline_clean_mask_gain_dB']) = NaN;
    end

    % 自适应方法 vs Calibrated CLEAN gap
    if summary.has_calibrated && ~isempty(idx_val)
        summary.([method_prefix, '_vs_calibrated_pdr_gap_dB']) = ...
            metrics_table.target_pdr_dB(idx_cal) - ...
            metrics_table.target_pdr_dB(idx_val);
    else
        summary.([method_prefix, '_vs_calibrated_pdr_gap_dB']) = NaN;
    end
end

summary.phase_rmse_active_deg = phase_metrics.rmse_active_deg;
summary.phase_rmse_all_deg = phase_metrics.rmse_all_deg;
summary.phase_max_abs_deg = phase_metrics.max_abs_deg;
end
