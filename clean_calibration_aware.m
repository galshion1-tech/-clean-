function [Y_res, clean_data, info] = clean_calibration_aware(cfg, Y, mf_data, s, mf, paths, phi_hat)
%CLEAN_CALIBRATION_AWARE 校准感知 CLEAN，在匹配滤波域抑制直达波/强多径。
%
% 用法:
%   [Y_res, clean_data, info] = ...
%       clean_calibration_aware(cfg, Y, mf_data, s, mf, paths, phi_hat);
%
%   若不使用 paths，需要显式传入 []:
%       clean_calibration_aware(cfg, Y, mf_data, s, mf, [], phi_hat);
%
% 输入:
%   cfg      sonar_config() 配置结构体。
%   Y        M x N_mf 匹配滤波输出矩阵，来自 matched_filter_bank()。
%   mf_data  matched_filter_bank() 输出结构体，至少包含 t_mf。
%   s        N_tx x 1 发射 LFM 信号，来自 generate_lfm()。
%   mf       N_tx x 1 匹配滤波器，来自 generate_lfm()。
%   paths    generate_virtual_sources() 输出的路径结构体数组。
%            若提供，则只在直达波/强多径干扰窗口内搜索。
%   phi_hat  M x 1 估计阵列通道相位误差，rad。
%
% 输出:
%   Y_res       校准感知 CLEAN 后的匹配滤波残差。
%   clean_data  CLEAN 结果结构体。
%   info        算法信息结构体。
%
% 物理模型:
%   匹配滤波输出中的单条路径近似为:
%       Y_p(tau) = alpha_p * a_true(theta_p) * g(tau - tau_p)
%
%   传统 CLEAN 使用:
%       a_nom(theta)
%
%   本方法使用校准后的导向矢量:
%       a_cal(theta) = diag(exp(1j*phi_hat)) * a_nom(theta)
%
%   因此 replica 更接近真实失配阵列产生的直达波/多径结构。

%% 1. 输入检查
if nargin < 6
    paths = [];
end

assert(nargin >= 7 && ~isempty(phi_hat), ...
    'clean_calibration_aware 需要输入 phi_hat。');

required_fields = { ...
    'M', 'fs', 'N_tx', ...
    'max_iter', 'loop_gain', 'residual_threshold_rel', ...
    'angle_grid', ...
    'direct_search_half_window', ...
    'multipath_search_half_window', ...
    'array_pos', 'array_center', 'array_axis', 'broadside_axis', ...
    'lambda'};

for k = 1:numel(required_fields)
    field_name = required_fields{k};
    assert(isfield(cfg, field_name), 'cfg.%s 是必需字段。', field_name);
end

assert(isfield(mf_data, 't_mf'), 'mf_data.t_mf 是必需字段。');

M = cfg.M;
fs = cfg.fs;
N_tx = cfg.N_tx;

assert(size(Y, 1) == M, 'Y 的行数必须等于 cfg.M。');
assert(all(isfinite(Y(:))), 'Y 中存在非有限数。');

t_mf = mf_data.t_mf(:);
N_mf = numel(t_mf);

assert(size(Y, 2) == N_mf, 'Y 的列数必须等于 numel(mf_data.t_mf)。');

s = s(:);
mf = mf(:);

assert(numel(s) == N_tx, 's 的长度必须等于 cfg.N_tx。');
assert(numel(mf) == N_tx, 'mf 的长度必须等于 cfg.N_tx。');
assert(all(isfinite(s)), 's 中存在非有限数。');
assert(all(isfinite(mf)), 'mf 中存在非有限数。');

phi_hat = phi_hat(:);
assert(numel(phi_hat) == M, 'phi_hat 的长度必须等于 cfg.M。');
assert(isreal(phi_hat) && all(isfinite(phi_hat)), ...
    'phi_hat 必须是有限实数。');

assert(cfg.max_iter > 0 && mod(cfg.max_iter, 1) == 0, ...
    'cfg.max_iter 必须是正整数。');
assert(cfg.loop_gain > 0 && cfg.loop_gain <= 1, ...
    'cfg.loop_gain 必须位于 (0, 1]。');
assert(cfg.residual_threshold_rel > 0, ...
    'cfg.residual_threshold_rel 必须为正数。');

angle_grid = cfg.angle_grid(:).';
assert(~isempty(angle_grid), 'cfg.angle_grid 不能为空。');
assert(isreal(angle_grid) && all(isfinite(angle_grid)), ...
    'cfg.angle_grid 必须是有限实数。');

%% 2. 构造匹配滤波点扩散函数 g(tau)
% 与 baseline CLEAN 保持相同物理模型:
%   g = conv(s, mf, 'full')
%   g 的峰值位置定义为 tau = 0。
pc_kernel = conv(s, mf, 'full');
pc_kernel = pc_kernel(:);

L_pc = numel(pc_kernel);
pc_sample_axis = (0:L_pc-1).';
pc_group_delay_samples = N_tx - 1;
pc_delay_axis = (pc_sample_axis - pc_group_delay_samples) / fs;

[~, pc_peak_idx] = min(abs(pc_delay_axis));
pc_peak_value = pc_kernel(pc_peak_idx);

assert(abs(pc_peak_value) > eps, ...
    '匹配滤波点扩散函数峰值过小，请检查 s 和 mf。');

%% 3. 构造搜索时延索引
[search_idx, search_mask, search_info] = build_clean_search_indices(cfg, t_mf, paths);

assert(~isempty(search_idx), ...
    'CLEAN 搜索索引为空，请检查路径时延或搜索窗口设置。');

%% 4. 预计算校准后的导向矢量角度字典
% 这是本方法区别于 baseline CLEAN 的核心:
%   baseline: A_grid = steering_vector(angle_grid, cfg)
%   proposed: A_grid = steering_vector(angle_grid, cfg, phi_hat)
A_grid = steering_vector(angle_grid, cfg, phi_hat);
K_angle = numel(angle_grid);

% steering_vector 每个阵元幅度为 1，因此每列范数平方约为 M。
A_norm2 = sum(abs(A_grid).^2, 1).';
A_norm2 = max(A_norm2, eps);

%% 5. 初始化 CLEAN
R = Y;
Y_model = complex(zeros(size(Y)));

max_iter = cfg.max_iter;
loop_gain = cfg.loop_gain;

initial_mask_energy = energy_sum(R(:, search_idx));
initial_total_energy = energy_sum(R);

if initial_mask_energy <= eps
    warning('初始搜索区域能量接近 0，校准感知 CLEAN 不执行迭代。');
    Y_res = R;

    clean_data = make_empty_clean_data(Y_model, R, search_idx, search_mask, ...
        search_info, initial_mask_energy, initial_total_energy);
    clean_data.phi_hat = phi_hat;
    clean_data.phi_hat_deg = phi_hat * 180/pi;

    info = make_info(cfg, t_mf, pc_delay_axis, pc_kernel, pc_peak_value, ...
        0, initial_mask_energy, initial_total_energy, ...
        'zero_initial_search_energy', phi_hat);

    info.final_mask_energy = clean_data.final_mask_energy;
    info.final_mask_energy_rel = clean_data.final_mask_energy_rel;
    info.final_total_energy = clean_data.final_total_energy;
    info.final_total_energy_rel = clean_data.final_total_energy_rel;
    info.search_delay_min = min(t_mf(search_idx));
    info.search_delay_max = max(t_mf(search_idx));
    info.num_search_delay_bins = numel(search_idx);
    info.num_angle_bins = K_angle;
    return;
end

iter_delay = zeros(max_iter, 1);
iter_theta = zeros(max_iter, 1);
iter_theta_deg = zeros(max_iter, 1);
iter_alpha_hat = complex(zeros(max_iter, 1));
iter_alpha_sub = complex(zeros(max_iter, 1));
iter_peak_idx = zeros(max_iter, 1);
iter_angle_idx = zeros(max_iter, 1);
iter_score = zeros(max_iter, 1);
iter_mask_energy = zeros(max_iter, 1);
iter_mask_energy_rel = zeros(max_iter, 1);
iter_total_energy = zeros(max_iter, 1);
iter_total_energy_rel = zeros(max_iter, 1);

stop_reason = 'max_iter_reached';
num_iter_done = 0;

%% 6. 校准感知 CLEAN 迭代
for it = 1:max_iter
    current_mask_energy = energy_sum(R(:, search_idx));
    current_total_energy = energy_sum(R);

    rel_mask_energy = current_mask_energy / max(initial_mask_energy, eps);
    rel_total_energy = current_total_energy / max(initial_total_energy, eps);

    if rel_mask_energy < cfg.residual_threshold_rel
        stop_reason = 'residual_threshold_reached';
        break;
    end

    % 6.1 在搜索区域内做 delay-angle 匹配搜索。
    R_search = R(:, search_idx);

    % alpha_hat(theta, tau) = a_cal^H(theta)y(tau) / (||a_cal||^2 g(0))
    alpha_grid = bsxfun(@rdivide, ...
        A_grid' * R_search, ...
        A_norm2 * pc_peak_value);

    score_grid = abs(alpha_grid).^2;

    [best_score, linear_idx] = max(score_grid(:));

    if ~isfinite(best_score) || best_score <= eps
        stop_reason = 'no_valid_peak';
        break;
    end

    [best_angle_idx, best_delay_local_idx] = ind2sub(size(score_grid), linear_idx);

    best_delay_idx = search_idx(best_delay_local_idx);
    tau_hat = t_mf(best_delay_idx);
    theta_hat = angle_grid(best_angle_idx);
    alpha_hat = alpha_grid(best_angle_idx, best_delay_local_idx);

    % 6.2 构造校准后的匹配滤波域 replica。
    % g_shifted(k) = g(t_mf(k) - tau_hat)
    g_shifted = interp1(pc_delay_axis, pc_kernel, ...
        t_mf - tau_hat, 'linear', 0);

    a_cal = A_grid(:, best_angle_idx);
    alpha_sub = loop_gain * alpha_hat;

    Y_replica = a_cal * (alpha_sub * g_shifted.');

    % 6.3 从残差中剥离校准 replica。
    R = R - Y_replica;
    Y_model = Y_model + Y_replica;

    % 6.4 记录迭代信息。
    num_iter_done = it;

    iter_delay(it) = tau_hat;
    iter_theta(it) = theta_hat;
    iter_theta_deg(it) = theta_hat * 180/pi;
    iter_alpha_hat(it) = alpha_hat;
    iter_alpha_sub(it) = alpha_sub;
    iter_peak_idx(it) = best_delay_idx;
    iter_angle_idx(it) = best_angle_idx;
    iter_score(it) = best_score;
    iter_mask_energy(it) = current_mask_energy;
    iter_mask_energy_rel(it) = rel_mask_energy;
    iter_total_energy(it) = current_total_energy;
    iter_total_energy_rel(it) = rel_total_energy;
end

Y_res = R;

%% 7. 截断迭代记录
valid_idx = 1:num_iter_done;

iter_table = table();
if num_iter_done > 0
    iter_table.iter = valid_idx(:);
    iter_table.delay = iter_delay(valid_idx);
    iter_table.theta = iter_theta(valid_idx);
    iter_table.theta_deg = iter_theta_deg(valid_idx);
    iter_table.alpha_hat = iter_alpha_hat(valid_idx);
    iter_table.alpha_subtracted = iter_alpha_sub(valid_idx);
    iter_table.peak_idx = iter_peak_idx(valid_idx);
    iter_table.angle_idx = iter_angle_idx(valid_idx);
    iter_table.score = iter_score(valid_idx);
    iter_table.mask_energy = iter_mask_energy(valid_idx);
    iter_table.mask_energy_rel = iter_mask_energy_rel(valid_idx);
    iter_table.total_energy = iter_total_energy(valid_idx);
    iter_table.total_energy_rel = iter_total_energy_rel(valid_idx);
end

%% 8. 输出结构体
clean_data = struct();

clean_data.Y_original = Y;
clean_data.Y_residual = Y_res;
clean_data.Y_model = Y_model;

clean_data.t_mf = t_mf;
clean_data.search_idx = search_idx;
clean_data.search_mask = search_mask;
clean_data.search_info = search_info;

clean_data.pc_kernel = pc_kernel;
clean_data.pc_delay_axis = pc_delay_axis;
clean_data.pc_peak_value = pc_peak_value;

clean_data.angle_grid = angle_grid;
clean_data.iter_table = iter_table;

clean_data.num_iter = num_iter_done;
clean_data.stop_reason = stop_reason;

clean_data.phi_hat = phi_hat;
clean_data.phi_hat_deg = phi_hat * 180/pi;

clean_data.final_mask_energy = energy_sum(Y_res(:, search_idx));
clean_data.final_mask_energy_rel = clean_data.final_mask_energy / max(initial_mask_energy, eps);
clean_data.final_total_energy = energy_sum(Y_res);
clean_data.final_total_energy_rel = clean_data.final_total_energy / max(initial_total_energy, eps);

clean_data.reconstruction_energy = energy_sum(Y_model);

%% 9. 信息结构体
info = make_info(cfg, t_mf, pc_delay_axis, pc_kernel, pc_peak_value, ...
    num_iter_done, initial_mask_energy, initial_total_energy, ...
    stop_reason, phi_hat);

info.final_mask_energy = clean_data.final_mask_energy;
info.final_mask_energy_rel = clean_data.final_mask_energy_rel;
info.final_total_energy = clean_data.final_total_energy;
info.final_total_energy_rel = clean_data.final_total_energy_rel;

info.search_delay_min = min(t_mf(search_idx));
info.search_delay_max = max(t_mf(search_idx));
info.num_search_delay_bins = numel(search_idx);
info.num_angle_bins = K_angle;
end

%% ========================================================================
% 局部辅助函数
% ========================================================================

function [search_idx, search_mask, search_info] = build_clean_search_indices(cfg, t_mf, paths)
%BUILD_CLEAN_SEARCH_INDICES 构造 CLEAN 的时延搜索区域。
%
% 若提供 paths，则只在干扰路径附近搜索:
%   direct 使用 cfg.direct_search_half_window；
%   surface/bottom/multipath 使用 cfg.multipath_search_half_window。
%
% 这样可以避免校准感知 CLEAN 主动剥离目标路径。
% 若未提供 paths，则默认使用整个 t_mf 轴搜索。

N_mf = numel(t_mf);
search_mask = false(N_mf, 1);

search_info = struct();
search_info.mode = '';
search_info.windows = [];

if nargin < 3 || isempty(paths)
    search_mask(:) = true;
    search_info.mode = 'full_delay_axis_without_paths';
else
    search_info.mode = 'interference_windows_from_paths';

    win_list = struct('type', {}, 'delay', {}, 'half_window', {}, ...
        'idx_start', {}, 'idx_end', {});

    for p = 1:numel(paths)
        if ~isfield(paths(p), 'is_interference') || ~paths(p).is_interference
            continue;
        end

        tau0 = paths(p).delay;

        if strcmp(paths(p).type, 'direct')
            half_window = cfg.direct_search_half_window;
        else
            half_window = cfg.multipath_search_half_window;
        end

        local_mask = abs(t_mf - tau0) <= half_window;
        search_mask = search_mask | local_mask;

        idx_local = find(local_mask);
        if ~isempty(idx_local)
            w = struct();
            w.type = paths(p).type;
            w.delay = tau0;
            w.half_window = half_window;
            w.idx_start = idx_local(1);
            w.idx_end = idx_local(end);
            win_list(end+1) = w; %#ok<AGROW>
        end
    end

    search_info.windows = win_list;

    if ~any(search_mask)
        warning('未能根据 paths 构造干扰搜索窗口，退化为全时延搜索。');
        search_mask(:) = true;
        search_info.mode = 'fallback_full_delay_axis';
    end
end

search_idx = find(search_mask);
end

function clean_data = make_empty_clean_data(Y_model, R, search_idx, search_mask, ...
    search_info, initial_mask_energy, initial_total_energy)
%MAKE_EMPTY_CLEAN_DATA 构造空 CLEAN 输出。

clean_data = struct();
clean_data.Y_original = R;
clean_data.Y_residual = R;
clean_data.Y_model = Y_model;
clean_data.search_idx = search_idx;
clean_data.search_mask = search_mask;
clean_data.search_info = search_info;
clean_data.iter_table = table();
clean_data.num_iter = 0;
clean_data.stop_reason = 'zero_initial_search_energy';

clean_data.final_mask_energy = energy_sum(R(:, search_idx));
clean_data.final_mask_energy_rel = clean_data.final_mask_energy / max(initial_mask_energy, eps);
clean_data.final_total_energy = energy_sum(R);
clean_data.final_total_energy_rel = clean_data.final_total_energy / max(initial_total_energy, eps);
clean_data.reconstruction_energy = energy_sum(Y_model);
end

function info = make_info(cfg, t_mf, pc_delay_axis, pc_kernel, pc_peak_value, ...
    num_iter_done, initial_mask_energy, initial_total_energy, stop_reason, phi_hat)
%MAKE_INFO 构造校准感知 CLEAN 信息结构体。

info = struct();

info.model = 'calibration_aware_clean_in_matched_filter_domain';
info.algorithm = 'CLEAN_with_estimated_phase_calibrated_steering_vector';
info.uses_phase_calibration = true;

info.physical_model = ...
    'Y_p(tau)=alpha_p*a_cal(theta_p)*g(tau-tau_p)';
info.important_note = ...
    'calibration-aware CLEAN uses phi_hat to update steering vector; phi_hat estimates equivalent channel phase error, not element position';

info.num_iter = num_iter_done;
info.max_iter = cfg.max_iter;
info.loop_gain = cfg.loop_gain;
info.residual_threshold_rel = cfg.residual_threshold_rel;
info.stop_reason = stop_reason;

info.initial_mask_energy = initial_mask_energy;
info.initial_total_energy = initial_total_energy;

info.t_mf_start = t_mf(1);
info.t_mf_end = t_mf(end);

info.pc_delay_start = pc_delay_axis(1);
info.pc_delay_end = pc_delay_axis(end);
info.pc_peak_value = pc_peak_value;
info.pc_peak_abs = abs(pc_peak_value);
info.pc_peak_phase = angle(pc_peak_value);
info.pc_kernel_energy = energy_sum(pc_kernel);

info.phi_hat = phi_hat;
info.phi_hat_deg = phi_hat * 180/pi;
info.phi_hat_rms_deg = sqrt(mean(phi_hat.^2)) * 180/pi;
info.phi_hat_max_abs_deg = max(abs(phi_hat)) * 180/pi;
end

function e = energy_sum(X)
%ENERGY_SUM 计算复矩阵总能量，兼容旧版 MATLAB。

e = sum(abs(X(:)).^2);
end
