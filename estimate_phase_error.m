function [phi_hat, calib_data, info] = estimate_phase_error(cfg, Y, mf_data, s, mf, paths)
%ESTIMATE_PHASE_ERROR 利用直达波局部匹配滤波结构估计阵列相位误差。
%
% 用法:
%   [phi_hat, calib_data, info] = estimate_phase_error(cfg, Y, mf_data, s, mf, paths);
%
% 输入:
%   cfg      sonar_config() 配置结构体。
%   Y        M x N_mf 匹配滤波输出矩阵，来自 matched_filter_bank()。
%   mf_data  matched_filter_bank() 输出结构体，至少包含 t_mf。
%   s        N_tx x 1 发射 LFM 信号。
%   mf       N_tx x 1 匹配滤波器。
%   paths    generate_virtual_sources() 输出的路径结构体数组。
%
% 输出:
%   phi_hat     M x 1 估计阵列通道相位误差，rad。
%               第 cfg.ref_sensor 个阵元固定为 0。
%   calib_data  校准中间量结构体。
%   info        校准信息结构体。
%
% 定位说明:
%   本模块估计的是等效阵列通道相位误差，不是阵元物理位置扰动。
%   本模块利用 paths 中的直达波理论时延和 DOA，属于几何辅助的
%   direct-path-aided phase calibration，不是完全盲估计。
%
% 物理模型:
%   在直达波附近:
%       Y(tau) ≈ alpha_d * diag(exp(1j*phi)) * a0(theta_d) * g(tau - tau_d)
%
%   先通过 g(tau - tau_d) 做局部最小二乘，估计:
%       b_hat ≈ alpha_d * diag(exp(1j*phi)) * a0(theta_d)
%
%   再除以名义导向矢量:
%       z = b_hat ./ a0(theta_d) ≈ alpha_d * exp(1j*phi)
%
%   最后用参考阵元消去共同复幅度 alpha_d:
%       phi_hat(m) = angle(z_m / z_ref)

%% 1. 输入检查
required_fields = { ...
    'M', 'fs', 'B', 'N_tx', 'ref_sensor', ...
    'phase_est_clip', ...
    'direct_window', ...
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
ref_sensor = cfg.ref_sensor;

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
assert(ref_sensor >= 1 && ref_sensor <= M, ...
    'cfg.ref_sensor 必须位于 [1, M]。');
assert(~isempty(paths), 'paths 不能为空。当前相位估计模块需要 direct path 信息。');

%% 2. 找到直达波路径
path_types = {paths.type};
direct_idx = find(strcmp(path_types, 'direct'), 1);

assert(~isempty(direct_idx), ...
    'paths 中未找到 type = ''direct'' 的直达波路径。');

tau_d = paths(direct_idx).delay;
theta_d = paths(direct_idx).theta;
theta_d_deg = paths(direct_idx).theta_deg;

%% 3. 构造匹配滤波点扩散函数 g(tau)
% 与 clean_baseline.m 保持一致:
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

%% 4. 构造直达波局部估计窗口
% 默认使用 cfg.direct_window；若配置了 cfg.calib_half_window，则优先使用。
% 当直达波与其他路径很近时，自适应缩小窗口，避免直接覆盖相邻路径峰值。
half_window_nominal = cfg.direct_window;

if isfield(cfg, 'calib_half_window')
    half_window_nominal = cfg.calib_half_window;
end

half_window = half_window_nominal;

other_delays = [paths.delay].';
other_delays(direct_idx) = [];

if ~isempty(other_delays)
    min_sep = min(abs(other_delays - tau_d));

    if isfinite(min_sep) && min_sep > 0
        half_window = min(half_window, 0.45 * min_sep);
    end
end

% 至少保留若干采样点，避免窗口过窄导致局部 LS 不稳定。
min_half_window = max(5 / fs, 0.5 / max(get_bandwidth(cfg), eps));
half_window = max(half_window, min_half_window);

window_idx = find(abs(t_mf - tau_d) <= half_window);

assert(~isempty(window_idx), ...
    '直达波校准窗口为空，请检查 direct_window 或 t_mf。');

t_win = t_mf(window_idx);
Y_win = Y(:, window_idx);

%% 5. 构造直达波局部压缩脉冲模板
% g_win(k) = g(t_win(k) - tau_d)
g_win = interp1(pc_delay_axis, pc_kernel, ...
    t_win - tau_d, 'linear', 0);

g_win = g_win(:);

denom = sum(abs(g_win).^2);
assert(denom > eps, ...
    '直达波窗口内的匹配滤波模板能量过小。');

%% 6. 局部最小二乘估计阵列复向量 b_hat
% 若 Y_win ≈ b * g_win^T，则:
%   b_hat = Y_win * conj(g_win) / ||g_win||^2
b_hat = Y_win * conj(g_win) / denom;

%% 7. 除掉名义导向矢量，得到 z ≈ alpha_d * exp(j*phi)
a0 = steering_vector(theta_d, cfg);
a0 = a0(:);

assert(numel(a0) == M, '名义导向矢量维度错误。');

z = b_hat ./ safe_divisor(a0);

z_ref = z(ref_sensor);
assert(abs(z_ref) > eps, ...
    '参考阵元 z_ref 幅度过小，无法进行相对相位校准。');

%% 8. 相对相位误差估计
phi_hat = angle(z ./ z_ref);
phi_hat = wrap_to_pi_local(phi_hat);

% 参考阵元固定为 0。
phi_hat(ref_sensor) = 0;

%% 9. 估计值限幅
if isfield(cfg, 'phase_est_clip') && isfinite(cfg.phase_est_clip)
    phi_hat = min(max(phi_hat, -cfg.phase_est_clip), cfg.phase_est_clip);
    phi_hat(ref_sensor) = 0;
end

%% 10. 输出中间量
calib_data = struct();

calib_data.method = 'direct_path_local_ls_relative_phase';
calib_data.direct_path_index = direct_idx;
calib_data.tau_d = tau_d;
calib_data.theta_d = theta_d;
calib_data.theta_d_deg = theta_d_deg;

calib_data.window_idx = window_idx;
calib_data.window_delay = t_win;
calib_data.half_window = half_window;
calib_data.half_window_nominal = half_window_nominal;

calib_data.pc_kernel = pc_kernel;
calib_data.pc_delay_axis = pc_delay_axis;
calib_data.pc_peak_value = pc_peak_value;

calib_data.g_win = g_win;
calib_data.template_energy = denom;

calib_data.b_hat = b_hat;
calib_data.a0 = a0;
calib_data.z = z;
calib_data.z_ref = z_ref;

calib_data.phi_hat = phi_hat;
calib_data.phi_hat_deg = phi_hat * 180/pi;

calib_data.ref_sensor = ref_sensor;

% 用估计的 phi_hat 重构直达波局部阵列向量，便于调试校准效果。
a_cal = steering_vector(theta_d, cfg, phi_hat);
calib_data.a_cal_direct = a_cal;
calib_data.alpha_common_hat = z_ref;

b_nom_fit = z_ref * a0;
b_cal_fit = z_ref * a_cal;

calib_data.fit_error_nominal = norm(b_hat - b_nom_fit) / max(norm(b_hat), eps);
calib_data.fit_error_calibrated = norm(b_hat - b_cal_fit) / max(norm(b_hat), eps);

%% 11. 输出信息
active_idx = setdiff((1:M).', ref_sensor);

info = struct();

info.model = 'direct_path_phase_error_estimation';
info.method = calib_data.method;
info.note = 'estimated equivalent channel phase error, not physical element-position perturbation';
info.prior = 'direct path delay and DOA are provided by geometry/path model; not blind estimation';

info.M = M;
info.ref_sensor = ref_sensor;

info.tau_d = tau_d;
info.theta_d = theta_d;
info.theta_d_deg = theta_d_deg;

info.window_num_samples = numel(window_idx);
info.window_half_width = half_window;
info.window_start = t_win(1);
info.window_end = t_win(end);

info.phase_est_clip = cfg.phase_est_clip;
info.phi_hat = phi_hat;
info.phi_hat_deg = phi_hat * 180/pi;

info.phi_hat_rms_deg = sqrt(mean(phi_hat.^2)) * 180/pi;
info.phi_hat_active_rms_deg = sqrt(mean(phi_hat(active_idx).^2)) * 180/pi;
info.phi_hat_max_abs_deg = max(abs(phi_hat)) * 180/pi;

info.fit_error_nominal = calib_data.fit_error_nominal;
info.fit_error_calibrated = calib_data.fit_error_calibrated;
end

%% ========================================================================
% 局部辅助函数
% ========================================================================

function y = wrap_to_pi_local(x)
%WRAP_TO_PI_LOCAL 将相位包裹到 [-pi, pi)。

y = mod(x + pi, 2*pi) - pi;
end

function y = safe_divisor(x)
%SAFE_DIVISOR 为复数逐元素除法提供安全分母。
%
% 对于正常导向矢量，abs(x) = 1，因此不会改变结果。
% 该函数仅防止极端异常输入导致除零。

mag = abs(x);
y = x;

bad_idx = mag < eps;
if any(bad_idx)
    y(bad_idx) = eps;
end
end

function B = get_bandwidth(cfg)
%GET_BANDWIDTH 读取带宽参数。

if isfield(cfg, 'B')
    B = cfg.B;
else
    B = NaN;
end
end
