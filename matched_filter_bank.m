function [Y, mf_data, info] = matched_filter_bank(cfg, X, mf, paths)
%MATCHED_FILTER_BANK 对阵列接收数据逐阵元进行匹配滤波。
%
% 用法:
%   [Y, mf_data, info] = matched_filter_bank(cfg, X, mf);
%   [Y, mf_data, info] = matched_filter_bank(cfg, X, mf, paths);
%
% 输入:
%   cfg    sonar_config() 与 generate_virtual_sources() 更新后的配置结构体。
%   X      M x N_rx 阵列接收数据矩阵。
%   mf     N_tx x 1 匹配滤波器，通常来自 generate_lfm(cfg)。
%          当前工程中 mf 已按信号能量归一化，使单路径峰值近似等于路径复幅度。
%   paths  可选，路径结构体数组，用于标注各路径理论时延位置。
%
% 输出:
%   Y        M x N_mf 匹配滤波输出矩阵。
%   mf_data  匹配滤波结果结构体。
%   info     匹配滤波信息结构体。
%
% 时间轴约定:
%   X 的列对应 cfg.t_rx。
%   使用 conv(x, mf, 'full') 后，输出长度为:
%       N_mf = N_rx + N_tx - 1
%
%   full 卷积会引入 N_tx-1 个采样点的峰值偏移。
%   本模块定义:
%       t_mf(k) = t_rx_start + (k-1 - (N_tx-1))/fs
%   使匹配滤波峰值位置直接对应物理传播时延 tau。

%% 1. 输入检查
if nargin < 4
    paths = [];
end

required_fields = {'M', 'fs', 'N_rx', 't_rx', 'N_tx', 't_tx'};

for k = 1:numel(required_fields)
    field_name = required_fields{k};
    assert(isfield(cfg, field_name), 'cfg.%s 是必需字段。', field_name);
end

M = cfg.M;
fs = cfg.fs;
N_rx = cfg.N_rx;
N_tx = cfg.N_tx;

assert(size(X, 1) == M, 'X 的行数必须等于 cfg.M。');
assert(size(X, 2) == N_rx, 'X 的列数必须等于 cfg.N_rx。');

mf = mf(:);
assert(numel(mf) == N_tx, 'mf 的长度必须等于 cfg.N_tx。');
assert(all(isfinite(X(:))), 'X 中存在非有限数。');
assert(all(isfinite(mf(:))), 'mf 中存在非有限数。');

assert(numel(cfg.t_rx) == N_rx, 'cfg.t_rx 的长度必须等于 cfg.N_rx。');
assert(numel(cfg.t_tx) == N_tx, 'cfg.t_tx 的长度必须等于 cfg.N_tx。');

dt_rx = diff(cfg.t_rx(:));
if ~isempty(dt_rx)
    assert(max(abs(dt_rx - 1/fs)) < 1e-9/fs, ...
        'cfg.t_rx 必须按 1/fs 均匀采样。');
end

%% 2. 逐阵元匹配滤波
N_mf = N_rx + N_tx - 1;
Y = complex(zeros(M, N_mf));

for m = 1:M
    y_m = conv(X(m, :).', mf, 'full');
    Y(m, :) = y_m.';
end

%% 3. 构造匹配滤波时延轴
sample_axis = (0:N_mf-1).';
mf_group_delay_samples = N_tx - 1;

t_rx_start = cfg.t_rx(1);
t_mf = t_rx_start + (sample_axis - mf_group_delay_samples) / fs;

%% 4. 阵列非相干时延能量谱
% 这里只做逐阵元匹配滤波后的非相干能量叠加，不做波束形成。
% 后续 CLEAN / DOA 搜索会再利用 steering_vector 进行空间匹配。
power_incoh = sum(abs(Y).^2, 1).';

power_incoh_norm = power_incoh / (max(power_incoh) + eps);
power_incoh_dB = 10*log10(power_incoh_norm + eps);

[global_peak_power, global_peak_idx] = max(power_incoh);
global_peak_delay = t_mf(global_peak_idx);

%% 5. 路径理论时延标注
path_info = struct();

if ~isempty(paths)
    P = numel(paths);

    path_delay = [paths.delay].';
    path_theta = [paths.theta].';
    path_theta_deg = [paths.theta_deg].';
    path_type = {paths.type}.';
    path_name = {paths.name}.';
    path_amp = [paths.amp].';
    path_power = [paths.power].';

    path_nearest_idx = zeros(P, 1);
    path_nearest_delay = zeros(P, 1);
    path_delay_error = zeros(P, 1);
    path_incoh_power = zeros(P, 1);
    path_sensor1_abs = zeros(P, 1);

    for p = 1:P
        [~, idx] = min(abs(t_mf - path_delay(p)));

        path_nearest_idx(p) = idx;
        path_nearest_delay(p) = t_mf(idx);
        path_delay_error(p) = t_mf(idx) - path_delay(p);
        path_incoh_power(p) = power_incoh(idx);

        % 注意: 这是“总匹配滤波输出”在该路径理论时延附近的
        % 第 1 阵元幅值，不是单独路径分量的幅值。
        path_sensor1_abs(p) = abs(Y(1, idx));
    end

    path_info.delay = path_delay;
    path_info.theta = path_theta;
    path_info.theta_deg = path_theta_deg;
    path_info.type = path_type;
    path_info.name = path_name;
    path_info.amp = path_amp;
    path_info.power = path_power;
    path_info.nearest_idx = path_nearest_idx;
    path_info.nearest_delay = path_nearest_delay;
    path_info.delay_error = path_delay_error;
    path_info.incoh_power = path_incoh_power;
    path_info.sensor1_abs = path_sensor1_abs;
    path_info.sensor1_abs_definition = ...
        'abs(Y(1, nearest_idx)); total matched-filter output, not isolated path amplitude';
end

%% 6. 输出结构体
mf_data = struct();

mf_data.Y = Y;
mf_data.t_mf = t_mf;
mf_data.delay_axis = t_mf;

mf_data.power_incoh = power_incoh;
mf_data.power_incoh_norm = power_incoh_norm;
mf_data.power_incoh_dB = power_incoh_dB;

mf_data.global_peak_idx = global_peak_idx;
mf_data.global_peak_delay = global_peak_delay;
mf_data.global_peak_power = global_peak_power;

mf_data.path_info = path_info;

%% 7. 信息结构体
info = struct();

info.model = 'per_sensor_matched_filter_bank';
info.operation = 'Y_m = conv(X_m, mf, full)';
info.delay_axis_definition = ...
    't_mf is absolute propagation delay after removing matched-filter group delay';

info.M = M;
info.N_rx = N_rx;
info.N_tx = N_tx;
info.N_mf = N_mf;
info.fs = fs;

info.t_rx_start = cfg.t_rx(1);
info.t_rx_end = cfg.t_rx(end);
info.t_mf_start = t_mf(1);
info.t_mf_end = t_mf(end);

info.mf_length = N_tx;
info.mf_group_delay_samples = mf_group_delay_samples;
info.mf_group_delay_time = mf_group_delay_samples / fs;

info.mf_energy = sum(abs(mf).^2);
info.mf_peak_abs = max(abs(mf));

% 若 mf 来自 generate_lfm，且 s 为单位平均功率、mf = conj(flip(s))/Es，
% 则 sum(abs(mf).^2) 约为 1/Es。
% 复白噪声经过匹配滤波后的理论方差为:
%   noise_power_after_mf = noise_power * sum(abs(mf).^2)
if isfield(cfg, 'noise_power')
    info.noise_power_before_mf = cfg.noise_power;
    info.noise_power_after_mf_theory = cfg.noise_power * info.mf_energy;
else
    info.noise_power_before_mf = NaN;
    info.noise_power_after_mf_theory = NaN;
end

info.global_peak_idx = global_peak_idx;
info.global_peak_delay = global_peak_delay;
info.global_peak_power = global_peak_power;

if ~isempty(paths)
    info.num_paths = numel(paths);
else
    info.num_paths = 0;
end
end
