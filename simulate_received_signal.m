function [X, components, info] = simulate_received_signal(cfg, s, paths, phi_true, mc_index)
%SIMULATE_RECEIVED_SIGNAL 生成阵列接收信号。
%
% 用法:
%   [X, components, info] = simulate_received_signal(cfg, s, paths);
%   [X, components, info] = simulate_received_signal(cfg, s, paths, phi_true);
%   [X, components, info] = simulate_received_signal(cfg, s, paths, phi_true, mc_index);
%
% 输入:
%   cfg       sonar_config() 与 generate_virtual_sources() 更新后的配置结构体。
%   s         发射复基带 LFM 信号，N_tx x 1。
%   paths     generate_virtual_sources() 输出的路径结构体数组。
%   phi_true  M x 1 阵列通道相位误差，rad。可选，默认全 0。
%   mc_index  Monte Carlo 序号，可选，默认 1，用于噪声随机种子。
%
% 输出:
%   X           M x N_rx 阵列接收数据矩阵。
%   components  信号分量结构体，包含无噪声信号、噪声和每条路径贡献。
%   info        接收信号仿真信息结构体。
%
% 信号模型:
%   X(t) = sum_p alpha_p * a_true(theta_p) * s(t - tau_p) + n(t)
%
% 重要约定:
%   paths(p).alpha 已包含阵列中心处的载频传播相位
%       exp(-1j*2*pi*fc*tau_p)
%   因此本模块不得再次乘 exp(-1j*2*pi*fc*tau_p)，否则相位会重复计算。

%% 1. 输入检查
if nargin < 4 || isempty(phi_true)
    phi_true = [];
end

if nargin < 5 || isempty(mc_index)
    mc_index = 1;
end

required_fields = { ...
    'M', 'fs', 'N_tx', 't_tx', ...
    'N_rx', 't_rx', ...
    'signal_model', ...
    'noise_power', ...
    'random_seed', ...
    'array_pos', 'array_center', 'array_axis', 'broadside_axis', ...
    'lambda'};

for k = 1:numel(required_fields)
    field_name = required_fields{k};
    assert(isfield(cfg, field_name), 'cfg.%s 是必需字段。', field_name);
end

assert(strcmp(cfg.signal_model, 'complex_baseband'), ...
    ['simulate_received_signal 当前仅支持 complex_baseband。', ...
     'paths.alpha 已经包含载频传播相位。']);

M = cfg.M;
N_tx = cfg.N_tx;
N_rx = cfg.N_rx;

s = s(:);
assert(numel(s) == N_tx, '发射信号 s 的长度必须等于 cfg.N_tx。');
assert(numel(cfg.t_tx) == N_tx, 'cfg.t_tx 的长度必须等于 cfg.N_tx。');
assert(numel(cfg.t_rx) == N_rx, 'cfg.t_rx 的长度必须等于 cfg.N_rx。');

assert(~isempty(paths), 'paths 不能为空。');
assert(cfg.noise_power >= 0, 'cfg.noise_power 必须非负。');
assert(mc_index >= 1 && mod(mc_index, 1) == 0, ...
    'mc_index 必须是正整数。');

%% 2. 相位误差处理
if isempty(phi_true)
    phi_true = zeros(M, 1);
else
    phi_true = phi_true(:);
    assert(numel(phi_true) == M, 'phi_true 的长度必须等于 cfg.M。');
    assert(isreal(phi_true) && all(isfinite(phi_true)), ...
        'phi_true 必须是有限实数。');
end

%% 3. 初始化接收矩阵
X_noiseless = complex(zeros(M, N_rx));
P = numel(paths);

% 默认保存每条路径分量，便于调试、画图和指标分解。
if isfield(cfg, 'store_path_components')
    store_path_components = cfg.store_path_components;
else
    store_path_components = true;
end

if store_path_components
    X_paths = complex(zeros(M, N_rx, P));
else
    X_paths = [];
end

%% 4. 遍历路径并叠加
t_tx = cfg.t_tx(:);
t_rx = cfg.t_rx(:);

path_peak_amp = zeros(P, 1);
path_energy = zeros(P, 1);

for p = 1:P
    check_path_fields(paths(p), p);

    tau = paths(p).delay;
    theta = paths(p).theta;
    alpha = paths(p).alpha;

    % 4.1 生成延迟复基带包络 s(t_rx - tau)。
    %
    % t_rx 是绝对接收时间轴，s 的定义时间为 t_tx = 0 ... T。
    % 使用线性插值处理非整数采样延迟；超出脉冲支撑区间处补 0。
    t_query = t_rx - tau;
    s_delayed = interp1(t_tx, s, t_query, 'linear', 0);

    % 4.2 生成真实阵列导向矢量，包含阵列相位失配。
    a_true = steering_vector(theta, cfg, phi_true);

    % 4.3 当前路径对阵列接收信号的贡献。
    %
    % a_true:    M x 1
    % s_delayed: N_rx x 1
    % X_p:       M x N_rx
    X_p = a_true * (alpha * s_delayed.');

    X_noiseless = X_noiseless + X_p;

    if store_path_components
        X_paths(:, :, p) = X_p;
    end

    path_peak_amp(p) = max(abs(X_p(:)));
    path_energy(p) = sum(abs(X_p(:)).^2);
end

%% 5. 生成复基带白噪声
noise_seed = make_noise_seed(cfg, mc_index);
noise_stream = RandStream('mt19937ar', 'Seed', noise_seed);

if cfg.noise_power > 0
    noise_sigma = sqrt(cfg.noise_power / 2);
    noise = noise_sigma * ...
        (randn(noise_stream, M, N_rx) + 1j * randn(noise_stream, M, N_rx));
else
    noise_sigma = 0;
    noise = complex(zeros(M, N_rx));
end

%% 6. 输出总接收信号
X = X_noiseless + noise;

%% 7. 分量输出
components = struct();

components.received = X;
components.noiseless = X_noiseless;
components.noise = noise;
components.path_signals = X_paths;

components.path_peak_amp = path_peak_amp;
components.path_energy = path_energy;

components.path_type = {paths.type}.';
components.path_name = {paths.name}.';
components.path_delay = [paths.delay].';
components.path_theta = [paths.theta].';
components.path_theta_deg = [paths.theta_deg].';
components.path_alpha = [paths.alpha].';
components.path_amp = [paths.amp].';
components.path_power = [paths.power].';
components.is_interference = [paths.is_interference].';
components.is_target = [paths.is_target].';

components.interference = [];
components.target = [];

%% 8. 信息输出
info = struct();

info.model = 'complex_baseband_array_received_signal';
info.equation = 'X(t)=sum_p alpha_p*a_true(theta_p)*s(t-tau_p)+n(t)';

info.M = M;
info.N_rx = N_rx;
info.N_tx = N_tx;
info.num_paths = P;

info.t_rx_start = cfg.t_rx(1);
info.t_rx_end = cfg.t_rx(end);
info.fs = cfg.fs;

info.noise_seed = noise_seed;
info.noise_power_setting = cfg.noise_power;
info.noise_sigma_per_real_dimension = noise_sigma;

info.phi_true = phi_true;
info.phase_error_used = any(abs(phi_true) > 0);

info.signal_power_mean_noiseless = mean(abs(X_noiseless(:)).^2);
info.noise_power_empirical = mean(abs(noise(:)).^2);
info.rx_power_mean_total = mean(abs(X(:)).^2);

%% 9. 能量分解
% window_energy_sinr_dB 是整个接收观测窗口上的能量 SINR。
% 它不等同于 cfg.SNR；cfg.SNR 是目标回波有效脉冲区域内的功率参考。
if store_path_components
    idx_interf = find([paths.is_interference]);
    idx_target = find([paths.is_target]);

    if ~isempty(idx_interf)
        X_interf = sum(X_paths(:, :, idx_interf), 3);
    else
        X_interf = complex(zeros(M, N_rx));
    end

    if ~isempty(idx_target)
        X_target = sum(X_paths(:, :, idx_target), 3);
    else
        X_target = complex(zeros(M, N_rx));
    end

    components.interference = X_interf;
    components.target = X_target;

    info.interference_energy = sum(abs(X_interf(:)).^2);
    info.target_energy = sum(abs(X_target(:)).^2);
    info.noise_energy = sum(abs(noise(:)).^2);

    info.window_energy_sinr_dB = 10*log10( ...
        info.target_energy / max(info.interference_energy + info.noise_energy, eps));
    info.input_sinr_energy_dB = info.window_energy_sinr_dB;
else
    info.interference_energy = NaN;
    info.target_energy = NaN;
    info.noise_energy = sum(abs(noise(:)).^2);
    info.window_energy_sinr_dB = NaN;
    info.input_sinr_energy_dB = NaN;
end
end

%% ========================================================================
% 局部辅助函数
% ========================================================================

function check_path_fields(path, p)
%CHECK_PATH_FIELDS 检查路径结构体字段。

required_path_fields = {'delay', 'theta', 'alpha', 'type', 'name', ...
    'amp', 'power', 'is_interference', 'is_target', 'theta_deg'};

for k = 1:numel(required_path_fields)
    field_name = required_path_fields{k};
    assert(isfield(path, field_name), ...
        'paths(%d).%s 是必需字段。', p, field_name);
end

assert(isfinite(path.delay) && path.delay >= 0, ...
    'paths(%d).delay 必须是非负有限数。', p);
assert(isfinite(path.theta) && isreal(path.theta), ...
    'paths(%d).theta 必须是有限实数。', p);
assert(isfinite(path.alpha), ...
    'paths(%d).alpha 必须是有限复数。', p);
end

function seed = make_noise_seed(cfg, mc_index)
%MAKE_NOISE_SEED 生成噪声随机种子。

if isfield(cfg, 'noise_seed_offset')
    noise_seed_offset = cfg.noise_seed_offset;
else
    noise_seed_offset = 20001;
end

seed_raw = cfg.random_seed + noise_seed_offset + 100000 * (mc_index - 1);
seed = mod(seed_raw, 2^32 - 1);
end
