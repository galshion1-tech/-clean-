function [phi_true, info] = generate_phase_error(cfg, mc_index)
%GENERATE_PHASE_ERROR 生成阵列通道相位误差。
%
% 用法:
%   [phi_true, info] = generate_phase_error(cfg);
%   [phi_true, info] = generate_phase_error(cfg, mc_index);
%
% 输入:
%   cfg       sonar_config() 生成的配置结构体。
%   mc_index  Monte Carlo 序号，可选，默认 mc_index = 1。
%
% 输出:
%   phi_true  M x 1 阵列通道相位误差，单位 rad。
%             第 cfg.ref_sensor 个阵元作为参考阵元，相位误差固定为 0。
%
%   info      相位误差生成信息。
%
% 物理约定:
%   这里建模的是接收阵列各通道的相对相位误差。
%   全局共同相位不可辨识，会被路径复幅度 alpha 吸收。
%   因此参考阵元相位误差固定为 0。

%% 1. 输入检查
if nargin < 2 || isempty(mc_index)
    mc_index = 1;
end

required_fields = { ...
    'M', ...
    'enable_mismatch', ...
    'phase_error_std', ...
    'phase_error_clip', ...
    'ref_sensor', ...
    'random_seed'};

for k = 1:numel(required_fields)
    field_name = required_fields{k};
    assert(isfield(cfg, field_name), 'cfg.%s 是必需字段。', field_name);
end

M = cfg.M;
ref_sensor = cfg.ref_sensor;

assert(M > 0 && mod(M, 1) == 0, 'cfg.M 必须是正整数。');
assert(ref_sensor >= 1 && ref_sensor <= M && mod(ref_sensor, 1) == 0, ...
    'cfg.ref_sensor 必须是 [1, M] 范围内的整数。');

assert(cfg.phase_error_std >= 0, 'cfg.phase_error_std 必须非负。');
assert(cfg.phase_error_clip >= 0 || isinf(cfg.phase_error_clip), ...
    'cfg.phase_error_clip 必须非负，或为 Inf 表示不限幅。');

assert(mc_index >= 1 && mod(mc_index, 1) == 0, ...
    'mc_index 必须是正整数。');

%% 2. 初始化
phi_true = zeros(M, 1);

% 如果不开启阵列失配，或标准差为 0，则所有通道相位误差为 0。
if ~cfg.enable_mismatch || cfg.phase_error_std == 0
    info = make_info(cfg, mc_index, NaN, phi_true, 'disabled_or_zero');
    return;
end

%% 3. 构造局部随机流
% 不直接调用 rng，避免污染全局随机数状态。
% mc_index 用于 Monte Carlo 实验，使每次误差不同但可复现。
if isfield(cfg, 'phase_seed_offset')
    phase_seed_offset = cfg.phase_seed_offset;
else
    phase_seed_offset = 10001;
end

seed_raw = cfg.random_seed + phase_seed_offset + 100000 * (mc_index - 1);
seed = mod(seed_raw, 2^32 - 1);

stream = RandStream('mt19937ar', 'Seed', seed);

%% 4. 生成相对相位误差
% 参考阵元固定为 0，其余阵元独立服从 N(0, sigma_phi^2)。
active_idx = setdiff((1:M).', ref_sensor);

phi_true(active_idx) = cfg.phase_error_std * randn(stream, numel(active_idx), 1);

%% 5. 相位误差限幅
% phase_error_clip = Inf 表示不限幅。
% 当前配置默认使用 30*pi/180，避免个别阵元误差过大。
if isfinite(cfg.phase_error_clip)
    phi_true = min(max(phi_true, -cfg.phase_error_clip), cfg.phase_error_clip);
end

% 再次强制参考阵元为 0，确保相位误差是相对参考阵元定义的。
phi_true(ref_sensor) = 0;

%% 6. 输出信息
info = make_info(cfg, mc_index, seed, phi_true, 'iid_gaussian_relative_phase');
end

%% ========================================================================
% 局部辅助函数
% ========================================================================

function info = make_info(cfg, mc_index, seed, phi_true, model_name)
%MAKE_INFO 构造相位误差信息结构体。

active_idx = setdiff((1:cfg.M).', cfg.ref_sensor);

info = struct();

info.model = model_name;
info.mc_index = mc_index;
info.seed = seed;

info.M = cfg.M;
info.ref_sensor = cfg.ref_sensor;
info.active_idx = active_idx;
info.enable_mismatch = cfg.enable_mismatch;

info.phase_error_std_rad = cfg.phase_error_std;
info.phase_error_std_deg = cfg.phase_error_std * 180/pi;

info.phase_error_clip_rad = cfg.phase_error_clip;
info.phase_error_clip_deg = cfg.phase_error_clip * 180/pi;

info.phi_true = phi_true;
info.phi_true_deg = phi_true * 180/pi;

info.rms_rad = sqrt(mean(phi_true.^2));
info.rms_deg = info.rms_rad * 180/pi;

info.active_rms_rad = sqrt(mean(phi_true(active_idx).^2));
info.active_rms_deg = info.active_rms_rad * 180/pi;

info.max_abs_rad = max(abs(phi_true));
info.max_abs_deg = info.max_abs_rad * 180/pi;

info.ref_sensor_value_rad = phi_true(cfg.ref_sensor);
info.ref_sensor_value_deg = phi_true(cfg.ref_sensor) * 180/pi;
end
