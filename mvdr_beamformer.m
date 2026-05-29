function [Y_res, clean_data, info] = mvdr_beamformer(cfg, X, Y, mf_data, s, mf, paths)
%MVDR_BEAMFORMER MVDR (Capon) 自适应波束形成器，在阵元域抑制直达波。
%
% 用法:
%   [Y_res, clean_data, info] = mvdr_beamformer(cfg, X, Y, mf_data, s, mf, paths);
%
% 输入:
%   cfg      sonar_config() 配置结构体。
%   X        M x N_rx 阵列接收数据矩阵，来自 simulate_received_signal()。
%   Y        M x N_mf 原始匹配滤波输出，来自 matched_filter_bank()。
%            仅用于计算残余能量比值的基准，不参与 MVDR 波束形成计算。
%   mf_data  matched_filter_bank() 输出结构体，至少包含 t_mf。
%   s        N_tx x 1 发射 LFM 信号，来自 generate_lfm()。
%   mf       N_tx x 1 匹配滤波器，来自 generate_lfm()。
%   paths    generate_virtual_sources() 输出的路径结构体数组。
%
% 输出:
%   Y_res       M x N_mf MVDR 处理后输出（阵元维度已做缩放复制，
%              保证非相干求和与真实 MVDR 输出功率一致）。
%   clean_data  结果结构体，与 clean_baseline 输出字段对齐。
%   info        算法信息结构体。
%
% 物理模型:
%   MVDR 在阵元域计算采样协方差 Rxx = X*X^H / N_rx，
%   对目标方向施加单位增益约束，最小化总输出功率:
%       w_mvdr = Rxx^{-1} * a(θ_target) / (a^H * Rxx^{-1} * a)
%
%   强直达波来自不同方向 → MVDR 自动形成零陷。
%   阵列相位误差被包含在 Rxx 中 → MVDR 自适应补偿，无需显式校准。
%
%   与 CLEAN 类方法的本质区别:
%       CLEAN:    在匹配滤波域迭代剥离最强分量（标量处理）
%       MVDR:     在阵元域通过空间滤波一次性抑制干扰（矢量处理）

%% 1. 输入检查
required_fields = { ...
    'M', 'fs', 'N_rx', 'N_tx', 't_rx', ...
    'array_pos', 'array_center', 'array_axis', 'broadside_axis', ...
    'lambda', 'direct_window', 'multipath_window'};

for k = 1:numel(required_fields)
    field_name = required_fields{k};
    assert(isfield(cfg, field_name), 'cfg.%s 是必需字段。', field_name);
end

assert(isfield(mf_data, 't_mf'), 'mf_data.t_mf 是必需字段。');

M = cfg.M;
fs = cfg.fs;
N_rx = cfg.N_rx;
N_tx = cfg.N_tx;

assert(size(X, 1) == M, 'X 的行数必须等于 cfg.M。');
assert(size(X, 2) == N_rx, 'X 的列数必须等于 cfg.N_rx。');
assert(all(isfinite(X(:))), 'X 中存在非有限数。');

assert(size(Y, 1) == M, 'Y 的行数必须等于 cfg.M。');
assert(all(isfinite(Y(:))), 'Y 中存在非有限数。');

s = s(:);
mf = mf(:);

assert(numel(s) == N_tx, 's 的长度必须等于 cfg.N_tx。');
assert(numel(mf) == N_tx, 'mf 的长度必须等于 cfg.N_tx。');
assert(~isempty(paths), 'paths 不能为空。');

t_mf = mf_data.t_mf(:);
N_mf = numel(t_mf);

assert(size(Y, 2) == N_mf, 'Y 的列数必须等于 numel(mf_data.t_mf)。');

%% 2. 获取目标 DOA
target_idx = find(strcmp({paths.type}, 'target'), 1);
assert(~isempty(target_idx), 'paths 中未找到 target 路径。');

theta_target = paths(target_idx).theta;

%% 3. 获取 MVDR 参数
if isfield(cfg, 'mvdr_diagonal_loading')
    lambda = cfg.mvdr_diagonal_loading;
else
    lambda = 1e-2;
end

assert(lambda >= 0, 'mvdr_diagonal_loading 必须非负。');

%% 4. 计算采样协方差矩阵
Rxx = (X * X') / N_rx;

% 对角加载，保证数值稳定性。
tr_Rxx = abs(sum(diag(Rxx)));
if tr_Rxx < eps
    tr_Rxx = 1;
end
Rxx_reg = Rxx + lambda * (tr_Rxx / M) * eye(M);

%% 5. 计算 MVDR 权向量
a_target = steering_vector(theta_target, cfg);
a_target = a_target(:);

assert(numel(a_target) == M, '导向矢量维度错误。');

% w = R^{-1} a / (a^H R^{-1} a)
R_inv_a = Rxx_reg \ a_target;
denom = a_target' * R_inv_a;

assert(abs(denom) > eps, ...
    'MVDR 分母接近零，请增大 diagonal_loading 或检查目标方向。');

w_mvdr = R_inv_a / denom;

%% 6. 波束形成
y_bf = w_mvdr' * X;   % 1 x N_rx
y_bf = y_bf(:).';     % 确保行向量

%% 7. 匹配滤波
y_mf = conv(y_bf(:), mf, 'full');
y_mf = y_mf(:).';

assert(numel(y_mf) == N_mf, ...
    'MVDR 匹配滤波输出长度与 matched_filter_bank 不一致。');

%% 8. 构造阵元域输出（缩放复制，保证非相干求和正确）
% compute_metrics 对 Y_res 做 non-coherent summation over elements:
%   power(tau) = sum_m |Y_res(m, tau_idx)|^2
%
% 每行放 y_mf / sqrt(M)，则 sum_m |y_mf/sqrt(M)|^2 = |y_mf|^2。
Y_res = repmat(y_mf / sqrt(M), M, 1);

%% 9. 构造干扰搜索区域（与 CLEAN 一致，用于计算残余能量指标）
[search_idx, search_mask, search_info] = build_search_indices(cfg, t_mf, paths);

%% 10. 计算残余能量指标
initial_mask_energy = energy_sum(Y(:, search_idx));
initial_total_energy = energy_sum(Y);

final_mask_energy = energy_sum(Y_res(:, search_idx));
final_total_energy = energy_sum(Y_res);

final_mask_energy_rel = final_mask_energy / max(initial_mask_energy, eps);
final_total_energy_rel = final_total_energy / max(initial_total_energy, eps);

%% 11. 输出 clean_data（与 CLEAN 接口对齐）
clean_data = struct();

clean_data.Y_original = Y;
clean_data.Y_residual = Y_res;
clean_data.Y_model = [];  % MVDR 没有重建的干扰模型

clean_data.t_mf = t_mf;
clean_data.search_idx = search_idx;
clean_data.search_mask = search_mask;
clean_data.search_info = search_info;

clean_data.iter_table = table();
clean_data.num_iter = 1;        % MVDR 是单次处理，非迭代
clean_data.stop_reason = 'mvdr_adaptive_beamforming';

clean_data.final_mask_energy = final_mask_energy;
clean_data.final_mask_energy_rel = final_mask_energy_rel;
clean_data.final_total_energy = final_total_energy;
clean_data.final_total_energy_rel = final_total_energy_rel;
clean_data.reconstruction_energy = NaN;

%% 12. 信息结构体
info = struct();

info.model = 'mvdr_adaptive_beamformer_in_element_time_domain';
info.algorithm = 'MVDR_Capon_minimum_variance_distortionless_response';
info.uses_phase_calibration = false;

info.physical_model = ...
    'w_mvdr = Rxx^{-1} a(theta_target) / (a^H Rxx^{-1} a); y_bf = w^H X; y_mf = conv(y_bf, mf)';

info.look_direction_theta = theta_target;
info.look_direction_theta_deg = theta_target * 180/pi;
info.diagonal_loading = lambda;
info.sample_covariance_trace = tr_Rxx;

info.num_iter = 1;
info.stop_reason = 'mvdr_adaptive_beamforming';

info.initial_mask_energy = initial_mask_energy;
info.initial_total_energy = initial_total_energy;
info.final_mask_energy = final_mask_energy;
info.final_mask_energy_rel = final_mask_energy_rel;
info.final_total_energy = final_total_energy;
info.final_total_energy_rel = final_total_energy_rel;

info.t_mf_start = t_mf(1);
info.t_mf_end = t_mf(end);

info.search_delay_min = min(t_mf(search_idx));
info.search_delay_max = max(t_mf(search_idx));
info.num_search_delay_bins = numel(search_idx);
end

%% ========================================================================
% 局部辅助函数
% ========================================================================

function [search_idx, search_mask, search_info] = build_search_indices(cfg, t_mf, paths)
%BUILD_SEARCH_INDICES 构造干扰路径附近的时延搜索区域。
%
% 与 clean_baseline.m 中 build_clean_search_indices 使用相同的逻辑:
%   direct 路径使用 cfg.direct_window；
%   其他多径使用 cfg.multipath_window。

N_mf = numel(t_mf);
search_mask = false(N_mf, 1);

search_info = struct();
search_info.mode = 'interference_windows_from_paths';
search_info.windows = struct('type', {}, 'delay', {}, 'half_window', {}, ...
    'idx_start', {}, 'idx_end', {});

for p = 1:numel(paths)
    if ~isfield(paths(p), 'is_interference') || ~paths(p).is_interference
        continue;
    end

    tau0 = paths(p).delay;

    if strcmp(paths(p).type, 'direct')
        half_window = cfg.direct_window;
    else
        half_window = cfg.multipath_window;
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
        search_info.windows(end+1) = w; %#ok<AGROW>
    end
end

if ~any(search_mask)
    warning('MVDR: 未能构造干扰搜索窗口，退化为全时延搜索。');
    search_mask(:) = true;
    search_info.mode = 'fallback_full_delay_axis';
end

search_idx = find(search_mask);
end

function e = energy_sum(X)
%ENERGY_SUM 计算复矩阵总能量。

e = sum(abs(X(:)).^2);
end
