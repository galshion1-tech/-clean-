function [Y_res, clean_data, info] = rls_canceller(cfg, X, Y, mf_data, s, mf, paths)
%RLS_CANCELLER RLS 自适应直达波对消器。
%
% 用法:
%   [Y_res, clean_data, info] = rls_canceller(cfg, X, Y, mf_data, s, mf, paths);
%
% 输入:
%   cfg      sonar_config() 配置结构体。
%   X        M x N_rx 阵列接收数据矩阵。
%   Y        M x N_mf 原始匹配滤波输出，仅用于残余能量基准计算。
%   mf_data  matched_filter_bank() 输出结构体，至少包含 t_mf。
%   s        N_tx x 1 发射 LFM 信号。
%   mf       N_tx x 1 匹配滤波器。
%   paths    generate_virtual_sources() 输出路径结构体。
%
% 输出:
%   Y_res       M x N_mf RLS 对消后匹配滤波输出。
%   clean_data  结果结构体。
%   info        算法信息结构体。
%
% 物理模型:
%   与 LMS 对消器结构相同，但使用 RLS 更新权重。
%   RLS 每步最小化指数加权的误差平方和，收敛远快于 LMS，
%   适合直达波持续时间短（T=20ms）、需要快速收敛的场景。
%
%   RLS 更新方程（复数）:
%       k = (lambda^{-1} P r) / (1 + lambda^{-1} r^H P r)
%       xi = d - w^H r
%       w = w + k xi^*
%       P = lambda^{-1} (P - k r^H P)

%% 1. 输入检查
required_fields = { ...
    'M', 'fs', 'N_rx', 'N_tx', 't_rx', ...
    'direct_window', 'multipath_window'};

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

%% 2. 获取直达波时延
direct_idx = find(strcmp({paths.type}, 'direct'), 1);
assert(~isempty(direct_idx), 'paths 中未找到 direct 路径。');

tau_d = paths(direct_idx).delay;

%% 3. 获取 RLS 参数
if isfield(cfg, 'rls_forgetting_factor')
    lam = cfg.rls_forgetting_factor;
else
    lam = 0.99;
end

if isfield(cfg, 'rls_delta')
    delta = cfg.rls_delta;
else
    delta = 100;
end

if isfield(cfg, 'rls_num_taps')
    L = cfg.rls_num_taps;
else
    L = 1;
end

assert(lam > 0 && lam <= 1, 'rls_forgetting_factor 必须位于 (0, 1]。');
assert(delta > 0, 'rls_delta 必须为正数。');
assert(L >= 1 && mod(L, 1) == 0, 'rls_num_taps 必须是正整数。');

%% 4. 构造参考信号
t_rx = cfg.t_rx(:);
t_tx = cfg.t_tx(:);

r_full = interp1(t_tx, s, t_rx - tau_d, 'linear', 0);
r_full = r_full(:);

r_active = abs(r_full) > 1e-6;
active_indices = find(r_active);

if isempty(active_indices)
    warning('RLS: 参考信号支撑区与接收时间窗无重叠，回退为不做处理。');
    Y_res = Y;
    [search_idx, search_mask, ~] = build_search_indices(cfg, t_mf, paths);
    init_energy = energy_sum(Y(:, search_idx));
    clean_data = make_empty_result(Y, init_energy, search_idx, search_mask, ...
        'no_reference_overlap');
    info = make_info(cfg, tau_d, lam, delta, L, 'no_reference_overlap');
    info.initial_mask_energy = init_energy;
    return;
end

%% 5. 逐阵元 RLS 对消
lam_inv = 1 / lam;
X_clean = complex(zeros(M, N_rx));
w_final = complex(zeros(M, L));

for m = 1:M
    x_m = X(m, :).';

    w = complex(zeros(L, 1));
    P = delta * eye(L);

    % 两次遍历：前向 + 后向，消除暂态方向性偏差。
    for pass = 1:2
        if pass == 1
            idx_iter = active_indices;
        else
            idx_iter = flipud(active_indices);
        end

        for k = 1:numel(idx_iter)
            n = idx_iter(k);

            % 构造参考向量。
            r_vec = complex(zeros(L, 1));
            for tap = 1:L
                n_tap = n - tap + 1;
                if n_tap >= 1 && n_tap <= N_rx
                    r_vec(tap) = r_full(n_tap);
                end
            end

            % RLS 增益向量。
            pi_vec = P * r_vec;
            denom = lam + r_vec' * pi_vec;

            if abs(denom) < eps
                continue;
            end

            k_vec = pi_vec / denom;

            % 先验误差。
            xi = x_m(n) - w' * r_vec;

            % 权重更新。
            w = w + k_vec * conj(xi);

            % 逆相关矩阵更新。
            P = lam_inv * (P - k_vec * (r_vec' * P));
        end
    end

    w_final(m, :) = w.';

    % 用最终权重对消。
    for n = 1:N_rx
        r_vec = build_rvec(n, r_full, L);
        X_clean(m, n) = x_m(n) - w' * r_vec;
    end
end

%% 6. 匹配滤波
[Y_clean, ~, ~] = matched_filter_bank(cfg, X_clean, mf, paths);

Y_res = Y_clean;

%% 7. 构造干扰搜索区域并计算残余能量
[search_idx, search_mask, search_info] = build_search_indices(cfg, t_mf, paths);

initial_mask_energy = energy_sum(Y(:, search_idx));
initial_total_energy = energy_sum(Y);

final_mask_energy = energy_sum(Y_res(:, search_idx));
final_total_energy = energy_sum(Y_res);

final_mask_energy_rel = final_mask_energy / max(initial_mask_energy, eps);
final_total_energy_rel = final_total_energy / max(initial_total_energy, eps);

%% 8. 输出 clean_data
clean_data = struct();

clean_data.Y_original = Y;
clean_data.Y_residual = Y_res;
clean_data.Y_model = [];

clean_data.t_mf = t_mf;
clean_data.search_idx = search_idx;
clean_data.search_mask = search_mask;
clean_data.search_info = search_info;

clean_data.iter_table = table();
clean_data.num_iter = 1;
clean_data.stop_reason = 'rls_adaptive_cancellation';

clean_data.final_mask_energy = final_mask_energy;
clean_data.final_mask_energy_rel = final_mask_energy_rel;
clean_data.final_total_energy = final_total_energy;
clean_data.final_total_energy_rel = final_total_energy_rel;
clean_data.reconstruction_energy = NaN;

%% 9. 信息结构体
info = struct();

info.model = 'rls_adaptive_interference_canceller';
info.algorithm = 'complex_RLS_per_element_direct_path_cancellation';
info.uses_phase_calibration = false;

info.physical_model = ...
    'x_clean,m(t) = x_m(t) - w_m^H * r_vec(t); r(t) = s(t-tau_d)';

info.direct_delay = tau_d;
info.forgetting_factor = lam;
info.delta = delta;
info.num_taps = L;

info.w_final = w_final;

info.num_iter = 1;
info.stop_reason = 'rls_adaptive_cancellation';

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

function r_vec = build_rvec(n, r_full, L)
%BUILD_RVEC 构造抽头延迟线参考向量。

r_vec = complex(zeros(L, 1));
N = numel(r_full);

for tap = 1:L
    n_tap = n - tap + 1;
    if n_tap >= 1 && n_tap <= N
        r_vec(tap) = r_full(n_tap);
    end
end
end

function [search_idx, search_mask, search_info] = build_search_indices(cfg, t_mf, paths)
%BUILD_SEARCH_INDICES 构造干扰路径附近的时延搜索区域。

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
    search_mask(:) = true;
    search_info.mode = 'fallback_full_delay_axis';
end

search_idx = find(search_mask);
end

function e = energy_sum(X)
%EENERGY_SUM 计算复矩阵总能量。

e = sum(abs(X(:)).^2);
end

function clean_data = make_empty_result(Y, init_energy, search_idx, search_mask, reason)
%MAKE_EMPTY_RESULT 构造空结果。

clean_data = struct();
clean_data.Y_original = Y;
clean_data.Y_residual = Y;
clean_data.Y_model = [];
clean_data.search_idx = search_idx;
clean_data.search_mask = search_mask;
clean_data.iter_table = table();
clean_data.num_iter = 0;
clean_data.stop_reason = reason;
clean_data.final_mask_energy = init_energy;
clean_data.final_mask_energy_rel = 1;
clean_data.reconstruction_energy = NaN;
end

function info = make_info(cfg, tau_d, lam, delta, L, stop_reason)
%MAKE_INFO 构造 RLS 信息结构体。

info = struct();
info.model = 'rls_adaptive_interference_canceller';
info.algorithm = 'complex_RLS_per_element_direct_path_cancellation';
info.uses_phase_calibration = false;
info.direct_delay = tau_d;
info.forgetting_factor = lam;
info.delta = delta;
info.num_taps = L;
info.num_iter = 0;
info.stop_reason = stop_reason;
end
