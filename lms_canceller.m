function [Y_res, clean_data, info] = lms_canceller(cfg, X, Y, mf_data, s, mf, paths)
%LMS_CANCELLER LMS 自适应直达波对消器。
%
% 用法:
%   [Y_res, clean_data, info] = lms_canceller(cfg, X, Y, mf_data, s, mf, paths);
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
%   Y_res       M x N_mf LMS 对消后匹配滤波输出。
%   clean_data  结果结构体。
%   info        算法信息结构体。
%
% 物理模型:
%   逐阵元使用已知直达波波形 s(t-τ_d) 作为参考信号，
%   通过复 LMS 自适应估计各通道直达波复增益 w_m，
%   从接收信号中减去 w_m * s(t-τ_d) 得到对消后信号:
%       e_m(t) = x_m(t) - w_m * s(t-τ_d)
%
%   LMS 更新:
%       w_m(n+1) = w_m(n) + μ * conj(e_m(n)) * s(t_n - τ_d)
%
%   对消后的信号再做匹配滤波和波束形成。
%
%   与校准感知 CLEAN 的关系:
%       两者都利用了几何先验（τ_d 已知）；
%       LMS 隐式学习 α_d * e^{jφ_m} * a_nom,m(θ_d)，不做显式相位估计；
%       校准 CLEAN 显式估计 φ_m 然后用于导向矢量校准。

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

%% 3. 获取 LMS 参数
if isfield(cfg, 'lms_step_size')
    mu = cfg.lms_step_size;
else
    mu = 0.01;
end

if isfield(cfg, 'lms_num_taps')
    L = cfg.lms_num_taps;
else
    L = 1;
end

if isfield(cfg, 'lms_num_passes')
    num_passes = cfg.lms_num_passes;
else
    num_passes = 2;
end

assert(mu > 0 && mu < 1, 'lms_step_size 必须位于 (0, 1)。');
assert(L >= 1 && mod(L, 1) == 0, 'lms_num_taps 必须是正整数。');
assert(num_passes >= 1 && mod(num_passes, 1) == 0, ...
    'lms_num_passes 必须是正整数。');

%% 4. 构造参考信号
t_rx = cfg.t_rx(:);
t_tx = cfg.t_tx(:);

% 参考信号 r(t) = s(t - τ_d)，线性插值实现非整数延迟。
r_full = interp1(t_tx, s, t_rx - tau_d, 'linear', 0);
r_full = r_full(:);

% 判断参考信号活跃区域（非零支撑区）。
r_active = abs(r_full) > 1e-6;
active_indices = find(r_active);

if isempty(active_indices)
    warning('LMS: 参考信号支撑区与接收时间窗无重叠，回退为不做处理。');
    Y_res = Y;
    clean_data = make_empty_result(Y, t_mf, paths, cfg, 'no_reference_overlap');
    info = make_info(cfg, tau_d, mu, L, num_passes, 'no_reference_overlap');
    return;
end

%% 5. 逐阵元 LMS 对消
X_clean = complex(zeros(M, N_rx));
w_final = complex(zeros(M, L));

lms_pass1_mse = zeros(1, M);  %#ok<NASGU>

for m = 1:M
    x_m = X(m, :).';

    w = complex(zeros(L, 1));
    w_final_m = w;
    pass1_err = NaN;

    for pass = 1:num_passes
        % 前向/后向交替，减少暂态效应。
        if mod(pass, 2) == 1
            idx_iter = active_indices;
        else
            idx_iter = flipud(active_indices);
        end

        for k = 1:numel(idx_iter)
            n = idx_iter(k);

            % 构造抽头延迟线参考向量。
            r_vec = complex(zeros(L, 1));
            for tap = 1:L
                n_tap = n - tap + 1;
                if n_tap >= 1 && n_tap <= N_rx
                    r_vec(tap) = r_full(n_tap);
                end
            end

            % 滤波输出。
            y_hat = w' * r_vec;

            % 误差。
            e = x_m(n) - y_hat;

            % LMS 权重更新。
            w = w + mu * conj(e) * r_vec;

            if pass == 1 && k == numel(idx_iter)
                w_final_m = w;
            end
        end

        if pass == 1
            pass1_err = mean(abs(x_m(active_indices) - ...
                arrayfun(@(nn) w_final_m' * build_rvec(nn, r_full, L), ...
                active_indices)).^2);
        end
    end

    w_final(m, :) = w.';
    lms_pass1_mse(m) = pass1_err;

    % 用最终权重对消整个接收信号。
    for n = 1:N_rx
        r_vec = build_rvec(n, r_full, L);
        X_clean(m, n) = x_m(n) - w' * r_vec;
    end
end

%% 6. 匹配滤波
mf_data_clean = struct();
mf_data_clean.t_mf = t_mf;
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
clean_data.num_iter = num_passes;
clean_data.stop_reason = 'lms_adaptive_cancellation';

clean_data.final_mask_energy = final_mask_energy;
clean_data.final_mask_energy_rel = final_mask_energy_rel;
clean_data.final_total_energy = final_total_energy;
clean_data.final_total_energy_rel = final_total_energy_rel;
clean_data.reconstruction_energy = NaN;

%% 9. 信息结构体
info = struct();

info.model = 'lms_adaptive_interference_canceller';
info.algorithm = 'complex_LMS_per_element_direct_path_cancellation';
info.uses_phase_calibration = false;

info.physical_model = ...
    'x_clean,m(t) = x_m(t) - w_m^H * r_vec(t); r(t) = s(t-tau_d)';

info.direct_delay = tau_d;
info.step_size = mu;
info.num_taps = L;
info.num_passes = num_passes;

info.w_final = w_final;
info.lms_pass1_mse = lms_pass1_mse;

info.num_iter = num_passes;
info.stop_reason = 'lms_adaptive_cancellation';

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
%ENERGY_SUM 计算复矩阵总能量。

e = sum(abs(X(:)).^2);
end

function clean_data = make_empty_result(Y, t_mf, paths, cfg, reason)
%MAKE_EMPTY_RESULT 构造空结果。

[search_idx, search_mask, search_info] = build_search_indices(cfg, t_mf, paths);

initial_mask_energy = energy_sum(Y(:, search_idx));
initial_total_energy = energy_sum(Y);

clean_data = struct();
clean_data.Y_original = Y;
clean_data.Y_residual = Y;
clean_data.Y_model = [];
clean_data.t_mf = t_mf;
clean_data.search_idx = search_idx;
clean_data.search_mask = search_mask;
clean_data.search_info = search_info;
clean_data.iter_table = table();
clean_data.num_iter = 0;
clean_data.stop_reason = reason;
clean_data.final_mask_energy = initial_mask_energy;
clean_data.final_mask_energy_rel = 1;
clean_data.final_total_energy = initial_total_energy;
clean_data.final_total_energy_rel = 1;
clean_data.reconstruction_energy = NaN;
end

function info = make_info(cfg, tau_d, mu, L, num_passes, stop_reason)
%MAKE_INFO 构造 LMS 信息结构体。

info = struct();
info.model = 'lms_adaptive_interference_canceller';
info.algorithm = 'complex_LMS_per_element_direct_path_cancellation';
info.uses_phase_calibration = false;
info.direct_delay = tau_d;
info.step_size = mu;
info.num_taps = L;
info.num_passes = num_passes;
info.num_iter = 0;
info.stop_reason = stop_reason;
end
