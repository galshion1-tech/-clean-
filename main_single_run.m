function result = main_single_run(cfg, mc_index, overrides, opts)
%MAIN_SINGLE_RUN 运行一次完整的声呐直达波抑制仿真实验。
%
% 用法:
%   result = main_single_run();
%   result = main_single_run(cfg);
%   result = main_single_run(cfg, mc_index);
%   result = main_single_run(cfg, mc_index, overrides);
%   result = main_single_run(cfg, mc_index, overrides, opts);

%% 0. 添加项目子目录到 MATLAB 路径
root_dir = fileparts(mfilename('fullpath'));
addpath(fullfile(root_dir, 'lib'));
addpath(fullfile(root_dir, 'adaptive'));
%
% 输入:
%   cfg        sonar_config() 输出配置。若为空，则自动调用 sonar_config()。
%   mc_index   Monte Carlo 序号，默认 1。
%   overrides  参数覆盖结构体，用于参数扫描。
%              例如:
%                 overrides.phase_error_std = 10*pi/180;
%                 overrides.SNR = 0;
%                 overrides.INR = 30;
%
%   opts       运行选项结构体:
%              opts.run_baseline      默认 true
%              opts.run_calibrated    默认 true
%              opts.run_oracle        默认 true
%              opts.store_signals     默认 false
%              opts.store_clean_data  默认 false
%              opts.verbose           默认 true
%
% 输出:
%   result     单次实验结果结构体。
%
% 说明:
%   该函数是参数扫描的基础入口。后续 run_phase_error_sweep.m、
%   run_snr_sweep.m 等脚本应调用本函数，而不是重复写完整流程。

%% 1. 输入默认值
if nargin < 1 || isempty(cfg)
    cfg = sonar_config();
end

if nargin < 2 || isempty(mc_index)
    mc_index = 1;
end

if nargin < 3 || isempty(overrides)
    overrides = struct();
end

if nargin < 4 || isempty(opts)
    opts = struct();
end

opts = fill_default_opts(opts);

assert(mc_index >= 1 && mod(mc_index, 1) == 0, ...
    'mc_index 必须是正整数。');

%% 2. 应用参数覆盖，并刷新派生参数
cfg = apply_overrides(cfg, overrides);
cfg = refresh_dependent_config(cfg, overrides);

t_start = tic;

%% 3. 生成发射 LFM 和匹配滤波器
[s, mf, wav_info] = generate_lfm(cfg);

%% 4. 生成虚源路径，并更新接收时间轴
[paths, cfg] = generate_virtual_sources(cfg);

%% 5. 生成真实阵列相位误差
[phi_true, phase_info] = generate_phase_error(cfg, mc_index);

%% 6. 生成阵列接收信号
[X, components, rx_info] = simulate_received_signal(cfg, s, paths, phi_true, mc_index);

%% 7. 匹配滤波
[Y, mf_data, mf_info] = matched_filter_bank(cfg, X, mf, paths);

%% 8. Baseline CLEAN
clean_base = [];
info_base = [];

if opts.run_baseline
    [~, clean_base, info_base] = ...
        clean_baseline(cfg, Y, mf_data, s, mf, paths);
end

%% 9. 相位误差估计
phi_hat = [];
calib_data = [];
est_info = [];

if opts.run_calibrated && cfg.enable_calibration
    if cfg.enable_phase_error_estimation
        [phi_hat, calib_data, est_info] = ...
            estimate_phase_error(cfg, Y, mf_data, s, mf, paths);
    else
        phi_hat = zeros(cfg.M, 1);
        est_info = struct();
        est_info.method = 'disabled_phase_error_estimation_zero_phi';
        est_info.phi_hat = phi_hat;
        est_info.phi_hat_deg = phi_hat * 180/pi;
    end
end

%% 10. Calibration-aware CLEAN
clean_cal = [];
info_cal = [];

if opts.run_calibrated && cfg.enable_calibration
    assert(~isempty(phi_hat), ...
        '校准感知 CLEAN 需要 phi_hat，但当前 phi_hat 为空。');

    [~, clean_cal, info_cal] = ...
        clean_calibration_aware(cfg, Y, mf_data, s, mf, paths, phi_hat);
end

%% 11. Oracle CLEAN
clean_oracle = [];
info_oracle = [];

if opts.run_oracle
    [~, clean_oracle, info_oracle] = ...
        clean_calibration_aware(cfg, Y, mf_data, s, mf, paths, phi_true);
end

%% 12. MVDR 自适应波束形成
clean_mvdr = [];
info_mvdr = [];

if opts.run_mvdr && cfg.enable_mvdr
    [~, clean_mvdr, info_mvdr] = ...
        mvdr_beamformer(cfg, X, Y, mf_data, s, mf, paths);
end

%% 13. LMS 自适应对消
clean_lms = [];
info_lms = [];

if opts.run_lms && cfg.enable_lms
    [~, clean_lms, info_lms] = ...
        lms_canceller(cfg, X, Y, mf_data, s, mf, paths);
end

%% 14. RLS 自适应对消
clean_rls = [];
info_rls = [];

if opts.run_rls && cfg.enable_rls
    [~, clean_rls, info_rls] = ...
        rls_canceller(cfg, X, Y, mf_data, s, mf, paths);
end

%% 15. 统一指标计算
metrics = compute_metrics(cfg, Y, mf_data, paths, ...
    clean_base, clean_cal, clean_oracle, clean_mvdr, clean_lms, clean_rls, ...
    phi_true, phi_hat);

runtime_sec = toc(t_start);

%% 16. 整理输出
result = struct();

result.ok = true;
result.mc_index = mc_index;
result.runtime_sec = runtime_sec;

result.cfg = cfg;
result.overrides = overrides;
result.opts = opts;

result.paths = paths;

result.phi_true = phi_true;
result.phi_hat = phi_hat;

result.metrics = metrics;

result.info = struct();
result.info.wav = wav_info;
result.info.phase = phase_info;
result.info.rx = rx_info;
result.info.mf = mf_info;
result.info.baseline = info_base;
result.info.estimation = est_info;
result.info.calibrated = info_cal;
result.info.oracle = info_oracle;
result.info.mvdr = info_mvdr;
result.info.lms = info_lms;
result.info.rls = info_rls;

if opts.store_clean_data
    result.clean = struct();
    result.clean.baseline = clean_base;
    result.clean.calibrated = clean_cal;
    result.clean.oracle = clean_oracle;
    result.clean.mvdr = clean_mvdr;
    result.clean.lms = clean_lms;
    result.clean.rls = clean_rls;
else
    result.clean = struct();
    result.clean.baseline = slim_clean_data(clean_base);
    result.clean.calibrated = slim_clean_data(clean_cal);
    result.clean.oracle = slim_clean_data(clean_oracle);
    result.clean.mvdr = slim_clean_data(clean_mvdr);
    result.clean.lms = slim_clean_data(clean_lms);
    result.clean.rls = slim_clean_data(clean_rls);
end

if opts.store_signals
    result.signals = struct();
    result.signals.s = s;
    result.signals.mf = mf;
    result.signals.X = X;
    result.signals.Y = Y;
    result.signals.components = components;
    result.signals.mf_data = mf_data;
else
    result.signals = [];
end

%% 17. 命令行摘要
if opts.verbose
    print_single_run_summary(result);
end
end

%% ========================================================================
% 局部辅助函数
% ========================================================================

function opts = fill_default_opts(opts)
%FILL_DEFAULT_OPTS 补齐运行选项默认值。

if ~isfield(opts, 'run_baseline')
    opts.run_baseline = true;
end

if ~isfield(opts, 'run_calibrated')
    opts.run_calibrated = true;
end

if ~isfield(opts, 'run_oracle')
    opts.run_oracle = true;
end

if ~isfield(opts, 'run_mvdr')
    opts.run_mvdr = true;
end

if ~isfield(opts, 'run_lms')
    opts.run_lms = true;
end

if ~isfield(opts, 'run_rls')
    opts.run_rls = true;
end

if ~isfield(opts, 'store_signals')
    opts.store_signals = false;
end

if ~isfield(opts, 'store_clean_data')
    opts.store_clean_data = false;
end

if ~isfield(opts, 'verbose')
    opts.verbose = true;
end
end

function cfg = apply_overrides(cfg, overrides)
%APPLY_OVERRIDES 将 overrides 中的字段覆盖到 cfg 中。

assert(isstruct(overrides), 'overrides 必须是结构体。');

field_names = fieldnames(overrides);

for k = 1:numel(field_names)
    field_name = field_names{k};
    cfg.(field_name) = overrides.(field_name);
end
end

function cfg = refresh_dependent_config(cfg, overrides)
%REFRESH_DEPENDENT_CONFIG 覆盖参数后刷新派生参数。
%
% 该函数主要用于参数扫描，例如改变 SNR、INR、phase_error_std、
% max_order、enable_multipath、target_window 等。
%
% 推荐参数扫描时优先覆盖:
%   SNR, INR, phase_error_std, max_order, enable_multipath
%
% 不建议直接覆盖 direct_amp/noise_power，除非非常清楚功率标定关系。

%% 1. 功率相关派生量
if ~isfield(cfg, 'target_power')
    cfg.target_power = 1;
end

cfg.target_amp = sqrt(cfg.target_power);

if isfield(cfg, 'INR')
    cfg.direct_to_target_dB = cfg.INR;
else
    cfg.INR = cfg.direct_to_target_dB;
end

cfg.direct_amp = cfg.target_amp * sqrt(10^(cfg.INR/10));

if isfield(cfg, 'SNR')
    cfg.noise_power = cfg.target_power / 10^(cfg.SNR/10);
end

cfg.target_echo_amp = cfg.target_amp;

%% 1.1 传播损耗模型默认字段
% 兼容旧版 cfg：hybrid_physical_loss 模式真正使用这些字段；
% received_power 模式下保留它们不会影响原有算法。
if ~isfield(cfg, 'spreading_exponent')
    cfg.spreading_exponent = 1.0;
end

if ~isfield(cfg, 'enable_absorption_loss')
    cfg.enable_absorption_loss = true;
end

if ~isfield(cfg, 'absorption_dB_per_km')
    cfg.absorption_dB_per_km = 1.0;
end

%% 2. 波形与阵列相关派生量
cfg.lambda = cfg.c / cfg.fc;

cfg.N_tx = round(cfg.T * cfg.fs);
cfg.t_tx = (0:cfg.N_tx-1).' / cfg.fs;

cfg.N = cfg.N_tx;
cfg.t = cfg.t_tx;

if ~isfield(overrides, 'd')
    cfg.d = cfg.c / (2 * cfg.fc);
end

cfg.array_axis = cfg.array_axis(:).';
cfg.broadside_axis = cfg.broadside_axis(:).';

assert(norm(cfg.array_axis) > 0, 'array_axis 不能为零向量。');
assert(norm(cfg.broadside_axis) > 0, 'broadside_axis 不能为零向量。');

cfg.array_axis = cfg.array_axis / norm(cfg.array_axis);
cfg.broadside_axis = cfg.broadside_axis / norm(cfg.broadside_axis);

assert(abs(dot(cfg.array_axis, cfg.broadside_axis)) < 1e-10, ...
    'array_axis 与 broadside_axis 应近似正交。');

if ~isfield(overrides, 'array_center')
    cfg.array_center = cfg.rx_pos;
end

m = (0:cfg.M-1).' - (cfg.M-1)/2;
cfg.array_pos = cfg.array_center + (m * cfg.d) * cfg.array_axis;

%% 3. 搜索网格相关派生量
cfg.angle_grid = cfg.angle_min:cfg.angle_step:cfg.angle_max;
cfg.delay_step = 1 / cfg.fs;

%% 4. 指标窗口相关派生量
% 若用户没有显式覆盖 target_guard_window，则随 target_window 自动更新。
if ~isfield(overrides, 'target_guard_window')
    cfg.target_guard_window = 2 * cfg.target_window;
end

%% 5. 初始接收时间轴
% generate_virtual_sources() 后还会根据所有路径重新更新 t_rx。
cfg.tau_direct_nominal = norm(cfg.rx_pos - cfg.tx_pos) / cfg.c;
cfg.tau_target_nominal = ...
    (norm(cfg.target_pos - cfg.tx_pos) + norm(cfg.target_pos - cfg.rx_pos)) / cfg.c;

cfg.t_rx_start = max(0, cfg.tau_direct_nominal - cfg.rx_margin_before);
cfg.t_rx_end = cfg.tau_target_nominal + cfg.T + cfg.rx_margin_after;

cfg.N_rx = ceil((cfg.t_rx_end - cfg.t_rx_start) * cfg.fs) + 1;
cfg.t_rx = cfg.t_rx_start + (0:cfg.N_rx-1).' / cfg.fs;

%% 6. 基本合法性检查
switch cfg.signal_model
    case 'complex_baseband'
        assert(cfg.fs >= cfg.B, '复基带主线要求 fs >= B。');
    case 'real_passband'
        assert(cfg.fs >= 2 * (cfg.fc + cfg.B/2), ...
            '实数带通信号要求 fs >= 2*(fc+B/2)。');
    otherwise
        error('signal_model 必须是 ''complex_baseband'' 或 ''real_passband''。');
end

assert(cfg.ref_sensor >= 1 && cfg.ref_sensor <= cfg.M, ...
    'ref_sensor 必须位于 [1, M]。');
assert(cfg.target_window > 0, 'target_window 必须为正数。');
assert(cfg.target_guard_window > 0, 'target_guard_window 必须为正数。');
assert(cfg.direct_window > 0, 'direct_window 必须为正数。');
assert(cfg.interference_window > 0, 'interference_window 必须为正数。');
assert(cfg.multipath_window > 0, 'multipath_window 必须为正数。');
assert(cfg.spreading_exponent > 0, 'spreading_exponent 必须为正数。');
assert(cfg.absorption_dB_per_km >= 0, 'absorption_dB_per_km 必须非负。');
end

function clean_slim = slim_clean_data(clean_data)
%SLIM_CLEAN_DATA 去掉大矩阵，只保留扫描实验需要的小字段。

if isempty(clean_data)
    clean_slim = [];
    return;
end

clean_slim = struct();

copy_fields = { ...
    'num_iter', ...
    'stop_reason', ...
    'final_mask_energy', ...
    'final_mask_energy_rel', ...
    'final_total_energy', ...
    'final_total_energy_rel', ...
    'reconstruction_energy', ...
    'phi_hat', ...
    'phi_hat_deg', ...
    'iter_table'};

for k = 1:numel(copy_fields)
    f = copy_fields{k};
    if isfield(clean_data, f)
        clean_slim.(f) = clean_data.(f);
    end
end
end

function print_single_run_summary(result)
%PRINT_SINGLE_RUN_SUMMARY 打印单次实验摘要。

fprintf('\nMAIN_SINGLE_RUN summary:\n');
fprintf('-------------------------------------------------------------\n');
fprintf('mc_index       = %d\n', result.mc_index);
fprintf('runtime        = %.3f s\n', result.runtime_sec);
fprintf('SNR            = %.3f dB\n', result.cfg.SNR);
fprintf('INR            = %.3f dB\n', result.cfg.INR);
fprintf('phase std      = %.3f deg\n', result.cfg.phase_error_std * 180/pi);

summary = result.metrics.summary;

fprintf('\nCLEAN residual:\n');
fprintf('Baseline MaskRel   = %.6g\n', summary.baseline_clean_mask_rel);
fprintf('Calibrated MaskRel = %.6g\n', summary.calibrated_clean_mask_rel);
fprintf('Oracle MaskRel     = %.6g\n', summary.oracle_clean_mask_rel);

fprintf('\nImprovement:\n');
fprintf('Calibrated vs Baseline clean-mask gain = %.3f dB\n', ...
    summary.calibrated_vs_baseline_clean_mask_gain_dB);
fprintf('Calibrated vs Baseline PDR gain        = %.3f dB\n', ...
    summary.calibrated_vs_baseline_pdr_gain_dB);

fprintf('\nPDR:\n');
fprintf('Baseline PDR   = %.3f dB\n', summary.baseline_pdr_dB);
fprintf('Calibrated PDR = %.3f dB\n', summary.calibrated_pdr_dB);
fprintf('Oracle PDR     = %.3f dB\n', summary.oracle_pdr_dB);

fprintf('\nPhase estimation:\n');
fprintf('phase RMSE active = %.6f deg\n', summary.phase_rmse_active_deg);
fprintf('-------------------------------------------------------------\n');
end
