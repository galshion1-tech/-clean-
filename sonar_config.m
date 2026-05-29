function cfg = sonar_config()
%SONAR_CONFIG 声呐仿真主配置文件。
%
% 用法:
%   cfg = sonar_config();
%
% 单位约定:
%   距离: m
%   时间: s
%   频率/带宽: Hz
%   角度/相位: rad

%% 1. 可复现实验设置
cfg.random_seed = 20260512;       % 随机种子，用于噪声、阵列相位误差和 Monte Carlo 实验复现。

%% 2. 坐标系约定
% z = 0 表示海面，z = H 表示海底，z 轴正方向向下。
cfg.z_positive_down = true;

%% 3. 声呐波形与传播参数
cfg.c = 1500;                     % 声速 c，m/s。
cfg.fc = 20e3;                    % 载频 fc，Hz。
cfg.fs = 100e3;                   % 采样率 fs，Hz。
cfg.B = 5e3;                      % 发射信号带宽 B，Hz。
cfg.T = 20e-3;                    % 脉冲宽度 T，s。
cfg.signal_model = 'complex_baseband'; % 信号模型: 'complex_baseband' 或 'real_passband'。

%% 4. 接收阵列参数
cfg.M = 16;                       % 阵元数 M。
cfg.d = cfg.c / (2 * cfg.fc);     % 阵元间距 d，默认半波长间距。
cfg.array_axis = [0, 1, 0];       % ULA 排列方向。默认沿 y 轴排列。

%% 5. 海洋环境与几何位置
cfg.H = 100;                      % 海深 H，m。

cfg.tx_pos = [0, 0, 20];          % 发射机位置 [x, y, z]，m。

cfg.rx_pos = [1000, 0, 20];       % 接收阵列中心位置 [x, y, z]，m。
cfg.array_center = cfg.rx_pos;    % 阵列中心位置，默认与 rx_pos 一致。

cfg.target_pos = [600, 300, 50];  % 目标位置 [x, y, z]，m。

%% 6. 功率标定、SNR 与 INR 定义
% 会议版仿真采用“接收端功率归一化”作为默认约定，便于做算法对比:
%   目标回波接收功率 target_power = 1
%   直达波功率 / 目标回波功率 = 10^(INR/10)
%   目标回波功率 / 噪声功率 = 10^(SNR/10)
%
% power_control_mode = 'received_power':
%   在接收端直接设定目标、直达波和噪声功率。
%   gamma_surface/gamma_bottom 只描述多径相对强弱，避免与传播损耗重复标定。
%
% power_control_mode = 'physical_spreading':
%   后续路径生成模块根据传播距离、反射系数和扩展损耗计算所有路径幅度。
%
% power_control_mode = 'hybrid_physical_loss':
%   直达波和目标强度仍由 direct_amp/target_amp 控制；
%   海面/海底多径在直达波幅度基础上叠加反射、相对扩展和吸收损耗。
cfg.SNR = 10;                     % 目标回波 SNR，dB。
cfg.INR = 30;                     % 直达波相对目标回波的功率比，dB。
cfg.direct_to_target_dB = cfg.INR; % INR 定义的显式别名，避免与传统干噪比混淆。
cfg.power_control_mode = 'hybrid_physical_loss'; % 功率控制模式，默认保留稳定的接收端功率归一化。

cfg.target_power = 1;             % 目标回波参考功率。
cfg.target_amp = sqrt(cfg.target_power);
cfg.direct_amp = cfg.target_amp * sqrt(10^(cfg.INR/10));
cfg.noise_power = cfg.target_power / 10^(cfg.SNR/10);

cfg.target_echo_amp = cfg.target_amp; % 兼容旧脚本的字段名。

%% 7. 虚源法与多径参数
cfg.enable_multipath = true;      % 是否开启多径。
cfg.max_order = 1;                % 虚源反射阶数。调试默认 1；复杂多径实验可设为 2。

cfg.gamma_direct = 1.0;           % 直达路径系数。
cfg.gamma_surface = -0.8;         % 海面反射系数，负号表示相位反转。
cfg.gamma_bottom = 0.5;           % 海底反射系数。

cfg.enable_spreading_loss = true; % 是否考虑扩展损耗。
cfg.spreading_type = 'spherical'; % 扩展损耗类型: 'none'、'cylindrical'、'spherical'。
cfg.spreading_exponent = 1.0;     % hybrid 模式下的相对扩展损耗指数: 1 球面，0.5 近似柱面。

cfg.enable_absorption_loss = true; % hybrid 模式下是否考虑吸收损耗。
cfg.absorption_dB_per_km = 1.0;    % 吸收损耗，单位 dB/km，用于相对直达波的额外路程。

%% 8. 阵列失配与校准设置
cfg.enable_mismatch = true;       % 是否开启阵列失配。
cfg.phase_error_std = 5*pi/180;   % 阵列相位误差标准差 σ_phi。
cfg.phase_error_clip = 30*pi/180; % 真实相位误差限幅。

cfg.enable_calibration = true;    % 是否开启校准。
cfg.ref_sensor = 1;               % 相位校准参考阵元，默认第 1 个阵元相位为 0。
cfg.phase_est_clip = 30*pi/180;   % 相位误差估计值限幅。

cfg.enable_oracle_calibration = false;     % 是否使用真实相位误差进行 Oracle 校准。
cfg.enable_one_shot_calibration = true;    % 是否开启一次校准。
cfg.enable_iterative_calibration = false;  % 是否开启迭代校准。
cfg.enable_phase_error_estimation = true;  % 是否估计阵列相位误差。
cfg.enable_use_true_doa = false;           % 是否使用真实 DOA。

%% 9. 目标与消融实验开关
cfg.enable_target = true;         % 是否加入目标回波。

%% 10. CLEAN 算法参数
cfg.max_iter = 100;               % CLEAN 最大迭代次数。
cfg.loop_gain = 0.1;              % CLEAN loop gain。
cfg.residual_threshold_rel = 1e-3; % 相对残差停止阈值。

cfg.residual_threshold = cfg.residual_threshold_rel; % 兼容旧脚本；新代码优先使用 residual_threshold_rel。

%% 11. CLEAN 搜索网格
cfg.angle_min = -60*pi/180;       % 最小搜索角 θ_min。
cfg.angle_max =  60*pi/180;       % 最大搜索角 θ_max。
cfg.angle_step = 0.5*pi/180;      % 搜索角步长 Δθ。
cfg.angle_grid = cfg.angle_min:cfg.angle_step:cfg.angle_max;

cfg.delay_step = 1 / cfg.fs;      % 时延搜索步长，默认一个采样间隔。
cfg.direct_search_half_window = 5e-3;    % 直达波搜索半窗宽。
cfg.multipath_search_half_window = 20e-3; % 多径搜索半窗宽。

%% 12. 对比算法参数: MVDR
cfg.enable_mvdr = true;              % 是否运行 MVDR 波束形成器。
cfg.mvdr_diagonal_loading = 1e-2;    % MVDR 对角加载系数，保证协方差矩阵可逆。

%% 13. 对比算法参数: LMS
cfg.enable_lms = true;               % 是否运行 LMS 自适应对消器。
cfg.lms_step_size = 0.05;            % LMS 步长 μ。
cfg.lms_num_taps = 1;                % LMS 抽头数，1 为单复增益对消。
cfg.lms_num_passes = 2;              % LMS 遍历次数，≥2 时前向+后向消除暂态。

%% 14. 对比算法参数: RLS
cfg.enable_rls = true;               % 是否运行 RLS 自适应对消器。
cfg.rls_forgetting_factor = 0.99;    % RLS 遗忘因子 λ。
cfg.rls_delta = 100;                 % RLS 初始 P 矩阵缩放。
cfg.rls_num_taps = 1;                % RLS 抽头数，1 为单复增益对消。

%% 16. 角度定义
% θ 表示水平 DOA，参考方向为阵列 broadside。
% 默认 ULA 沿 y 轴排列，因此 broadside 方向为 +x 轴。
% 对应导向矢量常用形式:
%   a_m(θ) = exp(-j * 2π/λ * y_m * sin(θ))
cfg.angle_definition = 'broadside_horizontal';
cfg.broadside_axis = [1, 0, 0];

%% 17. 指标计算窗口
% 以下窗口均为“半窗宽”，使用时通常取 tau0 ± window。
% direct_window 用于单独直达波峰值附近指标。
% interference_window 用于直达波及其近邻强多径区域整体指标。
% target_window 用于目标峰值保持、定位误差等指标。
cfg.direct_window = 0.5e-3;       % 直达波指标半窗宽。
cfg.interference_window = 2e-3;   % 强干扰区域指标半窗宽。
cfg.multipath_window = 2e-3;      % 多径指标半窗宽。
cfg.target_window = 2e-3;         % 目标指标半窗宽。
cfg.target_guard_window = 2 * cfg.target_window; % PDR 计算中的目标保护半窗宽。

%% 18. Monte Carlo 与输出目录
cfg.num_mc = 1;                   % Monte Carlo 次数。调试设 1，正式曲线可设 50 或 100。
cfg.result_dir = './results';     % 数值结果保存目录。
cfg.figure_dir = './figures';     % 图片结果保存目录。

%% 19. 发射波形派生参数
cfg.lambda = cfg.c / cfg.fc;      % 波长 λ，m。
cfg.N_tx = round(cfg.T * cfg.fs); % 发射脉冲采样点数。
cfg.t_tx = (0:cfg.N_tx-1).' / cfg.fs; % 发射波形时间轴，仅覆盖一个脉冲。

cfg.N = cfg.N_tx;                 % 兼容旧脚本。
cfg.t = cfg.t_tx;                 % 兼容旧脚本；接收信号仿真请使用 t_rx。

%% 20. 阵元坐标
m = (0:cfg.M-1).' - (cfg.M-1)/2;
cfg.array_pos = cfg.array_center + (m * cfg.d) * cfg.array_axis;

%% 21. 接收观测时间轴
% 当前接收窗口由直达波和目标双基地路径的名义时延确定。
% 后续完成 generate_virtual_sources.m 后，建议根据所有路径时延重新确定:
%   all_tau = [paths.delay];
%   t_rx_start = min(all_tau) - rx_margin_before;
%   t_rx_end   = max(all_tau) + T + rx_margin_after;
cfg.tau_direct_nominal = norm(cfg.rx_pos - cfg.tx_pos) / cfg.c;
cfg.tau_target_nominal = ...
    (norm(cfg.target_pos - cfg.tx_pos) + norm(cfg.target_pos - cfg.rx_pos)) / cfg.c;

cfg.rx_margin_before = 10e-3;     % 最早到达路径前的保护时间。
cfg.rx_margin_after = 50e-3;      % 最晚回波后的保护时间。

cfg.t_rx_start = max(0, cfg.tau_direct_nominal - cfg.rx_margin_before);
cfg.t_rx_end = cfg.tau_target_nominal + cfg.T + cfg.rx_margin_after;

cfg.N_rx = ceil((cfg.t_rx_end - cfg.t_rx_start) * cfg.fs) + 1;
cfg.t_rx = cfg.t_rx_start + (0:cfg.N_rx-1).' / cfg.fs;

validate_config(cfg);
end

function validate_config(cfg)
%VALIDATE_CONFIG 检查配置参数的基本合法性。

mustBePositive(cfg.c);
mustBePositive(cfg.fc);
mustBePositive(cfg.fs);
mustBePositive(cfg.B);
mustBePositive(cfg.T);
mustBeInteger(cfg.M);
mustBePositive(cfg.M);
mustBePositive(cfg.d);
mustBePositive(cfg.H);
mustBeNonnegative(cfg.max_order);
mustBeInteger(cfg.max_order);
mustBePositive(cfg.spreading_exponent);
mustBeNonnegative(cfg.absorption_dB_per_km);
mustBeNonnegative(cfg.phase_error_std);
mustBeNonnegative(cfg.phase_error_clip);
mustBeNonnegative(cfg.phase_est_clip);
mustBeInteger(cfg.max_iter);
mustBePositive(cfg.max_iter);
mustBePositive(cfg.loop_gain);
mustBeLessThanOrEqual(cfg.loop_gain, 1);
mustBePositive(cfg.residual_threshold_rel);
mustBeInteger(cfg.ref_sensor);
mustBePositive(cfg.ref_sensor);
mustBeInteger(cfg.num_mc);
mustBePositive(cfg.num_mc);

assert(numel(cfg.tx_pos) == 3, 'tx_pos 必须是 1x3 位置向量。');
assert(numel(cfg.rx_pos) == 3, 'rx_pos 必须是 1x3 位置向量。');
assert(numel(cfg.target_pos) == 3, 'target_pos 必须是 1x3 位置向量。');

assert(numel(cfg.array_axis) == 3, 'array_axis 必须是 1x3 方向向量。');
assert(norm(cfg.array_axis) > 0, 'array_axis 不能为零向量。');
assert(abs(norm(cfg.array_axis) - 1) < 1e-12, 'array_axis 必须是单位向量。');

assert(numel(cfg.broadside_axis) == 3, 'broadside_axis 必须是 1x3 方向向量。');
assert(norm(cfg.broadside_axis) > 0, 'broadside_axis 不能为零向量。');
assert(abs(norm(cfg.broadside_axis) - 1) < 1e-12, ...
    'broadside_axis 必须是单位向量。');
assert(abs(dot(cfg.array_axis, cfg.broadside_axis)) < 1e-12, ...
    '对于 ULA，array_axis 与 broadside_axis 应当正交。');

assert(cfg.ref_sensor >= 1 && cfg.ref_sensor <= cfg.M, ...
    'ref_sensor 必须位于 [1, M] 范围内。');

assert(cfg.tx_pos(3) >= 0 && cfg.tx_pos(3) <= cfg.H, ...
    'tx_pos 的 z 坐标必须位于海面和海底之间。');
assert(cfg.rx_pos(3) >= 0 && cfg.rx_pos(3) <= cfg.H, ...
    'rx_pos 的 z 坐标必须位于海面和海底之间。');
assert(cfg.target_pos(3) >= 0 && cfg.target_pos(3) <= cfg.H, ...
    'target_pos 的 z 坐标必须位于海面和海底之间。');

switch cfg.signal_model
    case 'complex_baseband'
        assert(cfg.fs >= cfg.B, ...
            '复基带仿真要求 fs >= B。');
    case 'real_passband'
        assert(cfg.fs >= 2 * (cfg.fc + cfg.B/2), ...
            '实数带通信号仿真要求 fs >= 2*(fc+B/2)。');
    otherwise
        error('signal_model 必须是 ''complex_baseband'' 或 ''real_passband''。');
end

valid_spreading_types = {'none', 'cylindrical', 'spherical'};
assert(any(strcmp(cfg.spreading_type, valid_spreading_types)), ...
    'spreading_type 必须是 ''none''、''cylindrical'' 或 ''spherical''。');

valid_power_modes = {'received_power', 'physical_spreading', 'hybrid_physical_loss'};
assert(any(strcmp(cfg.power_control_mode, valid_power_modes)), ...
    ['power_control_mode 必须是 ''received_power''、''physical_spreading'' ', ...
     '或 ''hybrid_physical_loss''。']);

valid_angle_definitions = {'broadside_horizontal'};
assert(any(strcmp(cfg.angle_definition, valid_angle_definitions)), ...
    'angle_definition 必须是 ''broadside_horizontal''。');

assert(cfg.t_rx_end > cfg.t_rx_start, 't_rx_end 必须大于 t_rx_start。');
assert(cfg.angle_max > cfg.angle_min, 'angle_max 必须大于 angle_min。');
assert(cfg.angle_step > 0, 'angle_step 必须为正数。');
assert(cfg.delay_step > 0, 'delay_step 必须为正数。');
assert(cfg.direct_window > 0, 'direct_window 必须为正数。');
assert(cfg.interference_window > 0, 'interference_window 必须为正数。');
assert(cfg.multipath_window > 0, 'multipath_window 必须为正数。');
assert(cfg.target_window > 0, 'target_window 必须为正数。');
assert(cfg.target_guard_window > 0, 'target_guard_window 必须为正数。');
end
