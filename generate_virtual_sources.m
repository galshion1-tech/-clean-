function [paths, cfg] = generate_virtual_sources(cfg)
%GENERATE_VIRTUAL_SOURCES 生成直达波、多径虚源和目标双基地路径。
%
% 用法:
%   [paths, cfg] = generate_virtual_sources(cfg);
%
% 输入:
%   cfg   sonar_config() 生成的声呐仿真配置结构体。
%
% 输出:
%   paths 路径结构体数组，包含每条路径的时延、DOA、幅度等信息。
%   cfg   更新后的配置结构体，其中 t_rx 会根据所有路径时延重新计算。
%
% hybrid_physical_loss 模式下，直达波和目标幅度仍由 direct_amp/target_amp
% 控制；海面/海底多径幅度由直达波幅度、反射系数、相对扩展损耗和
% 吸收损耗共同决定。paths 中会记录各损耗因子，便于调试。
%
% DOA 约定:
%   theta 表示声波传播方向 source -> receiver 相对于阵列 broadside 的水平角。
%   后续 steering_vector.m 应使用:
%       a_m(theta) = exp(-1j * 2*pi/lambda * y_m * sin(theta))

%% 1. 基本字段检查
required_fields = { ...
    'c', 'fc', 'fs', 'T', ...
    'tx_pos', 'rx_pos', 'target_pos', 'H', ...
    'enable_multipath', 'max_order', ...
    'gamma_direct', 'gamma_surface', 'gamma_bottom', ...
    'power_control_mode', 'direct_amp', 'target_amp', ...
    'enable_spreading_loss', 'spreading_type', 'spreading_exponent', ...
    'enable_absorption_loss', 'absorption_dB_per_km', ...
    'array_axis', 'broadside_axis', ...
    'rx_margin_before', 'rx_margin_after'};

for k = 1:numel(required_fields)
    field_name = required_fields{k};
    assert(isfield(cfg, field_name), 'cfg.%s 是必需字段。', field_name);
end

tx = cfg.tx_pos(:).';
rx = cfg.rx_pos(:).';
target = cfg.target_pos(:).';

assert(numel(tx) == 3, 'cfg.tx_pos 必须是 1x3 向量。');
assert(numel(rx) == 3, 'cfg.rx_pos 必须是 1x3 向量。');
assert(numel(target) == 3, 'cfg.target_pos 必须是 1x3 向量。');

%% 2. 初始化路径列表
paths = repmat(empty_path(), 0, 1);
path_id = 0;

%% 3. 直达波路径
path_id = path_id + 1;

direct_distance = norm(rx - tx);
direct_delay = direct_distance / cfg.c;
direct_theta = compute_horizontal_doa(tx, rx, cfg);

direct_gamma = cfg.gamma_direct;
[direct_alpha, direct_loss_info] = compute_interference_alpha( ...
    cfg, direct_gamma, '', direct_distance, direct_delay, direct_distance);

paths(end+1, 1) = make_path( ...
    path_id, ...
    'direct', ...
    '直达波路径', ...
    '', ...
    0, ...
    tx, ...
    tx, ...
    direct_distance, ...
    direct_delay, ...
    direct_theta, ...
    direct_gamma, ...
    direct_alpha, ...
    true, ...
    false, ...
    0, ...
    0, ...
    NaN, ...
    NaN, ...
    direct_loss_info);

%% 4. 海面/海底虚源多径
if cfg.enable_multipath && cfg.max_order > 0
    sequences = generate_reflection_sequences(cfg.max_order);

    for i = 1:numel(sequences)
        seq = sequences{i};

        image_pos = tx;
        gamma = cfg.gamma_direct;
        num_surface = 0;
        num_bottom = 0;

        for q = 1:numel(seq)
            switch seq(q)
                case 'S'
                    image_pos = reflect_surface(image_pos);
                    gamma = gamma * cfg.gamma_surface;
                    num_surface = num_surface + 1;

                case 'B'
                    image_pos = reflect_bottom(image_pos, cfg.H);
                    gamma = gamma * cfg.gamma_bottom;
                    num_bottom = num_bottom + 1;

                otherwise
                    error('未知反射序列符号: %s', seq(q));
            end
        end

        distance = norm(rx - image_pos);
        delay = distance / cfg.c;
        theta = compute_horizontal_doa(image_pos, rx, cfg);
        [alpha, loss_info] = compute_interference_alpha( ...
            cfg, gamma, seq, distance, delay, direct_distance);

        path_id = path_id + 1;
        [path_type, path_name] = multipath_name(seq);

        paths(end+1, 1) = make_path( ...
            path_id, ...
            path_type, ...
            path_name, ...
            seq, ...
            numel(seq), ...
            image_pos, ...
            image_pos, ...
            distance, ...
            delay, ...
            theta, ...
            gamma, ...
            alpha, ...
            true, ...
            false, ...
            num_surface, ...
            num_bottom, ...
            NaN, ...
            NaN, ...
            loss_info);
    end
end

%% 5. 目标双基地回波路径
enable_target = ~isfield(cfg, 'enable_target') || cfg.enable_target;

if enable_target
    tx_to_target = norm(target - tx);
    target_to_rx = norm(rx - target);

    target_distance = tx_to_target + target_to_rx;
    target_delay = target_distance / cfg.c;

    % 接收阵列看到的目标 DOA 由 target -> receiver 的传播方向决定。
    target_theta = compute_horizontal_doa(target, rx, cfg);

    target_gamma = 1.0;
    [target_alpha, target_loss_info] = compute_target_alpha(cfg, target_gamma, ...
        tx_to_target, target_to_rx, target_delay);

    path_id = path_id + 1;

    paths(end+1, 1) = make_path( ...
        path_id, ...
        'target', ...
        '目标双基地回波', ...
        'T', ...
        0, ...
        target, ...
        target, ...
        target_distance, ...
        target_delay, ...
        target_theta, ...
        target_gamma, ...
        target_alpha, ...
        false, ...
        true, ...
        0, ...
        0, ...
        tx_to_target, ...
        target_to_rx, ...
        target_loss_info);
end

%% 6. 按时延排序并重编号
[~, idx] = sort([paths.delay]);
paths = paths(idx);

for k = 1:numel(paths)
    paths(k).id = k;
end

%% 7. 记录路径数量并检查 DOA 搜索范围
cfg.num_paths = numel(paths);
cfg.num_interference_paths = sum([paths.is_interference]);
cfg.num_target_paths = sum([paths.is_target]);

if isfield(cfg, 'angle_min') && isfield(cfg, 'angle_max')
    theta_all = [paths.theta];

    if any(theta_all < cfg.angle_min) || any(theta_all > cfg.angle_max)
        warning('部分路径 DOA 超出 CLEAN 搜索角度范围，请检查 angle_grid 设置。');
    end
end

%% 8. 根据所有路径时延更新接收观测时间轴
all_tau = [paths.delay];

cfg.t_rx_start = max(0, min(all_tau) - cfg.rx_margin_before);
cfg.t_rx_end = max(all_tau) + cfg.T + cfg.rx_margin_after;

cfg.N_rx = ceil((cfg.t_rx_end - cfg.t_rx_start) * cfg.fs) + 1;
cfg.t_rx = cfg.t_rx_start + (0:cfg.N_rx-1).' / cfg.fs;
end

%% ========================================================================
% 局部辅助函数
% ========================================================================

function p = empty_path()
%EMPTY_PATH 返回包含所有预期字段的空路径结构体。

p = struct();
p.id = [];
p.type = '';
p.name = '';
p.sequence = '';
p.order = [];
p.image_pos = [];
p.source_pos_for_doa = [];
p.distance = [];
p.delay = [];
p.theta = [];
p.theta_deg = [];
p.gamma = [];
p.alpha = [];
p.amp = [];
p.power = [];
p.is_interference = [];
p.is_target = [];
p.num_surface_reflect = [];
p.num_bottom_reflect = [];
p.tx_to_target = NaN;
p.target_to_rx = NaN;
p.excess_distance = 0;
p.refl_coeff = 1;
p.refl_abs = 1;
p.refl_phase = 0;
p.spreading_factor = 1;
p.absorption_factor = 1;
p.loss_factor = 1;
p.power_control_mode = '';
end

function p = make_path( ...
    id, type, name, sequence, order, image_pos, source_pos_for_doa, ...
    distance, delay, theta, gamma, alpha, is_interference, is_target, ...
    num_surface, num_bottom, tx_to_target, target_to_rx, loss_info)
%MAKE_PATH 构造一条路径记录。

if nargin < 19 || isempty(loss_info)
    loss_info = default_loss_info();
end

p = empty_path();

p.id = id;
p.type = type;
p.name = name;
p.sequence = sequence;
p.order = order;
p.image_pos = image_pos;
p.source_pos_for_doa = source_pos_for_doa;
p.distance = distance;
p.delay = delay;
p.theta = theta;
p.theta_deg = theta * 180/pi;
p.gamma = gamma;
p.alpha = alpha;
p.amp = abs(alpha);
p.power = abs(alpha)^2;
p.is_interference = is_interference;
p.is_target = is_target;
p.num_surface_reflect = num_surface;
p.num_bottom_reflect = num_bottom;
p.tx_to_target = tx_to_target;
p.target_to_rx = target_to_rx;

p.excess_distance = loss_info.excess_distance;
p.refl_coeff = loss_info.refl_coeff;
p.refl_abs = loss_info.refl_abs;
p.refl_phase = loss_info.refl_phase;
p.spreading_factor = loss_info.spreading_factor;
p.absorption_factor = loss_info.absorption_factor;
p.loss_factor = loss_info.loss_factor;
p.power_control_mode = loss_info.power_control_mode;
end

function sequences = generate_reflection_sequences(max_order)
%GENERATE_REFLECTION_SEQUENCES 生成简化虚源法反射序列。
%
% 对平行海面/海底边界，会议版先使用交替反射序列:
%   order 1: S, B
%   order 2: SB, BS
%   order 3: SBS, BSB
%
% 其中 S 表示 sea surface，B 表示 bottom。

sequences = {};

if max_order <= 0
    return;
end

for order = 1:max_order
    if order == 1
        sequences{end+1} = 'S'; %#ok<AGROW>
        sequences{end+1} = 'B'; %#ok<AGROW>
    else
        seq1 = repmat('S', 1, order);
        seq2 = repmat('B', 1, order);

        for k = 2:order
            if seq1(k-1) == 'S'
                seq1(k) = 'B';
            else
                seq1(k) = 'S';
            end

            if seq2(k-1) == 'B'
                seq2(k) = 'S';
            else
                seq2(k) = 'B';
            end
        end

        sequences{end+1} = seq1; %#ok<AGROW>
        sequences{end+1} = seq2; %#ok<AGROW>
    end
end
end

function p_img = reflect_surface(p)
%REFLECT_SURFACE 关于海面 z = 0 做镜像。

p_img = p;
p_img(3) = -p(3);
end

function p_img = reflect_bottom(p, H)
%REFLECT_BOTTOM 关于海底 z = H 做镜像。

p_img = p;
p_img(3) = 2*H - p(3);
end

function theta = compute_horizontal_doa(source_pos, rx_pos, cfg)
%COMPUTE_HORIZONTAL_DOA 计算相对于阵列 broadside 的水平 DOA。
%
% 约定:
%   theta 使用声波传播方向 source -> receiver。
%   prop_vec = rx_pos - source_pos。
%
% 默认配置:
%   array_axis     = [0, 1, 0]
%   broadside_axis = [1, 0, 0]
%
% 因此 theta = 0 表示沿 +x 方向入射，即 broadside 入射。

prop_vec = rx_pos(:).' - source_pos(:).';
prop_vec(3) = 0;

if norm(prop_vec) < eps
    theta = 0;
    return;
end

array_axis = cfg.array_axis(:).';
broadside_axis = cfg.broadside_axis(:).';

array_axis = array_axis / norm(array_axis);
broadside_axis = broadside_axis / norm(broadside_axis);

x_comp = dot(prop_vec, broadside_axis);
y_comp = dot(prop_vec, array_axis);

theta = atan2(y_comp, x_comp);
end

function [path_type, path_name] = multipath_name(seq)
%MULTIPATH_NAME 根据反射序列生成路径类型和名称。

if numel(seq) == 1 && seq == 'S'
    path_type = 'surface';
    path_name = '海面一次反射路径';
elseif numel(seq) == 1 && seq == 'B'
    path_type = 'bottom';
    path_name = '海底一次反射路径';
else
    path_type = 'multipath';
    path_name = ['多径路径 ', seq];
end
end

function [alpha, loss_info] = compute_interference_alpha( ...
    cfg, gamma, sequence, distance, delay, direct_distance)
%COMPUTE_INTERFERENCE_ALPHA 计算直达波/多径干扰路径复幅度。
%
% received_power 模式:
%   直达波强度由 cfg.direct_amp 控制；
%   多径相对强弱由 gamma 控制。
%
% physical_spreading 模式:
%   保留旧的全物理扩展损耗逻辑。
%
% hybrid_physical_loss 模式:
%   直达波仍由 cfg.direct_amp 控制；多径在直达波幅度基础上叠加
%   反射系数、相对扩展损耗和吸收损耗。
%
% alpha 中包含阵列中心处的载频传播相位 exp(-j*2*pi*fc*tau)。

if nargin < 6 || isempty(direct_distance)
    direct_distance = distance;
end

carrier_phase = exp(-1j * 2*pi*cfg.fc * delay);
distance = max(distance, eps);
direct_distance = max(direct_distance, eps);

switch cfg.power_control_mode
    case 'received_power'
        refl_coeff = gamma;
        spreading_factor = 1;
        absorption_factor = 1;
        amp_complex = cfg.direct_amp * refl_coeff;

    case 'physical_spreading'
        refl_coeff = gamma;
        spreading_factor = spreading_loss(distance, cfg);
        absorption_factor = 1;
        amp_complex = refl_coeff * spreading_factor;

    case 'hybrid_physical_loss'
        refl_coeff = reflection_product(sequence, cfg);

        if cfg.enable_spreading_loss
            spreading_factor = (direct_distance / distance)^cfg.spreading_exponent;
        else
            spreading_factor = 1;
        end

        if cfg.enable_absorption_loss
            excess_km = max(distance - direct_distance, 0) / 1000;
            absorption_factor = 10^(-cfg.absorption_dB_per_km * excess_km / 20);
        else
            absorption_factor = 1;
        end

        amp_complex = cfg.direct_amp * refl_coeff * ...
            spreading_factor * absorption_factor;

    otherwise
        error('不支持的 power_control_mode: %s', cfg.power_control_mode);
end

alpha = amp_complex * carrier_phase;
loss_info = make_loss_info(cfg.power_control_mode, direct_distance, ...
    distance, refl_coeff, spreading_factor, absorption_factor);
end

function [alpha, loss_info] = compute_target_alpha( ...
    cfg, gamma, tx_to_target, target_to_rx, delay)
%COMPUTE_TARGET_ALPHA 计算目标双基地回波复幅度。
%
% hybrid_physical_loss 模式下目标幅度仍由 cfg.target_amp 控制，方便
% target_intensity_sweep 中用 Ae 明确控制目标强度。

carrier_phase = exp(-1j * 2*pi*cfg.fc * delay);
target_distance = tx_to_target + target_to_rx;

switch cfg.power_control_mode
    case {'received_power', 'hybrid_physical_loss'}
        refl_coeff = gamma;
        spreading_factor = 1;
        absorption_factor = 1;
        amp_complex = cfg.target_amp * refl_coeff;

    case 'physical_spreading'
        % 双基地目标回波按 tx -> target 和 target -> rx 两段传播损耗处理。
        loss_tx = spreading_loss(tx_to_target, cfg);
        loss_rx = spreading_loss(target_to_rx, cfg);
        refl_coeff = gamma;
        spreading_factor = loss_tx * loss_rx;
        absorption_factor = 1;
        amp_complex = refl_coeff * spreading_factor;

    otherwise
        error('不支持的 power_control_mode: %s', cfg.power_control_mode);
end

alpha = amp_complex * carrier_phase;
loss_info = make_loss_info(cfg.power_control_mode, NaN, target_distance, ...
    refl_coeff, spreading_factor, absorption_factor);
end

function refl_coeff = reflection_product(sequence, cfg)
%REFLECTION_PRODUCT 根据反射序列计算总反射系数。
%
% 这里继续复用 cfg.gamma_surface / cfg.gamma_bottom，避免引入新的
% 反射系数字段后发生重复标定。sequence 为空时返回直达路径系数。

refl_coeff = cfg.gamma_direct;

for k = 1:numel(sequence)
    switch sequence(k)
        case 'S'
            refl_coeff = refl_coeff * cfg.gamma_surface;

        case 'B'
            refl_coeff = refl_coeff * cfg.gamma_bottom;

        otherwise
            error('未知反射序列符号: %s', sequence(k));
    end
end
end

function loss_info = default_loss_info()
%DEFAULT_LOSS_INFO 构造默认损耗信息字段。

loss_info = struct();
loss_info.excess_distance = 0;
loss_info.refl_coeff = 1;
loss_info.refl_abs = 1;
loss_info.refl_phase = 0;
loss_info.spreading_factor = 1;
loss_info.absorption_factor = 1;
loss_info.loss_factor = 1;
loss_info.power_control_mode = '';
end

function loss_info = make_loss_info(power_control_mode, direct_distance, ...
    distance, refl_coeff, spreading_factor, absorption_factor)
%MAKE_LOSS_INFO 汇总路径损耗调试信息。

loss_info = default_loss_info();

if nargin < 2 || ~isfinite(direct_distance)
    direct_distance = NaN;
end

if nargin < 3 || ~isfinite(distance)
    distance = NaN;
end

if nargin < 4 || isempty(refl_coeff)
    refl_coeff = 1;
end

if nargin < 5 || isempty(spreading_factor)
    spreading_factor = 1;
end

if nargin < 6 || isempty(absorption_factor)
    absorption_factor = 1;
end

if isfinite(direct_distance) && isfinite(distance)
    loss_info.excess_distance = max(distance - direct_distance, 0);
else
    loss_info.excess_distance = 0;
end

loss_info.refl_coeff = refl_coeff;
loss_info.refl_abs = abs(refl_coeff);
loss_info.refl_phase = angle(refl_coeff);
loss_info.spreading_factor = spreading_factor;
loss_info.absorption_factor = absorption_factor;
loss_info.loss_factor = abs(refl_coeff) * ...
    abs(spreading_factor) * abs(absorption_factor);
loss_info.power_control_mode = power_control_mode;
end

function loss = spreading_loss(distance, cfg)
%SPREADING_LOSS 简化扩展损耗幅度模型。
%
% 该函数主要服务于 physical_spreading 模式。
% 会议版主实验建议使用 power_control_mode = 'received_power'。

if ~cfg.enable_spreading_loss
    loss = 1;
    return;
end

distance = max(distance, eps);

switch cfg.spreading_type
    case 'none'
        loss = 1;

    case 'cylindrical'
        loss = 1 / sqrt(distance);

    case 'spherical'
        loss = 1 / distance;

    otherwise
        error('不支持的 spreading_type: %s', cfg.spreading_type);
end
end
