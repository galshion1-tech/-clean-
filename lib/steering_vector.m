function [a, info] = steering_vector(theta, cfg, phase_error)
%STEERING_VECTOR 生成 ULA 水平 DOA 导向矢量。
%
% 用法:
%   [a, info] = steering_vector(theta, cfg);
%   [a, info] = steering_vector(theta, cfg, phase_error);
%
% 输入:
%   theta       水平 DOA，rad。可以是标量，也可以是向量。
%               约定为声波传播方向 source -> receiver 相对于阵列 broadside 的角度。
%
%   cfg         sonar_config() 生成的配置结构体。
%
%   phase_error 可选，M x 1 阵列通道相位误差，rad。
%               若提供，则生成带相位失配的真实导向矢量:
%                   a_true = diag(exp(1j*phase_error)) * a_nominal
%
% 输出:
%   a           M x K 导向矢量矩阵。
%               M 为阵元数，K 为 theta 的个数。
%               每一列对应一个 DOA。
%
%   info        导向矢量元信息。
%
% 符号约定:
%   theta 使用声波传播方向 source -> receiver。
%   默认 ULA 沿 y 轴排列，broadside 为 +x 方向。
%   因此导向矢量采用:
%       a_m(theta) = exp(-1j * k * y_m * sin(theta))

%% 1. 输入检查
if nargin < 3
    phase_error = [];
end

required_fields = {'M', 'array_pos', 'array_center', ...
    'array_axis', 'broadside_axis'};

for kf = 1:numel(required_fields)
    field_name = required_fields{kf};
    assert(isfield(cfg, field_name), 'cfg.%s 是必需字段。', field_name);
end

theta = theta(:).';             % 统一为 1 x K 行向量。
assert(isreal(theta) && all(isfinite(theta)), ...
    'theta 必须是有限实数。');

if isfield(cfg, 'lambda')
    lambda = cfg.lambda;
else
    assert(isfield(cfg, 'c') && isfield(cfg, 'fc'), ...
        'cfg.lambda 不存在时，必须提供 cfg.c 和 cfg.fc。');
    lambda = cfg.c / cfg.fc;
end

assert(lambda > 0, 'lambda 必须为正数。');

M = cfg.M;
K_theta = numel(theta);

assert(size(cfg.array_pos, 1) == M, ...
    'cfg.array_pos 的行数必须等于 cfg.M。');
assert(size(cfg.array_pos, 2) == 3, ...
    'cfg.array_pos 必须是 M x 3 矩阵。');
assert(numel(cfg.array_center) == 3, ...
    'cfg.array_center 必须是 1 x 3 向量。');
assert(numel(cfg.array_axis) == 3, ...
    'cfg.array_axis 必须是 1 x 3 向量。');
assert(numel(cfg.broadside_axis) == 3, ...
    'cfg.broadside_axis 必须是 1 x 3 向量。');

%% 2. 阵列轴与 broadside 轴归一化
array_axis = cfg.array_axis(:);
broadside_axis = cfg.broadside_axis(:);

assert(norm(array_axis) > 0, 'cfg.array_axis 不能为零向量。');
assert(norm(broadside_axis) > 0, 'cfg.broadside_axis 不能为零向量。');

array_axis = array_axis / norm(array_axis);
broadside_axis = broadside_axis / norm(broadside_axis);

assert(abs(dot(array_axis, broadside_axis)) < 1e-12, ...
    'array_axis 与 broadside_axis 应当正交。');

%% 3. 计算阵元相对坐标
% pos_rel: M x 3，每个阵元相对于阵列中心的位置。
% element_coord: M x 1，阵元在 array_axis 方向上的投影坐标。
%
% 对默认 ULA 而言，element_coord 就是相对于阵列中心的 y 坐标。
array_center = cfg.array_center(:).';
pos_rel = cfg.array_pos - repmat(array_center, M, 1);

element_coord = pos_rel * array_axis;

%% 4. 生成名义导向矢量
k0 = 2*pi / lambda;

% phase_matrix 为 M x K，每一列对应一个 theta。
% 注意负号来自当前 DOA 约定: theta 是 source -> receiver 的传播方向。
phase_matrix = -1j * k0 * element_coord * sin(theta);
a = exp(phase_matrix);

%% 5. 可选加入阵列通道相位误差
if ~isempty(phase_error)
    phase_error = phase_error(:);

    assert(numel(phase_error) == M, ...
        'phase_error 的长度必须等于 cfg.M。');
    assert(isreal(phase_error) && all(isfinite(phase_error)), ...
        'phase_error 必须是有限实数。');

    % 带相位失配的真实导向矢量:
    %   a_true = diag(exp(1j*phi)) * a_nominal
    phase_gain = exp(1j * phase_error);

    % 使用 bsxfun 兼容旧版 MATLAB。
    a = bsxfun(@times, phase_gain, a);
end

%% 6. 输出信息
info = struct();
info.M = M;
info.num_angles = K_theta;
info.theta = theta;
info.theta_deg = theta * 180/pi;
info.lambda = lambda;
info.k0 = k0;
info.array_axis = array_axis(:).';
info.broadside_axis = broadside_axis(:).';
info.pos_rel = pos_rel;
info.element_coord = element_coord;
info.phase_error_applied = ~isempty(phase_error);
info.normalization = 'elementwise_unit_magnitude';
info.sign_convention = 'exp(-1j*k*y_rel*sin(theta))';
info.doa_convention = 'theta is propagation direction source-to-receiver relative to broadside';
end
