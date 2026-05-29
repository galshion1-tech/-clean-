function [s, mf, info] = generate_lfm(cfg)
%GENERATE_LFM 生成 LFM 发射信号、匹配滤波器和波形信息。
%
% 用法:
%   [s, mf, info] = generate_lfm(cfg);
%
% 输入:
%   cfg.fs            采样率，Hz
%   cfg.B             LFM 带宽，Hz
%   cfg.T             脉冲宽度，s
%   cfg.N_tx          发射脉冲采样点数
%   cfg.t_tx          发射波形时间轴
%   cfg.signal_model  'complex_baseband' 或 'real_passband'
%   cfg.fc            载频，Hz；仅 real_passband 生成时直接使用
%   cfg.c             声速，m/s；用于计算路径长度分辨率
%
% 输出:
%   s     发射信号，列向量
%   mf    匹配滤波器，列向量，已按信号能量归一化
%   info  波形元信息结构体

%% 1. 基本字段检查
required_fields = {'fs', 'B', 'T', 'N_tx', 't_tx', ...
    'signal_model', 'fc', 'c'};

for k = 1:numel(required_fields)
    field_name = required_fields{k};
    assert(isfield(cfg, field_name), 'cfg.%s 是必需字段。', field_name);
end

fs = cfg.fs;
B = cfg.B;
T = cfg.T;
N = cfg.N_tx;
t = cfg.t_tx(:);

assert(numel(t) == N, 'cfg.t_tx 的长度必须等于 cfg.N_tx。');
assert(B > 0, 'cfg.B 必须为正数。');
assert(T > 0, 'cfg.T 必须为正数。');
assert(fs > 0, 'cfg.fs 必须为正数。');

%% 2. 采样率检查
switch cfg.signal_model
    case 'complex_baseband'
        assert(fs >= B, '复基带 LFM 要求 fs >= B。');

    case 'real_passband'
        assert(fs >= 2 * (cfg.fc + B/2), ...
            '实数带通 LFM 要求 fs >= 2*(fc+B/2)。');

    otherwise
        error('signal_model 必须是 ''complex_baseband'' 或 ''real_passband''。');
end

%% 3. 生成复基带 LFM
% LFM 调频率:
%   K = B / T
%
% 使用居中时间轴 t0 = t - T/2 后，复基带瞬时频率为:
%   f_inst(t) = K * t0
%
% 因此频率从 -B/2 扫到 +B/2。
K = B / T;
t0 = t - T/2;
s_bb = exp(1j * pi * K * t0.^2);

%% 4. 根据信号模型选择输出
switch cfg.signal_model
    case 'complex_baseband'
        % 复基带模型下不要乘 exp(j*2*pi*fc*t)。
        % fc 只用于波长、导向矢量和传播相位等模块。
        s = s_bb;

    case 'real_passband'
        % 实数带通信号主要用于对比或可视化。
        s = real(s_bb .* exp(1j * 2*pi*cfg.fc*t));
end

%% 5. 发射信号功率归一化
% 归一化到脉冲内单位平均功率:
%   mean(abs(s).^2) = 1
%
% 这样后续路径幅度 direct_amp、target_amp 等更容易解释。
s = s / sqrt(mean(abs(s).^2));

%% 6. 匹配滤波器
% 常规匹配滤波器:
%   h[n] = conj(s[N-1-n])
%
% 这里除以信号能量 Es，使无噪声单路径情况下匹配滤波峰值约等于
% 路径复幅度 alpha，而不是 alpha * Es。
Es = sum(abs(s).^2);
mf = conj(flipud(s)) / Es;

%% 7. 波形元信息
info = struct();
info.fs = fs;
info.B = B;
info.T = T;
info.N = N;
info.K = K;
info.Es = Es;
info.time_bandwidth_product = B * T;
info.delay_resolution = 1 / B;
info.path_length_resolution = cfg.c / B;
info.signal_model = cfg.signal_model;
info.normalization = 'unit_average_power';
info.matched_filter_normalization = 'unit_peak_gain';
end
