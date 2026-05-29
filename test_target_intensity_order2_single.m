clear; clc; close all;

%TEST_TARGET_INTENSITY_ORDER2_SINGLE
% 二阶多径 + hybrid physical loss 下的目标强度单点测试。
%
% 这个脚本对应 run_target_intensity_sweep.m 的单次 sanity check:
%   固定直达波功率
%   固定噪声功率
%   通过 Ae 改变目标回波功率
%   开启二阶多径 max_order = 2
%
% Ae 定义:
%   Ae = 10*log10(P_target / P_direct)

%% 1. 基础配置
cfg = sonar_config();

%% 2. 单点目标强度设置
% 可改成 -33, -30, -27, -24, -21, -18, -15 等。
Ae_dB = -30;

% 固定直达波功率和噪声功率，和 run_target_intensity_sweep.m 保持一致。
direct_power_ref = cfg.target_power * 10^(cfg.INR/10);
noise_power_ref = cfg.noise_power;

target_power = direct_power_ref * 10^(Ae_dB/10);
INR_dB = -Ae_dB;
SNR_dB = 10 * log10(target_power / noise_power_ref);

%% 3. 覆盖参数：二阶多径 + 物理传播损耗环境
overrides = struct();
overrides.power_control_mode = 'hybrid_physical_loss';
overrides.max_order = 2;
overrides.target_power = target_power;
overrides.INR = INR_dB;
overrides.SNR = SNR_dB;

%% 4. 单次运行选项
opts = struct();
opts.verbose = true;
opts.store_signals = false;
opts.store_clean_data = false;
opts.run_baseline = true;
opts.run_calibrated = true;
opts.run_oracle = true;

%% 5. 运行一次完整实验
mc_index = 1;
result = main_single_run(cfg, mc_index, overrides, opts);

%% 6. 打印路径损耗表
fprintf('\nOrder-2 hybrid physical loss path table:\n');
fprintf('------------------------------------------------------------------------------------------\n');
fprintf('%-3s %-12s %-6s %-12s %-10s %-10s %-10s %-10s %-10s\n', ...
    'ID', 'Type', 'Seq', 'Distance', 'Amp', 'ReflAbs', 'Spread', 'Absorp', 'Loss');
fprintf('------------------------------------------------------------------------------------------\n');

for k = 1:numel(result.paths)
    fprintf('%-3d %-12s %-6s %-12.3f %-10.4g %-10.4g %-10.5f %-10.5f %-10.5f\n', ...
        result.paths(k).id, ...
        result.paths(k).type, ...
        result.paths(k).sequence, ...
        result.paths(k).distance, ...
        result.paths(k).amp, ...
        result.paths(k).refl_abs, ...
        result.paths(k).spreading_factor, ...
        result.paths(k).absorption_factor, ...
        result.paths(k).loss_factor);
end

fprintf('------------------------------------------------------------------------------------------\n');

%% 7. 打印目标强度和核心指标
fprintf('\nTarget intensity single-point check:\n');
fprintf('Ae                         = %.3f dB\n', Ae_dB);
fprintf('target_power               = %.6g\n', target_power);
fprintf('INR                        = %.3f dB\n', INR_dB);
fprintf('SNR                        = %.3f dB\n', SNR_dB);
fprintf('power_control_mode         = %s\n', result.cfg.power_control_mode);
fprintf('max_order                  = %d\n', result.cfg.max_order);
fprintf('num_paths                  = %d\n', numel(result.paths));

fprintf('\nMetrics table:\n');
disp(result.metrics.table);

fprintf('\nSummary:\n');
disp(result.metrics.summary);

%% 8. 基本一致性检查
assert(result.ok == true, 'main_single_run 未正常完成。');
assert(result.cfg.max_order == 2, 'max_order 没有成功覆盖为 2。');
assert(strcmp(result.cfg.power_control_mode, 'hybrid_physical_loss'), ...
    'power_control_mode 没有成功覆盖为 hybrid_physical_loss。');
assert(numel(result.paths) >= 6, ...
    '二阶多径场景下应至少包含 direct, S, B, SB, BS, target。');
assert(isfinite(result.metrics.summary.calibrated_vs_baseline_pdr_gain_dB), ...
    'PDR gain 不是有限数。');

fprintf('\nOrder-2 target intensity single-point check passed.\n');
