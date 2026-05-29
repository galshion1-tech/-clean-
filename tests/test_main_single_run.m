% 添加项目目录到 MATLAB 路径
root_dir = fileparts(fileparts(mfilename('fullpath')));
addpath(root_dir);
addpath(fullfile(root_dir, 'lib'));
addpath(fullfile(root_dir, 'adaptive'));

cfg = sonar_config();

overrides = struct();
overrides.power_control_mode = 'hybrid_physical_loss';
overrides.max_order = 2;

opts = struct();
opts.verbose = true;
opts.run_oracle = true;

result = main_single_run(cfg, 1, overrides, opts);
