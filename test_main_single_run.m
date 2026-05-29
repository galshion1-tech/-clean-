cfg = sonar_config();

overrides = struct();
overrides.power_control_mode = 'hybrid_physical_loss';
overrides.max_order = 2;

opts = struct();
opts.verbose = true;
opts.run_oracle = true;

result = main_single_run(cfg, 1, overrides, opts);
