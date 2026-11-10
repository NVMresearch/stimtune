% Timing options
on_durations  = 1:24;                            % Range of fixed on durations for each condition (seconds)
min_gap_range = 2:2:32;                          % Minimum gap duration (seconds)
max_gap_range = 2:2:32;                          % Maximum gap duration (seconds)

% Specify design
design = struct();
design.nrep          = 20; % number of random trial orders per design
design.seed          = 1; % seed for random number generator
design.nbf           = 1; % number of haemodynamic basis functions
design.maxtime       = 10*60; % 10 mins in seconds
design.durations_on  = on_durations;
design.max_gap_range = max_gap_range;
design.min_gap_range = min_gap_range;
design.TR            = 3.22; % volume repetition time in seconds
design.strict        = true; % whether to be strict about the time limit

% Add the StimTune toolbox to the MATLAB path
toolbox_dir = fullfile(pwd,'../toolbox');
addpath(toolbox_dir);

% Select where to store results
results_dir = fullfile(pwd,'stimtune_results');
mkdir(results_dir);
%% Evaluate a single condition design (vs baseline)

% Specify factor structure
factors = struct();
factors(1).name = 'Factor 1';
factors(1).nlevels = 1;

design.factors = factors;

[T_factorial,T_designs,eff,opt,con] = optimise_design(design);

fn = fullfile(results_dir,'opt_design_single_condition.mat');
save(fn,'T_factorial','T_designs','eff','opt','con');

report_optimise_design_results(T_factorial,T_designs,eff,opt,con);
%% Evaluate a 1-factor design with 2 levels

% Specify factor structure
factors = struct();
factors(1).name = 'Factor 1';
factors(1).nlevels = 2;

design.factors = factors;

[T_factorial,T_designs,eff,opt,con] = optimise_design(design);

fn = fullfile(results_dir,'opt_design_1factor.mat');
save(fn,'T_factorial','T_designs','eff','opt','con');

report_optimise_design_results(T_factorial,T_designs,eff,opt,con);

%% Evaluate a 2x2 design 

% Specify factor structure
factors = struct();
factors(1).name = 'Factor 1';
factors(1).nlevels = 2;
factors(2).name = 'Factor 2';
factors(2).nlevels = 2;

design.factors = factors;

[T_factorial,T_designs,eff,opt,con]= optimise_design(design);

fn = fullfile(results_dir,'opt_design_2by2.mat');
save(fn,'T_factorial','T_designs','eff','opt','con');

report_optimise_design_results(T_factorial,T_designs,eff,opt,con);
%% Now try a 2x2 design but with 3 HRF basis functions

% Specify factor structure
factors = struct();
factors(1).name = 'Factor 1';
factors(1).nlevels = 2;
factors(2).name = 'Factor 2';
factors(2).nlevels = 2;

design_3bf         = design;
design_3bf.factors = factors;
design_3bf.nbf     = 3;

[T_factorial,T_designs,eff,opt,con]= optimise_design(design_3bf);

fn = fullfile(results_dir,'opt_design_2by2_3bf.mat');
save(fn,'T_factorial','T_designs','eff','opt','con');

report_optimise_design_results(T_factorial,T_designs,eff,opt,con);
