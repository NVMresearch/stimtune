function [T_factorial,T_designs,onset_tables,eff,opt,con] = optimise_design(settings)
% Optimizes an experimental design for fMRI.
%
% An experimental design is specified by blocks of a certain duration and 
% jittered inter-trial intervals (ITIs), sampled from a uniform 
% distribution between min_gap_range and max_gap_range. This tool performs 
% a grid search to select the block length and range of ITIs that maximize 
% the efficiency for each contrast.
%
% Returns:
% T_factorial  - table with the optimal settings for each contrast
% T_designs    - table with the max efficiency (over contrasts) for each 
%                candidate design
% onset_tables - tables with names, onsets, durations
% eff          - raw efficiencies arranged in a matrix 
%                [designs x contrasts x repetitions]
% opt{c}       - metadata for the optimal design for each contrast c
% con(c)       - definition of each contrast c
% _________________________________________________________________________

% Unpack inputs
on_durations    = settings.durations_on;
iti_dist        = settings.iti_dist;
maxtime         = settings.maxtime;
nrep            = settings.nrep;
seed            = settings.seed;
factors         = settings.factors;
nbf             = settings.nbf;
TR              = settings.TR;

rng(seed);

% Generate table of designs to evaluate. A "design" consists of a duration
% and parameters for the ITI distribution
% -------------------------------------------------------------------------

% Build a cell array param_ranges{i} is the vector of candidate values for
% parameter i
param_ranges = [on_durations; iti_dist.param_ranges(:)];
n_params     = length(param_ranges);
n_iti_params = n_params-1;

% Generate all permutations (requires MATLAB 2023a)
% Each row is a design and each column a parameter (duration followed by
% ITI distribution parameters)
param_names = ['Duration', cellstr(compose("ITI_Param_%d",1:n_iti_params))];
params = combinations(param_ranges{:});
params.Properties.VariableNames = param_names;
T_designs = table(params);

% Remove invalid designs
% -------------------------------------------------------------------------
to_delete = [];

switch iti_dist.distribution
    case 'uniform'
        to_delete = T_designs.params.ITI_Param_1 > T_designs.params.ITI_Param_2;        
end

T_designs(to_delete, :) = [];
ndesigns = size(T_designs,1);

% Generate list of experimental conditions (cells in the factorial design)
% -------------------------------------------------------------------------
[names,nlevels] = gen_experimental_conditions(factors); 

nconditions = length(names);

% Prepare contrasts
% -------------------------------------------------------------------------
con = spm_make_contrasts(nlevels);

% Remove any contrasts if undefined (e.g. for a single condition design)
to_remove = arrayfun(@(contrast) isempty(contrast.c), con);
con(to_remove) = [];

% Adapt the contrasts for the number of HRF basis functions
for i = 1:length(con)
    con(i).c = kron(con(i).c, eye(nbf));
end

ncontrasts = length(con);

% Main loop
% -------------------------------------------------------------------------
eff = zeros(ndesigns, ncontrasts, nrep);
max_eff = repmat(-999,1,ncontrasts);
opt = cell(1,ncontrasts);

for i = 1:ndesigns

    if mod(i,10) == 0
        fprintf('Evaluating design %d of %d\n',i,ndesigns);
    end

    % Get design parameters
    params = T_designs{i,'params'};
    duration = params.Duration;
    iti_sampler = get_iti_distribution(iti_dist.distribution, params);
    
    % Repeatedly randomize trial order and assess efficiency
    for rep = 1:nrep

        % Make a table of trials
        all_trials = make_schedule(duration, nconditions, iti_sampler, maxtime);

        % Skip if design was not feasible
        if isempty(all_trials)
            warning('Design was not feasible');
            continue
        end

        % Randomize the trial order [condition, duration, ITI]
        ntrials    = size(all_trials,1);
        idx_rnd    = randperm(ntrials);
        trials_rnd = all_trials(idx_rnd,:);

        % Get onsets and durations in SPM format
        [onsets,durations] = generate_spm_onsets(trials_rnd);

        % Calculate required number of scans
        nscans = ceil(maxtime / TR);

        % Generate design matrix
        pmod = [];
        [X,Xf,Xnames,U] = generate_spm(names, onsets, durations, pmod, nscans, TR, nbf);
    
        % Assess design's efficiency for each contrast
        for c = 1:length(con)

            % Get contrast vector or matrix
            cvec = con(c).c';

            % Evaluate efficiency
            eff(i,c,rep) = compute_glm_efficiency(X,cvec);

            % Store the best design for each contrast
            if eff(i,c,rep) > max_eff(c)

                % Details of the design
                design = struct();
                design.params = params;
                design.eff    = eff(i,c,rep);
                
                % Details of the specific best repetition of this design
                meta = struct();
                meta.names      = names;
                meta.onsets     = onsets;
                meta.durations  = durations;
                meta.U          = U;
                meta.design     = design;
                meta.settings   = settings;
                opt{c} = meta;
                
                % Overall best efficiency
                max_eff(c) = eff(i,c,rep);
            end
        end
    end
end
fprintf('Evaluated %d designs\n',ndesigns);

% Build a table with the optimal design parameters for each contrast
% -------------------------------------------------------------------------
contrast      = {con.name}';
params        = cell2mat(cellfun(@(x)x.design.params,opt,'UniformOutput',false)');
efficiency    = cellfun(@(x)x.design.eff,opt)';

T_factorial = table(contrast, params, efficiency);

% Build a table with all individual repetitions of designs that were
% assessed, sorted by total efficiency across contrasts
% -------------------------------------------------------------------------

% variable eff has dimension [designs x contrasts x repetitions]. Take the
% max over the third dimension (repetitions)
eff_max = max(eff, [], 3);

% temporary table for efficiencies
T_eff = table();

% Add columns with the max efficiency of each contrast
for i = 1:ncontrasts
    fieldname = sprintf('eff_%s',contrast{i});
    T_eff.(fieldname) = eff_max(:,i);
end

% Total efficiency over contrasts for each design
eff_total = sum(eff_max,2);
T_eff.eff_total = eff_total;

% Append to T_designs table
T_designs.efficiency = T_eff;

% Sort designs by eff_total
[~, sortOrder] = sort(eff_total, 'descend');
T_designs = T_designs(sortOrder,:);

% Sort raw eff matrix by the same
eff = eff(sortOrder,:,:);

% Generate and return a set of repetitions of the best design per
% contrast. We don't store these for every design (above) to save memory
% -------------------------------------------------------------------------
onset_tables = cell(1,length(con));

repetition = (1:nrep)';

for c = 1:length(con)
    % Get the parameters of the winning design
    params      = opt{c}.design.params;
    duration    = params.Duration;
    iti_sampler = get_iti_distribution(iti_dist.distribution, params);

    % Get contrast vector or matrix
    cvec = con(c).c';

    % Output variables
    all_names       = cell(nrep,1);
    all_onsets      = cell(nrep,1);
    all_durations   = cell(nrep,1);
    all_itis        = cell(nrep,1);
    all_efficiency  = zeros(nrep,1);    
    
    for rep = 1:nrep

        % Make a table of trials
        all_trials = make_schedule(duration, nconditions, iti_sampler, maxtime);

        % Skip if design was not feasible
        if isempty(all_trials)
            warning('Design was not feasible');
            continue
        end

        % Randomize the trial order [condition, duration, ITI]
        ntrials    = size(all_trials,1);
        idx_rnd    = randperm(ntrials);
        trials_rnd = all_trials(idx_rnd,:);

        % Get onsets and durations in SPM format
        [onsets,durations] = generate_spm_onsets(trials_rnd);

        % Generate design matrix
        pmod = [];
        [X,Xf,Xnames,U] = generate_spm(names, onsets, durations, pmod, nscans, TR, nbf);

        % Evaluate efficiency
        all_efficiency(rep)  = compute_glm_efficiency(X,cvec);
        
        % Store names, onsets & durations
        all_names{rep}       = names;
        all_onsets{rep}      = onsets;
        all_durations{rep}   = durations;  
        all_itis{rep}        = trials_rnd.ITI;
    end
    onset_tables{c} = table(repetition,all_names,all_onsets,all_durations,all_itis,all_efficiency,...
        'VariableNames',{'Repetition','names','onsets','durations','ITIs','efficiency'});
end

% Sort repetitions by efficiency
for c = 1:numel(onset_tables)
    onset_tables{c} = sortrows(onset_tables{c}, 'efficiency', 'descend');
end

% -------------------------------------------------------------------------
function iti_sampler = get_iti_distribution(distribution, params)
% Returns a handle to a function for sampling a single ITI
% distribution - name of the distribution
% params       - array of parameters

switch distribution
    case 'uniform'
        % Uniform distribution bounded at a,b
        uniformrnd  = @(a,b) a + (b - a) * rand(1,1);
        iti_sampler = @()uniformrnd(params.ITI_Param_1,params.ITI_Param_2);
    case 'gamma'
        % Gamma random distribution with shape and scale parameters a,b
        % lower bounded by minval
        try
            minval = params.ITI_Param_3;
        catch
            minval = 0;
        end
        gammrnd     = @(a,b) spm_gamrnd(a,b);
        iti_sampler = @() minval + gammrnd(params.ITI_Param_1,params.ITI_Param_2); 
    otherwise
        error('Unknown distribution');
end

% -------------------------------------------------------------------------
function [names,nlevels] = gen_experimental_conditions(factors)
% Generates a list of experimental conditions (cells of the
% factorial design)
% 
% names   - names of the conditions
% nlevels - number of levels per factor.  nlevels(i) is the number of 
%           levels for factor i.
nfactors = length(factors);

% Number of levels per factor
nlevels = [factors.nlevels];

% Create inputs to MATLAB function ndgrid
grids = cell(1, nfactors);
for i = 1:nfactors
    grids{i} = 1:nlevels(nfactors - i + 1);
end

% Generate full grid
[grids{:}] = ndgrid(grids{:});

% Create a "conditions" matrix, where each row is a condition (i.e., a cell
% of the factorial design). conditions(i,j)=k indicates that condition i has
% level k of factor j. Importantly, factors that change more slowly over 
% levels appear first, for compatibility with spm_make_contrasts.m
conditions = zeros(numel(grids{1}), nfactors);
for i = 1:nfactors
    conditions(:, i) = grids{nfactors - i + 1}(:);
end
nconditions = size(conditions,1);

% Generate names
names = cell(nconditions,1);
for i = 1:nconditions
    str = sprintf('Condition %d: ',i);
    for j = 1:nfactors
        str = [str sprintf('L%d of %s. ',conditions(i,j),factors(j).name)];
    end
    names{i} = str;
end

% -------------------------------------------------------------------------
function [onsets,durations] = generate_spm_onsets(trials)
% Converts from the trials matrix format used here: 
% [condition, duration, iti] to standard SPM format for onsets and
% durations

% Unpack
cond = trials.Condition;
dur  = trials.Duration;
iti  = trials.ITI;

nconditions = max(cond);
ntrials     = size(trials,1);

% Calculate onsets
ons = cumsum([0; (dur+iti)]);
ons = ons(1:ntrials);

% Convert to standard SPM names, onsets, durations format
onsets    = cell(1,nconditions);
durations = cell(1,nconditions);
for c = 1:nconditions
    is_relevant_trial = (cond == c);
    onsets{c}    = ons(is_relevant_trial);
    durations{c} = dur(is_relevant_trial);
end
