function [T_factorial,T_designs,eff,opt,con] = optimise_design(settings)
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
% eff          - raw efficiencies arranged in a matrix 
%                [designs x contrasts x repetitions]
% opt{c}       - metadata for the optimal design for each contrast c
% con(c)       - definition of each contrast c
% _________________________________________________________________________

% Unpack inputs
on_durations    = settings.durations_on;
max_gap_range   = settings.max_gap_range;
min_gap_range   = settings.min_gap_range;
maxtime         = settings.maxtime;
nrep            = settings.nrep;
seed            = settings.seed;
factors         = settings.factors;
nbf             = settings.nbf;
TR              = settings.TR;
strict          = settings.strict;

rng(seed);

% Generate table of designs to evaluate. A "design" consists of a duration,
% min_gap and max_gap.
% -------------------------------------------------------------------------
[duration,min_gap,max_gap] = ndgrid(on_durations,min_gap_range,max_gap_range);

% Vectors -> table
duration = duration(:);
min_gap  = min_gap(:);
max_gap  = max_gap(:);
designs  = table(duration,min_gap,max_gap);

% Remove designs with max_gap > min_gap
to_delete = designs.min_gap > designs.max_gap;
designs(to_delete, :) = [];

ndesigns = size(designs,1);

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
    duration = designs(i,:).duration;
    min_gap  = designs(i,:).min_gap;
    max_gap  = designs(i,:).max_gap; 

    % Make a matrix of trials for this design
    all_trials = make_trial_matrix(duration,min_gap,max_gap,nconditions,settings);

    % Calculate required number of scans
    if strict
        nscans = ceil(maxtime / TR);
    else    
        nscans = ceil( max(maxtime,t) / TR);
    end
    
    % Repeatedly randomize trial order and assess efficiency
    for rep = 1:nrep
        
        % Randomize the trial order [condition, duration, ITI]
        ntrials    = size(all_trials,1);
        idx_rnd    = randperm(ntrials);
        trials_rnd = all_trials(idx_rnd,:);

        % Get onsets and durations in SPM format
        [onsets,durations] = generate_spm_onsets(trials_rnd);

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
                design.duration = duration;
                design.min_gap  = min_gap;
                design.max_gap  = max_gap; 
                design.eff      = eff(i,c,rep);
                
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
stim_duration = cellfun(@(x)x.design.duration,opt)';
min_gap       = cellfun(@(x)x.design.min_gap,opt)';
max_gap       = cellfun(@(x)x.design.max_gap,opt)';
efficiency    = cellfun(@(x)x.design.eff,opt)';
T_factorial  = table(contrast, stim_duration, min_gap, max_gap, efficiency);

% Build a table with all individual repetitions of designs that were
% assessed, sorted by total efficiency across contrasts
% -------------------------------------------------------------------------

% variable eff has dimension [designs x contrasts x repetitions]. Take the
% max over the third dimension (repetitions)
eff_max = max(eff, [], 3);

% Create output table
T_designs = designs;

% Add columns with the max efficiency of each contrast
for i = 1:ncontrasts
    fieldname = sprintf('eff_%s',contrast{i});
    T_designs.(fieldname) = eff_max(:,i);
end

% Total efficiency over contrasts for each design
eff_total = sum(eff_max,2);
T_designs.eff_total = eff_total;

% Now sort designs by eff_total
[~, sortOrder] = sort(eff_total, 'descend');
T_designs = T_designs(sortOrder,:);

% Sort raw eff matrix by the same
eff = eff(sortOrder,:,:);

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
function all_trials = make_trial_matrix(duration,min_gap,max_gap,nconditions,settings)
% Generates a list of trials
%
% duration    - [1 x 1] duration of trials
% min_gap     - [1 x 1] minimum gap between trials
% max_gap     - [1 x 1] maximum gap between trials
% nconditions - [1 x 1] number of cells in the factorial design
% settings    - (array) settings
% 
% Returns:
% all_trials - [ntrials x 3] with columns: condition, duration, ITI

% Unpack settings
maxtime = settings.maxtime;
strict  = settings.strict;

% Select ITIs to use
nitis = 6;
itis  = linspace(min_gap,max_gap, nitis);

% Prepare trials
t = 0;        
iti_idx = 1;
all_trials = [];    
while t < maxtime                
    % Make a table with one trial of each condition. Each row
    % contains: [condition number, duration, ITI]
    cond   = (1:nconditions)';
    dur    = repmat(duration,nconditions,1);
    iti    = repmat(itis(iti_idx),nconditions,1);            
    trials = [cond, dur, iti];

    % Advance the choice of ITI for the next batch of trials
    if iti_idx == nitis
        iti_idx = 1;
    else
        iti_idx = iti_idx + 1;
    end            

    % Total time for these trials
    dt = sum(dur + iti);

    % If logging these trials will exceed our max time, abort
    if strict && ((t + dt) > maxtime)
        break
    end

    % Log the trials       
    all_trials = [all_trials; trials];

    % Calculate total time
    t = t + dt;
end


% -------------------------------------------------------------------------
function [onsets,durations] = generate_spm_onsets(trials)
% Converts from the trials matrix format used here: 
% [condition, duration, iti] to standard SPM format for onsets and
% durations

% Unpack
cond = trials(:,1);
dur  = trials(:,2);
iti  = trials(:,3);

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
