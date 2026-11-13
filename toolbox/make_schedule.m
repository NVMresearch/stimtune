function [schedule,total_time] = make_schedule(trial_duration, n_conditions, iti_sampler, max_scan_time)
% MAKE_SCHEDULE Build a trial schedule with ITIs balanced across conditions.
%
% Inputs:
%   trial_duration - duration of each trial (secs, positive scalar)
%   n_conditions   - number of experimental conditions (positive integer)
%   iti_sampler     - function handle returning one nonnegative ITI (secs)
%   max_scan_time  - maximum total scan time (secs), including the final ITI
%
% Outputs:
%   schedule - one row per trial; columns are:
%              1) condition number
%              2) trial duration (secs)
%              3) post-trial ITI duration (secs)
% 
% This function was initially generated using Copilot and manually validated

% ---- Input validation ----
validateattributes(trial_duration, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'positive'});
validateattributes(n_conditions, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'integer', 'positive'});
validateattributes(max_scan_time, {'numeric'}, ...
    {'real', 'finite', 'scalar', 'nonnegative'});
assert(isa(iti_sampler, 'function_handle'), ...
    'sample_iti must be a function handle.');

% ---- Sample one ITI per set of trials ----
% Each "set" has one trial per condition, and all trials in a set share the
% same ITI. This guarantees equal trial numbers and identical ITI
% distributions across conditions.
set_itis   = zeros(0, 1);
total_time = 0;

% Stop early if even a set with zero ITI would not fit
while total_time + n_conditions * trial_duration <= max_scan_time
    iti = iti_sampler();
    validateattributes(iti, {'numeric'}, ...
        {'real', 'finite', 'scalar', 'nonnegative'});

    % Time taken by this set: every condition's trial plus its ITI
    set_time = n_conditions * (trial_duration + iti);

    % Stop at the first set that does not fit
    if total_time + set_time > max_scan_time
        break
    end

    set_itis(end + 1, 1) = iti;
    total_time = total_time + set_time;
end

% No complete set fits: return an empty schedule
if isempty(set_itis)
    schedule = zeros(0, 3);
    return
end

% ---- Build the schedule matrix ----
n_sets = numel(set_itis);

conditions     = repmat((1:n_conditions)', n_sets, 1);   % 1..N repeated per set
post_trial_iti = repelem(set_itis, n_conditions);        % same ITI for each condition in a set
durations      = repmat(trial_duration, numel(conditions), 1);

schedule = table(conditions, durations, post_trial_iti, ...
    'VariableNames',{'Condition','Duration','ITI'});
end