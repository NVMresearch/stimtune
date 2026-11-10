function report_optimise_design_results(T_factorial,T_designs,eff,opt,con)

min_gap_range = opt{1}.settings.min_gap_range;
max_gap_range = opt{1}.settings.max_gap_range;

% Display table of optimal designs
% -------------------------------------------------------------------------
disp(T_factorial);

% Plot efficiencies across landscape of jitters for the optimal duration
% -------------------------------------------------------------------------

% Maximum efficiencies over repetitions [design x contrast]
eff_max = max(eff,[],3);

% Create figure
figure;
nrows = 1;
ncols = length(con);

% Loop over contrasts
for i = 1:ncols
    subplot(nrows,ncols,i);

    % Get the optimal stimulus duration
    opt_duration = T_factorial{i,"stim_duration"};

    % Output variable
    efficiency_matrix = zeros(length(min_gap_range), length(max_gap_range));

    % Find all designs with the optimal duration
    relevant_design_idx = find(T_designs.duration==opt_duration);

    % Loop over relevant designs
    for j = 1:size(relevant_design_idx,1)
        % Get this design's gap range
        min_gap = T_designs{relevant_design_idx(j),'min_gap'};
        max_gap = T_designs{relevant_design_idx(j),'max_gap'};

        % Get this design's max efficiency over repetitions
        efficiency = eff_max(relevant_design_idx(j),i);

        % Place in matrix
        min_gap_idx = find(min_gap_range == min_gap);
        max_gap_idx = find(max_gap_range == max_gap);
        efficiency_matrix(min_gap_idx,max_gap_idx) = efficiency;
    end

    % Plot efficiency_matrix
    imagesc(max_gap_range, min_gap_range, efficiency_matrix);
    set(gca, 'YDir', 'normal');
    xticks(max_gap_range);
    yticks(min_gap_range);
    
    % Decorate
    colorbar
    colormap gray;
    axis square
    str = sprintf('Efficiency of %s\n For optimal stim duration: %2.2fs',...
        con(i).name, opt_duration);
    title(str);
    xlabel('Max jitter (s)');
    ylabel('Min jitter (s)');

    % Identify row and column index of maximum value in efficiency_matrix
    [maxEfficiency, linearIdx] = max(efficiency_matrix(:));
    [maxRowIdx, maxColIdx] = ind2sub(size(efficiency_matrix), linearIdx);
    hold on;
    plot(max_gap_range(maxColIdx), min_gap_range(maxRowIdx), 'r*', ...
        'MarkerSize', 10, 'LineWidth', 1.5);
    hold off;

end

% Plot effect of duration on efficiency for optimal jitter range
% -------------------------------------------------------------------------

% Create figure
figure;
nrows = 1;
ncols = length(con);
ax = [];

% Loop over contrasts
for i = 1:ncols
    ax(i)=subplot(nrows,ncols,i);

    % Get the optimal ITI durations
    opt_min_gap = T_factorial{i,"min_gap"};
    opt_max_gap = T_factorial{i,"max_gap"};

    % Find all designs with the optimal durations
    relevant_design_idx = find(T_designs.min_gap==opt_min_gap & ...
        T_designs.max_gap==opt_max_gap);

    % Extract efficiencies across stimulus durations
    efficiency = eff_max(relevant_design_idx,i);
    durations  = T_designs.duration(relevant_design_idx);

    % Plot efficiency by duration
    plot(durations, efficiency, '.','MarkerSize',20);
    xlabel('Stimulus duration (s)');
    ylabel('Maximum efficiency');
    str = sprintf('Efficiency of %s\nFor optimal gap range [%2.2fs-%2.2fs]',...
        con(i).name, opt_min_gap, opt_max_gap);    
    title(str);
    grid on;
    axis square
end
linkaxes(ax);

% Plot best repetition of best design
% -------------------------------------------------------------------------
figure;
for i = 1:ncols

    subplot(1,ncols,i);

    u = [opt{i}.U.u];
    u = u(33:end,:);
    dt = opt{i}.U.dt;

    t = (0:(length(u)-1)) .* dt;
    bar(t,u);

    str = sprintf('Design for contrast %s', ...
        con(i).name);
    title(str);

    xlabel('Scan time (seconds)');
    axis square
end