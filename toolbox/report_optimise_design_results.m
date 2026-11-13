function report_optimise_design_results(T_factorial,T_designs,eff,opt,con)

iti_param1_range = opt{1}.settings.iti_dist.param_ranges{1};
iti_param2_range = opt{1}.settings.iti_dist.param_ranges{2};

% Display table of optimal designs
% -------------------------------------------------------------------------
disp(T_factorial);

% Plot efficiencies across landscape of jitters for the optimal duration
% -------------------------------------------------------------------------

% Maximum efficiencies per contrast over repetitions [design x contrast]
eff_max = max(eff,[],3);

% Create figure
figure;
nrows = 1;
ncols = length(con);

% Loop over contrasts
for i = 1:ncols
    subplot(nrows,ncols,i);

    % Get the optimal stimulus duration
    opt_duration = T_factorial{i,"params"}.Duration;

    % Output variable
    efficiency_matrix = zeros(length(iti_param1_range), length(iti_param2_range));

    % Find all designs with the optimal duration
    relevant_design_idx = find(T_designs.params.Duration==opt_duration);

    % Loop over relevant designs
    for j = 1:size(relevant_design_idx,1)
        % Get this design's ITI parameter range
        iti_param1 = T_designs.params.ITI_Param_1(relevant_design_idx(j));
        iti_param2 = T_designs.params.ITI_Param_2(relevant_design_idx(j));

        % Get this design's max efficiency over repetitions
        efficiency = eff_max(relevant_design_idx(j),i);

        % Place in matrix
        param1_idx = find(iti_param1_range == iti_param1);
        param2_idx = find(iti_param2_range == iti_param2);
        efficiency_matrix(param1_idx,param2_idx) = efficiency;
    end

    % Plot efficiency_matrix
    imagesc(iti_param2_range, iti_param1_range, efficiency_matrix);
    set(gca, 'YDir', 'normal');
    xticks(iti_param2_range);
    yticks(iti_param1_range);
    
    % Decorate
    colorbar
    colormap gray;
    axis square
    str = sprintf('Efficiency of %s\n For optimal stim duration: %2.2fs',...
        con(i).name, opt_duration);
    title(str);
    xlabel('ITI parameter 2');
    ylabel('ITI parameter 1');

    % Identify row and column index of maximum value in efficiency_matrix
    [maxEfficiency, linearIdx] = max(efficiency_matrix(:));
    [maxRowIdx, maxColIdx] = ind2sub(size(efficiency_matrix), linearIdx);
    hold on;
    plot(iti_param2_range(maxColIdx), iti_param1_range(maxRowIdx), 'r*', ...
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
    opt_iti_param_1 = T_factorial.params.ITI_Param_1(i);
    opt_iti_param_2 = T_factorial.params.ITI_Param_2(i);

    % Find all designs with the optimal durations
    relevant_design_idx = find(T_designs.params.ITI_Param_1==opt_iti_param_1 & ...
        T_designs.params.ITI_Param_2==opt_iti_param_2);

    % Extract efficiencies across stimulus durations
    efficiency = eff_max(relevant_design_idx,i);
    durations  = T_designs.params.Duration(relevant_design_idx);

    % Plot efficiency by duration
    plot(durations, efficiency, '.','MarkerSize',20);
    xlabel('Stimulus duration (s)');
    ylabel('Maximum efficiency');
    str = sprintf('Efficiency of %s\nFor optimal ITI params [%2.2fs-%2.2fs]',...
        con(i).name, opt_iti_param_1, opt_iti_param_2);    
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