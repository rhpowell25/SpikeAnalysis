function [merged_fract_contam, ISI_modes] = ...
    InterstimulusInterval_Morn_v_Noon(xds_morn, xds_noon, unit_name, Plot_Figs, Save_File)

%% Find the unit of interest
[N_morn] = Find_Unit(xds_morn, unit_name);
[N_noon] = Find_Unit(xds_noon, unit_name);

% Extracting the spikes of the designated unit
spikes_morn = xds_morn.spikes{N_morn};
spikes_noon = xds_noon.spikes{N_noon};

%% Catch possible sources of error
% If there is no unit of that name
if isempty(N_morn) || isempty(N_noon)
    fprintf('%s does not exist \n', unit_name);
    merged_fract_contam = NaN;
    ISI_modes = NaN(1,3);
    return
end

% If there are less than 1000 spikes
spike_limit = 1000;
if length(spikes_morn) < spike_limit || length(spikes_noon) < spike_limit
    disp('Too few spikes!')
    merged_fract_contam = NaN;
    ISI_modes = NaN(1,3);
    return
end

%% Loop through the number of units selected

% Define the output variables
ISI_modes = zeros(1,3);

 % Font specifications
label_font_size = 20;
title_font_size = 15;
legend_font_size = 15;
font_name = 'Arial';

%% Some variable extraction & definitions

% Extracting the spikes of the designated unit
spikes_morn = xds_morn.spikes{N_morn};
spikes_noon = xds_noon.spikes{N_noon};

% Calculate the interstimulus intervals in the morning
ISI_morn = zeros(length(spikes_morn) - 1, 1);
for ii = 1:length(ISI_morn)
    ISI_morn(ii) = spikes_morn(ii+1) - spikes_morn(ii);
end

% Calculate the interstimulus intervals in the afternoon
ISI_noon = zeros(length(spikes_noon) - 1, 1);
for ii = 1:length(ISI_noon)
    ISI_noon(ii) = spikes_noon(ii+1) - spikes_noon(ii);
end

% Merge the morning & afternoon interstimulus intervals
ISI_merged = cat(1, ISI_morn, ISI_noon);

% Convert the time to milliseconds
ISI_morn = ISI_morn*1000;
ISI_noon = ISI_noon*1000;
ISI_merged = ISI_merged*1000;

% Define the edges of the histogram (0 - 100 ms)
ISI_edges = (0:1:100);

ISI_modes(1,1) = mode(ISI_morn);
ISI_modes(1,2) = mode(ISI_noon);
ISI_modes(1,3) = mode(ISI_merged);

%% Calculate the number of refractory period violations

% Define the refractory period cutoff (2 ms)
refract_min = 2;
% Find the number of refractory period violations
merged_refract_violations = length(find(ISI_merged < refract_min));

% Find the percent of spikes this corresponds to
merged_violation_ratio = merged_refract_violations / (length(spikes_morn) + length(spikes_noon));

%% Calculate the fractional level of contamination

% Define the censored period (from Central's settings: 1 ms)
censored_period = 1;

% Define the number of spikes
morn_spikes = length(spikes_morn);
noon_spikes = length(spikes_noon);
total_spikes = morn_spikes + noon_spikes;
% Define the length of the experiment (in ms)
morn_experiment_length = xds_morn.time_frame(end)*1000;
noon_experiment_length = xds_noon.time_frame(end)*1000;
total_experiment_length = morn_experiment_length + noon_experiment_length;
% Find the average firing rate of the experiment
%avg_fr = (length(spikes_morn) + length(spikes_noon)) / (xds_morn.time_frame(end) + xds_noon.time_frame(end));

syms Fract_Contam

% The Sample Numbers From The Publication
%merged_refract_violations = 20;
%refract_min = 3;
%censored_period = 1;
%total_spikes = 10000;
%experiment_length = 1000000;

merged_eqn = merged_refract_violations == 2*(refract_min - censored_period)*(total_spikes^2)*... 
    (1 - Fract_Contam)*(Fract_Contam/total_experiment_length);

fract_contam = solve(merged_eqn, Fract_Contam);

%% Plotting the histograms
if isequal(Plot_Figs, 1)

    % Merged histogram
    figure
    hold on
    histogram(ISI_morn, ISI_edges, 'EdgeColor', 'k', 'FaceColor', [0.9290, 0.6940, 0.1250])
    histogram(ISI_noon, ISI_edges, 'EdgeColor', 'k', 'FaceColor', [.5 0 .5])
    % Set the title
    Fig_Title = sprintf('Interstimulus Intervals - %s (Merged)', ... 
        char(xds_morn.unit_names(N_morn)));
    title(Fig_Title, 'FontSize', title_font_size)
    xlabel('ISI (mSec)', 'FontSize', label_font_size)
    ylabel('Counts', 'FontSize', label_font_size)

    % Collect the current axis limits
    y_limits = ylim;
    x_limits = xlim;

    % Annotation of the fractional contamination
    if isreal(fract_contam)
        legend_dims = [0.52 0.35 0.44 0.44];
        fract_contam_string = strcat('f_{1}^{p} =', {' '}, char(round(min(fract_contam), 3)));
        legend_string = {char(fract_contam_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
            'FitBoxToText', 'on', 'EdgeColor','none', ... 
            'verticalalignment', 'top', 'horizontalalignment', 'center');
        ann_legend.FontSize = legend_font_size;
        ann_legend.FontName = font_name;
    end
    if ~isreal(fract_contam)
        legend_dims = [0.52 0.35 0.44 0.44];
        fract_contam_string = strcat('f_{1}^{p} >', {' '}, '0.5');
        legend_string = {char(fract_contam_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
            'FitBoxToText', 'on', 'EdgeColor','none', ... 
            'verticalalignment', 'top', 'horizontalalignment', 'center');
        ann_legend.FontSize = legend_font_size;
        ann_legend.FontName = font_name;
    end

    % Plot dummy points for the legend
    dummy_yellow = scatter(-1, -1, 's', 'filled', 'Color', [0.9290, 0.6940, 0.1250]);
    dummy_purple = scatter(-1.5, -1.5, 's', 'filled', 'Color', [.5 0 .5]);
    % Legend
    legend([dummy_yellow, dummy_purple], ... 
        {'Morning', 'Afternoon'}, 'FontSize', legend_font_size, 'Location', 'northeast')
    legend boxoff

    % Reset the axis limits
    xlim([x_limits(1),x_limits(2)])
    ylim([y_limits(1),y_limits(2) + 1])

    % Only label every other tick
    figure_axes = gca;
    x_labels = string(figure_axes.XAxis.TickLabels);
    y_labels = string(figure_axes.YAxis.TickLabels);
    x_labels(2:2:end) = NaN;
    y_labels(2:2:end) = NaN;
    figure_axes.XAxis.TickLabels = x_labels;
    figure_axes.YAxis.TickLabels = y_labels;
    % Set ticks to outside
    set(figure_axes,'TickDir','out');
    % Remove the top and right tick marks
    set(figure_axes,'box','off')
    % Set The Font
    set(figure_axes,'fontname', font_name);

    %% Save the file if selected
    Save_Figs(Fig_Title, Save_File)

end

%% Print the percentages & fractional contamination

% Merged
fprintf("%0.1f%% of the spikes in %s have an ISI less than %0.1f ms. \n", ...
merged_violation_ratio*100, string(xds_noon.unit_names{N_noon}), refract_min);

% Merged
if ~isreal(fract_contam)
    fprintf("The fractional level of contamination in %s > 0.5 \n", ...
        string(xds_noon.unit_names{N_noon}));
    merged_fract_contam = 1;
end
if isreal(fract_contam)
    fprintf("The fractional level of contamination is %0.2f in %s \n", ...
        min(fract_contam), string(xds_noon.unit_names{N_noon}));
    merged_fract_contam = double(min(fract_contam));
end




