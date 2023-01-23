function [euc_dis, nonlin_p_val] = NonLinearEnergy_Morn_v_Noon(xds_morn, xds_noon, unit_name, Plot_Figs, Save_Figs)

%% Load the excel file
if ~ischar(unit_name)

    [xds_output] = Find_Excel(xds_morn);

    %% Find the unit of interest

    unit = xds_output.unit_names(unit_name);

    %% Identify the index of the unit
    N_morn = find(strcmp(xds_morn.unit_names, unit));
    N_noon = find(strcmp(xds_noon.unit_names, unit));

else
    N_morn = find(strcmp(xds_morn.unit_names, unit_name));
    N_noon = find(strcmp(xds_noon.unit_names, unit_name));
end

%% If The Unit Doesn't Exist

if isempty(N_morn) || isempty(N_noon)
    fprintf('%s does not exist \n', unit_name);
    euc_dis = NaN;
    nonlin_p_val = NaN;
    return
end

%% Some variable extraction & definitions

% Extracting the spike waveforms of the designated unit
waveforms_morn = xds_morn.spike_waveforms{N_morn};
waveforms_noon = xds_noon.spike_waveforms{N_noon};

if isfield(xds_morn, 'nonlin_waveforms') == 1
    % Extracting the nonlinear waveforms of the all units
    nonlin_waveforms_morn = xds_morn.nonlin_waveforms;
    nonlin_waveforms_noon = xds_noon.nonlin_waveforms;
else
    nonlin_waveforms_morn = zeros(height(waveforms_morn), width(waveforms_morn) - 2);
    nonlin_waveforms_noon = zeros(height(waveforms_noon), width(waveforms_noon) - 2);
    % Morning
    for jj = 1:height(waveforms_morn)
        for ii = 1:width(waveforms_morn) - 2
            nonlin_waveforms_morn(jj,ii) = waveforms_morn(jj,ii+1)^2 - ... 
                waveforms_morn(jj,ii)*waveforms_morn(jj,ii+2);
        end
    end
    % Afternoon
    for jj = 1:height(waveforms_noon)
        for ii = 1:width(waveforms_noon) - 2
            nonlin_waveforms_noon(jj,ii) = waveforms_noon(jj,ii+1)^2 - ... 
                waveforms_noon(jj,ii)*waveforms_noon(jj,ii+2);
        end
    end
end

% Extracting the spike times of the designated unit
morn_spikes = xds_morn.spikes{N_morn};
noon_spikes = xds_noon.spikes{N_noon};

% Font specifications
label_font_size = 18;
title_font_size = 15;
legend_font_size = 15;
font_name = 'Arial';

%% Calculating the means & standard deviations

% Nonlinear waveform means
if isfield(xds_morn, 'nonlin_waveforms') == 1
    morn_nonlin_mean = mean(nonlin_waveforms_morn{N_morn},1);
    noon_nonlin_mean = mean(nonlin_waveforms_noon{N_noon},1);
    morn_perspike_nonlin = mean(nonlin_waveforms_morn{N_morn}, 2);
    noon_perspike_nonlin = mean(nonlin_waveforms_noon{N_noon}, 2);
else
    morn_nonlin_mean = mean(nonlin_waveforms_morn,1);
    noon_nonlin_mean = mean(nonlin_waveforms_noon,1);
    morn_perspike_nonlin = mean(nonlin_waveforms_morn, 2);
    noon_perspike_nonlin = mean(nonlin_waveforms_noon, 2);
end

% Calculate the Euclidean distance
euc_dis = sqrt(sum((morn_nonlin_mean - noon_nonlin_mean) .^2));

%% Find the change in mean nonlinear energy between the two sessions 
% Finding the mean & difference between morning & afternoon
nonlin_morn = mean(morn_perspike_nonlin);
nonlin_noon = mean(noon_perspike_nonlin);
std_nonlin_morn = std(morn_perspike_nonlin);
std_nonlin_noon = std(noon_perspike_nonlin);

% Finding the difference
nonlin_difference = abs(nonlin_morn - nonlin_noon);

%% Find the average nonlinear energy across all units
if isfield(xds_morn, 'nonlin_waveforms') == 1
    % The means of each nonlinear waveform
    all_units_nonlin_morn = struct([]);
    all_units_nonlin_noon = struct([]);
    for ii = 1:length(xds_morn.unit_names)
        all_units_nonlin_morn{ii,1} = mean(nonlin_waveforms_morn{ii});
    end
    for ii = 1:length(xds_noon.unit_names)
        all_units_nonlin_noon{ii,1} = mean(nonlin_waveforms_noon{ii});
    end

    % The mean of every point in the mean of the nonlinear waveforms
    all_units_nonlin_means_morn = zeros(length(xds_morn.unit_names),1);
    all_units_nonlin_means_noon = zeros(length(xds_noon.unit_names),1);
    for ii = 1:length(xds_morn.unit_names)
        all_units_nonlin_means_morn(ii) = mean((all_units_nonlin_morn{ii,1}));
    end
    for ii = 1:length(xds_noon.unit_names)
        all_units_nonlin_means_noon(ii) = mean((all_units_nonlin_noon{ii,1}));
    end

    % Find the differences between morning & afternoon
    all_units_nonlin_difference = abs(all_units_nonlin_means_morn - all_units_nonlin_means_noon);
end

%% Determine if the unit is well sorted 
% Pick a random 1000 threshold crossings 1000 times
loop_size = 1000;
p_val_idx = zeros(loop_size,1);
for ii = 1:loop_size
    % Morning
    rand_morn_nonlin_idx = randperm(size(morn_perspike_nonlin,1), 50);
    rand_morn_perspike_nonlin = morn_perspike_nonlin(rand_morn_nonlin_idx);
    % Afternoon
    rand_noon_nonlin_idx = randperm(size(noon_perspike_nonlin,1), 50);
    rand_noon_perspike_nonlin = noon_perspike_nonlin(rand_noon_nonlin_idx);
    [~, nonlin_p_val] = ttest2(rand_morn_perspike_nonlin, rand_noon_perspike_nonlin);
    p_val_idx(ii) = nonlin_p_val;
end
nonlin_p_val = mean(p_val_idx);

%% Determine if the change in nonlinear energy is more or less than the average change
if isfield(xds_morn, 'nonlin_waveforms') == 1
    mean_nonlin_difference = mean(all_units_nonlin_difference);

    if nonlin_difference < mean_nonlin_difference
        fprintf("The change in mean nonlinear energy of %s is less than the mean of all units \n", ...
            string(xds_morn.unit_names{N_morn}));
    end

    if nonlin_difference > mean_nonlin_difference
        fprintf("The change in mean nonlinear energy of %s is greater than the mean of all units \n", ...
            string(xds_morn.unit_names{N_morn}));
    end
end

%% Plotting

if isequal(Plot_Figs, 1)

    %% Plot the morning & afternoon histograms
    figure
    hold on
    histogram(morn_perspike_nonlin, 'EdgeColor', 'k', 'FaceColor', [0.9290, 0.6940, 0.1250])
    histogram(noon_perspike_nonlin, 'EdgeColor', 'k', 'FaceColor', [.5 0 .5])
    title(sprintf('Nonlinear Energy Histogram - %s', ... 
        char(xds_morn.unit_names(N_morn))), 'FontSize', title_font_size)
    xlabel('Mean Nonlinear Energy', 'FontSize', label_font_size)
    ylabel('Counts', 'FontSize', label_font_size)

    % Collect the current axis limits
    y_limits = ylim;
    x_limits = xlim;

    % Plot dummy points for the legend
    dummy_yellow = scatter(-1, -1, 's', 'filled', 'Color', [0.9290, 0.6940, 0.1250]);
    dummy_purple = scatter(-1.5, -1.5, 's', 'filled', 'Color', [.5 0 .5]);
    % Legend
    legend([dummy_yellow, dummy_purple], ... 
        {'Morning', 'Afternoon'}, 'FontSize', legend_font_size, 'Location', 'northeast')
    legend boxoff

    % Annotation of the p_value
    if round(nonlin_p_val, 3) > 0
        legend_dims = [0.525 0.325 0.44 0.44];
        p_value_string = strcat('p =', {' '}, mat2str(round(nonlin_p_val, 3)));
        legend_string = {char(p_value_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
            'FitBoxToText', 'on', 'verticalalignment', 'top', ...
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = legend_font_size;
        ann_legend.FontName = font_name;
    end

        if isequal(round(nonlin_p_val, 3), 0)
            legend_dims = [0.525 0.325 0.44 0.44];
            p_value_string = strcat('p <', {' '}, '0.001');
            legend_string = {char(p_value_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = legend_font_size;
            ann_legend.FontName = font_name;
        end

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
    set(figure_axes,'box','off');
    % Set The Font
    set(figure_axes,'fontname', font_name);

    %% Calculate the per-spike nonlinear energy

    % Convert time into minutes
    spike_min_morn = morn_spikes / 60;
    spike_min_noon = noon_spikes / 60;

    % Find the axis limits
    morn_ymax = nonlin_morn + 4*std_nonlin_morn;
    morn_ymin = nonlin_morn - 4*std_nonlin_morn;
    noon_ymax = nonlin_noon + 4*std_nonlin_noon;
    noon_ymin = nonlin_noon - 4*std_nonlin_noon;
    morn_xmax = max(spike_min_morn);
    noon_xmax = max(spike_min_noon);
    axis_limits = cat(1, morn_ymax, morn_ymin, noon_ymax, noon_ymin);
    max_YLim = max(axis_limits);
    min_YLim = min(axis_limits);

    % Size of the markers
    sz = 2;

    nonlin_timeline = figure;
    nonlin_timeline.Position = [200 200 750 350];
    % Morning
    morn_nonlin_plot = subplot(121); % (Top Left Plot)
    hold on
    scatter(spike_min_morn, morn_perspike_nonlin, sz, 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'k');
    % Lines indicating the morning mean
    line([0 max(spike_min_morn)],[nonlin_morn nonlin_morn], ... 
        'Color', 'red', 'LineWidth',1);
    % Lines indicating the morning standard deviation
    line([0 max(spike_min_morn)], [nonlin_morn + std_nonlin_morn, nonlin_morn + std_nonlin_morn], ...
        'Color', 'red', 'LineStyle','--', 'LineWidth',1);
    line([0 max(spike_min_morn)], [nonlin_morn - std_nonlin_morn, nonlin_morn - std_nonlin_morn], ...
        'Color', 'red', 'LineStyle','--', 'LineWidth',1);

    % Setting the morning axis limits
    ylim([min_YLim, max_YLim])
    xlim([0, morn_xmax + 0.1])

    % Set the ticks for morning
    figure_axes = gca;
    x_labels = string(figure_axes.XAxis.TickLabels);
    morn_end_tick = find(x_labels == num2str(round(max(spike_min_morn))));
    noon_ticks = str2double(x_labels(morn_end_tick + 1:end)) - str2double(x_labels(morn_end_tick));
    x_labels(morn_end_tick + 1:end) = num2str(noon_ticks);
    y_labels = string(figure_axes.YAxis.TickLabels);
    % Only label every other tick
    x_labels(1:2:end) = NaN;
    y_labels(2:2:end) = NaN;
    figure_axes.XAxis.TickLabels = x_labels;
    figure_axes.YAxis.TickLabels = y_labels;
    % Set ticks to outside
    set(figure_axes,'TickDir','out');
    % Remove the top and right tick marks
    set(figure_axes,'box','off');
    % Set The Font
    set(figure_axes,'fontname', font_name);

    % Bring the morning plot to the right
    set(morn_nonlin_plot, 'Units', 'normalized')
    set(morn_nonlin_plot, 'Position', [0.1, 0.65, 0.4, 0.2]);

    xlabel('Morning', 'FontSize', label_font_size)

    % Afternoon
    noon_nonlin_plot = subplot(122); % (Top Right Plot)
    scatter(spike_min_noon, noon_perspike_nonlin, sz, 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'k');
    % Lines indicating the afternoon mean
    line([min(spike_min_noon) max(spike_min_noon)],[nonlin_noon nonlin_noon], ... 
        'Color', 'red', 'LineWidth',1);
    % Lines indicating the afternoon standard deviation
    line([min(spike_min_noon) max(spike_min_noon)], [nonlin_noon + std_nonlin_noon, nonlin_noon + std_nonlin_noon], ...
        'Color', 'red', 'LineStyle','--', 'LineWidth',1);
    line([min(spike_min_noon) max(spike_min_noon)], [nonlin_noon - std_nonlin_noon, nonlin_noon - std_nonlin_noon], ...
        'Color', 'red', 'LineStyle','--', 'LineWidth',1);

    % Setting the afternoon axis limits
    ylim([min_YLim, max_YLim])
    xlim([0, noon_xmax + 0.1])

    % Set the ticks for afternoon
    figure_axes = gca;
    x_labels = string(figure_axes.XAxis.TickLabels);
    morn_end_tick = find(x_labels == num2str(round(max(spike_min_morn))));
    noon_ticks = str2double(x_labels(morn_end_tick + 1:end)) - str2double(x_labels(morn_end_tick));
    x_labels(morn_end_tick + 1:end) = num2str(noon_ticks);
    set(figure_axes, 'YTickLabel', []);
    % Only label every other tick
    x_labels(1:2:end) = NaN;
    figure_axes.XAxis.TickLabels = x_labels;
    % Set ticks to outside
    set(figure_axes,'TickDir','out');
    % Remove the top and right tick marks
    set(figure_axes,'box','off');
    % Set The Font
    set(figure_axes,'fontname', font_name);

    % Bring the afternoon plot to the left
    set(noon_nonlin_plot, 'Units', 'normalized')
    set(noon_nonlin_plot, 'Position', [0.525, 0.65, .4, .2]);
    
    xlabel('Afternoon', 'FontSize', label_font_size)

    % Annotation of the p_value
    if round(nonlin_p_val, 3) > 0
        legend_dims = [0.65 0.525 0.44 0.44];
        p_value_string = strcat('p =', {' '}, mat2str(round(nonlin_p_val, 3)));
        legend_string = {char(p_value_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = legend_font_size;
        ann_legend.FontName = font_name;
    end

    if isequal(round(nonlin_p_val, 3), 0)
        legend_dims = [0.65 0.525 0.44 0.44];
        p_value_string = strcat('p <', {' '}, '0.001');
        legend_string = {char(p_value_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = legend_font_size;
        ann_legend.FontName = font_name;
    end

    % Set the common title
    sgt = sgtitle(sprintf('Nonlinear Energy - %s', ...
        char(xds_morn.unit_names(N_morn))), 'FontSize', (title_font_size + 5));
    % Set the common xlabel
    common_x_label = axes(nonlin_timeline, 'visible', 'off');
    common_x_label.XLabel.Visible = 'on';
    xlabel(common_x_label, 'Time (min.)', 'FontSize', label_font_size);
    common_x_label.Position(2) = common_x_label.Position(2) + 0.5;

end

%% Define the save directory & save the figures
if ~isequal(Save_Figs, 0)
    save_dir = 'C:\Users\rhpow\Desktop\';
    for ii = 1:length(findobj('type','figure'))
        fig_info = get(gca,'title');
        fig_title = get(fig_info, 'string');
        if isempty(fig_title)
            fig_info = sgt;
            fig_title = get(fig_info, 'string');
        end
        fig_title = strrep(fig_title, ':', '');
        fig_title = strrep(fig_title, 'vs.', 'vs');
        fig_title = strrep(fig_title, 'mg.', 'mg');
        fig_title = strrep(fig_title, 'kg.', 'kg');
        fig_title = strrep(fig_title, '.', '_');
        fig_title = strrep(fig_title, '/', '_');
        %title '';
        if strcmp(Save_Figs, 'All')
            saveas(gcf, fullfile(save_dir, char(fig_title)), 'png')
            saveas(gcf, fullfile(save_dir, char(fig_title)), 'pdf')
            saveas(gcf, fullfile(save_dir, char(fig_title)), 'fig')
        else
            saveas(gcf, fullfile(save_dir, char(fig_title)), Save_Figs)
        end
        close gcf
    end
end


