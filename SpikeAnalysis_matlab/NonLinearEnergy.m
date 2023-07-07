function [euc_dis, nonlin_p_val] = NonLinearEnergy(xds, unit_name, Plot_Figs, Save_Figs)

%% Find the unit of interest
[N] = Find_Unit(xds, unit_name);

% Extracting the spikes of the designated unit
spikes = xds.spikes{N};

%% Catch possible sources of error
% If there is no unit of that name
if isempty(N)
    fprintf('%s does not exist \n', unit_name);
    euc_dis = NaN;
    nonlin_p_val = NaN;
    return
end

% If there are less than 1000 spikes
spike_limit = 1000;
if length(spikes) < spike_limit
    disp('Too few spikes!')
    euc_dis = NaN;
    nonlin_p_val = NaN;
    return
end

%% Some variable extraction & definitions

% Extracting the spike waveforms of the designated unit
waveforms = xds.spike_waveforms{N};

% Extracting the nonlinear waveforms of the all units
nonlin_waveforms = xds.nonlin_waveforms{N};

% Extracting the spike times of the designated unit
spike_times = xds.spikes{N};

% Font specifications
label_font_size = 18;
title_font_size = 15;
legend_font_size = 15;
font_name = 'Arial';

%% Calculating the means & standard deviations

first_quarter_crossings = length(find(spike_times <= xds.time_frame(end)/4));
last_quarter_crossings = length(find(spike_times >= 3*xds.time_frame(end)/4));

% Nonlinear waveform means
nonlin_mean = mean(nonlin_waveforms,1);
start_nonlin_mean = mean(nonlin_waveforms(1:first_quarter_crossings,:), 1);
end_nonlin_mean = mean(nonlin_waveforms(last_quarter_crossings:end,:), 1);

perspike_nonlin = mean(nonlin_waveforms, 2);
start_perspike_nonlin = mean(nonlin_waveforms(1:first_quarter_crossings,:), 2);
end_perspike_nonlin = mean(nonlin_waveforms(last_quarter_crossings:end,:), 2);

% Calculate the Euclidean distance
euc_dis = sqrt(sum((start_nonlin_mean - end_nonlin_mean) .^2));

% Nonlinear waveform standard deviations
standard_nonlin = std(nonlin_waveforms(1:height(nonlin_waveforms),:));

%% Find the change in mean nonlinear energy
nonlin_mean_mean = mean(perspike_nonlin);
std_nonlin = std(perspike_nonlin);

%% Determine if the unit is well sorted 
% Pick a random 1000 threshold crossings 1000 times
loop_size = 1000;
p_val_idx = zeros(loop_size,1);
for ii = 1:loop_size
    % First quarter
    rand_start_nonlin_idx = randperm(size(start_perspike_nonlin,1), 50);
    rand_start_perspike_nonlin = start_perspike_nonlin(rand_start_nonlin_idx);
    % Last quarter
    rand_end_nonlin_idx = randperm(size(end_perspike_nonlin,1), 50);
    rand_end_perspike_nonlin = end_perspike_nonlin(rand_end_nonlin_idx);
    [~, nonlin_p_val] = ttest2(rand_start_perspike_nonlin, rand_end_perspike_nonlin);
    p_val_idx(ii) = nonlin_p_val;
end
nonlin_p_val = mean(p_val_idx);

%% Plotting
if isequal(Plot_Figs, 1)

    %% Plot random 100 waveforms

    rand_nonlin_idx = randperm(size(waveforms,1),100);
    rand_nonlin = nonlin_waveforms(rand_nonlin_idx,:);

    figure
    hold on
    plot(rand_nonlin', 'k')

    % Set the title
    nonlin_title = strcat('Nonlinear Energy -', {' '}, char(xds.unit_names(N)));
    if contains(xds.meta.rawFileName, 'Pre')
        nonlin_title = strcat(nonlin_title, {' '}, '(Morning)');
    end
    if contains(xds.meta.rawFileName, 'Post')
        nonlin_title = strcat(nonlin_title, {' '}, '(Afternoon)');
    end
    title(nonlin_title, 'FontSize', title_font_size)

    % Axis labels
    xlabel('Time', 'FontSize', label_font_size)
    ylabel('µV^2', 'FontSize', label_font_size)

    % Annotation of the n-count
    legend_dims = [0.555 0.425 0.44 0.44];
    n_count_string = strcat('n =', {' '}, mat2str(length(rand_nonlin)));
    legend_string = {char(n_count_string)};
    ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
        'FitBoxToText', 'on', 'verticalalignment', 'top', ...
        'EdgeColor','none', 'horizontalalignment', 'center');
    ann_legend.FontSize = legend_font_size;
    ann_legend.FontName = font_name;

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

    %% Plotting the mean & standard deviations

    figure
    hold on
    plot(nonlin_mean,'k')
    plot(mean(transpose(nonlin_waveforms(1:length(nonlin_waveforms),:) ...
        + standard_nonlin), 2),'r')
    plot(mean(transpose(nonlin_waveforms(1:length(nonlin_waveforms),:) ... 
        - standard_nonlin), 2),'r')

    % Set the title
    avg_title = strcat('Average Nonlinear Energy -', {' '}, char(xds.unit_names(N)));
    if contains(xds.meta.rawFileName, 'Pre')
        avg_title = strcat(avg_title, {' '}, '(Morning)');
    end
    if contains(xds.meta.rawFileName, 'Post')
        avg_title = strcat(avg_title, {' '}, '(Afternoon)');
    end
    title(avg_title, 'FontSize', title_font_size)

    % Axis Labels
    xlabel('Time', 'FontSize', label_font_size)
    ylabel('µV^2', 'FontSize', label_font_size)

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

    %% Plot the histogram
    figure
    hold on
    histogram(perspike_nonlin, 'EdgeColor', 'k', 'FaceColor', [1, 1, 1])

    % Set the title
    hist_title = strcat('Nonlinear Energy Histogram -', {' '}, char(xds.unit_names(N)));
    if contains(xds.meta.rawFileName, 'Pre')
        hist_title = strcat(hist_title, {' '}, '(Morning)');
    end
    if contains(xds.meta.rawFileName, 'Post')
        hist_title = strcat(hist_title, {' '}, '(Afternoon)');
    end
    title(hist_title, 'FontSize', title_font_size)

    % Axis Labels
    xlabel('Mean Nonlinear Energy', 'FontSize', label_font_size)
    ylabel('Counts', 'FontSize', label_font_size)

    % Collect the current axis limits
    y_limits = ylim;
    x_limits = xlim;

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
    spike_min = spike_times / 60;

    % Find the axis limits
    y_max = nonlin_mean_mean + 4*std_nonlin;
    y_min = nonlin_mean_mean - 4*std_nonlin;
    x_max = max(spike_min);
    axis_limits = cat(1, y_max, y_min);
    max_YLim = max(axis_limits);
    min_YLim = min(axis_limits);

    % Size of the markers
    sz = 2;

    nonlin_timeline = figure;
    nonlin_timeline.Position = [200 200 700 200];
    hold on
    scatter(spike_min, perspike_nonlin, sz, 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'k');

    % Set the title
    mean_mean_title = strcat('Nonlinear Energy -', {' '}, char(xds.unit_names(N)));
    if contains(xds.meta.rawFileName, 'Pre')
        mean_mean_title = strcat(mean_mean_title, {' '}, '(Morning)');
    end
    if contains(xds.meta.rawFileName, 'Post')
        mean_mean_title = strcat(mean_mean_title, {' '}, '(Afternoon)');
    end
    title(mean_mean_title, 'FontSize', title_font_size)

    % Lines indicating the morning mean
    line([0 max(spike_min)],[nonlin_mean_mean nonlin_mean_mean], ... 
        'Color', 'red', 'LineWidth',1);
    % Lines indicating the morning standard deviation
    line([0 max(spike_min)], [nonlin_mean_mean + std_nonlin, nonlin_mean_mean + std_nonlin], ...
        'Color', 'red', 'LineStyle','--', 'LineWidth',1);
    line([0 max(spike_min)], [nonlin_mean_mean - std_nonlin, nonlin_mean_mean - std_nonlin], ...
        'Color', 'red', 'LineStyle','--', 'LineWidth',1);

    % Setting the morning axis limits
    ylim([min_YLim, max_YLim])
    xlim([0, x_max + 0.1])

    % Set the ticks for morning
    figure_axes = gca;
    x_labels = string(figure_axes.XAxis.TickLabels);
    morn_end_tick = find(x_labels == num2str(round(max(spike_min))));
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

    % Set the x-label
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


