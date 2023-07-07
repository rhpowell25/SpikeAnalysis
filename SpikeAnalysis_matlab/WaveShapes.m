function [peaktopeak_amp, euc_dis, wave_p_val] = WaveShapes(xds, unit_name, Plot_Figs, Save_Figs)

%% Find the unit of interest
[N] = Find_Unit(xds, unit_name);

% Extracting the spikes of the designated unit
spikes = xds.spikes{N};

%% Catch possible sources of error
% If there is no unit of that name
if isempty(N)
    fprintf('%s does not exist \n', unit_name);
    peaktopeak_amp = NaN;
    euc_dis = NaN;
    wave_p_val = NaN;
    return
end

% If there are less than 1000 spikes
spike_limit = 1000;
if length(spikes) < spike_limit
    disp('Too few spikes!')
    peaktopeak_amp = NaN;
    euc_dis = NaN;
    wave_p_val = NaN;
    return
end

%% Some variable extraction & definitions

waveform_length = 1.6; % Time (ms.) 

% Extracting the waveforms of the designated unit
unit_waveforms = xds.spike_waveforms{N};

% Extracting the spike times of the designated unit
spike_times = xds.spikes{N};

% Font specifications
label_font_size = 18;
title_font_size = 15;
legend_font_size = 15;
font_name = 'Arial';
figure_width = 550;
figure_height = 300;

%% Calculating the means & standard deviations

first_quarter_crossings = length(find(spike_times <= xds.time_frame(end)/4));
last_quarter_crossings = length(find(spike_times >= 3*xds.time_frame(end)/4));

% Calculating the means
mean_amp = mean(unit_waveforms, 1);
start_amp_mean = mean(unit_waveforms(1:first_quarter_crossings,:), 1);
end_amp_mean = mean(unit_waveforms(last_quarter_crossings:end,:), 1);

perspike_amp = mean(unit_waveforms, 2);
start_perspike_amp = mean(unit_waveforms(1:first_quarter_crossings,:), 2);
end_perspike_amp = mean(unit_waveforms(last_quarter_crossings:end,:), 2);

% Calculate the Euclidean distance
euc_dis = sqrt(sum((start_amp_mean - end_amp_mean) .^2));

% Calculating their standard deviations
standard_dev = std(unit_waveforms);

%% Find the change in mean amplitude
% Finding the mean and difference between morning and afternoon
amp_mean_mean = mean(perspike_amp);
std_amp = std(perspike_amp);

%% Calculate the peak to peak amplitudes
waveform_crest = max(mean_amp);
waveform_trough = min(mean_amp);

% Peak to peak amplitudes
peaktopeak_amp = abs(waveform_crest) + abs(waveform_trough);

% Length of each bin
bin_size = waveform_length/length(mean_amp); % Time (ms.)
spike_time = (-(10*bin_size):bin_size:waveform_length - bin_size - (10*bin_size));

%% Determine if the unit is stable
% Pick a random 50 threshold crossings 1000 times
p_val_idx = zeros(1000,1);
for ii = 1:1000
    % First quarter
    try
        rand_start_idx = randperm(size(start_perspike_amp,1), 50);
    catch
        disp('Too few spikes!')
        euc_dis = NaN;
        wave_p_val = NaN;
        return
    end
    rand_start_perspike = start_perspike_amp(rand_start_idx);
    % Last quarter
    try
        rand_end_idx = randperm(size(end_perspike_amp,1), 50);
    catch
        disp('Too few spikes!')
        euc_dis = NaN;
        wave_p_val = NaN;
        return
    end
    rand_end_perspike = end_perspike_amp(rand_end_idx);
    [~, wave_p_val] = ttest2(rand_start_perspike, rand_end_perspike);
    p_val_idx(ii) = wave_p_val;
end
wave_p_val = mean(p_val_idx);

%% Plotting 

if isequal(Plot_Figs, 1)

    %% Plotting random 100 waverforms

    rand_wave_idx = randperm(size(unit_waveforms,1),100);
    rand_waves = unit_waveforms(rand_wave_idx,:);

    rand_figure = figure;
    rand_figure.Position = [300 300 figure_width figure_height];
    hold on
    plot(spike_time, rand_waves', 'k')

    % Set the title
    rand_title = strcat('Individual Waveforms -', {' '}, char(xds.unit_names(N)));
    if contains(xds.meta.rawFileName, 'Pre')
        rand_title = strcat(rand_title, {' '}, '(Morning)');
    end
    if contains(xds.meta.rawFileName, 'Post')
        rand_title = strcat(rand_title, {' '}, '(Afternoon)');
    end
    title(rand_title, 'FontSize', title_font_size)

    % Axis Labels
    xlabel('Time (ms)', 'FontSize', label_font_size)
    ylabel('Amplitude (µV)', 'FontSize', label_font_size)

    % Annotation of the n-count
    legend_dims = [0.555 0.425 0.44 0.44];
    n_count_string = strcat('n =', {' '}, mat2str(length(rand_waves)));
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
    set(figure_axes,'fontname',font_name);

    %% Plotting the mean and standard deviations

    avg_figure = figure;
    avg_figure.Position = [300 300 figure_width figure_height];
    hold on
    plot(spike_time, mean_amp,'k')
    plot(spike_time, mean_amp + standard_dev,'r')
    plot(spike_time, mean_amp - standard_dev,'r')

    % Set the title
    avg_title = strcat('Average Waveform -', {' '}, char(xds.unit_names(N)));
    if contains(xds.meta.rawFileName, 'Pre')
        avg_title = strcat(avg_title, {' '}, '(Morning)');
    end
    if contains(xds.meta.rawFileName, 'Post')
        avg_title = strcat(avg_title, {' '}, '(Afternoon)');
    end
    title(avg_title, 'FontSize', title_font_size)

    % Axis Labels
    xlabel('Time (ms)', 'FontSize', label_font_size)
    ylabel('Amplitude (µV)', 'FontSize', label_font_size)
    
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
    histogram(perspike_amp, 'EdgeColor', 'k', 'FaceColor', [1, 1, 1])

    % Set the title
    hist_title = strcat('Average Amplitude Histogram -', {' '}, char(xds.unit_names(N)));
    if contains(xds.meta.rawFileName, 'Pre')
        hist_title = strcat(hist_title, {' '}, '(Morning)');
    end
    if contains(xds.meta.rawFileName, 'Post')
        hist_title = strcat(hist_title, {' '}, '(Afternoon)');
    end
    title(hist_title, 'FontSize', title_font_size)

    % Axis Labels
    xlabel('Mean Waveform Amplitude', 'FontSize', label_font_size)
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

    %% Calculate the per-spike waveform amplitude

    % Convert time into minutes
    spike_min = spike_times / 60;

    % Find the axis limits
    y_max = amp_mean_mean + 4*std_amp;
    y_min = amp_mean_mean - 4*std_amp;
    x_max = max(spike_min);
    axis_limits = cat(1, y_max, y_min);
    max_YLim = max(axis_limits);
    min_YLim = min(axis_limits);

    % Size of the markers
    sz = 2;

    amp_timeline = figure;
    amp_timeline.Position = [200 200 700 200];
    hold on
    scatter(spike_min, perspike_amp, sz, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'k');

    % Set the title
    mean_mean_title = strcat('Mean Waveform Amplitude -', {' '}, char(xds.unit_names(N)));
    if contains(xds.meta.rawFileName, 'Pre')
        mean_mean_title = strcat(mean_mean_title, {' '}, '(Morning)');
    end
    if contains(xds.meta.rawFileName, 'Post')
        mean_mean_title = strcat(mean_mean_title, {' '}, '(Afternoon)');
    end
    title(mean_mean_title, 'FontSize', title_font_size)

    % Line indicating the mean
    line([0 max(spike_min)],[amp_mean_mean amp_mean_mean], ... 
        'Color', 'red', 'LineWidth',1);
    % Lines indicating the standard deviation
    line([0 max(spike_min)], [amp_mean_mean + std_amp, amp_mean_mean + std_amp], ...
        'Color', 'red', 'LineStyle','--', 'LineWidth',1);
    line([0 max(spike_min)], [amp_mean_mean - std_amp, amp_mean_mean - std_amp], ...
        'Color', 'red', 'LineStyle','--', 'LineWidth',1);

    % Setting the afternoon axis limits
    ylim([min_YLim, max_YLim])
    xlim([0, x_max + 0.1])

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

end

%% Print the peak to peak amplitudes
fprintf("The average peak to peak amplitude of %s is %0.1f µV \n", ...
string(xds.unit_names{N}), peaktopeak_amp);

%% Figure Saving
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
        title '';
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



