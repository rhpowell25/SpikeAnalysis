function [euc_dis, wave_p_val] = WaveShapes_Morn_v_Noon(xds_morn, xds_noon, unit_name, Plot_Figs, Save_Figs)

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
    wave_p_val = NaN;
    return
end

%% Some variable extraction & definitions

% Extracting the waveforms of the designated unit
waveforms_morn = xds_morn.spike_waveforms;
waveforms_noon = xds_noon.spike_waveforms;

% Extracting the spike times of the designated unit
spike_times_morn = xds_morn.spikes{N_morn};
spike_times_noon = xds_noon.spikes{N_noon};

% Font specifications
label_font_size = 18;
title_font_size = 15;
legend_font_size = 15;
font_name = 'Arial';

%% Calculating the means and standard deviations

% Calculating the means
morn_amp_mean = mean(waveforms_morn{N_morn}, 1);
noon_amp_mean = mean(waveforms_noon{N_noon}, 1);
morn_perspike_amp = mean(waveforms_morn{N_morn}, 2);
noon_perspike_amp = mean(waveforms_noon{N_noon}, 2);

% Calculate the Euclidean distance
euc_dis = sqrt(sum((morn_amp_mean - noon_amp_mean) .^2));

%% Find the change in mean amplitude between the two sessions
% Finding the mean and difference between morning and afternoon
amp_morn = mean(morn_perspike_amp);
amp_noon = mean(noon_perspike_amp);
std_amp_morn = std(morn_perspike_amp);
std_amp_noon = std(noon_perspike_amp);

% Finding the difference
amp_difference = abs(amp_morn - amp_noon);

%% Calculate the peak to peak amplitudes
% Morning
crest_morn = max(morn_amp_mean);
trough_morn = min(morn_amp_mean);
% Afternoon
crest_noon = max(noon_amp_mean);
trough_noon = min(noon_amp_mean);

% Peak to peak amplitudes
peaktopeak_morn = abs(crest_morn) + abs(trough_morn);
peaktopeak_noon = abs(crest_noon) + abs(trough_noon);

%% Finding the average amplitude across all units
% The means of each waveform
all_units_amp_morn = struct([]);
all_units_amp_noon = struct([]);
for ii = 1:length(xds_morn.unit_names)
    all_units_amp_morn{ii,1} = mean(waveforms_morn{ii});
end
for ii = 1:length(xds_noon.unit_names)
    all_units_amp_noon{ii,1} = mean(waveforms_noon{ii});
end

% The mean of every point in the mean of the waveforms
all_units_amp_means_morn = zeros(length(xds_morn.unit_names),1);
all_units_amp_means_noon = zeros(length(xds_noon.unit_names),1);
for ii = 1:length(xds_morn.unit_names)
    all_units_amp_means_morn(ii) = mean((all_units_amp_morn{ii,1}));
end
for ii = 1:length(xds_noon.unit_names)
    all_units_amp_means_noon(ii) = mean((all_units_amp_noon{ii,1}));
end

% Find the differences between morning and afternoon
all_units_amp_difference = abs(all_units_amp_means_morn - all_units_amp_means_noon);

%% Determine if the unit is stable
% Pick a random 50 threshold crossings 1000 times
p_val_idx = zeros(1000,1);
for ii = 1:1000
    % Morning
    try
        rand_morn_idx = randperm(size(morn_perspike_amp,1), 50);
    catch
        disp('Too few spikes!')
        euc_dis = NaN;
        wave_p_val = NaN;
        return
    end
    rand_morn_perspike = morn_perspike_amp(rand_morn_idx);
    % Afternoon
    try
        rand_noon_idx = randperm(size(noon_perspike_amp,1), 50);
    catch
        disp('Too few spikes!')
        euc_dis = NaN;
        wave_p_val = NaN;
        return
    end
    rand_noon_perspike = noon_perspike_amp(rand_noon_idx);
    [~, wave_p_val] = ttest2(rand_morn_perspike, rand_noon_perspike);
    p_val_idx(ii) = wave_p_val;
end
wave_p_val = mean(p_val_idx);

%% Determine if the change in mean amplitude is more or less than the average change

mean_amp_difference = mean(all_units_amp_difference);

if amp_difference < mean_amp_difference
    fprintf("The change in waveform amplitude of %s is less than the mean of all units \n", ...
        string(xds_morn.unit_names{N_morn}));
end

if amp_difference > mean_amp_difference
        fprintf("The change in waveform amplitude of %s is greater than the mean of all units \n", ...
        string(xds_morn.unit_names{N_morn}));
end

%% Plotting 

if isequal(Plot_Figs, 1)

    %% Plot the morning and afternoon histograms
    figure
    hold on
    histogram(morn_perspike_amp, 'EdgeColor', 'k', 'FaceColor', [0.9290, 0.6940, 0.1250])
    histogram(noon_perspike_amp, 'EdgeColor', 'k', 'FaceColor', [.5 0 .5])
    title(sprintf('Average Amplitude Histogram - %s', ... 
        char(xds_morn.unit_names(N_morn))), 'FontSize', title_font_size)
    xlabel('Mean Waveform Amplitude', 'FontSize', label_font_size)
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

    % Reset the axis limits
    xlim([x_limits(1),x_limits(2)])
    ylim([y_limits(1),y_limits(2) + 1])

    % Annotation of the p-value
        if round(wave_p_val, 3) > 0
            legend_dims = [0.525 0.325 0.44 0.44];
            p_value_string = strcat('p =', {' '}, mat2str(round(wave_p_val, 3)));
            legend_string = {char(p_value_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = legend_font_size;
            ann_legend.FontName = font_name;
        end

        if isequal(round(wave_p_val, 3), 0)
            legend_dims = [0.525 0.325 0.44 0.44];
            p_value_string = strcat('p <', {' '}, '0.001');
            legend_string = {char(p_value_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = legend_font_size;
            ann_legend.FontName = font_name;
        end

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
    spike_min_morn = spike_times_morn / 60;
    spike_min_noon = spike_times_noon / 60;

    % Find the axis limits
    morn_ymax = amp_morn + 4*std_amp_morn;
    morn_ymin = amp_morn - 4*std_amp_morn;
    noon_ymax = amp_noon + 4*std_amp_noon;
    noon_ymin = amp_noon - 4*std_amp_noon;
    morn_xmax = max(spike_min_morn);
    noon_xmax = max(spike_min_noon);
    axis_limits = cat(1, morn_ymax, morn_ymin, noon_ymax, noon_ymin);
    max_YLim = max(axis_limits);
    min_YLim = min(axis_limits);

    % Size of the markers
    sz = 2;

    amp_timeline = figure;
    amp_timeline.Position = [200 200 750 350];
    % Morning
    morn_amp_plot = subplot(121); % (Top Left Plot)
    hold on
    scatter(spike_min_morn, morn_perspike_amp, sz, 'MarkerEdgeColor', 'k', ...
        'MarkerFaceColor', 'k');
    % Line indicating the morning mean
    line([0 max(spike_min_morn)],[amp_morn amp_morn], ... 
        'Color', 'red', 'LineWidth',1);
    % Lines indicating the morning standard deviation
    line([0 max(spike_min_morn)], [amp_morn + std_amp_morn, amp_morn + std_amp_morn], ...
        'Color', 'red', 'LineStyle','--', 'LineWidth',1);
    line([0 max(spike_min_morn)], [amp_morn - std_amp_morn, amp_morn - std_amp_morn], ...
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
    set(morn_amp_plot, 'Units', 'normalized')
    set(morn_amp_plot, 'Position', [0.1, 0.65, 0.4, 0.2]);

    xlabel('Morning', 'FontSize', label_font_size)

    % Afternoon
    noon_amp_plot = subplot(122); % (Top Right Plot)
    scatter(spike_min_noon, noon_perspike_amp, sz, 'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'k');
    % Lines indicating the afternoon mean
    line([min(spike_min_noon) max(spike_min_noon)],[amp_noon amp_noon], ... 
        'Color', 'red', 'LineWidth',1);
    % Lines indicating the afternoon standard deviation
    line([min(spike_min_noon) max(spike_min_noon)], [amp_noon + std_amp_noon, amp_noon + std_amp_noon], ...
        'Color', 'red', 'LineStyle','--', 'LineWidth',1);
    line([min(spike_min_noon) max(spike_min_noon)], [amp_noon - std_amp_noon, amp_noon - std_amp_noon], ...
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
    set(noon_amp_plot, 'Units', 'normalized')
    set(noon_amp_plot, 'Position', [0.525, 0.65, .4, .2]);
    
    xlabel('Afternoon', 'FontSize', label_font_size)

    % Annotation of the p_value
    if round(wave_p_val, 3) > 0
        legend_dims = [0.65 0.525 0.44 0.44];
        p_value_string = strcat('p =', {' '}, mat2str(round(wave_p_val, 3)));
        legend_string = {char(p_value_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = legend_font_size;
        ann_legend.FontName = font_name;
    end

    if isequal(round(wave_p_val, 3), 0)
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
    sgt = sgtitle(sprintf('Mean Waveform Amplitude - %s', ...
        char(xds_morn.unit_names(N_morn))), 'FontSize', (title_font_size + 5));
    % Set the common x-label
    common_x_label = axes(amp_timeline, 'visible', 'off');
    common_x_label.XLabel.Visible = 'on';
    xlabel(common_x_label, 'Time (min.)', 'FontSize', label_font_size);
    common_x_label.Position(2) = common_x_label.Position(2) + 0.5;

end

%% Print the peak to peak amplitudes
% Morning
fprintf("The average peak to peak amplitude of %s is %0.1f µV in the morning \n", ...
string(xds_morn.unit_names{N_morn}), peaktopeak_morn);
% Afternoon
fprintf("The average peak to peak amplitude of %s is %0.1f µV in the afternoon \n", ...
string(xds_noon.unit_names{N_noon}), peaktopeak_noon);

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
        %title '';
        saveas(gcf, fullfile(save_dir, char(fig_title)), Save_Figs) 
        close gcf
    end
end



