function [peaktopeak_amp, spike_width, repol_time, wave_p_val] = ...
    WaveShapes_Morn_v_Noon(xds_morn, xds_noon, unit_name, Plot_Figs, Save_File)

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
    peaktopeak_amp = NaN;
    spike_width = NaN;
    repol_time = NaN;
    wave_p_val = NaN;
    return
end

% If there are less than 1000 spikes
spike_limit = 1000;
if length(spikes_morn) < spike_limit || length(spikes_noon) < spike_limit
    disp('Too few spikes!')
    peaktopeak_amp = NaN;
    spike_width = NaN;
    repol_time = NaN;
    wave_p_val = NaN;
    return
end

%% Some variable extraction & definitions

waveform_length = 1.6; % Time (ms.) 

% Extracting the waveforms of the designated unit
waveforms_morn = xds_morn.spike_waveforms{N_morn};
waveforms_noon = xds_noon.spike_waveforms{N_noon};

% Extracting the spike times of the designated unit
spike_times_morn = xds_morn.spikes{N_morn};
spike_times_noon = xds_noon.spikes{N_noon};

% Font & plotting specifications
[Plot_Params] = Plot_Specs;

%% Calculating the means and standard deviations

% Calculating the means
%morn_amp_mean = mean(waveforms_morn, 1);
%noon_amp_mean = mean(waveforms_noon, 1);
mean_amp = mean([waveforms_morn; waveforms_noon], 1);
morn_perspike_amp = mean(waveforms_morn, 2);
noon_perspike_amp = mean(waveforms_noon, 2);

% Calculate the Euclidean distance
%euc_dis = sqrt(sum((morn_amp_mean - noon_amp_mean) .^2));

%% Find the change in mean amplitude between the two sessions
% Finding the mean and difference between morning & afternoon
amp_morn = mean(morn_perspike_amp);
amp_noon = mean(noon_perspike_amp);
std_amp_morn = std(morn_perspike_amp);
std_amp_noon = std(noon_perspike_amp);

%% Calculate the peak to peak amplitude & spike width
waveform_crest = max(mean_amp);
crest_idx = find(mean_amp == waveform_crest);
waveform_trough = min(mean_amp);
trough_idx = find(mean_amp == waveform_trough);

% Peak to peak amplitudes
peaktopeak_amp = abs(waveform_crest) + abs(waveform_trough);

% Length of each bin
bin_size = waveform_length/length(mean_amp); % Time (ms.)

% Spike width
spike_width = (crest_idx - trough_idx)*bin_size; % Time (ms.)

% Repolarization time
detrend_amp = detrend(mean_amp(crest_idx:end),1);
deriv_amp = gradient(detrend_amp) ./ gradient(crest_idx:length(mean_amp));
inflection_idx = round(interp1(deriv_amp, (crest_idx:length(mean_amp)), 0));
repol_time = (inflection_idx - crest_idx)*bin_size;

%% Determine if the unit is stable
% Pick a random 50 threshold crossings 1000 times
p_val_idx = zeros(1000,1);
for ii = 1:1000
    % Morning
    rand_morn_idx = randperm(size(morn_perspike_amp,1), 50);
    rand_morn_perspike = morn_perspike_amp(rand_morn_idx);
    % Afternoon
    rand_noon_idx = randperm(size(noon_perspike_amp,1), 50);
    rand_noon_perspike = noon_perspike_amp(rand_noon_idx);
    [~, wave_p_val] = ttest2(rand_morn_perspike, rand_noon_perspike);
    p_val_idx(ii) = wave_p_val;
end
wave_p_val = mean(p_val_idx);

%% Plotting 

if isequal(Plot_Figs, 1)

    %% Plot the morning and afternoon histograms
    figure
    hold on
    histogram(morn_perspike_amp, 'EdgeColor', 'k', 'FaceColor', [0.9290, 0.6940, 0.1250])
    histogram(noon_perspike_amp, 'EdgeColor', 'k', 'FaceColor', [.5 0 .5])
    Fig_Title = sprintf('Average Amplitude Histogram - %s', ... 
        char(xds_morn.unit_names(N_morn)));
    title(Fig_Title, 'FontSize', Plot_Params.title_font_size)
    xlabel('Mean Waveform Amplitude', 'FontSize', Plot_Params.label_font_size)
    ylabel('Counts', 'FontSize', Plot_Params.label_font_size)

    % Collect the current axis limits
    y_limits = ylim;
    x_limits = xlim;

    % Plot dummy points for the legend
    dummy_yellow = scatter(-1, -1, 's', 'filled', 'Color', [0.9290, 0.6940, 0.1250]);
    dummy_purple = scatter(-1.5, -1.5, 's', 'filled', 'Color', [.5 0 .5]);
    % Legend
    legend([dummy_yellow, dummy_purple], ... 
        {'Morning', 'Afternoon'}, 'FontSize', Plot_Params.legend_size, 'Location', 'northeast')
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
            ann_legend.FontSize = Plot_Params.legend_size;
            ann_legend.FontName = Plot_Params.font_name;
        end

        if isequal(round(wave_p_val, 3), 0)
            legend_dims = [0.525 0.325 0.44 0.44];
            p_value_string = strcat('p <', {' '}, '0.001');
            legend_string = {char(p_value_string)};
            ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
                'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
                'EdgeColor','none', 'horizontalalignment', 'center');
            ann_legend.FontSize = Plot_Params.legend_size;
            ann_legend.FontName = Plot_Params.font_name;
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
    set(figure_axes,'fontname', Plot_Params.font_name);

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
    set(figure_axes,'fontname', Plot_Params.font_name);

    % Bring the morning plot to the right
    set(morn_amp_plot, 'Units', 'normalized')
    set(morn_amp_plot, 'Position', [0.1, 0.65, 0.4, 0.2]);

    xlabel('Morning', 'FontSize', Plot_Params.label_font_size)

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
    set(figure_axes,'fontname', Plot_Params.font_name);

    % Bring the afternoon plot to the left
    set(noon_amp_plot, 'Units', 'normalized')
    set(noon_amp_plot, 'Position', [0.525, 0.65, .4, .2]);
    
    xlabel('Afternoon', 'FontSize', Plot_Params.label_font_size)

    % Annotation of the p_value
    if round(wave_p_val, 3) > 0
        legend_dims = [0.65 0.525 0.44 0.44];
        p_value_string = strcat('p =', {' '}, mat2str(round(wave_p_val, 3)));
        legend_string = {char(p_value_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = Plot_Params.legend_size;
        ann_legend.FontName = Plot_Params.font_name;
    end

    if isequal(round(wave_p_val, 3), 0)
        legend_dims = [0.65 0.525 0.44 0.44];
        p_value_string = strcat('p <', {' '}, '0.001');
        legend_string = {char(p_value_string)};
        ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
            'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
            'EdgeColor','none', 'horizontalalignment', 'center');
        ann_legend.FontSize = Plot_Params.legend_size;
        ann_legend.FontName = Plot_Params.font_name;
    end

    % Set the common title
    sgtitle(sprintf('Mean Waveform Amplitude - %s', ...
        char(xds_morn.unit_names(N_morn))), 'FontSize', (Plot_Params.title_font_size + 5));
    % Set the common x-label
    common_x_label = axes(amp_timeline, 'visible', 'off');
    common_x_label.XLabel.Visible = 'on';
    xlabel(common_x_label, 'Time (min.)', 'FontSize', Plot_Params.label_font_size);
    common_x_label.Position(2) = common_x_label.Position(2) + 0.5;

    %% Save the file if selected
    Save_Figs(Fig_Title, Save_File)

end


