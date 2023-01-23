function Unit_Summary(xds_morn, xds_noon, unit_name, Save_Figs)

%% Load the excel file
if ~ischar(unit_name)

    [xds_output] = Find_Excel(xds_morn);

    %% Find the unit of interest

    unit = xds_output.unit_names(unit_name);

    %% Identify the index of the unit
    N = find(strcmp(xds_morn.unit_names, unit));

end

if ischar(unit_name) && ~strcmp(unit_name, 'All')
    N = find(strcmp(xds_morn.unit_names, unit_name));
end

%% Define the figure

% Font specifications
label_font_size = 12;
title_font_size = 15;
legend_font_size = 8;
font_name = 'Arial';
figure_width = 900;
figure_height = 700;

axis_expansion = 75;

% Define the annotation heights & widths
first_row_height = 0.45;
second_row_height = 0.16;
first_column_width = 0.09;
third_column_width = 0.65;

unit_sum_fig = figure;
unit_sum_fig.Position = [200 50 figure_width figure_height];
hold on

% Set the common title
fig_title = strcat('Unit Summary -', {' '}, char(xds_morn.unit_names(N)));
sgtitle(fig_title, 'FontSize', (title_font_size + 5));

%% Plot random 100 morning waverforms

% Extracting the waveforms of the designated unit
morn_unit_waveforms = xds_morn.spike_waveforms{N};

morn_rand_wave_idx = randperm(size(morn_unit_waveforms,1),100);
morn_rand_waves = morn_unit_waveforms(morn_rand_wave_idx,:);

subplot(3, 3, 1)
hold on
plot(morn_rand_waves', 'k')

% Set the title
title('Random Waveforms', 'FontSize', title_font_size)

% Axis Labels
ylabel({'Morning'; 'Amplitude (µV)'}, 'FontSize', label_font_size)

% Annotation of the n-count
legend_dims = [first_column_width first_row_height 0.44 0.44];
n_count_string = strcat('n =', {' '}, mat2str(length(morn_rand_waves)));
legend_string = {char(n_count_string)};
ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
    'FitBoxToText', 'on', 'verticalalignment', 'top', ...
    'EdgeColor','none', 'horizontalalignment', 'center');
ann_legend.FontSize = legend_font_size;
ann_legend.FontName = font_name;

% Collect the current axis limits
y_limits = ylim;
x_limits = xlim;

% Reset the axis limits
xlim([x_limits(1), x_limits(2)])
ylim([y_limits(1) - axis_expansion, y_limits(2) + axis_expansion])

% Only label every other tick
subplot_axes = gca;
x_labels = string(subplot_axes.XAxis.TickLabels);
y_labels = string(subplot_axes.YAxis.TickLabels);
x_labels(1:end) = NaN;
y_labels(2:2:end) = NaN;
subplot_axes.XAxis.TickLabels = x_labels;
subplot_axes.YAxis.TickLabels = y_labels;
% Set ticks to outside
set(subplot_axes,'TickDir','out');
% Remove the top and right tick marks
set(subplot_axes,'box','off');
% Set The Font
set(subplot_axes,'fontname',font_name);

%% Plot random 100 afternoon waverforms

% Extracting the waveforms of the designated unit
noon_unit_waveforms = xds_noon.spike_waveforms{N};

noon_rand_wave_idx = randperm(size(noon_unit_waveforms,1),100);
noon_rand_waves = noon_unit_waveforms(noon_rand_wave_idx,:);

subplot(3, 3, 4)
hold on
plot(noon_rand_waves', 'k')

% Axis Labels
ylabel({'Afternoon'; 'Amplitude (µV)'}, 'FontSize', label_font_size)
xlabel('Time', 'FontSize', label_font_size)

% Annotation of the n-count
legend_dims = [first_column_width second_row_height 0.44 0.44];
n_count_string = strcat('n =', {' '}, mat2str(length(morn_rand_waves)));
legend_string = {char(n_count_string)};
ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ...
    'FitBoxToText', 'on', 'verticalalignment', 'top', ...
    'EdgeColor','none', 'horizontalalignment', 'center');
ann_legend.FontSize = legend_font_size;
ann_legend.FontName = font_name;

% Reset the axis limits
xlim([x_limits(1), x_limits(2)])
ylim([y_limits(1) - axis_expansion, y_limits(2) + axis_expansion])

% Only label every other tick
subplot_axes = gca;
x_labels = string(subplot_axes.XAxis.TickLabels);
y_labels = string(subplot_axes.YAxis.TickLabels);
x_labels(2:2:end) = NaN;
y_labels(2:2:end) = NaN;
subplot_axes.XAxis.TickLabels = x_labels;
subplot_axes.YAxis.TickLabels = y_labels;
% Set ticks to outside
set(subplot_axes,'TickDir','out');
% Remove the top and right tick marks
set(subplot_axes,'box','off');
% Set The Font
set(subplot_axes,'fontname',font_name);

%% Plotting the morning mean & standard deviations

% Calculating the means
morn_amp_mean = mean(morn_unit_waveforms, 1);

% Calculating their standard deviations
morn_standard_dev = std(morn_unit_waveforms);

subplot(3, 3, 2)
hold on

% Plotting
plot(morn_amp_mean,'k')
plot(morn_amp_mean + morn_standard_dev,'r')
plot(morn_amp_mean - morn_standard_dev,'r')

% Set the title
title('Average Waveforms', 'FontSize', title_font_size)

% Reset the axis limits
xlim([x_limits(1), x_limits(2)])
ylim([y_limits(1) - axis_expansion, y_limits(2) + axis_expansion])

% Only label every other tick
subplot_axes = gca;
x_labels = string(subplot_axes.XAxis.TickLabels);
y_labels = string(subplot_axes.YAxis.TickLabels);
x_labels(1:end) = NaN;
y_labels(1:end) = NaN;
subplot_axes.XAxis.TickLabels = x_labels;
subplot_axes.YAxis.TickLabels = y_labels;
% Set ticks to outside
set(subplot_axes,'TickDir','out');
% Remove the top and right tick marks
set(subplot_axes,'box','off');
% Set The Font
set(subplot_axes,'fontname', font_name);

%% Plotting the afternoon mean & standard deviations

% Calculating the means
noon_amp_mean = mean(noon_unit_waveforms, 1);

subplot(3, 3, 5)
hold on

% Plotting
plot(noon_amp_mean,'k')
plot(noon_amp_mean + morn_standard_dev,'r')
plot(noon_amp_mean - morn_standard_dev,'r')

% Reset the axis limits
xlim([x_limits(1), x_limits(2)])
ylim([y_limits(1) - axis_expansion, y_limits(2) + axis_expansion])

% Axis Labels
xlabel('Time', 'FontSize', label_font_size)

% Only label every other tick
subplot_axes = gca;
x_labels = string(subplot_axes.XAxis.TickLabels);
y_labels = string(subplot_axes.YAxis.TickLabels);
x_labels(2:2:end) = NaN;
y_labels(1:end) = NaN;
subplot_axes.XAxis.TickLabels = x_labels;
subplot_axes.YAxis.TickLabels = y_labels;
% Set ticks to outside
set(subplot_axes,'TickDir','out');
% Remove the top and right tick marks
set(subplot_axes,'box','off');
% Set The Font
set(subplot_axes,'fontname', font_name);

%% Plotting the morning histogram

% Extracting the spikes of the designated unit
morn_spikes = xds_morn.spikes{N};
morn_num_spikes = length(morn_spikes);

% Calculate the interstimulus intervals
morn_ISI = diff(morn_spikes);

% Convert the time to milliseconds
morn_ISI = morn_ISI*1000;

% Define the edges of the histogram (0 - 100 ms)
ISI_edges = (0:1:100);

% Define the refractory period cutoff (2 ms)
refract_min = 2;
% Define the censored period (from Central's settings: 1 ms)
censored_period = 1;
% Find the number of refractory period violations
morn_refract_violations = length(find(morn_ISI < refract_min));

% Define the length of the experiment (in ms)
morn_experiment_length = xds_morn.time_frame(end)*1000;

syms Fract_Contam

% The Sample Numbers From The Publication
%merged_refract_violations = 20;
%refract_min = 3;
%censored_period = 1;
%total_spikes = 10000;
%experiment_length = 1000000;

refract_eqn = morn_refract_violations == 2*(refract_min - censored_period)*(morn_num_spikes^2)*... 
    (1 - Fract_Contam)*(Fract_Contam/morn_experiment_length);

morn_fract_contam = solve(refract_eqn, Fract_Contam);

subplot(3, 3, 3)
hold on

% Set the title
title('Interstimulus Intervals', 'FontSize', title_font_size)

% Plot the histogram
histogram(morn_ISI, ISI_edges, 'EdgeColor', 'k', 'FaceColor', [0.9290, 0.6940, 0.1250])

% Axis Labels
ylabel('Counts', 'FontSize', label_font_size)

% Collect the current axis limits
y_limits = ylim;
x_limits = xlim;

% Annotation of the fractional contamination
if isreal(morn_fract_contam)
    legend_dims = [third_column_width first_row_height 0.44 0.44];
    fract_contam_string = strcat('f_{1}^{p} =', {' '}, char(round(min(morn_fract_contam), 2)));
    legend_string = {char(fract_contam_string)};
    ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
        'FitBoxToText', 'on', 'EdgeColor','none', ... 
        'verticalalignment', 'top', 'horizontalalignment', 'center');
    ann_legend.FontSize = legend_font_size;
    ann_legend.FontName = font_name;
end
if ~isreal(morn_fract_contam)
    legend_dims = [third_column_width first_row_height 0.44 0.44];
    fract_contam_string = strcat('f_{1}^{p} >', {' '}, '0.5');
    legend_string = {char(fract_contam_string)};
    ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
        'FitBoxToText', 'on', 'EdgeColor','none', ... 
        'verticalalignment', 'top', 'horizontalalignment', 'center');
    ann_legend.FontSize = legend_font_size;
    ann_legend.FontName = font_name;
end

% Reset the axis limits
xlim([x_limits(1),x_limits(2)])
ylim([y_limits(1),y_limits(2) + 1])

% Only label every other tick
subplot_axes = gca;
x_labels = string(subplot_axes.XAxis.TickLabels);
x_labels(2:2:end) = NaN;
subplot_axes.XAxis.TickLabels = x_labels;
y_labels = string(subplot_axes.YAxis.TickLabels);
y_labels(2:2:end) = NaN;
subplot_axes.YAxis.TickLabels = y_labels;
% Set ticks to outside
set(subplot_axes,'TickDir','out');
% Remove the top and right tick marks
set(subplot_axes,'box','off')
% Set The Font
set(subplot_axes,'fontname', font_name);

%% Plotting the afternoon histogram

% Extracting the spikes of the designated unit
noon_spikes = xds_noon.spikes{N};
noon_num_spikes = length(noon_spikes);

% Calculate the interstimulus intervals
noon_ISI = diff(noon_spikes);

% Convert the time to milliseconds
noon_ISI = noon_ISI*1000;

% Find the number of refractory period violations
noon_refract_violations = length(find(noon_ISI < refract_min));

% Define the length of the experiment (in ms)
noon_experiment_length = xds_noon.time_frame(end)*1000;

syms Fract_Contam

% The Sample Numbers From The Publication
%merged_refract_violations = 20;
%refract_min = 3;
%censored_period = 1;
%total_spikes = 10000;
%experiment_length = 1000000;

refract_eqn = noon_refract_violations == 2*(refract_min - censored_period)*(noon_num_spikes^2)*... 
    (1 - Fract_Contam)*(Fract_Contam/noon_experiment_length);

noon_fract_contam = solve(refract_eqn, Fract_Contam);

subplot(3, 3, 6)
hold on

% Plot the histogram
histogram(noon_ISI, ISI_edges, 'EdgeColor', 'k', 'FaceColor', [.5 0 .5])

% Axis Labels
ylabel('Counts', 'FontSize', label_font_size)
xlabel('ISI (ms)', 'FontSize', label_font_size)

% Collect the current axis limits
y_limits = ylim;
x_limits = xlim;

% Annotation of the fractional contamination
if isreal(noon_fract_contam)
    legend_dims = [third_column_width second_row_height 0.44 0.44];
    fract_contam_string = strcat('f_{1}^{p} =', {' '}, char(round(min(noon_fract_contam), 2)));
    legend_string = {char(fract_contam_string)};
    ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
        'FitBoxToText', 'on', 'EdgeColor','none', ... 
        'verticalalignment', 'top', 'horizontalalignment', 'center');
    ann_legend.FontSize = legend_font_size;
    ann_legend.FontName = font_name;
end
if ~isreal(noon_fract_contam)
    legend_dims = [third_column_width second_row_height 0.44 0.44];
    fract_contam_string = strcat('f_{1}^{p} >', {' '}, '0.5');
    legend_string = {char(fract_contam_string)};
    ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
        'FitBoxToText', 'on', 'EdgeColor','none', ... 
        'verticalalignment', 'top', 'horizontalalignment', 'center');
    ann_legend.FontSize = legend_font_size;
    ann_legend.FontName = font_name;
end

% Reset the axis limits
xlim([x_limits(1),x_limits(2)])
ylim([y_limits(1),y_limits(2) + 1])

% Only label every other tick
subplot_axes = gca;
x_labels = string(subplot_axes.XAxis.TickLabels);
x_labels(2:2:end) = NaN;
subplot_axes.XAxis.TickLabels = x_labels;
y_labels = string(subplot_axes.YAxis.TickLabels);
y_labels(2:2:end) = NaN;
subplot_axes.YAxis.TickLabels = y_labels;
% Set ticks to outside
set(subplot_axes,'TickDir','out');
% Remove the top and right tick marks
set(subplot_axes,'box','off')
% Set The Font
set(subplot_axes,'fontname', font_name);

%% Calculate the per-spike nonlinear energy

% Extracting the nonlinear waveforms of the all units
nonlin_waveforms_morn = xds_morn.nonlin_waveforms;
nonlin_waveforms_noon = xds_noon.nonlin_waveforms;

% Nonlinear waveform means
morn_perspike_nonlin = mean(nonlin_waveforms_morn{N}, 2);
noon_perspike_nonlin = mean(nonlin_waveforms_noon{N}, 2);

% Finding the mean & difference between morning & afternoon
nonlin_morn = mean(morn_perspike_nonlin);
nonlin_noon = mean(noon_perspike_nonlin);
std_nonlin_morn = std(morn_perspike_nonlin);
std_nonlin_noon = std(noon_perspike_nonlin);

% Convert time into minutes
spike_min_morn = morn_spikes / 60;
spike_min_noon = noon_spikes / 60;

% Find the axis limits
morn_ymax = nonlin_morn + 5*std_nonlin_morn;
morn_ymin = nonlin_morn - 5*std_nonlin_morn;
noon_ymax = nonlin_noon + 5*std_nonlin_noon;
noon_ymin = nonlin_noon - 5*std_nonlin_noon;
morn_xmax = max(spike_min_morn);
noon_xmax = max(spike_min_noon);
axis_limits = cat(1, morn_ymax, morn_ymin, noon_ymax, noon_ymin);
max_YLim = max(axis_limits);
min_YLim = min(axis_limits);

% Size of the markers
sz = 2;

% Determine if the unit is well sorted 
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

subplot(3, 3, [7 7.5])
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
subplot_axes = gca;
x_labels = string(subplot_axes.XAxis.TickLabels);
morn_end_tick = find(x_labels == num2str(round(max(spike_min_morn))));
noon_ticks = str2double(x_labels(morn_end_tick + 1:end)) - str2double(x_labels(morn_end_tick));
x_labels(morn_end_tick + 1:end) = num2str(noon_ticks);
y_labels = string(subplot_axes.YAxis.TickLabels);
% Only label every other tick
x_labels(1:2:end) = NaN;
y_labels(2:2:end) = NaN;
subplot_axes.XAxis.TickLabels = x_labels;
subplot_axes.YAxis.TickLabels = y_labels;
% Set ticks to outside
set(subplot_axes,'TickDir','out');
% Remove the top and right tick marks
set(subplot_axes,'box','off');
% Set The Font
set(subplot_axes,'fontname', font_name);

xlabel('Morning', 'FontSize', label_font_size + 5)

subplot(3, 3, [8.5 9])
hold on

% Afternoon
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

xlabel('Afternoon', 'FontSize', label_font_size + 5)

% Set the ticks for afternoon
subplot_axes = gca;
x_labels = string(subplot_axes.XAxis.TickLabels);
morn_end_tick = find(x_labels == num2str(round(max(spike_min_morn))));
noon_ticks = str2double(x_labels(morn_end_tick + 1:end)) - str2double(x_labels(morn_end_tick));
x_labels(morn_end_tick + 1:end) = num2str(noon_ticks);
set(subplot_axes, 'YTickLabel', []);
% Only label every other tick
x_labels(1:2:end) = NaN;
subplot_axes.XAxis.TickLabels = x_labels;
% Set ticks to outside
set(subplot_axes,'TickDir','out');
% Remove the top and right tick marks
set(subplot_axes,'box','off');
% Set The Font
set(subplot_axes,'fontname', font_name);

% Annotation of the p_value
if round(nonlin_p_val, 3) > 0
    legend_dims = [third_column_width 0.1 0.44 0.44];
    p_value_string = strcat('p =', {' '}, mat2str(round(nonlin_p_val, 3)));
    legend_string = {char(p_value_string)};
    ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
        'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
        'EdgeColor','none', 'horizontalalignment', 'center');
    ann_legend.FontSize = legend_font_size;
    ann_legend.FontName = font_name;
end

if isequal(round(nonlin_p_val, 3), 0)
    legend_dims = [third_column_width 0.1 0.44 0.44];
    p_value_string = strcat('p <', {' '}, '0.001');
    legend_string = {char(p_value_string)};
    ann_legend = annotation('textbox', legend_dims, 'String', legend_string, ... 
        'FitBoxToText', 'on', 'verticalalignment', 'top', ... 
        'EdgeColor','none', 'horizontalalignment', 'center');
    ann_legend.FontSize = legend_font_size;
    ann_legend.FontName = font_name;
end

% Set the common x label
fig_handle = axes(unit_sum_fig,'visible','off'); 
fig_handle.XLabel.Visible='on';
xlabel(fig_handle,'Time (min.)', 'fontsize', 20);
fig_handle.Position = [0.13 0.08 0.775 0.815];

%% Define the save directory & save the figures
if ~isequal(Save_Figs, 0)
    save_dir = 'C:\Users\rhpow\Desktop\';
    for ii = 1:length(findobj('type','figure'))
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


