function [fract_contam, ISI_modes] = InterstimulusInterval(xds, unit_name, Plot_Figs, Save_File)

%% Find the unit of interest
[N] = Find_Unit(xds, unit_name);

% Extracting the spikes of the designated unit
spikes = xds.spikes{N};

%% Catch possible sources of error
% If there is no unit of that name
if isempty(N)
    fprintf('%s does not exist \n', unit_name);
    fract_contam = NaN;
    ISI_modes = NaN;
    return
end

% If there are less than 1000 spikes
spike_limit = 1000;
if length(spikes) < spike_limit
    disp('Too few spikes!')
    fract_contam = NaN;
    ISI_modes = NaN;
    return
end

%% Basic settings, some variable extractions, & definitions

 % Font specifications
label_font_size = 25;
title_font_size = 15;
legend_font_size = 15;
font_name = 'Arial';

%% Some variable extraction & definitions

% Calculate the interstimulus intervals
ISI = diff(spikes);

% Convert the time to milliseconds
ISI = ISI*1000;

% Define the edges of the histogram (0 - 100 ms)
ISI_edges = (0:1:100);

ISI_modes = mode(ISI);

%% Calculate the number of refractory period violations

% Define the refractory period cutoff (2 ms)
refract_min = 2;
% Find the number of refractory period violations
refract_violations = length(find(ISI < refract_min));

% Find the percent of spikes this corresponds to
violation_ratio = refract_violations / length(spikes);

%% Calculate the fractional level of contamination

% Define the censored period (from Central's settings: 1 ms)
censored_period = 1;

% Define the number of spikes
spikes = length(spikes);

% Define the length of the experiment (in ms)
experiment_length = xds.time_frame(end)*1000;

% Find the average firing rate of the experiment
%avg_fr = (length(spikes_morn) + length(spikes_noon)) / (xds_morn.time_frame(end) + xds_noon.time_frame(end));

syms Fract_Contam

% The Sample Numbers From The Publication
%merged_refract_violations = 20;
%refract_min = 3;
%censored_period = 1;
%total_spikes = 10000;
%experiment_length = 1000000;

refract_eqn = refract_violations == 2*(refract_min - censored_period)*(spikes^2)*... 
    (1 - Fract_Contam)*(Fract_Contam/experiment_length);

fract_contam = solve(refract_eqn, Fract_Contam);

%% Plotting the histograms
if isequal(Plot_Figs, 1)

    % Histogram
    figure
    hold on

    % Set the title
    Fig_Title = strcat('Interstimulus Intervals -', {' '}, char(xds.unit_names(N)));
    if contains(xds.meta.rawFileName, 'Pre')
        Fig_Title = strcat(Fig_Title, {' '}, '(Morning)');
        hist_color = [0.9290, 0.6940, 0.1250];
    elseif contains(xds.meta.rawFileName, 'Post')
        Fig_Title = strcat(Fig_Title, {' '}, '(Afternoon)');
        hist_color = [.5 0 .5];
    else
        hist_color = 'k';
    end
    title(Fig_Title, 'FontSize', title_font_size)

    % Plot the histogram
    histogram(ISI, ISI_edges, 'EdgeColor', 'k', 'FaceColor', hist_color)

    % Axis Labels
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

    % Reset the axis limits
    xlim([x_limits(1),x_limits(2)])
    ylim([y_limits(1),y_limits(2) + 1])

    % Only label every other tick
    figure_axes = gca;
    x_labels = string(figure_axes.XAxis.TickLabels);
    x_labels(2:2:end) = NaN;
    figure_axes.XAxis.TickLabels = x_labels;
    y_labels = string(figure_axes.YAxis.TickLabels);
    y_labels(2:2:end) = NaN;
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
fprintf("%0.1f%% of the spikes in %s have an ISI less than %0.1f ms. \n", ...
violation_ratio*100, string(xds.unit_names{N}), refract_min);

if isreal(fract_contam)
    fprintf("The fractional level of contamination is %0.2f in %s \n", ...
        min(fract_contam), string(xds.unit_names{N}));
    fract_contam = double(min(fract_contam));
end

if ~isreal(fract_contam)
    fprintf("The fractional level of contamination in %s > 0.5 \n", ...
        string(xds.unit_names{N}));
    fract_contam = 1;
end









