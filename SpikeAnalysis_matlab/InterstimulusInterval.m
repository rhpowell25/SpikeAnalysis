function [fract_contam, ISI_modes] = InterstimulusInterval(xds, unit_name, Plot_Figs, Save_Figs)

%% Load the excel file
if ~ischar(unit_name)

    [xds_output] = Find_Excel(xds);

    %% Find the unit of interest

    unit = xds_output.unit_names(unit_name);

    %% Identify the index of the unit
    N = find(strcmp(xds.unit_names, unit));

end

if ischar(unit_name) && ~strcmp(unit_name, 'All')
    N = find(strcmp(xds.unit_names, unit_name));
end

%% Basic settings, some variable extractions, & definitions

% Define the output variables
ISI_modes = NaN(length(N),1);
fract_contam = NaN(length(N),1);

 % Font specifications
label_font_size = 25;
title_font_size = 15;
legend_font_size = 15;
font_name = 'Arial';

if ~isequal(Save_Figs, 0)
    % Do you want a save title or blank title (1 = save_title, 0 = blank)
    Fig_Save_Title = 0;
end
   
for uu = 1:length(N)

    %% If The Unit Doesn't Exist

    if isempty(N(uu))
        fprintf('%s does not exist \n', unit_name);
        ISI_modes(uu,:) = NaN;
        fract_contam(uu) = NaN;
        continue
    end

    %% Some variable extraction & definitions

    % Extracting the spikes of the designated unit
    spikes = xds.spikes{N(uu)};

    % Calculate the interstimulus intervals
    ISI = diff(spikes);

    % Convert the time to milliseconds
    ISI = ISI*1000;

    % Define the edges of the histogram (0 - 100 ms)
    ISI_edges = (0:1:100);

    ISI_modes(uu,1) = mode(ISI);

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
        hist_title = strcat('Interstimulus Intervals -', {' '}, char(xds.unit_names(N)));
        if contains(xds.meta.rawFileName, 'Pre')
            hist_title = strcat(hist_title, {' '}, '(Morning)');
            hist_color = [0.9290, 0.6940, 0.1250];
        elseif contains(xds.meta.rawFileName, 'Post')
            hist_title = strcat(hist_title, {' '}, '(Afternoon)');
            hist_color = [.5 0 .5];
        else
            hist_color = 'k';
        end
        title(hist_title, 'FontSize', title_font_size)

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

    end

    %% Print the percentages & fractional contamination
    fprintf("%0.1f%% of the spikes in %s have an ISI less than %0.1f ms. \n", ...
    violation_ratio*100, string(xds.unit_names{N(uu)}), refract_min);

    if isreal(fract_contam)
        fprintf("The fractional level of contamination is %0.2f in %s \n", ...
            min(fract_contam), string(xds.unit_names{N(uu)}));
        fract_contam = double(min(fract_contam));
    end

    if ~isreal(fract_contam)
        fprintf("The fractional level of contamination in %s > 0.5 \n", ...
            string(xds.unit_names{N(uu)}));
        fract_contam = 1;
    end

end % End of the unit loop

%% Define the save directory & save the figures
if ~isequal(Save_Figs, 0)
    save_dir = 'C:\Users\rhpow\Desktop\';
    for ii = 1:length(findobj('type','figure'))
        fig_info = get(gca,'title');
        fig_title = get(fig_info, 'string');
        fig_title = strrep(fig_title, ':', '');
        fig_title = strrep(fig_title, 'vs.', 'vs');
        fig_title = strrep(fig_title, 'mg.', 'mg');
        fig_title = strrep(fig_title, 'kg.', 'kg');
        fig_title = strrep(fig_title, '.', '_');
        fig_title = strrep(fig_title, '/', '_');
        if isequal(Fig_Save_Title, 0)
            title '';
        end
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







