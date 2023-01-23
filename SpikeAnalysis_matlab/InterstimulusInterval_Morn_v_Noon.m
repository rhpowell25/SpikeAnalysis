function [merged_fract_contam, ISI_modes] = InterstimulusInterval_Morn_v_Noon(xds_morn, xds_noon, unit_name, Plot_Figs, Save_Figs)

%% Load the excel file
if ~ischar(unit_name)

    [xds_output] = Find_Excel(xds_morn);

    %% Find the unit of interest

    unit = xds_output.unit_names(unit_name);

    %% Identify the index of the unit
    N_morn = find(strcmp(xds_morn.unit_names, unit));
    N_noon = find(strcmp(xds_noon.unit_names, unit));

end

if ischar(unit_name) && ~strcmp(unit_name, 'All')
    N_morn = find(strcmp(xds_morn.unit_names, unit_name));
    N_noon = find(strcmp(xds_noon.unit_names, unit_name));
end

%% Loop through the number of units selected

% Define the output variables
ISI_modes = NaN(length(N_morn),3);
merged_fract_contam = NaN(length(N_morn),1);

 % Font specifications
label_font_size = 20;
title_font_size = 15;
legend_font_size = 15;
font_name = 'Arial';

if ~isequal(Save_Figs, 0)
    % Do you want a save title or blank title (1 = save_title, 0 = blank)
    Fig_Save_Title = 0;
end
   
for uu = 1:length(N_morn)

    %% If The Unit Doesn't Exist

    if isempty(N_morn(uu)) || isempty(N_noon(uu))
        fprintf('%s does not exist \n', unit_name);
        ISI_modes(uu,:) = NaN;
        merged_fract_contam(uu) = NaN;
        continue
    end

    %% Some variable extraction & definitions

    % Extracting the spikes of the designated unit
    spikes_morn = xds_morn.spikes{N_morn(uu)};
    spikes_noon = xds_noon.spikes{N_noon(uu)};

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

    ISI_modes(uu,1) = mode(ISI_morn);
    ISI_modes(uu,2) = mode(ISI_noon);
    ISI_modes(uu,3) = mode(ISI_merged);

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
        title(sprintf('Interstimulus Intervals - %s (Merged)', ... 
            char(xds_morn.unit_names(N_morn(uu)))), 'FontSize', title_font_size)
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

    end

    %% Print the percentages & fractional contamination

    % Merged
    fprintf("%0.1f%% of the spikes in %s have an ISI less than %0.1f ms. \n", ...
    merged_violation_ratio*100, string(xds_noon.unit_names{N_noon(uu)}), refract_min);

    % Merged
    if ~isreal(fract_contam)
        fprintf("The fractional level of contamination in %s > 0.5 \n", ...
            string(xds_noon.unit_names{N_noon(uu)}));
        merged_fract_contam(uu) = 1;
    end
    if isreal(fract_contam)
        fprintf("The fractional level of contamination is %0.2f in %s \n", ...
            min(fract_contam), string(xds_noon.unit_names{N_noon(uu)}));
        merged_fract_contam(uu) = double(min(fract_contam));
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







