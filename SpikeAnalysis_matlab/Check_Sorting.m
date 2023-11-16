function Check_Sorting(xds_sorted, xds_unsorted, unit_name, Save_File)

%% What do you want plotted?
% 'PC1', 'PC2', 'PC3', 'Nonlin', or 'Time'
x_plot = 'PC1';
y_plot = 'Nonlin';

% Do you want to check your sorting or k-means clustering? ('Sort' or 'K_Means')
Sort_Check = 'Sort';

%% Basic settings, some variable extractions, & definitions

% Monkey
Monkey = xds_sorted.meta.monkey;

% Date
file_name = xds_sorted.meta.rawFileName;
xtra_info = extractAfter(file_name, '_');
Date = erase(file_name, strcat('_', xtra_info));

% Find the unit of interest
[N] = Find_Unit(xds_unsorted, unit_name);
if ~ischar(unit_name)
    % Identify the index of the unit
    current_unit = strcat(xds_unsorted.unit_names(N), '_');
else
    current_unit = strcat(unit_name, '_');
end

% Load unsorted waveforms
unsorted_waveforms = xds_unsorted.spike_waveforms{1,N};

% Load the unsorted spikes
unsorted_spikes = xds_unsorted.spikes{1,N};

% Find the associated sorted units
sorted_unit_idxs = find(contains(xds_sorted.unit_names, current_unit));

% Calculate the principle components
if contains(x_plot, 'PC') || contains(y_plot, 'PC')
    [~, Principal_Comps, ~] = pca(unsorted_waveforms);
end

% Calculate the nonlinear energy
if contains(x_plot, 'Nonlin') || contains(y_plot, 'Nonlin')
    unsorted_nonlin = mean(xds_unsorted.nonlin_waveforms{1,N}, 2);
end

%% Define the figure

% X Axis
if strcmp(x_plot, 'PC1')
    x_axis = Principal_Comps(:,1);
    x_label = 'Principal Component 1';
elseif strcmp(x_plot, 'PC2')
    x_axis = Principal_Comps(:,2);
    x_label = 'Principal Component 2';
elseif strcmp(x_plot, 'PC3')
    x_axis = Principal_Comps(:,3);
    x_label = 'Principal Component 3';
elseif strcmp(x_plot, 'Nonlin')
    x_axis = unsorted_nonlin;
    x_label = 'Nonlinear Energy';
elseif strcmp(x_plot, 'Time')
    x_axis = unsorted_spikes;
    x_label = 'Time';
end
% Y Axis
if strcmp(y_plot, 'PC1')
    y_axis = Principal_Comps(:,1);
    y_label = 'Principal Component 1';
elseif strcmp(y_plot, 'PC2')
    y_axis = Principal_Comps(:,2);
    y_label = 'Principal Component 2';
elseif strcmp(y_plot, 'PC3')
    y_axis = Principal_Comps(:,3);
    y_label = 'Principal Component 3';
elseif strcmp(y_plot, 'Nonlin')
    y_axis = unsorted_nonlin;
    y_label = 'Nonlinear Energy';
elseif strcmp(y_plot, 'Time')
    y_axis = unsorted_spikes;
    y_label = 'Time';
end

% Font specifications
label_font_size = 20;
title_font_size = 15;
% Figure size
figure_width = 800;
figure_height = 600;

marker_size = 2;

% Generate the figure
sort_figure = figure;
sort_figure.Position = [150 150 figure_width figure_height];
hold on

% Set the background as black
set(gca,'Color','k')

% Set the common title
if strcmp(Sort_Check, 'Sort')
    Fig_Title = strcat('Spike Sorting -', {' '}, char(xds_unsorted.unit_names(N)));
elseif strcmp(Sort_Check, 'K_Means')
    Fig_Title = strcat('K-Means Clusters -', {' '}, char(xds_unsorted.unit_names(N)));
end
sgtitle(Fig_Title, 'FontSize', (title_font_size + 5));
Fig_Title = strcat(Date, '_', Monkey, '_', Fig_Title);


% Axis Labels
ylabel(y_label, 'FontSize', label_font_size)
xlabel(x_label, 'FontSize', label_font_size)

% Only label every other tick
figure_axes = gca;
x_labels = string(figure_axes.XAxis.TickLabels);
y_labels = string(figure_axes.YAxis.TickLabels);
x_labels(1:end) = NaN;
y_labels(1:end) = NaN;
figure_axes.XAxis.TickLabels = x_labels;
figure_axes.YAxis.TickLabels = y_labels;
% Set ticks to outside
set(figure_axes,'TickDir','out');
% Remove the top and right tick marks
set(figure_axes,'box','off');

%% Plotting

% Plot the unsorted threshold crossings
scatter(x_axis, y_axis, ...
    marker_size, 'MarkerEdgeColor', [1 1 1], 'MarkerFaceColor', [1 1 1])

if strcmp(Sort_Check, 'K_Means')
    % How many clusters do you want?
    num_clusters = length(sorted_unit_idxs);
    [k_means_idx, ~] = kmeans(cat(2, x_axis, y_axis), num_clusters, 'Replicates', 5);
    
    for ii = 1:num_clusters
        % Define the plotting colors
        if ii == 1
            plot_color = 'y';
        elseif ii == 2
            plot_color = 'g';
        elseif ii == 3
            plot_color = 'c';
        elseif ii == 4
            plot_color = 'r';
        end

        scatter(x_axis(k_means_idx == ii), y_axis(k_means_idx == ii), ...
            marker_size, 'MarkerEdgeColor', plot_color, 'MarkerFaceColor', plot_color)
        scatter(x_axis(k_means_idx == ii), y_axis(k_means_idx == ii), ...
            marker_size, 'MarkerEdgeColor', plot_color, 'MarkerFaceColor', plot_color)
    end
end

if strcmp(Sort_Check, 'Sort')
    % Loop through each sorted unit
    for ii = 1:length(sorted_unit_idxs)
       
        % Define the plotting colors
        if ii == 1
            plot_color = 'y';
        elseif ii == 2
            plot_color = 'g';
        elseif ii == 3
            plot_color = 'c';
        elseif ii == 4
            plot_color = 'r';
        end
    
        % Load the sorted spikes
        sorted_spikes = xds_sorted.spikes{1,sorted_unit_idxs(ii)};
    
        % Find the indices of the sorted spikes
        sorted_spike_idxs = find(ismember(unsorted_spikes, sorted_spikes));
        
        % Plot the unsorted threshold crossings
        scatter(x_axis(sorted_spike_idxs), y_axis(sorted_spike_idxs), ...
            marker_size, 'MarkerEdgeColor', plot_color, 'MarkerFaceColor', plot_color)
    
    end
end

%% Save the file if selected
Save_Figs(Fig_Title, Save_File)
