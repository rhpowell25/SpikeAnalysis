function Map_File_Figs(xds, Plot_Specs, Save_File)

%% Find & load the map file
Monkey = xds.meta.monkey;

% Date
file_name = xds.meta.rawFileName;
xtra_info = extractAfter(file_name, '_');
Date = erase(file_name, strcat('_', xtra_info));

% Find the map file
map_dir = strcat('C:\Users\rhpow\Documents\Work\Northwestern\Monkey_Data\', Monkey, '\');
dir_files = dir(strcat(map_dir, '*.cmp'));

% Name the map file
map_name = dir_files.name(1:end-4);
disp(map_name);

% Load the map file
map_file = readtable(strcat(map_dir, map_name, '.cmp'), 'Filetype', 'text', 'VariableNamingRule', 'Preserve');

%% Find the rows, columbs, & unit names

unit_columns = zeros(length(map_file.(1)), 1);
unit_rows = zeros(length(map_file.(1)), 1);
unit_names = struct([]);
cc = 1;
for uu = 1:length(map_file.(1))
    if isnan(map_file.(1)(uu))
        unit_columns(uu,:) = [];
        unit_rows(uu,:) = [];
        continue
    end
    unit_columns(cc,1) = map_file.(1)(uu);
    unit_rows(cc,1) = map_file.(2)(uu);
    unit_names{cc,1} = map_file.(5)(uu);
    cc = cc + 1;
end

% Find the width & height of the map file
unit_rows = unit_rows + 1;
unit_columns = unit_columns + 1;
map_height = max(unit_rows);
map_width = max(unit_columns);

%% Plot the map file

% Figure size
figure_width = 600;
figure_height = 600;

% Generate the figure
map_figure = figure;
map_figure.Position = [150 150 figure_width figure_height];

% Set the title
Fig_Title = strcat(Date, '_', Monkey, '_', Plot_Specs);
if contains(xds.meta.rawFileName, 'Pre')
    Fig_Title = strcat(Fig_Title, '_morn');
end
if contains(xds.meta.rawFileName, 'Post')
    Fig_Title = strcat(Fig_Title, '_noon');
end
sgtitle(Fig_Title, 'Interpreter', 'None')

% If you're making the colormap
if strcmp(Plot_Specs, 'Colormap')
    color_map = turbo;
    total_crossings = zeros(map_height*map_width,1);
    peak_firing_rate = zeros(map_height*map_width,1);
    for ii = 1:length(total_crossings)
        current_unit = strcat('elec', num2str(ii), '_');
        % Find the associated sorted units
        sorted_unit_idxs = find(contains(xds.unit_names, current_unit));
        if isempty(sorted_unit_idxs)
            total_crossings(ii) = NaN;
            peak_firing_rate(ii) = NaN;
            continue
        end
        temp_crossings = zeros(length(sorted_unit_idxs),1);
        temp_ISI_FR = [];
        for jj = 1:length(temp_crossings)
            temp_crossings(jj) = length(xds.spikes{sorted_unit_idxs(jj)});
            temp_ISI_FR = cat(1, temp_ISI_FR, (1./diff(xds.spikes{sorted_unit_idxs(jj)})));
        end
        peak_firing_rate(ii) = prctile(temp_ISI_FR, 90);
        total_crossings(ii) = sum(temp_crossings);
    end
    % Normalize the total spikes
    norm_crossings = round(total_crossings / max(total_crossings) * length(color_map));
    norm_crossings(isnan(norm_crossings)) = 1;
    norm_crossings(norm_crossings == 0) = 1;

end

% Loop through each unit
for uu = 1:length(unit_names)

    % Find the corresponding row, column, & unit
    current_row = unit_rows(uu);
    current_column = unit_columns(uu);
    current_unit = strcat(unit_names{uu}, '_');

    p = (current_row*map_width) - (map_width - current_column);

    % Define the area within the figure
    subplot(map_height, map_width, p)
    hold on
    set(gca,'Color','k')

    % Find the associated sorted units
    sorted_unit_idxs = find(contains(xds.unit_names, current_unit));
    if isempty(sorted_unit_idxs)

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

        continue
    end

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

        if strcmp(Plot_Specs, 'Wave')
            % Extracting the waveforms of the designated unit
            unit_waveforms = xds.spike_waveforms{sorted_unit_idxs(ii)};
            % Calculating the mean
            amp_mean = mean(unit_waveforms, 1);
            
            % Plot the waveform
            plot(amp_mean, plot_color, 'LineWidth',2)
        end

        if strcmp(Plot_Specs, 'ISI')
            % Extracting the spikes of the designated unit
            spikes = xds.spikes{sorted_unit_idxs(ii)};
        
            % Calculate the interstimulus intervals
            ISI = diff(spikes);
        
            % Convert the time to milliseconds
            ISI = ISI*1000;

            % Define the edges of the histogram (0 - 100 ms)
            ISI_edges = (0:1:100);

            % Plot the histogram
            histogram(ISI, ISI_edges, 'EdgeColor', plot_color, 'FaceColor', plot_color)
        end

        if strcmp(Plot_Specs, 'Nonlin')

            % Extracting the spike times of the designated unit
            spike_times = xds.spikes{sorted_unit_idxs(ii)};

            % Extracting the nonlinear waveforms of the all units
            nonlin_waveforms = xds.nonlin_waveforms{sorted_unit_idxs(ii)};

            perspike_nonlin = mean(nonlin_waveforms, 2);

            % Convert time into minutes
            spike_min = spike_times / 60;

            % Plot the scatter plot
            scatter(spike_min, perspike_nonlin, 1, 'MarkerEdgeColor', plot_color,...
                'MarkerFaceColor', plot_color);
        end

        if strcmp(Plot_Specs, 'Colormap')
            
            set(gca,'Color', color_map(norm_crossings(p), :))

            % Annotation of the peak firing rate
            if ~isnan(peak_firing_rate(p))
                firing_rate_string = mat2str(round(peak_firing_rate(p), 1));
                firing_rate_string = {char(firing_rate_string)};
                ann_fire_rate = text(0.5*max(xlim), 0.5*max(ylim), firing_rate_string, ...
                    'horizontalalignment', 'center', 'verticalalignment', 'middle');
                set(ann_fire_rate, 'Color', [1, 1, 1])
            end

        end

    end

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

    %% Save the file if selected
    Save_Figs(Fig_Title, Save_File)

end


