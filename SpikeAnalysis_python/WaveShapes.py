#%% Import basic packages
import os
import numpy as np
import pandas as pd
import scipy.stats as stats
from Plot_Specs import Font_Specs
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

class Wave_Shapes():
    
    def __init__(self, xds_morn, xds_noon, unit_name):
    
        #%% Find the meta info to load the output excel table
        if isinstance(unit_name, int):
            # Date
            file_name = xds_morn.file_name
            split_info = file_name.split("_", 1)
            trial_date = split_info[0]
            
            # Task
            if xds_morn._lab_data__meta['task_name'] == 'multi_gadget':
                trial_task = 'PG'
            elif xds_morn._lab_data__meta['task_name'] == 'WS':
                trial_task = 'WS'
                
            # Monkey
            monkey_name = xds_morn._lab_data__meta['monkey_name']
            
            # File path
            base_excel_dir = 'C:/Users/rhpow/Documents/Work/Northwestern/Excel Data/'
            # List of all files
            dir_list = os.listdir(base_excel_dir)
            
            # Selected file index
            file_substring = trial_date + '_' + monkey_name + '_' + trial_task
            select_file = list(filter(lambda x: file_substring in x, dir_list))[0]
            
            output_xds = pd.read_excel(base_excel_dir + select_file)
            
            #%% Find the unit of interest
            try:
                unit = output_xds.unit_names[unit_name]
            
                ## Identify the index of the unit
                self.N_morn = xds_morn.unit_names.index(unit)
                self.N_noon = xds_noon.unit_names.index(unit)
            except KeyError:
                ## If the unit doesn't exist
                print(unit + ' does not exist')
                self.euc_dis = float('NaN')
                self.wave_p_val = float('NaN')
                return
            
        elif isinstance(unit_name, str):
            try:
                self.N_morn = xds_morn.unit_names.index(unit_name)
                self.N_noon = xds_noon.unit_names.index(unit_name)
            except ValueError:
                ## If the unit doesn't exist
                print(unit_name + ' does not exist')
                self.euc_dis = float('NaN')
                self.wave_p_val = float('NaN')
                return
            
        #%% Some variable extraction & definitions
        
        # Name of the designated unit
        self.unit_name = xds_morn.unit_names[self.N_morn]
        
        # Extracting the waveforms of the designated unit
        self.waveforms_morn = xds_morn.spike_waveforms
        self.waveforms_noon = xds_noon.spike_waveforms
        
        # Extracting the spike times of the designated unit
        self.spike_times_morn = xds_morn.spikes[self.N_morn]
        self.spike_times_noon = xds_noon.spikes[self.N_noon]
        
        # Font & plotting specifications
        self.font_specs = Font_Specs()

        #%% Calculating the means & standard deviations
        
        # Calculating the means
        self.morn_amp_mean = np.asarray(self.waveforms_morn[self.N_morn]).mean(0)
        self.noon_amp_mean = np.asarray(self.waveforms_noon[self.N_noon]).mean(0)
        self.morn_perspike_amp = np.asarray(self.waveforms_morn[self.N_morn]).mean(1)
        self.noon_perspike_amp = np.asarray(self.waveforms_noon[self.N_noon]).mean(1)
        
        # Calculating the standard deviations
        self.standard_morn = np.asarray(self.waveforms_morn[self.N_morn]).std(0)
        self.standard_noon = np.asarray(self.waveforms_noon[self.N_noon]).std(0)
        
        # Calculate the euclidean distance
        self.euc_dis = np.linalg.norm(self.morn_amp_mean - self.noon_amp_mean)
        
        #%% Find the change in mean amplitude between the two sessions
        # Finding the mean & difference between morning & afternoon
        self.amp_morn = np.mean(self.morn_perspike_amp)
        self.amp_noon = np.mean(self.noon_perspike_amp)
        self.std_amp_morn = np.std(self.morn_perspike_amp)
        self.std_amp_noon = np.std(self.noon_perspike_amp)
        
        # Finding the difference
        amp_difference = np.abs(self.amp_morn - self.amp_noon)
        
        #%% Calculate the peak to peak amplitude
        # Morning
        crest_morn = np.max(self.morn_amp_mean)
        trough_morn = np.min(self.morn_amp_mean)
        # Afternoon
        crest_noon = np.max(self.noon_amp_mean)
        trough_noon = np.min(self.noon_amp_mean)
        
        # Peak to peak amplitudes
        peaktopeak_morn = np.abs(crest_morn) + np.abs(trough_morn)
        peaktopeak_noon = np.abs(crest_noon) + np.abs(trough_noon)
        
        #%% Finding the average amplitudes across all units
        # The means of each waveform
        all_units_amp_morn = []
        all_units_amp_noon = []
        for ii in range(len(xds_morn.unit_names)):
            all_units_amp_morn.append(np.asmatrix(self.waveforms_morn[ii]).mean(0))
        for ii in range(len(xds_noon.unit_names)):
            all_units_amp_noon.append(np.asmatrix(self.waveforms_noon[ii]).mean(0))
        
        # The mean of every point in the mean of the waveforms
        all_units_amp_means_morn = []
        all_units_amp_means_noon = []
        for ii in range(len(xds_morn.unit_names)):
            all_units_amp_means_morn.append(np.mean(self.waveforms_morn[ii]))
        for ii in range(len(xds_noon.unit_names)):
            all_units_amp_means_noon.append(np.mean(self.waveforms_noon[ii]))
        
        # Find the differences between morning & afternoon
        all_units_amp_difference = np.abs(np.array(all_units_amp_means_morn) - np.array(all_units_amp_means_noon))
        
        #%% Determine if the unit is well sorted
        # Pick a random 50 units 1000 times
        p_val_idx = []
        for ii in range(0, 1000):
            # Morning
            rand_morn_idx = np.random.randint(len(self.morn_perspike_amp), size = 50)
            rand_morn_perspike = self.morn_perspike_amp[rand_morn_idx]
            # Afternoon
            rand_noon_idx = np.random.randint(len(self.noon_perspike_amp), size = 50)
            rand_noon_perspike = self.noon_perspike_amp[rand_noon_idx]
            wave_p_val = stats.ttest_ind(rand_morn_perspike, rand_noon_perspike)[1]
            p_val_idx.append(wave_p_val)
        self.wave_p_val = np.mean(p_val_idx)
        
        #%% Determine if the change in mean amplitude is more or less than the average change
        mean_amp_difference = np.mean(all_units_amp_difference)
        
        if amp_difference < mean_amp_difference:
            print('The change in waveform amplitude of ' + self.unit_name + ' is less than the mean of all units')
            
        if amp_difference > mean_amp_difference:
            print('The change in waveform amplitude of ' + self.unit_name + ' is greater than the mean of all units')
            
        #%% Print the peak to peak amplitude
        # Morning
        print('The average peak to peak amplitude of ' + self.unit_name + ' is ' + str(round(peaktopeak_morn, 3)) + ' µV in the morning')
        # Afternoon
        print('The average peak to peak amplitude of ' + self.unit_name + ' is ' + str(round(peaktopeak_noon, 3)) + ' µV in the afternoon')
             
    #%% Plotting the random 100 waveforms
    
    def WaveShapes_Rand(self, Save_Figs):
        
        # End the function if the unit does not exist
        if np.isnan(self.euc_dis) or np.isnan(self.wave_p_val):
            print('Unit does not exist')
            return
    
        # Morning
        rand_wave_idx_morn = np.random.randint(len(self.waveforms_morn[self.N_morn]), size = 100)
        rand_waves_morn = self.waveforms_morn[self.N_morn][rand_wave_idx_morn].T
        
        # Afternoon
        rand_wave_idx_noon = np.random.randint(len(self.waveforms_noon[self.N_noon]), size = 100)
        rand_waves_noon = self.waveforms_noon[self.N_noon][rand_wave_idx_noon].T
        
        # Morning
        morn_fig, figure_axes = plt.subplots()
        plt.plot(rand_waves_morn, 'k')
        morn_title_string = 'Individual Waveforms - ' + self.unit_name + ' (Morning)'
        plt.title(morn_title_string, fontname = self.font_specs.font_name, fontsize = self.font_specs.title_font_size, fontweight = 'bold')
        # Axis Labels
        plt.xlabel('Time', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        plt.ylabel('Amplitude (µV)', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        
        # Annotation of the n-count
        plt.text(0.85, 0.9, 'n = ' + str(len(rand_wave_idx_morn)), verticalalignment = 'center', horizontalalignment = 'center', 
                 transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
        
        # Only label every other tick
        for label in figure_axes.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        for label in figure_axes.yaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        # Set The Font
        plt.xticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5)
        plt.yticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5)
    
        # Afternoon
        noon_fig, figure_axes = plt.subplots()
        plt.plot(rand_waves_noon, 'k')
        noon_title_string = 'Individual Waveforms - ' + self.unit_name + ' (Afternoon)'
        plt.title(noon_title_string, fontname = self.font_specs.font_name, fontsize = self.font_specs.title_font_size, fontweight = 'bold')
        # Axis Labels
        plt.xlabel('Time', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        plt.ylabel('Amplitude (µV)', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        
        # Annotation of the n-count
        plt.text(0.85, 0.9, 'n = ' + str(len(rand_wave_idx_noon)), verticalalignment = 'center', horizontalalignment = 'center', 
                 transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
        
        # Only label every other tick
        for label in figure_axes.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        for label in figure_axes.yaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        # Set The Font
        plt.xticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5)
        plt.yticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5) 
        
        # Figure Saving
        if Save_Figs != 0:
            save_dir = 'C:/Users/rhpow/Desktop/'
            for ii in range(2):
                if ii == 0:
                    fig_title = noon_title_string
                if ii == 1:
                    fig_title = morn_title_string
                plt.savefig(save_dir + fig_title + '.' + Save_Figs)
                plt.close()

    #%% Plotting the mean & standard deviations
    def WaveShapes_Mean(self, Save_Figs):
        
        # End the function if the unit does not exist
        if np.isnan(self.euc_dis) or np.isnan(self.wave_p_val):
            print('Unit does not exist')
            return
    
        # Morning
        morn_fig, figure_axes = plt.subplots()
        plt.plot(self.morn_amp_mean.T, 'k')
        plt.plot(self.morn_amp_mean.T + self.standard_morn, 'r')
        plt.plot(self.morn_amp_mean.T - self.standard_morn, 'r')
        morn_title_string = 'Average Waveform - ' + self.unit_name + ' (Morning)'
        plt.title(morn_title_string, fontname = self.font_specs.font_name, fontsize = self.font_specs.title_font_size, fontweight = 'bold')
        # Axis Labels
        plt.xlabel('Time', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        plt.ylabel('Amplitude (µV)', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        
        # Only label every other tick
        for label in figure_axes.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        for label in figure_axes.yaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        # Set The Font
        plt.xticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5)
        plt.yticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5)
        
        # Afternoon
        noon_fig, figure_axes = plt.subplots()
        plt.plot(self.noon_amp_mean.T, 'k')
        plt.plot(self.noon_amp_mean.T + self.standard_noon, 'r')
        plt.plot(self.noon_amp_mean.T - self.standard_noon, 'r')
        noon_title_string = 'Average Waveform - ' + self.unit_name + ' (Afternoon)'
        plt.title(noon_title_string, fontname = self.font_specs.font_name, fontsize = self.font_specs.title_font_size, fontweight = 'bold')
        # Axis Labels
        plt.xlabel('Time', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        plt.ylabel('Amplitude (µV)', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        
        # Only label every other tick
        for label in figure_axes.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        for label in figure_axes.yaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        # Set The Font
        plt.xticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5)
        plt.yticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5)
        
        # Figure Saving
        if Save_Figs != 0:
            save_dir = 'C:/Users/rhpow/Desktop/'
            for ii in range(2):
                if ii == 0:
                    fig_title = noon_title_string
                if ii == 1:
                    fig_title = morn_title_string
                plt.savefig(save_dir + fig_title + '.' + Save_Figs)
                plt.close()

    #%% Plotting the morning & afternoon histograms
    # (Number of bins determined by Freedman-Diaconis rule)
    def WaveShapes_Hist(self, Save_Figs):
        
        # End the function if the unit does not exist
        if np.isnan(self.euc_dis) or np.isnan(self.wave_p_val):
            print('Unit does not exist')
            return
    
        his_fig, figure_axes = plt.subplots()
        plt.hist(self.morn_perspike_amp, bins = 'fd', alpha = self.font_specs.hist_transparency, 
                 edgecolor = 'k', color = (0.9290, 0.6940, 0.1250), label = 'Morning')
        plt.hist(self.noon_perspike_amp, bins = 'fd', alpha = self.font_specs.hist_transparency, 
                 edgecolor = 'k', color = (0.5, 0, 0.5), label = 'Afternoon')
        title_string = 'Average Amplitude Histogram - ' + self.unit_name
        plt.title(title_string, fontname = self.font_specs.font_name, fontsize = self.font_specs.title_font_size, fontweight = 'bold')
        plt.xlabel('Mean Waveform Amplitude', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        plt.ylabel('Counts', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        
        # Collect the current axis limits
        y_limits = figure_axes.get_ylim()
        x_limits = figure_axes.get_xlim()
        
        # Reset the axis limits
        plt.xlim(x_limits)
        plt.ylim(y_limits[0], y_limits[1] + 1)
        
        legend_font = fm.FontProperties(family = self.font_specs.font_name, size = self.font_specs.legend_font_size)
        plt.legend(prop = legend_font)
        plt.legend(frameon = False)
        
        # Annotation of the p-value
        if round(self.wave_p_val, 3) > 0:
            plt.text(0.85, 0.825, 'p = ' + str(round(self.wave_p_val, 3)), 
                     verticalalignment = 'center', horizontalalignment = 'center', 
                     transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
        if round(self.wave_p_val, 3) == 0:
            plt.text(0.85, 0.825, 'p < 0.001', verticalalignment = 'center', horizontalalignment = 'center', 
                     transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
            
        # Only label every other tick
        for label in figure_axes.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        for label in figure_axes.yaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        # Set The Font
        plt.xticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5)
        plt.yticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5)
        
        # Figure Saving
        if Save_Figs != 0:
            save_dir = 'C:/Users/rhpow/Desktop/'
            fig_title = title_string
            fig_title = str.replace(fig_title, ':', '')
            fig_title = str.replace(fig_title, 'vs.', 'vs')
            fig_title = str.replace(fig_title, 'mg.', 'mg')
            fig_title = str.replace(fig_title, 'kg.', 'kg')
            fig_title = str.replace(fig_title, '.', '_')
            fig_title = str.replace(fig_title, '/', '_')
            plt.savefig(save_dir + fig_title + '.' + Save_Figs)
            plt.close()
    
    #%% Plotting the per-spike waveform amplitude
    def WaveShapes_Scatter(self, Save_Figs):
        
        # End the function if the unit does not exist
        if np.isnan(self.euc_dis) or np.isnan(self.wave_p_val):
            print('Unit does not exist')
            return
        
        # Convert the time into minutes
        spike_min_morn = self.spike_times_morn / 60
        spike_min_noon = self.spike_times_noon / 60
        
        # Find the axis limits
        morn_ymax = self.amp_morn + 4*self.std_amp_morn
        morn_ymin = self.amp_morn - 4*self.std_amp_morn
        noon_ymax = self.amp_noon + 4*self.std_amp_noon
        noon_ymin = self.amp_noon - 4*self.std_amp_noon
        morn_xmax = np.max(spike_min_morn)
        noon_xmax = np.max(spike_min_noon)
        axis_limits = np.array([morn_ymax, morn_ymin, noon_ymax, noon_ymin])
        max_YLim = np.max(axis_limits)
        min_YLim = np.min(axis_limits)
        
        # Size of the markers
        sz = 2
        
        amp_timeline = plt.figure(figsize = (11, 4))
        
        # Morning
        morn_amp_plot = plt.subplot(2,2,1)
        morn_amp_plot.scatter(spike_min_morn, self.morn_perspike_amp, s = sz, edgecolors = 'k', )
        # Line indicating the morning mean
        morn_amp_plot.plot([0, np.max(spike_min_morn)], [self.amp_morn, self.amp_morn], color = 'r', linewidth = 1)
        # Lines indicating the morning standard deviation
        morn_amp_plot.plot([0, np.max(spike_min_morn)], [self.amp_morn + self.std_amp_morn, self.amp_morn + self.std_amp_morn], 
                 color = 'r', linewidth = 1, linestyle = 'dashed', dashes = (5,5))
        morn_amp_plot.plot([0, np.max(spike_min_morn)], [self.amp_morn - self.std_amp_morn, self.amp_morn - self.std_amp_morn], 
                 color = 'r', linewidth = 1, linestyle = 'dashed', dashes = (5,5))
        
        # Setting the morning axis limits
        plt.ylim([min_YLim, max_YLim])
        plt.xlim([0, morn_xmax + 0.1])
        
        # Only label every other tick (Morning)
        morn_figure_axes = plt.gca()
        for label in morn_figure_axes.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        # Set The Font
        plt.xticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5)
        plt.yticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5)
        
        plt.xlabel('Morning', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        
        # Afternoon
        noon_amp_plot = plt.subplot(2,2,2)
        noon_amp_plot.scatter(spike_min_noon, self.noon_perspike_amp, s = sz, edgecolors = 'k', )
        # Line indicating the afternoon mean
        noon_amp_plot.plot([0, np.max(spike_min_noon)], [self.amp_noon, self.amp_noon], color = 'r', linewidth = 1)
        # Lines indicating the afternoon standard deviation
        noon_amp_plot.plot([0, np.max(spike_min_noon)], [self.amp_noon + self.std_amp_noon, self.amp_noon + self.std_amp_noon], 
                 color = 'r', linewidth = 1, linestyle = 'dashed', dashes = (5,5))
        noon_amp_plot.plot([0, np.max(spike_min_noon)], [self.amp_noon - self.std_amp_noon, self.amp_noon - self.std_amp_noon], 
                 color = 'r', linewidth = 1, linestyle = 'dashed', dashes = (5,5))
        
        # Setting the afternoon axis limits
        plt.ylim([min_YLim, max_YLim])
        plt.xlim([0, noon_xmax + 0.1])
        
        # Only label every other tick (Afternoon)
        noon_figure_axes = plt.gca()
        for label in noon_figure_axes.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        for label in noon_figure_axes.yaxis.get_ticklabels()[::1]:
            label.set_visible(False)
        # Set The Font
        plt.xticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5)
        plt.yticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5)
        
        plt.xlabel('Afternoon', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        
        # Annotation of the p-value
        if round(self.wave_p_val, 3) > 0:
            plt.text(0.6, 1.1, 'p = ' + str(round(self.wave_p_val, 3)), 
                     verticalalignment = 'center', horizontalalignment = 'center', 
                     transform = noon_figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size + 2)
        if round(self.wave_p_val, 3) == 0:
            plt.text(1, 1, 'p < 0.001', verticalalignment = 'center', horizontalalignment = 'center', 
                     transform = noon_figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
            
        # Squeeze the subplots closer
        amp_timeline.tight_layout()
        # Set the common title
        title_string = 'Mean Waveform Amplitude - ' + self.unit_name
        plt.suptitle(title_string, fontname = self.font_specs.font_name, 
                     fontsize = self.font_specs.title_font_size + 5, fontweight = 'bold')
        # Set the common x-label
        amp_timeline.text(0.47, 0.325, 'Time (min.)', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size + 2)

        # Figure Saving
        if Save_Figs != 0:
            save_dir = 'C:/Users/rhpow/Desktop/'
            fig_title = title_string
            fig_title = str.replace(fig_title, ':', '')
            fig_title = str.replace(fig_title, 'vs.', 'vs')
            fig_title = str.replace(fig_title, 'mg.', 'mg')
            fig_title = str.replace(fig_title, 'kg.', 'kg')
            fig_title = str.replace(fig_title, '.', '_')
            fig_title = str.replace(fig_title, '/', '_')
            plt.savefig(save_dir + fig_title + '.' + Save_Figs)
            plt.close()




















