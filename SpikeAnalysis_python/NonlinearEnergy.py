# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
from Plot_Specs import Font_Specs
import scipy.stats as stats
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

class Nonlinear_Energy():
    def __init__(self, xds_morn, xds_noon, unit_name):
        
        self.xds_morn = xds_morn
        self.xds_noon = xds_noon
        
        #%% Find the meta info to load the output excel table
        if isinstance(unit_name, int):
            # Date
            file_name = self.xds_morn.file_name
            split_info = file_name.split("_", 1)
            trial_date = split_info[0]
            
            # Task
            if self.xds_morn._lab_data__meta['task_name'] == 'multi_gadget':
                trial_task = 'PG'
            elif self.xds_morn._lab_data__meta['task_name'] == 'WS':
                trial_task = 'WS'
                
            # Monkey
            monkey_name = self.xds_morn._lab_data__meta['monkey_name']
            
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
                self.N_morn = self.xds_morn.unit_names.index(unit)
                self.N_noon = self.xds_noon.unit_names.index(unit)
            except KeyError:
                ## If the unit doesn't exist
                print(unit + ' does not exist')
                self.euc_dis = float('NaN')
                self.nonlin_p_val = float('NaN')
                return
            
        elif isinstance(unit_name, str):
            try:
                self.N_morn = self.xds_morn.unit_names.index(unit_name)
                self.N_noon = self.xds_noon.unit_names.index(unit_name)
            except ValueError:
                ## If the unit doesn't exist
                print(unit_name + ' does not exist')
                self.euc_dis = float('NaN')
                self.nonlin_p_val = float('NaN')
                return
            
        #%% Some variable extraction & definitions
        
        # Name of the designated unit
        self.unit_name = self.xds_morn.unit_names[self.N_morn]
        
        # Extracting the spike waveforms of the designated unit
        self.waveforms_morn = self.xds_morn.spike_waveforms[self.N_morn]
        self.waveforms_noon = self.xds_noon.spike_waveforms[self.N_noon]
        
        if hasattr(self.xds_morn, 'nonlin_waveforms'):
            # Extracting the nonlinear waveforms of all units
            self.nonlin_waveforms_morn = self.xds_morn.nonlin_waveforms
            self.nonlin_waveforms_noon = self.xds_noon.nonlin_waveforms
        else:
            self.nonlin_waveforms_morn = np.zeros([len(self.waveforms_morn), len(self.waveforms_morn[0]) - 2])
            self.nonlin_waveforms_noon = np.zeros([len(self.waveforms_noon), len(self.waveforms_noon[0]) - 2])
            # Morning
            for jj in range(len(self.waveforms_morn)):
                for ii in range(len(self.waveforms_morn[0]) - 2):
                    self.nonlin_waveforms_morn[jj,ii] = self.waveforms_morn[jj,ii+1]**2 - self.waveforms_morn[jj,ii]*self.waveforms_morn[jj,ii+2]
            # Afternoon
            for jj in range(len(self.waveforms_noon)):
                for ii in range(len(self.waveforms_noon[0]) - 2):
                    self.nonlin_waveforms_noon[jj,ii] = self.waveforms_noon[jj,ii+1]**2 - self.waveforms_noon[jj,ii]*self.waveforms_noon[jj,ii+2]
            
        # Extracting the spike times of the designated unit
        self.spike_times_morn = self.xds_morn.spikes[self.N_morn]
        self.spike_times_noon = self.xds_noon.spikes[self.N_noon]
        
        # Font & plotting specifications
        self.font_specs = Font_Specs()
        
        #%% Calculating the means & standard deviations
        
        # Nonlinear waveform means
        if hasattr(self.xds_morn, 'nonlin_waveforms'):
            self.morn_nonlin_mean = np.asarray(self.nonlin_waveforms_morn[self.N_morn]).mean(0)
            self.noon_nonlin_mean = np.asarray(self.nonlin_waveforms_noon[self.N_noon]).mean(0)
            self.morn_perspike_nonlin = np.asarray(self.nonlin_waveforms_morn[self.N_morn]).mean(1)
            self.noon_perspike_nonlin = np.asarray(self.nonlin_waveforms_noon[self.N_noon]).mean(1)
        else:
            self.morn_nonlin_mean = np.asarray(self.nonlin_waveforms_morn).mean(0)
            self.noon_nonlin_mean = np.asarray(self.nonlin_waveforms_noon).mean(0)
            self.morn_perspike_nonlin = np.asarray(self.nonlin_waveforms_morn).mean(1)
            self.noon_perspike_nonlin = np.asarray(self.nonlin_waveforms_noon).mean(1)
        
        # Nonlinear waveform standard deviations
        if hasattr(self.xds_morn, 'nonlin_waveforms'):
            self.standard_nonlin_morn = np.asarray(self.nonlin_waveforms_morn[self.N_morn]).std(0)
            self.standard_nonlin_noon = np.asarray(self.nonlin_waveforms_noon[self.N_noon]).std(0)
        else:
            self.standard_nonlin_morn = np.asarray(self.nonlin_waveforms_morn).std(0)
            self.standard_nonlin_noon = np.asarray(self.nonlin_waveforms_noon).std(0)
            
        # Calculate the Euclidean distance
        self.euc_dis = np.linalg.norm(self.morn_nonlin_mean - self.noon_nonlin_mean)
        
        #%% Find the change in mean nonlinear energy between the two sessions
        # Finding the mean & difference between morning & afternoon
        self.nonlin_morn = np.mean(self.morn_perspike_nonlin)
        self.nonlin_noon = np.mean(self.noon_perspike_nonlin)
        self.std_nonlin_morn = np.std(self.morn_perspike_nonlin)
        self.std_nonlin_noon = np.std(self.noon_perspike_nonlin)
        
        # Finding the difference
        nonlin_difference = np.abs(self.nonlin_morn - self.nonlin_noon)
        
        #%% Find the average nonlinear energy across all units
        if hasattr(self.xds_morn, 'nonlin_waveforms'):
            # The means of each nonlinear waveform
            all_units_nonlin_morn = []
            all_units_nonlin_noon = []
            for ii in range(len(self.xds_morn.unit_names)):
                all_units_nonlin_morn.append(np.asmatrix(self.nonlin_waveforms_morn[ii].mean(0)))
            for ii in range(len(self.xds_noon.unit_names)):
                all_units_nonlin_noon.append(np.asmatrix(self.nonlin_waveforms_noon[ii].mean(0)))
        
            # The means of every point in the mean of the nonlinear waveforms
            all_units_nonlin_means_morn = np.zeros(len(self.xds_morn.unit_names))
            all_units_nonlin_means_noon = np.zeros(len(self.xds_noon.unit_names))
            for ii in range(len(self.xds_morn.unit_names)):
                all_units_nonlin_means_morn[ii] = np.mean(all_units_nonlin_morn[ii])
            for ii in range(len(self.xds_noon.unit_names)):
                all_units_nonlin_means_noon[ii] = np.mean(all_units_nonlin_noon[ii])
                
            # Find the differences between morning & afternoon
            all_units_nonlin_difference = np.abs(all_units_nonlin_means_morn - all_units_nonlin_means_noon)
            
        #%% Determine if the unit is well sorted
        # Pick a random 50 units 1000 times
        p_val_idx = []
        for ii in range(0, 1000):
            # Morning
            rand_morn_idx = np.random.randint(len(self.morn_perspike_nonlin), size = 50)
            rand_morn_perspike = self.morn_perspike_nonlin[rand_morn_idx]
            # Afternoon
            rand_noon_idx = np.random.randint(len(self.noon_perspike_nonlin), size = 50)
            rand_noon_perspike = self.noon_perspike_nonlin[rand_noon_idx]
            nonlin_p_val = stats.ttest_ind(rand_morn_perspike, rand_noon_perspike)[1]
            p_val_idx.append(nonlin_p_val)
        self.nonlin_p_val = np.mean(p_val_idx)
        
        #%% Determine if the change in nonlinear energy is more or less than the average change
        
        if hasattr(self.xds_morn, 'nonlin_waveforms'):
            mean_nonlin_difference = np.mean(all_units_nonlin_difference)
            
            if nonlin_difference < mean_nonlin_difference:
                print('The change in nonlinear energy of ' + self.unit_name + ' is less than the mean of all units')
                
            if nonlin_difference > mean_nonlin_difference:
                print('The change in nonlinear energy of ' + self.unit_name + ' is greater than the mean of all units')
            
            
    #%% Plotting the random 100 waveforms
    
    def Nonlin_Rand(self, Save_Figs):
        
        # End the function if the unit does not exist
        if np.isnan(self.euc_dis) or np.isnan(self.nonlin_p_val):
            print('Unit does not exist')
            return
    
        # Morning
        rand_nonlin_idx_morn = np.random.randint(len(self.waveforms_morn[self.N_morn]), size = 100)
        if hasattr(self.xds_morn, 'nonlin_waveforms'):
            rand_nonlin_morn = self.nonlin_waveforms_morn[self.N_morn][rand_nonlin_idx_morn].T
        else:
            rand_nonlin_morn = self.nonlin_waveforms_morn[rand_nonlin_idx_morn].T
        
        # Afternoon
        rand_nonlin_idx_noon = np.random.randint(len(self.waveforms_noon[self.N_noon]), size = 100)
        if hasattr(self.xds_morn, 'nonlin_waveforms'):
            rand_nonlin_noon = self.nonlin_waveforms_noon[self.N_noon][rand_nonlin_idx_noon].T
        else:
            rand_nonlin_noon = self.nonlin_waveforms_noon[rand_nonlin_idx_noon].T
        
        # Morning
        morn_fig, figure_axes = plt.subplots()
        plt.plot(rand_nonlin_morn, 'k')
        morn_title_string = 'Nonlinear Energy - ' + self.unit_name + ' (Morning)'
        plt.title(morn_title_string, fontname = self.font_specs.font_name, fontsize = self.font_specs.title_font_size, fontweight = 'bold')
        # Axis labels
        plt.xlabel('Time', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        plt.ylabel('µV^2', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        
        # Annotation of the n-count
        plt.text(0.85, 0.9, 'n = ' + str(len(rand_nonlin_idx_morn)), verticalalignment = 'center', horizontalalignment = 'center', 
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
        plt.plot(rand_nonlin_noon, 'k')
        noon_title_string = 'Nonlinear Energy - ' + self.unit_name + ' (Afternoon)'
        plt.title(noon_title_string, fontname = self.font_specs.font_name, fontsize = self.font_specs.title_font_size, fontweight = 'bold')
        # Axis labels
        plt.xlabel('Time', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        plt.ylabel('µV^2', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
                
        # Annotation of the n-count
        plt.text(0.85, 0.9, 'n = ' + str(len(rand_nonlin_idx_noon)), verticalalignment = 'center', horizontalalignment = 'center', 
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
    def Nonlin_Mean(self, Save_Figs):
        
        # End the function if the unit does not exist
        if np.isnan(self.euc_dis) or np.isnan(self.nonlin_p_val):
            print('Unit does not exist')
            return
    
        # Morning
        morn_fig, figure_axes = plt.subplots()
        plt.plot(self.morn_nonlin_mean.T, 'k')
        plt.plot(self.morn_nonlin_mean.T + self.standard_nonlin_morn, 'r')
        plt.plot(self.morn_nonlin_mean.T - self.standard_nonlin_morn, 'r')
        morn_title_string = 'Average Nonlinear Energy - ' + self.unit_name + ' (Morning)'
        plt.title(morn_title_string, fontname = self.font_specs.font_name, fontsize = self.font_specs.title_font_size, fontweight = 'bold')
        # Axis Labels
        plt.xlabel('Time', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        plt.ylabel('µV^2', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        
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
        plt.plot(self.noon_nonlin_mean.T, 'k')
        plt.plot(self.noon_nonlin_mean.T + self.standard_nonlin_noon, 'r')
        plt.plot(self.noon_nonlin_mean.T - self.standard_nonlin_noon, 'r')
        noon_title_string = 'Average Nonlinear Energy - ' + self.unit_name + ' (Afternoon)'
        plt.title(noon_title_string, fontname = self.font_specs.font_name, fontsize = self.font_specs.title_font_size, fontweight = 'bold')
        # Axis Labels
        plt.xlabel('Time', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        plt.ylabel('µV^2', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        
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
    def Nonlin_Hist(self, Save_Figs):
    
        # End the function if the unit does not exist
        if np.isnan(self.euc_dis) or np.isnan(self.nonlin_p_val):
            print('Unit does not exist')
            return
    
        hist_fig, figure_axes = plt.subplots()
        plt.hist(self.morn_perspike_nonlin, bins = 'fd', alpha = self.font_specs.hist_transparency, 
                 edgecolor = 'k', color = (0.9290, 0.6940, 0.1250), label = 'Morning')
        plt.hist(self.noon_perspike_nonlin, bins = 'fd', alpha = self.font_specs.hist_transparency, 
                 edgecolor = 'k', color = (0.5, 0, 0.5), label = 'Afternoon')
        title_string = 'Nonlinear Energy Histogram - ' + self.unit_name
        plt.title(title_string, fontname = self.font_specs.font_name, fontsize = self.font_specs.title_font_size, fontweight = 'bold')
        plt.xlabel('Mean Nonlinear Energy', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
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
        if round(self.nonlin_p_val, 3) > 0:
            plt.text(0.85, 0.825, 'p = ' + str(round(self.nonlin_p_val, 3)), 
                     verticalalignment = 'center', horizontalalignment = 'center', 
                     transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
        if round(self.nonlin_p_val, 3) == 0:
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
    def Nonlin_Scatter(self, Save_Figs):
        
        # End the function if the unit does not exist
        if np.isnan(self.euc_dis) or np.isnan(self.nonlin_p_val):
            print('Unit does not exist')
            return

        # Convert the time into minutes
        spike_min_morn = self.spike_times_morn / 60
        spike_min_noon = self.spike_times_noon / 60
        
        # Find the axis limits
        morn_ymax = self.nonlin_morn + 4*self.std_nonlin_morn
        morn_ymin = self.nonlin_morn - 4*self.std_nonlin_morn
        noon_ymax = self.nonlin_noon + 4*self.std_nonlin_noon
        noon_ymin = self.nonlin_noon - 4*self.std_nonlin_noon
        morn_xmax = np.max(spike_min_morn)
        noon_xmax = np.max(spike_min_noon)
        axis_limits = np.array([morn_ymax, morn_ymin, noon_ymax, noon_ymin])
        max_YLim = np.max(axis_limits)
        min_YLim = np.min(axis_limits)
        
        # Size of the markers
        sz = 2
        
        nonlin_timeline = plt.figure(figsize = (11, 4))
                
        # Morning
        morn_nonlin_plot = plt.subplot(2,2,1)
        morn_nonlin_plot.scatter(spike_min_morn, self.morn_perspike_nonlin, s = sz, edgecolors = 'k', )
        # Line indicating the morning mean
        morn_nonlin_plot.plot([0, np.max(spike_min_morn)], [self.nonlin_morn, self.nonlin_morn], color = 'r', linewidth = 1)
        # Lines indicating the morning standard deviation
        morn_nonlin_plot.plot([0, np.max(spike_min_morn)], [self.nonlin_morn + self.std_nonlin_morn, self.nonlin_morn + self.std_nonlin_morn], 
                 color = 'r', linewidth = 1, linestyle = 'dashed', dashes = (5,5))
        morn_nonlin_plot.plot([0, np.max(spike_min_morn)], [self.nonlin_morn - self.std_nonlin_morn, self.nonlin_morn - self.std_nonlin_morn], 
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
        noon_nonlin_plot = plt.subplot(2,2,2)
        noon_nonlin_plot.scatter(spike_min_noon, self.noon_perspike_nonlin, s = sz, edgecolors = 'k', )
        # Line indicating the afternoon mean
        noon_nonlin_plot.plot([0, np.max(spike_min_noon)], [self.nonlin_noon, self.nonlin_noon], color = 'r', linewidth = 1)
        # Lines indicating the afternoon standard deviation
        noon_nonlin_plot.plot([0, np.max(spike_min_noon)], [self.nonlin_noon + self.std_nonlin_noon, self.nonlin_noon + self.std_nonlin_noon], 
                 color = 'r', linewidth = 1, linestyle = 'dashed', dashes = (5,5))
        noon_nonlin_plot.plot([0, np.max(spike_min_noon)], [self.nonlin_noon - self.std_nonlin_noon, self.nonlin_noon - self.std_nonlin_noon], 
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
        if round(self.nonlin_p_val, 3) > 0:
            plt.text(0.6, 1.1, 'p = ' + str(round(self.nonlin_p_val, 3)), 
                     verticalalignment = 'center', horizontalalignment = 'center', 
                     transform = noon_figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size + 2)
        if round(self.nonlin_p_val, 3) == 0:
            plt.text(1, 1, 'p < 0.001', verticalalignment = 'center', horizontalalignment = 'center', 
                     transform = noon_figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
            
        # Squeeze the subplots closer
        nonlin_timeline.tight_layout()
        # Set the common title
        title_string = 'Nonlinear Energy - ' + self.unit_name
        plt.suptitle(title_string, fontname = self.font_specs.font_name, 
                     fontsize = self.font_specs.title_font_size + 5, fontweight = 'bold')
        # Set the common x-label
        nonlin_timeline.text(0.47, 0.325, 'Time (min.)', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size + 2)

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
            
            
            
            
            
            
            
            
            
                
                
                
                
                
                
                
                
                
            