# -*- coding: utf-8 -*-

import os
import numpy as np
import pandas as pd
from Plot_Specs import Font_Specs
import scipy.stats as stats
from sympy import symbols, Eq, solve
import matplotlib.pyplot as plt
import matplotlib.font_manager as fm

class Interstimulus_Interval():
    def __init__(self, xds_morn, xds_noon, unit_name):
        
        self.xds_morn = xds_morn
        self.xds_noon = xds_noon
        self.unit_name = unit_name
        
        #%% Find the meta info to load the output excel table
        if isinstance(self.unit_name, int):
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
                unit = output_xds.unit_names[self.unit_name]
            
                ## Identify the index of the unit
                self.N_morn = self.xds_morn.unit_names.index(unit)
                self.N_noon = self.xds_noon.unit_names.index(unit)
            except KeyError:
                ## If the unit doesn't exist
                print(unit + ' does not exist')
                self.merged_fract_contam = float('NaN')
                self.ISI_modes = float('NaN')
                return
            
        elif isinstance(self.unit_name, str):
            try:
                self.N_morn = self.xds_morn.unit_names.index(self.unit_name)
                self.N_noon = self.xds_noon.unit_names.index(self.unit_name)
            except ValueError:
                ## If the unit doesn't exist
                print(self.unit_name + ' does not exist')
                self.merged_fract_contam = float('NaN')
                self.ISI_modes = float('NaN')
                return
            
        #%% Some variable extraction & definitions
        
        # Defind the output variables
        self.ISI_modes = np.empty((1,3,))
        self.ISI_modes[:] = np.nan
        
        # Extracting the spikes of the designated unit
        spikes_morn = self.xds_morn.spikes[self.N_morn]
        spikes_noon = self.xds_noon.spikes[self.N_noon]
        
        # Calculate the interstimulus intervals in the morning
        self.ISI_morn = np.zeros(len(spikes_morn) - 1)
        for ii in range(len(self.ISI_morn)):
            self.ISI_morn[ii] = spikes_morn[ii+1] - spikes_morn[ii]
            
        # Calculate the interstimulus intervals in the morning
        self.ISI_noon = np.zeros(len(spikes_noon) - 1)
        for ii in range(len(self.ISI_noon)):
            self.ISI_noon[ii] = spikes_noon[ii+1] - spikes_noon[ii]
        
        # Merge the morning & afternoon interstimulus intervals 
        ISI_merged = np.append(self.ISI_morn, self.ISI_noon)
        
        # Convert the time to milliseconds
        self.ISI_morn = self.ISI_morn*1000
        self.ISI_noon = self.ISI_noon*1000
        ISI_merged = ISI_merged*1000
        
        # Define the edges of the histogram (0 - 100 ms)
        self.ISI_edges = np.linspace(0, 100, 101)
        
        self.ISI_modes[0,0] = stats.mode(self.ISI_morn)[0]
        self.ISI_modes[0,1] = stats.mode(self.ISI_noon)[0]
        self.ISI_modes[0,2] = stats.mode(ISI_merged)[0]
        
        # Font & plotting specifications
        self.font_specs = Font_Specs()
        
        #%% Calculate the number of refractory period violations
        
        # Define thte refractory period cutoff (2 ms)
        refract_min = 2
        # Find the number of refractory period violations
        morn_refract_violations = np.sum(self.ISI_morn < refract_min)
        noon_refract_violations = np.sum(self.ISI_noon < refract_min)
        merged_refract_violations = np.sum(ISI_merged < refract_min)
        
        # Find the percent of spikes this corresponds to
        violation_ratio_morn = morn_refract_violations / len(spikes_morn)
        violation_ratio_noon = noon_refract_violations / len(spikes_noon)
        merged_violation_ratio = merged_refract_violations / (len(spikes_morn) + len(spikes_noon))
        
        #%% Calculate the fractional level of contamination
        
        # Define the censored period (from Central's settings: 1 ms)
        censored_period = 1
        
        # Define the number of spikes
        morn_spikes = len(spikes_morn)
        noon_spikes = len(spikes_noon)
        total_spikes = morn_spikes + noon_spikes
        # Define the length of the experiment (in ms.)
        morn_experiment_length = self.xds_morn.time_frame[-1]*1000
        noon_experiment_length = self.xds_noon.time_frame[-1]*1000
        total_experiment_length = morn_experiment_length + noon_experiment_length
        
        Fract_Contam = symbols('Fract_Contam')
        
        morn_eqn = Eq(2*(refract_min - censored_period)*(morn_spikes**2)*(1 - Fract_Contam)*
                      (Fract_Contam / morn_experiment_length), morn_refract_violations)
        
        noon_eqn = Eq(2*(refract_min - censored_period)*(noon_spikes**2)*(1 - Fract_Contam)*
                      (Fract_Contam / noon_experiment_length), noon_refract_violations)
        
        merged_eqn = Eq(2*(refract_min - censored_period)*(total_spikes**2)*(1 - Fract_Contam)*
                      (Fract_Contam / total_experiment_length), merged_refract_violations)
        
        self.morn_fract_contam = solve(morn_eqn, Fract_Contam)
        
        self.noon_fract_contam = solve(noon_eqn, Fract_Contam)
        
        fract_contam = solve(merged_eqn, Fract_Contam)
        
        #%% Print the percentages & fractional contamination
        # Morning
        print(str(violation_ratio_morn*100) + '% of the spikes in ' + self.xds_morn.unit_names[self.N_morn] + 
              ' have an ISI less than ' + str(refract_min) + ' ms. in the morning')
        # Afternoon
        print(str(violation_ratio_noon*100) + '% of the spikes in ' + self.xds_noon.unit_names[self.N_noon] + 
              ' have an ISI less than ' + str(refract_min) + ' ms. in the afternoon')
        # Merged
        print(str(merged_violation_ratio*100) + '% of the spikes in ' + self.xds_noon.unit_names[self.N_noon] + 
              ' have an ISI less than ' + str(refract_min) + ' ms.')
        
        # Morning
        if any(np.iscomplex(self.morn_fract_contam)):
            print('The fractional level of contamination in ' + self.xds_morn.unit_names[self.N_morn] + ' > 0.5 in the morning')

        if all(np.isreal(self.morn_fract_contam)):
            print('The fractional level of contamination is ' + str(min(self.morn_fract_contam)) + ' in ' + 
                  self.xds_morn.unit_names[self.N_morn] + ' in the morning')
            
        # Afternoon
        if any(np.iscomplex(self.noon_fract_contam)):
            print('The fractional level of contamination in ' + self.xds_noon.unit_names[self.N_noon] + ' > 0.5 in the afternoon')
        
        if all(np.isreal(self.noon_fract_contam)):
            print('The fractional level of contamination is ' + str(min(self.noon_fract_contam)) + ' in ' + 
                  self.xds_noon.unit_names[self.N_noon] + ' in the afternoon')
            
        # Merged
        if any(np.iscomplex(fract_contam)):
            print('The fractional level of contamination in ' + self.xds_noon.unit_names[self.N_noon] + ' > 0.5')
            self.merged_fract_contam = 1
            
        if all(np.isreal(fract_contam)):
            print('The fractional level of contamination is ' + str(min(fract_contam)) + ' in ' + 
                  self.xds_noon.unit_names[self.N_noon])
            self.merged_fract_contam = np.float64(min(fract_contam))
        
    #%% Plotting the histograms
    def ISI_Hist(self, Save_Figs):
        
        # End the function if the unit does not exist
        if np.isnan(self.merged_fract_contam) or np.any(np.isnan(self.ISI_modes)):
            print('Unit does not exist')
            return
        
        # Merged histogram
        merged_hist_fig, figure_axes = plt.subplots()
        plt.hist(self.ISI_morn, bins = self.ISI_edges, alpha = self.font_specs.hist_transparency, 
                 edgecolor = 'k', color = (0.9290, 0.6940, 0.1250), label = 'Morning')
        plt.hist(self.ISI_noon, bins = self.ISI_edges, alpha = self.font_specs.hist_transparency, 
                 edgecolor = 'k', color = (0.5, 0, 0.5), label = 'Afternoon')
        merged_title_string = 'Interstimulus Intervals - ' + self.unit_name + ' (Merged)'
        plt.title(merged_title_string, fontname = self.font_specs.font_name, fontsize = self.font_specs.title_font_size, fontweight = 'bold')
        plt.xlabel('ISI (mSec)', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        plt.ylabel('Counts', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        
        # Collect the current axis limits
        y_limits = figure_axes.get_ylim()
        x_limits = figure_axes.get_xlim()
        
        # Annotation of the fractional contamination
        if self.merged_fract_contam < 1:
            plt.text(0.85, 0.825, u'f\u1D56\u2081' + ' = ' + str(round(self.merged_fract_contam, 3)), 
                     verticalalignment = 'center', horizontalalignment = 'center', 
                     transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
        if self.merged_fract_contam == 1:
            plt.text(0.85, 0.825, u'f\u1D56\u2081' + ' < 0.5', verticalalignment = 'center', horizontalalignment = 'center', 
                     transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size) 
        
        legend_font = fm.FontProperties(family = self.font_specs.font_name, size = self.font_specs.legend_font_size)
        plt.legend(prop = legend_font)
        plt.legend(frameon = False)
        
        # Reset the axis limits
        plt.xlim(x_limits)
        plt.ylim(y_limits[0], y_limits[1] + 1)
        
        # Only label every other tick
        for label in figure_axes.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        for label in figure_axes.yaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        # Set The Font
        plt.xticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5)
        plt.yticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5)
        
        # Morning histogram
        morn_hist_fig, figure_axes = plt.subplots()
        plt.hist(self.ISI_morn, bins = self.ISI_edges, alpha = self.font_specs.hist_transparency, 
                 edgecolor = 'k', color = (0.9290, 0.6940, 0.1250), label = 'Morning')
        morn_title_string = 'Interstimulus Intervals - ' + self.unit_name + ' (Morning)'
        plt.title(morn_title_string, fontname = self.font_specs.font_name, fontsize = self.font_specs.title_font_size, fontweight = 'bold')
        plt.xlabel('ISI (mSec)', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        plt.ylabel('Counts', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        
        # Annotation of the fractional contamination
        if all(np.isreal(self.morn_fract_contam)):
            plt.text(0.85, 0.825, u'f\u1D56\u2081' + ' = ' + str(round(min(self.morn_fract_contam), 3)), 
                     verticalalignment = 'center', horizontalalignment = 'center', 
                     transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
        if any(np.iscomplex(self.morn_fract_contam)):
            plt.text(0.85, 0.825, u'f\u1D56\u2081' + ' < 0.5', verticalalignment = 'center', horizontalalignment = 'center', 
                     transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
        
        # Reset the axis limits
        plt.xlim(x_limits)
        plt.ylim(y_limits[0], y_limits[1] + 1)
        
        # Only label every other tick
        for label in figure_axes.xaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        for label in figure_axes.yaxis.get_ticklabels()[::2]:
            label.set_visible(False)
        # Set The Font
        plt.xticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5)
        plt.yticks(fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size - 5)
        
        # Afternoon histogram
        noon_hist_fig, figure_axes = plt.subplots()
        plt.hist(self.ISI_noon, bins = self.ISI_edges, alpha = self.font_specs.hist_transparency, 
                 edgecolor = 'k', color = (0.5, 0, 0.5), label = 'Afternoon')
        noon_title_string = 'Interstimulus Intervals - ' + self.unit_name + ' (Afternoon)'
        plt.title(noon_title_string, fontname = self.font_specs.font_name, fontsize = self.font_specs.title_font_size, fontweight = 'bold')
        plt.xlabel('ISI (mSec)', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        plt.ylabel('Counts', fontname = self.font_specs.font_name, fontsize = self.font_specs.label_font_size)
        
        # Annotation of the fractional contamination
        if all(np.isreal(self.noon_fract_contam)):
            plt.text(0.85, 0.825, u'f\u1D56\u2081' + ' = ' + str(round(min(self.noon_fract_contam), 3)), 
                     verticalalignment = 'center', horizontalalignment = 'center', 
                     transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
        if any(np.iscomplex(self.noon_fract_contam)):
            plt.text(0.85, 0.825, u'f\u1D56\u2081' + ' < 0.5', verticalalignment = 'center', horizontalalignment = 'center', 
                     transform = figure_axes.transAxes, fontname = self.font_specs.font_name, fontsize = self.font_specs.legend_font_size)
        
        # Reset the axis limits
        plt.xlim(x_limits)
        plt.ylim(y_limits[0], y_limits[1] + 1)
        
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
            for ii in range(3):
                if ii == 0:
                    fig_title = noon_title_string
                if ii == 1:
                    fig_title = morn_title_string
                if ii == 2:
                    fig_title = merged_title_string
                plt.savefig(save_dir + fig_title + '.' + Save_Figs)
                plt.close()
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
        
