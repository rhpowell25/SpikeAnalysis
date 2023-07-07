function [SNR_matrix] = SignalToNoise(xds, unit_name, wave_choice)

%% Find the unit of interest
[N] = Find_Unit(xds, unit_name);

%% Define the exported spreadsheet

SNR_matrix = cell(2 ,4);
SNR_matrix{1,1} = 'Signal:Noise';
SNR_matrix{1,2} = 'Resting Potential RMS';
SNR_matrix{1,3} = 'Average Crest';
SNR_matrix{1,4} = 'Average Trough';

%% Catch possible sources of error
% If there is no unit of that name
if isempty(N)
    fprintf('%s does not exist \n', unit_name);
    SNR_matrix = NaN(2,4);
    return
end

%% Extracting the waveforms of the designated unit
if strcmp(wave_choice, 'Wave')
    waves = xds.spike_waveforms{N};
elseif strcmp(wave_choice, 'Nonlin')
    waves = xds.nonlin_waveforms{N};
end

%% Find the first ten indices of each spike
rest_potential = zeros(height(waves),10);
for ii = 1:height(rest_potential)
    rest_potential(ii,:) = waves(ii,1:10);
end

%% Finding the root mean squared of the noise
rms_rest_pot = zeros(length(rest_potential),1);
for ii = 1:length(rest_potential)
    rms_rest_pot(ii) = rms(rest_potential(ii,:));
end

avg_rms_rest_pot = mean(rms_rest_pot);

%% Finding the mean of the waveforms
%Calculating the means
avg_waves = mean(waves);

%% Crests and troughs
% Amplitude of troughs and crests
avg_trough = min(avg_waves);
avg_crest = max(avg_waves);
% Peak to peak amplitude
peak_peak_amp = avg_crest - avg_trough;

SNR = abs(peak_peak_amp / avg_rms_rest_pot);

%% Putting the values in the SNR matrix
SNR_matrix{2,1} = SNR;
SNR_matrix{2,2} = avg_rms_rest_pot;
SNR_matrix{2,3} = avg_crest;
SNR_matrix{2,4} = avg_trough;





