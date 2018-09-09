%	OVERVIEW:
%       This basic demo will allow you to load an ECG file in matlab 
%       compatible wfdb format, detect the locations of the R peaks,
%       perform signal quality (SQI) analysis and plot the results.
%
%   OUTPUT:
%       A figure with the loaded ECG signal and detected peaks will be
%       generated
%
%   DEPENDENCIES & LIBRARIES:
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   REFERENCE: 
%       Vest et al. "An Open Source Benchmarked HRV Toolbox for Cardiovascular 
%       Waveform and Interval Analysis" Physiological Measurement (In Press), 2018. 
%	REPO:       
%       https://github.com/cliffordlab/PhysioNet-Cardiovascular-Signal-Toolbox
%   ORIGINAL SOURCE AND AUTHORS:     
%       Giulia Da Poian   
%	COPYRIGHT (C) 2018 
%   LICENSE:    
%       This software is offered freely and without warranty under 
%       the GNU (v3 or later) public license. See license file for
%       more information
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear all
clc
close all %% Added by Mahmoud Zeydabadinezhad,

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comments added by Mahmoud Zeydabadinezhad on 09/7/2018
% The function rdmat() was not in the path, so I to bypass that function.
% In fact, I used the MATLAB function "load" to read data in.
% Lines started with %%% are the ones commented out by me.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Where are the data, in this demo they are located in a subfolder
%%%InputFolder = [pwd filesep 'TestData' filesep 'mitdb-Arrhythmia']; % path to the folder where you data are located

newInputFolder = [strrep(mfilename('fullpath'),mfilename,'') 'TestData' filesep] % path to the new folder where you can find some .mat files.
%%%SigName = '200m';
newSigName = 'TestRawECG.mat'; % This is a random file that I picked to work on it.

% load the 'TestRawECG.mat' file
load([newInputFolder filesep newSigName]);
%%%[tm,sig,Fs] = rdmat([InputFolder filesep SigName]);
% the signal has two channels, from now on we will use just one 
%%%ecg = sig(:,1);
ecg = signal;
% plot the signal
figure(1)
%%%plot(tm,ecg);
plot(ecg)
xlabel('[s]');
ylabel('[mV]')


% Detection of the R-peaks using the jqrs.m function included in the
% toolbox, requires to set initialization parameters calling the
% InitializeHRVparams.m function

HRVparams = InitializeHRVparams('Demo');
% set the exact sampling frequency usign the one from the loaded signal
%%%HRVparams.Fs = Fs;
HRVparams.Fs = 360; % Set the sampling frequency to a number.

% call the function that perform peak detection
r_peaks = jqrs(ecg,HRVparams);

% plot the detected r_peaks on the top of the ecg signal
figure(1)
hold on;
%%%plot(r_peaks./HRVparams.Fs, ecg(r_peaks),'o');
plot(r_peaks, ecg(r_peaks),'o'); % There is no need to divide the signal by Fs.
legend('ecg signal', 'detected R peaks')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% From here on is the code added by Mahmoud Zeydabadinezhad on 09/7/2018 as
% part of the lab work for BMI500.
% What I'm doing here is to segment a quarter of a second before and after
% each peak and then put them into a matrix called peakSegmnt.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

nSample = 0.25*HRVparams.Fs; % This variable determines how many samples of signal are equal to the quarter of a second.
peakSegmnt = zeros(length(r_peaks)-2,2*nSample+1); % Memory allocation for matrix peakSegmnt. To avoid falling out of the range, we do not count the first and last peak.
figure, hold on % Depict the segmented peaks.
title('All segmented peaks');
xlabel('[s]');
xticks([0 nSample 2*nSample+1])
xticklabels({'-0.25','0','0.25'})
ylabel('[mV]')

for i=2:length(r_peaks)-2 % To avoid falling out of the range, we do not count the first and last peak.
    peakSegmnt(i,:) = ecg(r_peaks(i)-nSample:r_peaks(i)+nSample);
    plot(peakSegmnt(i,:));
end

meanPeak = mean(peakSegmnt, 1); % Average of all segmented paeks.
figure, plot(meanPeak); % Plot the average of all peaks
title('Mean of all peaks');
ylabel('[mV]')
xlabel('[s]');
xticks([0 nSample 2*nSample+1])
xticklabels({'-0.25','0','0.25'})

%% Here is the end of the script %%