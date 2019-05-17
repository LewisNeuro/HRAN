%% This is a demo script that shows how to run HRAN on simulated data
% The sections proceed as follows:
%   1) Simulate physiological noise voxel (e.g. ventricles)
%   2) Simulate voxel with neural activity and physiologial noise (e.g. V1)
%   3) Estimate frequencies from voxel with physiological noise
%   4) Remove physiological noise from voxel with neural activity

% Make sure chronux is downloaded! (See readme.md for instructions)
% Add path to chronux toolbox (used to compute multitaper spectra)
addpath(genpath('chronux'))

%% 1) Simulate physiological noise voxel
% First, we will create a simulate physiological noise time series which HRAN will 
% use to estimate the fundamental frequencies of cardiac and respiratory activity. 

% In practice, this physiological noise time series can be obtained by,
% for example, masking the fMRI data to physiologically noisy regions (e.g. ventricles)

% Initialize parameters
TR = .25; % repetition time (s)
T = 120; %total time (s)
time = 0:TR:T-TR; %time vector

% Cardiac and Respiratory Noise
cardiacFrequency = 1; % Hz
respiratoryFrequency = .3; % Hz
cardiacNoise = 10*cos(2*pi*cardiacFrequency.*time); 
respiratoryNoise = 20*sin(2*pi*respiratoryFrequency.*time);

% AR(1) background noise
alpha = .9; 
ARNoise = zeros(size(time));
ARNoise(1) = 5*randn(1);
for t=2:length(time)
    ARNoise(t) = alpha*ARNoise(t-1) + 5*randn(1);
end

% Physiological Noise
physiologicalNoise = cardiacNoise + respiratoryNoise + ARNoise;

%% 2) Simulate voxel with neural activity
% Now, we will simulate a voxel with oscillatory neural activity, and the
% physiological noise from our previous simulation.

% Neural Data
neuralFrequency = .1; % Hz
neuralSignal = 60*cos(2*pi*neuralFrequency.*time);
neuralData = neuralSignal + physiologicalNoise;

% Spectrograms of the simulated data
fig = figure('Position',[1 1 800 600]);
% Physio spectrogrm
subplot(2,2,1)
movingwin = [20 4];
params.Fs = 1/TR;
params.tapers = [2 3];
[spec_original, stime, sfreq] = mtspecgramc(detrend(physiologicalNoise), movingwin, params);
imagesc(stime, sfreq, 10*log10(spec_original)')
set(gca,'YDir','normal')
colormap('jet')
colorbar
ylabel('Freq (Hz)')
xlabel('Time (s)')
title('Physio Noise')

% Neural spectrogram
subplot(2,2,2)
[spec_original, stime, sfreq] = mtspecgramc(detrend(neuralData), movingwin, params);
imagesc(stime, sfreq, 10*log10(spec_original)')
set(gca,'YDir','normal')
colormap('jet')
colorbar
ylabel('Freq (Hz)')
xlabel('Time (s)')
title('Neural Data')

% Physio time series
subplot(2,2,3)
plot(time,physiologicalNoise,'r','LineWidth',2)
ylabel('Amplitude (a.u.)')
xlabel('Time (s)')

% Neural time series
subplot(2,2,4)
plot(time,neuralData,'k','LineWidth',2)
ylabel('Amplitude (a.u.)')
xlabel('Time (s)')

%% 3) Use HRAN to estimate frequencies 
% Using HRAN, we can estimate the fundamental cardiac and respiratory
% frequencies. First, we need to initialize the following parameters:

% Physiological Data
physiologicalData = physiologicalNoise'; % time series of physiological data

% Struct to store input parameters
inputParams = struct;

% TR, moving window length, percent overlap
inputParams.TR = TR; %TR of fMRI data (s)
inputParams.windowLength = 24; % length of moving window (s) (e.g. 24-30s)
inputParams.percentOverlap = .5; % percent overlap of windows (e.g. .5 or .75)

% Neural regressors of design matrix (time x # regressors)
% Estimated neural signal components to be used in GLM
% e.g. neural stimulus convolved with HRF
neuralZ = [cos(2*pi*neuralFrequency.*time)';sin(2*pi*neuralFrequency.*time)'];
inputParams.neuralZ = neuralZ;

% Physiological frequency range
inputParams.cardiacFreqRange = [40:80]; % cardiac fundamental freq range (bpm)
inputParams.respFreqRange = [8:24]; %resp fundamental freq range (bpm)

% Model orders used to estimate the cardiac and respiratory frequencies
% from physiological time series
inputParams.P_freq = 1; % AR order used to estimate physio frequencies (e.g. 1-4)
inputParams.R_freq = 1; % Resp order used to estimate physio frequencies (e.g. 1-4)
inputParams.C_freq = 1; % Cardiac order used to estimate physio frequencies (e.g. 1-3)
inputParams.N_freq = 0; % Number of neural regressors used to estimate physio freqs (e.g. 0 bc assuming no neural signal in this region)
inputParams.X_freq = 0; % Interaction order used to estimate physio frequencies (e.g. 0-1)

% Model orders used to create design matrix to regress physio noise from
% each voxel in the brain
inputParams.P_Z = 1; % AR order used to regress physio noise (e.g. 1-4)
inputParams.R_Z = 1; % Resp order used to regress physio noise (e.g. 1-4)
inputParams.C_Z = 1; % Cardiac order used to regress physio noise (e.g. 1-3)
inputParams.N_Z = size(neuralZ,2); % Number of neural regressors used to regress physio noise (e.g. number of columns in neuralZ)
inputParams.X_Z = 0; % Interaction order used to regress physio noise (e.g. 0-1)

% Toggle for waitbar on or off
waitbarBoolean = 1; % 1 - on, 0 - off

% Run function to estimate physiological frequencies
HRAN_estFreqs_outputParams = HRAN_estimatePhysFreqs(physiologicalData,inputParams,waitbarBoolean);

%% Plot estimated frequencies with simulated frequencies
fig = figure('Position',[1 1 800 600]);

% Physio spectrogrm
subplot(1,2,1)
movingwin = [20 4];
params.Fs = 1/TR;
params.tapers = [2 3];
[spec_original, stime, sfreq] = mtspecgramc(detrend(physiologicalNoise), movingwin, params);
hold on
imagesc(stime, sfreq, 10*log10(spec_original)')
scatter(HRAN_estFreqs_outputParams.t,HRAN_estFreqs_outputParams.w_hr_hat,'w','filled')
scatter(HRAN_estFreqs_outputParams.t,HRAN_estFreqs_outputParams.w_rr_hat,'w','filled')
hold off
set(gca,'YDir','normal')
colormap('jet')
colorbar
xlim([stime(1) stime(end)])
ylim([sfreq(1) sfreq(end)])
ylabel('Freq (Hz)')
xlabel('Time (s)')
title('Physio Noise')

% Simulated frequencies
hold on
subplot(1,2,2)
hold on
plt = plot(time,repmat(cardiacFrequency,length(time),1),'r','LineWidth',4);
plt.Color(4) = .5;
plt = plot(time,repmat(respiratoryFrequency,length(time),1),'b','LineWidth',4);
plt.Color(4) = .5;
scatter(HRAN_estFreqs_outputParams.t,HRAN_estFreqs_outputParams.w_hr_hat,'r')
scatter(HRAN_estFreqs_outputParams.t,HRAN_estFreqs_outputParams.w_rr_hat,'b')
hold off
xlabel('Time (s)')
ylabel('Physio Freq (Hz)')
legend('Sim Cardio','Sim Resp','Est Cardio','Est Resp')
ylim([sfreq(1) sfreq(end)])

%% 4) Use HRAN to remove physiological noise from separate voxel
% Now that physiological frequencies are known, can regress out
% physiological noise from across the brain (or in this demo, another
% brain ROI). All required parameters are stored in
% HRAN_estFreqs_outputParams, as described here: 

% Data to be de-noised (note: add functionality for 4d nifti)
data = neuralData';

% Parameters stored in HRAN_estFreqs_outputParams:
HRAN_estFreqs_outputParams.TR; % TR of fMRI data
HRAN_estFreqs_outputParams.Z_hat; % design matrices across time
HRAN_estFreqs_outputParams.W_freq; % physiological frequencies across time
HRAN_estFreqs_outputParams.windowIndices; % indices of each segment
HRAN_estFreqs_outputParams.windows; % windows applied to each segment
HRAN_estFreqs_outputParams.P_Z; % AR order
HRAN_estFreqs_outputParams.R_Z; % respiratory order
HRAN_estFreqs_outputParams.C_Z; % cardiac order
HRAN_estFreqs_outputParams.N_Z; % number of neural regressors
HRAN_estFreqs_outputParams.X_Z; % number of interaction terms
HRAN_estFreqs_outputParams.percentOverlap; % percent overlap of windows

% Run HRAN to remove physiological noise
[deNoisedData,HRAN_neural,HRAN_cardiac,HRAN_resp] = HRAN_regressPhysNoise(data,HRAN_estFreqs_outputParams);

%% Plot results

% Simulated Neural Data
fig = figure('Position',[1 1 800 600]);
subplot(5,2,[1 3])
[spec_original, stime, sfreq] = mtspecgramc(detrend(neuralData), movingwin, params);
imagesc(stime, sfreq, 10*log10(spec_original)')
set(gca,'YDir','normal')
colormap('jet')
colorbar
ylabel('Freq (Hz)')
xlabel('Time (s)')
title('Neural Data')

% De-noised neural data
subplot(5,2,[2 4])
[spec_original, stime, sfreq] = mtspecgramc(detrend(deNoisedData), movingwin, params);
imagesc(stime, sfreq, 10*log10(spec_original)')
set(gca,'YDir','normal')
colormap('jet')
colorbar
ylabel('Freq (Hz)')
xlabel('Time (s)')
title('De-noised Neural Data')

% Time series plots
subplot(5,2,[5 6])
hold on
plot(time, deNoisedData, 'g', 'LineWidth', 2)
plot(time, neuralData, 'k--', 'LineWidth', 2)
hold off
ylabel('De-noised data')
subplot(5,2,[7 8])
hold on
plot(time, HRAN_cardiac, 'r', 'LineWidth', 2)
plot(time, cardiacNoise, 'k--', 'LineWidth', 2)
ylabel('Cardiac')
hold off
subplot(5,2,[9 10])
hold on
plot(time, HRAN_resp, 'b', 'LineWidth', 2)
plot(time, respiratoryNoise, 'k--', 'LineWidth', 2)
hold off
ylabel('Resp')
xlabel('Time (s)')

%% 5) Examine residuals
% If the Econometrics toolbox is downloaded, can also request the residuals
% as an additional output. (Note: this will take longer)

% As before, but with residuals as fifth parameter
[deNoisedData,HRAN_neural,HRAN_cardiac,HRAN_resp,HRAN_residuals] = HRAN_regressPhysNoise(data,HRAN_estFreqs_outputParams);


fig = figure('Position',[1 1 300 480]);

% 1) Plot NCP
[P_r,f] = periodogram(HRAN_residuals,rectwin(length(HRAN_residuals)),length(HRAN_residuals),1/TR);
cum_P_r = cumsum(P_r)./sum(P_r);
[P_y,f] = periodogram(detrend(data),rectwin(length(HRAN_residuals)),length(HRAN_residuals),1/TR);
cum_P_y = cumsum(P_y)./sum(P_y);
[P_p,f] = periodogram(HRAN_cardiac + HRAN_resp,rectwin(length(HRAN_residuals)),length(HRAN_residuals),1/TR);
cum_P_p = cumsum(P_p)./sum(P_p);
[P_d,f] = periodogram(detrend(deNoisedData),rectwin(length(HRAN_residuals)),length(HRAN_residuals),1/TR);
cum_P_d = cumsum(P_d)./sum(P_d);
CI_bounds = 1.36*(1/sqrt(length(f)));
params.Fs = 1/TR;
params.tapers = [2 3];

subplot(2,1,1)
hold on
plot(f, linspace(0,1,length(f)),'k','LineWidth',2)
plot(f,cum_P_p,'LineWidth',2,'Color',[215,48,39]./256)
plot(f,cum_P_d,'LineWidth',2,'Color',[118,42,131]./256)
plot(f,cum_P_r,'LineWidth',2,'Color',[115,115,115]./256)
plot(f,cum_P_y,'LineWidth',2,'Color',[77,146,33]./256)
plot(f, f./f(end) + CI_bounds, 'k--','LineWidth',1)
plot(f, f./f(end) - CI_bounds, 'k--','LineWidth',1)
xlabel('Frequency (Hz)')
ylabel('NCP')
set(gca,'TickLength',[0 0])
ax = gca;
ax.FontSize = 10;
ax.FontName = 'Arial';
title('NCP')
ylim([0 1])
legend('Ideal White Noise','Physio','De-noised','Residuals','Data','Location','Southeast')
 
% 2) Plot QQ plot
subplot(2,1,2)
qq = qqplot(HRAN_residuals);
set(qq, 'Color', 'k', 'MarkerSize', 4, 'MarkerEdgeColor', [115,115,115]./256)
set(qq(2), 'Color', 'k', 'MarkerSize', 4, 'MarkerEdgeColor', 'k')
xlabel('Standard Normal Quantiles')
ylabel('Quantiles of Input Sample')
set(gca,'TickLength',[0 0])
ax = gca;
ax.FontSize = 10;
ax.FontName = 'Arial';
