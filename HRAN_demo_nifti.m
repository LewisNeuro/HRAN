%% This is a demo script to load a nifti, remove physio noise, display visual cortex data, and write data back to nifti
% Hopefully, it will help you to set up your particular analysis pipeline
% The sections proceed as follows:
%   1) Load nifti file and physiological noise time series
%   2) Estimate physiological frequencies
%   3) Regress physiological noise from each voxel in the nifti
%   4) Write de-noised nifti file to output to other programs

% Make sure chronux is downloaded! (See readme.md for instructions)
% Add path to chronux toolbox (used to compute multitaper spectra)
addpath(genpath('chronux'))

%% 1) Load nifti file and physiological noise time series

% First, we will load a nifti file. 

% The demo nifti has the following properties (Experiment C from Agrawal et. al. 2019):
% TR = .347s
% Dimensions: 76 (x) x 76 (y) x 10 (z) x 450 (t)
% Preprocessing: slice time corrected, motion corrected, smoothed with 5 mm FWHM Gaussian
% Functional task: presentation of visual stimulus at 0.13 Hz

% Load functional data and info
niftiFile = 'demo_data/demoData.nii.gz';
V = double(niftiread(niftiFile));
V_info = niftiinfo(niftiFile);

% Next, we will obtain a physiological time series
% We can obtain this physiological time series in many different ways.

% Here, we extracted a WM mask using FreeSurfer (Fischl 2012) and
% registered it to this functional run.
WM_mask = boolean(niftiread('demo_data/WhiteMatter_mask.nii'));

% Now, we can extract the appropriate physiological time series in MATLAB.
physiologicalData = squeeze(sum(WM_mask.*V,[1 2 3]))./sum(WM_mask,[1 2 3]);

% Alternatively, we could use any software to identify anatomical masks (eg
% CONN toolbox) and any software to extract the masks efficiently (eg FSL). 
% For example, in FSL the command to extract an ROI is: 
% fslmeants -i FUNC_RUN -o OUTPUT_FILE -m MASK_NIFTI

%% 2) Estimate physiological frequencies

% Now that we have the physiological time series, we will use it to
% estimate the physiological frequencies.

% 1) Initialize struct to store input parameters
inputParams = struct;

% 2) Set TR, moving window length, percent overlap
TR = V_info.PixelDimensions(4); %0.347s in the demo data
inputParams.TR = TR; %TR of fMRI data (s)
inputParams.windowLength = 30; % length of moving window (s) (e.g. 24-30s)
inputParams.percentOverlap = .75; % percent overlap of windows (e.g. .5 or .75)

% 3) Set neural regressors of design matrix (time x # regressors)
% In general, one can convolve the neural stimulus with HRF and input into
% design matrix (or use whatever design matrix they wish)

% Here, we use a cos and sin at the particular neural frequency (which
% approximate the convolution at steady state)
time = [0:inputParams.TR:size(V,4)*inputParams.TR-inputParams.TR];
neuralFrequency = .13; % Hz
neuralSignal = 60*cos(2*pi*neuralFrequency.*time);
neuralZ = [cos(2*pi*neuralFrequency.*time)';sin(2*pi*neuralFrequency.*time)'];
inputParams.neuralZ = neuralZ;

% 4) Set physiological frequency range
inputParams.cardiacFreqRange = [50:80]; % cardiac fundamental freq range (bpm)
inputParams.respFreqRange = [6:18]; %resp fundamental freq range (bpm)

% 5) Set model orders used to estimate the cardiac and respiratory
% frequencies from the physiological time series
inputParams.P_freq = 2; % AR order used to estimate physio frequencies (e.g. 1-4)
inputParams.R_freq = 1; % Resp order used to estimate physio frequencies (e.g. 1-4)
inputParams.C_freq = 1; % Cardiac order used to estimate physio frequencies (e.g. 1-3)
inputParams.N_freq = 0; % Number of neural regressors used to estimate physio freqs (e.g. 0 bc assuming no neural signal in this region)
inputParams.X_freq = 0; % Interaction order used to estimate physio frequencies (e.g. 0-1)

% 5) Set model orders used to create design matrix to regress physio noise
% from each voxel in the brain
inputParams.P_Z = 2; % AR order used to regress physio noise (e.g. 1-4)
inputParams.R_Z = 2; % Resp order used to regress physio noise (e.g. 1-4)
inputParams.C_Z = 3; % Cardiac order used to regress physio noise (e.g. 1-3)
inputParams.N_Z = size(neuralZ,2); % Number of neural regressors used to regress physio noise (e.g. number of columns in neuralZ)
inputParams.X_Z = 0; % Interaction order used to regress physio noise (e.g. 0-1)

% Toggle for waitbar on or off
waitbarBoolean = 1; % 1 - on, 0 - off

% Run function to estimate physiological frequencies
HRAN_estFreqs_outputParams = HRAN_estimatePhysFreqs(physiologicalData,inputParams,waitbarBoolean);

%% Plot estimated frequencies over spectrogram of physiological data

% It is important to visually inspect model performance. In general, one  
% can check if the frequency estimates track the observed physiological
% power bands in the physiological ROI spectrogram.

% Here, we demonstrate how to plot the spectrogram with frequency estimates
% overlaid. Also, we show how these estimates compare with HR and RR 
% obtained from an EKG and respiratory belt. 

% (Note: here, the EKG had some artifact, which made the HR obtained from
% the EKG more variable)

fig = figure('Position',[1 1 800 600]);

% Spectrogram of physiological data
subplot(1,2,1)
movingwin = [24 4];
params.Fs = 1/TR;
params.tapers = [2 3];
[spec_original, stime, sfreq] = mtspecgramc(detrend(physiologicalData), movingwin, params);
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
title('Physiological Noise ROI')

% Estimated frequencies from physiological data
physio = load('demo_data/physio.mat');
hold on
subplot(1,2,2)
hold on
plot(physio.T,physio.HR,'Color',[178,24,43]./256)
plot(physio.T,physio.RR,'Color',[33,102,172]./256)
scatter(HRAN_estFreqs_outputParams.t,HRAN_estFreqs_outputParams.w_hr_hat.*60,'MarkerEdgeColor',[178,24,43]./256,'MarkerFaceColor',[178,24,43]./256)
scatter(HRAN_estFreqs_outputParams.t,HRAN_estFreqs_outputParams.w_rr_hat.*60,'MarkerEdgeColor',[33,102,172]./256,'MarkerFaceColor',[33,102,172]./256)
hold off
xlabel('Time (s)')
ylabel('Frequency (bpm)')
lgd = legend('EKG Cardio','Belt Resp','HRAN Cardio','HRAN Resp','Location','best');
ylim([0 90])

% Save figure
%saveas(gcf,'demo_data/estimatedFrequencies_nifti.png')

% It is worth noting here that the actual EKG and respiratory traces were
% noisy! (Whereas HRAN's estimates of HR and RR were not...)


%% 3) Regress physiological noise from each voxel in the nifti
% Now that physiological frequencies are known, can regress out
% physiological noise from across the brain. All required parameters are stored in
% HRAN_estFreqs_outputParams, as described here: 

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

%% Now we will iterate through each voxel, and remove the physiological noise
% this could take 15 min - 1 hr depending on the computer, parallelization,
% etc.

% Store dimensions
[xDim,yDim,zDim,tDim] = size(V);

% Initialize arrays to store de noised data
deNoisedV = zeros(size(V));

% Load brain mask (registered to functional data) to speed up de-noising
brainMask = boolean(niftiread('demo_data/Brain_mask.nii'));

% Iterate through data
parfor z = 1:zDim % NOTE: this is a parfor loop to take advantage of parallelization, can replace with regular loop
    for y = 1:yDim
        for x = 1:xDim
            if brainMask(x,y,z) % NOTE: this line not required, but can speed up de-noising (by only de-noising voxels in the brain)
                %['Voxel: ' num2str(x) ',' num2str(y) ',' num2str(z)] % should comment out this line when running!
                deNoisedVoxel = HRAN_regressPhysNoise(squeeze(V(x,y,z,:)),HRAN_estFreqs_outputParams);
                deNoisedV(x,y,z,:) = deNoisedVoxel;
            end
        end
    end
end

%% Display visual cortex with and without physiological noise

% Load visual cortex anatomical mask (registered to functional run)
visualCortexMask = boolean(niftiread('demo_data/VisualCortex_mask.nii'));

% Extract a) original time series and b) de-noised time series
originalVisualCortex = squeeze(sum(visualCortexMask.*V,[1 2 3]))./sum(visualCortexMask,[1 2 3]);
deNoisedVisualCortex = squeeze(sum(visualCortexMask.*deNoisedV,[1 2 3]))./sum(visualCortexMask,[1 2 3]);

% Plot time series, power spectra, and spectrograms of visual cortex
fig = figure('Position',[1 1 800 600]);

% Time Series
subplot(2,3,[1 2 3])
hold on
plot(time, originalVisualCortex, 'k','LineWidth',2)
plt = plot(time, deNoisedVisualCortex,'Color',[152,78,163]./256,'LineWidth',2);
plt.Color(4) = .6;
hold off
xlabel('Time (s)')
ylabel('Amplitude (a.u.)')
legend('Original V1','De-noised V1')

% Power Spectra
subplot(2,3,4)
movingwin = [24 4];
params.Fs = 1/TR;
params.tapers = [2 3];
[P_orig,f] = mtspectrumc(detrend(originalVisualCortex),params);
[P_deNoised,f] = mtspectrumc(detrend(deNoisedVisualCortex),params);
hold on
plot(f,10*log10(P_orig),'k','LineWidth',2)
plt = plot(f,10*log10(P_deNoised),'Color',[152,78,163]./256,'LineWidth',2);
plt.Color(4) = .6;
xlabel('Frequency (Hz)')
ylabel('Power (dB)')
hold off

% Spectrograms
% Original
subplot(2,3,5)
movingwin = [24 4];
params.Fs = 1/TR;
params.tapers = [2 3];
[spec_original, stime, sfreq] = mtspecgramc(detrend(originalVisualCortex), movingwin, params);
imagesc(stime, sfreq, 10*log10(spec_original)')
set(gca,'YDir','normal')
colormap('jet')
%colorbar
caxis([-30 15])
xlim([stime(1) stime(end)])
ylim([sfreq(1) sfreq(end)])
ylabel('Freq (Hz)')
xlabel('Time (s)')
title('Original')

% De-noised
subplot(2,3,6)
movingwin = [24 4];
params.Fs = 1/TR;
params.tapers = [2 3];
[spec_original, stime, sfreq] = mtspecgramc(detrend(deNoisedVisualCortex), movingwin, params);
imagesc(stime, sfreq, 10*log10(spec_original)')
set(gca,'YDir','normal')
colormap('jet')
%colorbar
caxis([-30 15])
xlim([stime(1) stime(end)])
ylim([sfreq(1) sfreq(end)])
ylabel('Freq (Hz)')
xlabel('Time (s)')
title('De-noised')

% Save figure
%saveas(gcf,'demo_data/deNoisedVisualCortex_nifti.png')

% Note: we included aparcaseg.nii in case users want to look at other brain
% regions

%% 4) Write de-noised nifti file to output to other programs

% Write initial nifti file
fileName_denoised = 'demo_data/demoData_deNoised';
niftiwrite(deNoisedV,fileName_denoised);

% Make sure info is correct
deNoisedV_info = niftiinfo(fileName_denoised);
deNoisedV_info.Description = 'denoised';
deNoisedV_info.PixelDimensions = V_info.PixelDimensions;
deNoisedV_info.ImageSize = size(deNoisedV);
deNoisedV_info.BitsPerPixel = V_info.BitsPerPixel;
deNoisedV_info.SpaceUnits = V_info.SpaceUnits;
deNoisedV_info.TimeUnits = V_info.TimeUnits;
deNoisedV_info.MultiplicativeScaling = V_info.MultiplicativeScaling;
deNoisedV_info.SpatialDimension = V_info.SpatialDimension;
deNoisedV_info.TransformName = V_info.TransformName;
deNoisedV_info.Transform = V_info.Transform;
deNoisedV_info.raw = V_info.raw;

% Rewrite nifti file with correct info
niftiwrite(deNoisedV,fileName_denoised,deNoisedV_info,'Compressed',true);