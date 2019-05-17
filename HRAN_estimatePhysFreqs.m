function HRAN_estFreqs_outputParams = HRAN_estimatePhysFreqs(physiologicalData, inputParams, waitbarBoolean)
% HRAN_ESTIMATEPHYSFREQS estimates the fundamental cardiac and respiratory
% frequencies from a time series with physiological noise
%
%% Input (inputParams)
%
%   % Physiological data
%   physiologicalData: time series of physiological data (time x 1)
%
%   % TR, moving window length, percent overlap
%   .TR: TR of fMRI data
%   .windowLength: length of moving window (s) (e.g 24-30)
%   .percentOverlap: percent overlap of windows (e.g. .5 or .75)
%
%   % Neural regressors of design matrix
%   .neuralZ: neural regressors (time x # regressors)
%
%   % Physiological frequency range
%   .cardiacFreqRange: fundamental cardiac freq range (bpm) (e.g. [40:80])
%   .respFreqRange: fundamental respiratory freq range (bpm) (e.g. [8:24])
%
%   % Model orders to estimate physiological frequencies
%   .P_freq: AR order (e.g. 1-4)
%   .R_freq: Resp order (e.g. 1-4)
%   .C_freq: Cardiac order (e.g. 1-3)
%   .N_freq: Number of neural regressors (e.g. 0)
%   .X_freq: Interaction order (e.g. 0-1)
%
%   % Model orders to create design matrix and regress physio noise from
%   brain
%   .P_Z: AR order (e.g. 1-4)
%   .R_Z: Resp order (e.g. 1-4)
%   .C_Z: Cardiac order (e.g. 1-3)
%   .N_Z: Number of neural regressors (e.g. size(neuralZ,2))
%   .X_Z: Interaction order (e.g. 0)
%
%   % Boolean to display waitbar or not
%   waitbarBoolean: 1 (on), 0 (off)
%
%% Output (HRAN_estFreqs_outputParams)
%
%   % Estimated fundamental frequencies
%   .t: time points of estimated frequencies (s)
%   .w_hr_hat: estimated fundamental cardiac freq (bpm)
%   .w_rr_hat: estimated fundamental respiratory freq (bpm)
%   .W_freq: matrix of estimated frequencies (including harmonics) used to
%   initialize prior covariance (# regressors x number of windows)
%
%   % Design matrices used to regress physiological noise across brain
%   .Z_hat: estimated design matrix including drift terms, neural
%   regressors, physio noise (indices in one window x # regressors x number
%   of windows)
%
%   % Windowing parameters
%   .windowIndices = indices in data in which to apply windows
%   ([startIndex endIndex] x number of windows)
%   .windows = windows to be applied in each segment (indices in one window
%   x number of windows)


%% 1) Initialize parameters
TR = inputParams.TR;

% Parameters for windowing
windowLength = (floor(inputParams.windowLength/TR)); % length of window in indices

% Make sure windows sum to 1
if inputParams.percentOverlap == .5 || inputParams.percentOverlap == .75
    if ~mod(windowLength,2); windowLength = windowLength + 1; end %make window odd length (so that window symmetric around middle)
    win = hann(windowLength); % initialize window (default is hann)
else
    disp('Note: If percent overlap not .5 or .75, ensure that windows sum to 1')
end

% Make sure data is the right dimension
if size(physiologicalData,2) > size(physiologicalData,1)
    physiologicalData = physiologicalData';
end

% Compute appropriate number of windows
overlapWindowLength = floor(inputParams.percentOverlap*windowLength); % overlapping window length
nonOverlapWindowLength = floor((1-inputParams.percentOverlap)*windowLength); % non-overlapping window length
paddedData = [zeros(overlapWindowLength,1); physiologicalData; zeros(overlapWindowLength,1)]; % pad data for easy indexing
numWindows = (length(paddedData)-1-overlapWindowLength)/nonOverlapWindowLength; %number of windowed segments
if numWindows-floor(numWindows) > 0; numWindows = floor(numWindows) + 1;end % add entry to account for rounding
T = 0:TR:(length(physiologicalData)-1)*TR; %initialize time array

% Parameters for cyclic descent
w_rr_vector = inputParams.respFreqRange./60; % vector of RR
num_w_rr = length(w_rr_vector); % number of RR to iterate through
w_hr_vector = inputParams.cardiacFreqRange./60; %vector of HR
num_w_hr = length(w_hr_vector); % number of HR to iterate through
C_likelihood_vector = zeros(num_w_rr,num_w_hr); % concentrated likelihood
priorCov = 10^3.*eye(1+1+inputParams.N_freq+2*inputParams.C_freq+2*inputParams.R_freq+4*inputParams.X_freq); % prior covariance (for estimating phys freqs, set all equal)

% Number of beta parameters in output design matrix
% 1 (Mean) + 1 (Linear) + N_Z (neural regressors) + 2*C_Z (cardiac
% terms) + 2*R_Z (respiratory terms) + 4*X (interaction terms)
num_beta_parameters = 1+1+inputParams.N_Z+2*inputParams.C_Z+2*inputParams.R_Z+4*inputParams.X_Z;

% Parameters to output
t = zeros(1, numWindows); % centered time of each window
w_rr_hat = zeros(1,numWindows); % optimal w_rr per window
w_hr_hat = zeros(1,numWindows); % optimal w_hr per window
W_freq = zeros(num_beta_parameters,numWindows); % estimated frequencies and harmonics per segment
Z_hat = zeros(windowLength,num_beta_parameters,numWindows); % optimal design matrix per segment
windowIndices = zeros(2,numWindows); %data indices to iterate through per segment
windows = zeros(windowLength,numWindows); % windows to apply per segment

% Multiply physiological data by constant (because seems to be rounding
% error in MATLAB if not high enough amplitude)
physiologicalData = physiologicalData*1000;

%% 2) Iterate through physiological data time series in windows and estimate physiological frequencies
if waitbarBoolean
    wait = waitbar(0,'Estimating Frequencies...'); % initialize waitbar
end
for n = 1:numWindows
    
    %% Select appropriate indices of data
    startIndex = (n-1)*nonOverlapWindowLength+1;
    dataIndex = startIndex:startIndex+windowLength-1;
    
    % If at beginning segments...
    if startIndex < overlapWindowLength
        % Use first half window of data to predict frequencies (same amount of data as hann window)
        dataIndex = 1:floor(windowLength/2);
        % Store indices for Z, data (on which window will be applied)
        windowIndices(:,n) = [1 windowLength];
        % Store window to be applied on Z, data
        Z_Window = zeros(windowLength-startIndex-nonOverlapWindowLength,1);
        Z_Window = [win(windowLength+1-startIndex-nonOverlapWindowLength:end);Z_Window];
        windows(:,n) = Z_Window; %store applied window
        % Apply current window (in this case ones) to data
        appliedWindow = ones(length(dataIndex),1);
        %appliedWindow = Z_Window;
        Q_inv = eye(length(dataIndex)); % inverse covariance matrix
        
        % If at end segments...
    elseif startIndex+windowLength > length(physiologicalData)+overlapWindowLength+1
        % Use last windowLengthSeconds/2 of data to predict frequencies (same amount of data as hann window)
        dataIndex = length(physiologicalData)-floor(windowLength/2)+1:length(physiologicalData);
        % Store indices for Z, data (on which window will be applied)
        windowIndices(:,n) = [length(physiologicalData)-windowLength+1 length(physiologicalData)];
        % Store window to be applied on Z, data
        Z_Window = zeros(windowLength-(length(physiologicalData)-startIndex+overlapWindowLength+1),1);
        Z_Window = [Z_Window;win(1:length(physiologicalData)-startIndex+overlapWindowLength+1)];
        windows(:,n) = Z_Window; %store applied window
        % Apply current window (in this case ones) to data
        appliedWindow = ones(length(dataIndex),1);
        Q_inv = eye(length(dataIndex)); % inverse covariance matrix
        
        % For all other segments...
    else
        % Use windowLength data segments with hann window
        dataIndex = dataIndex-overlapWindowLength;
        % Store indices for Z, data (on which window will be applied)
        windowIndices(:,n) = [dataIndex(1) dataIndex(end)];
        % Store window to be applied on Z, data
        Z_Window = win;
        windows(:,n) = Z_Window; %store applied window
        % Apply current window (in this case ones) to data
        appliedWindow = win;%ones(length(dataIndex),1); %% BIG CHANGE
        Q_inv = eye(length(dataIndex)); % inverse covariance matrix
    end
    
    t(n) = T(floor(mean(windowIndices(:,n))));
    data = physiologicalData(dataIndex); % index through data
    
    if waitbarBoolean
        waitbar((n-1)/numWindows,wait) % increment the waitbar
    end
    
    %% Iterate through frequency grid and run cyclic descent algorithm
    for w_i_rr = 1:num_w_rr
        w_rr = w_rr_vector(w_i_rr); %set RR
        for w_j_hr = 1:num_w_hr
            w_hr = w_hr_vector(w_j_hr); %set HR
            
            % Initialize Design Matrix Z
            [Z, Z_illconditioned] = HRAN_createDesignMatrix(T(dataIndex),inputParams.neuralZ(dataIndex,:),w_hr,w_rr,inputParams.R_freq,inputParams.C_freq,inputParams.N_freq,inputParams.X_freq);
            
            % Run Cyclic Descent Algorithm on windowed segment
            if ~Z_illconditioned
                [~,~,~,C_likelihood,~] = HRAN_CyclicDescentAlgorithm(Z.*appliedWindow, inputParams.P_freq, data.*appliedWindow, Q_inv, priorCov);
            else
                C_likelihood = nan;
            end
            
            % Store concentrated likelihood for each frequency
            C_likelihood_vector(w_i_rr,w_j_hr) = C_likelihood;
        end
    end
    
    %% Optimize concentrated likelihood across all frequencies
    minC_likelihood = min(C_likelihood_vector(:)); %minimized concentrated negative log likelihood
    [w_rr_index,w_hr_index] = find(C_likelihood_vector==minC_likelihood); % w_hat_hr, w_hat_rr indices
    w_rr_hat(n) = w_rr_vector(w_rr_index); % w_rr
    w_hr_hat(n) = w_hr_vector(w_hr_index); % w_hr
    
    %% Store optimal Z_hat and W_freq
    [Z, ~] = HRAN_createDesignMatrix(T(windowIndices(1,n):windowIndices(2,n)),inputParams.neuralZ(windowIndices(1,n):windowIndices(2,n),:),w_hr_hat(n),w_rr_hat(n),inputParams.R_Z,inputParams.C_Z,inputParams.N_Z,inputParams.X_Z);
    Z_hat(:,:,n) = Z;
    
    W_f = HRAN_estimateAliasedFreq(TR,Z,w_hr_hat(n),w_rr_hat(n),inputParams.R_Z,inputParams.C_Z,inputParams.N_Z);
    W_freq(:,n) = W_f;
end

%% Output parameters for physio regression
HRAN_estFreqs_outputParams = struct;
HRAN_estFreqs_outputParams.t = t; % time point (s) of each window
HRAN_estFreqs_outputParams.w_rr_hat = w_rr_hat; % estimated respiratory frequency
HRAN_estFreqs_outputParams.w_hr_hat = w_hr_hat; % estimated cardiac frequency
HRAN_estFreqs_outputParams.Z_hat = Z_hat; % estimated design matrix
HRAN_estFreqs_outputParams.windowIndices = windowIndices; % indices of each window
HRAN_estFreqs_outputParams.windows = windows; % windows to apply
HRAN_estFreqs_outputParams.W_freq = W_freq; % frequencies of physio noise

% Also output redundant parameters so everything contained in HRAN_estFreqs_outputParams 
HRAN_estFreqs_outputParams.TR = TR; % TR of fMRI data
HRAN_estFreqs_outputParams.percentOverlap = inputParams.percentOverlap; % percent overlap
HRAN_estFreqs_outputParams.P_Z = inputParams.P_Z; % AR Order
HRAN_estFreqs_outputParams.R_Z = inputParams.R_Z; % Resp model order
HRAN_estFreqs_outputParams.C_Z = inputParams.C_Z; % Cardiac model order
HRAN_estFreqs_outputParams.N_Z = inputParams.N_Z; % Number of neural regressors
HRAN_estFreqs_outputParams.X_Z = inputParams.X_Z; % Interaction order

if waitbarBoolean
    close(wait);
end

end

