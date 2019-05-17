function BIC = HRAN_estimatePhysiologicalRegressorsBIC(data,TR,Z_hat,W_freq,windowIndices,windows,P,R,C,N,X)
%HRAN_REGRESSPHYSNOISE removes physiological noise from input data
%
%% Input (HRAN_estFreqs_outputParams)
%   % Data
%   data: data from which physio noise will be removed (time x 1)
%
%   % Windowing parameters
%   .TR: TR of fMRI data (s)
%   .percentOverlap: percent overlap of subsequent windows (e.g. .5 or .75)
%   .windowIndices: indices of each applied window
%   .windows: applied window to each segment
%
%   % Model orders
%   .P_Z: AR Order (e.g. 1-4)
%   .R_Z: Respiratory order (e.g. 1-4)
%   .C_Z: Cardiac order (e.g. 1-3)
%   .N_Z: Number of neural regressors
%   .X_Z: Interaction order (e.g. 0-1)
%
%   % Design Matrix
%   .Z_hat: Design matrix for each window (time x #regressors x #windows)
%
%   % Physio frequencies
%   .W_freq: physio freqs across time (to estimate prior covariance)
%
%% Output
%   deNoisedData: data minus estimated physiological noise
%   varargout{1} -- neuralSignal: estimated neural signal
%   varargout{2} -- cardiacNoise: estimated cardiac noise
%   varargout{3} -- respiratoryNoise: estimated respiratory noise
%   varargout{4} -- residual: estimated residual of model fit (note: time
%   intensive, requires econometrics toolbox in matlab!)

%% 1) Initialize Parameters
% Store for readability
detrendData = detrend(data); % detrend data (added back later)
numSegments = size(windowIndices,2); % number of segments to iterate through data
BIC = zeros(size(windowIndices,2),1); % estimate BIC
Q_inv = eye(windowIndices(2,1)-windowIndices(1,1)+1); % initialize AR covariance matrix as I

%% Iterate through data in segments
for n = 1:numSegments
    
    % Initialize parameters per segment
    Z = Z_hat(:,:,n).*windows(:,n); % window the design matrix
    dataSegment = detrendData(windowIndices(1,n):windowIndices(2,n)).*windows(:,n); % window the data
    W_f = W_freq(:,n);
    
    % Run cyclic descent algorithm on segment
    W = HRAN_estimatePriorCovariance(TR,dataSegment,Z,W_f,N);
    [Beta_hat,alpha_hat,sigmaSquared_hat,~,Q_inv] = HRAN_CyclicDescentAlgorithm(Z,P,dataSegment,Q_inv,W);
    
    % If noise removed correctly...
    if ~isnan(Beta_hat)
        
        % Compute BIC
        numParameters = 1 + 2*R + 2*C + P + N + 4*X;
        BIC(n) = length(dataSegment)*log10(sigmaSquared_hat) + numParameters*log10(length(dataSegment));
        
    else
        BIC(n) = nan;
    end
    
end

end

