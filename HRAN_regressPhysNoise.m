function [deNoisedData,varargout] = HRAN_regressPhysNoise(data,inputParams)
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
N = inputParams.N_Z;
C = inputParams.C_Z;
R = inputParams.R_Z;
P = inputParams.P_Z;
X = inputParams.X_Z;
windowIndices = inputParams.windowIndices;
detrendData = detrend(data); % detrend data (added back later)
numSegments = size(windowIndices,2); % number of segments to iterate through data
deNoisedData = detrendData; % initialize data
nout = max(nargout,1) - 1; % if additional outputs requested

% Initialize variable output arguments
% varargout{1} -- neural signal
if nout > 0
    neuralSignal = zeros(length(detrendData),1);
end
% varargout{2} -- estimated cardiac noise
if nout > 1
    cardiacNoise = zeros(length(detrendData),1);
end
% varargout{3} -- estimated respiratory noise
if nout > 2
    respiratoryNoise = zeros(length(detrendData),1);
end
% varargout{4} -- estimated residual (note: time intensive! requires
% econometrics toolbox)
if nout > 3
    trend = zeros(length(detrendData),1); % mean + linear trend
    ARestimate = zeros(length(detrendData),1); % AR estimate
    residual = zeros(length(detrendData),1); % residual
end

% Initialize taper factor based on window overlap
if inputParams.percentOverlap == .5 || inputParams.percentOverlap == 0
    taperFactor = 1;
elseif inputParams.percentOverlap == .75
    taperFactor = 2;
elseif inputParams.percentOverlap == .9
    taperFactor = 5;
end
Q_inv = eye(windowIndices(2,1)-windowIndices(1,1)+1); % initialize AR covariance matrix as I

%% Iterate through data in segments
for n = 1:numSegments
    
    % Initialize parameters per segment
    Z = inputParams.Z_hat(:,:,n).*inputParams.windows(:,n); % window the design matrix
    dataSegment = detrendData(windowIndices(1,n):windowIndices(2,n)).*inputParams.windows(:,n); % window the data
    W_f = inputParams.W_freq(:,n);
    
    % Run cyclic descent algorithm on segment
    W = HRAN_estimatePriorCovariance(inputParams.TR,dataSegment,Z,W_f,N);
    [Beta_hat,alpha_hat,sigmaSquared_hat,~,Q_inv] = HRAN_CyclicDescentAlgorithm(Z,P,dataSegment,Q_inv,W);
    
    % If noise removed correctly...
    if ~isnan(Beta_hat)
        
        % Remove physiological noise from data
        deNoisedData(windowIndices(1,n):windowIndices(2,n)) = deNoisedData(windowIndices(1,n):windowIndices(2,n)) - Z(:,1+2+N:end)*Beta_hat(1+2+N:end)./taperFactor;
        
        % If additional outputs requested...
        if nout > 0
            neuralSignal(windowIndices(1,n):windowIndices(2,n)) = neuralSignal(windowIndices(1,n):windowIndices(2,n)) + Z(:,1+1+1:1+1+N)*Beta_hat(1+1+1:1+1+N)./taperFactor;
        end
        % varargout{2} -- estimated cardiac noise
        if nout > 1
            cardiacNoise(windowIndices(1,n):windowIndices(2,n)) = cardiacNoise(windowIndices(1,n):windowIndices(2,n)) + Z(:,1+1+N+1:1+1+N+2*C)*Beta_hat(1+1+N+1:1+1+N+2*C)./taperFactor;
        end
        % varargout{3} -- estimated respiratory noise
        if nout > 2
            respiratoryNoise(windowIndices(1,n):windowIndices(2,n)) = respiratoryNoise(windowIndices(1,n):windowIndices(2,n)) + Z(:,1+1+N+2*C+1:1+1+N+2*C+2*R)*Beta_hat(1+1+N+2*C+1:1+1+N+2*C+2*R)./taperFactor;
        end
        % varargout{4} -- residual (requires estimation of trend and AR)
        if nout > 3
            trend(windowIndices(1,n):windowIndices(2,n)) = trend(windowIndices(1,n):windowIndices(2,n)) + Z(:,1:2)*Beta_hat(1:2)./taperFactor;
        end
        
        % If algorithm fails, do not remove any noise
    else
        
        deNoisedData(windowIndices(1,n):windowIndices(2,n)) = detrendData(windowIndices(1,n):windowIndices(2,n));
        
        % If additional outputs requested...
        if nout > 0
            neuralSignal(windowIndices(1,n):windowIndices(2,n)) = zeros(length(windowIndices(1,n):windowIndices(2,n)),1);
        end
        if nout > 1
            cardiacNoise(windowIndices(1,n):windowIndices(2,n)) = zeros(length(windowIndices(1,n):windowIndices(2,n)),1);
        end
        if nout > 2
            respiratoryNoise(windowIndices(1,n):windowIndices(2,n)) = zeros(length(windowIndices(1,n):windowIndices(2,n)),1);
        end
        if nout > 3
            trend(windowIndices(1,n):windowIndices(2,n)) = zeros(length(windowIndices(1,n):windowIndices(2,n)),1);
        end
    end
    
    % If residuals requested...
    if nout > 3
        if ~isnan(alpha_hat)
            
            % Obtain estimate of AR noise
            ARestimate_seg = dataSegment - Z*Beta_hat;
            ARmodel = arima('Constant',0,'AR',alpha_hat(2:end),'Variance',sigmaSquared_hat);

            % Residual of AR model fit
            r = infer(ARmodel,ARestimate_seg);
            ARestimate(windowIndices(1,n):windowIndices(2,n)) = ARestimate(windowIndices(1,n):windowIndices(2,n)) + (ARestimate_seg-r)./taperFactor;
        
        else
            ARestimate(windowIndices(1,n):windowIndices(2,n)) = ARestimate(windowIndices(1,n):windowIndices(2,n)) + dataSegment;
            
        end
    end
    
end

%% Assign output variables
deNoisedData = deNoisedData + (data-detrendData); % add back in original trend

if nout > 0
    varargout{1} = neuralSignal;
end
if nout > 1
    varargout{2} = cardiacNoise;
end
if nout > 2
    varargout{3} = respiratoryNoise;
end
if nout > 3
    residual = detrendData-(cardiacNoise+respiratoryNoise+trend+neuralSignal+ARestimate);
    varargout{4} = residual;
end

end

