function [Z_hat,varargout] = HRAN_recreateDesignMatrix(TR,w_rr_hat,w_hr_hat,dataIndices,predictedBoldResponse,R_out,C_out,N_out,X)
%% Recreate design matrix with different number of harmonics

% Initialize time
T = 0:TR:(length(predictedBoldResponse)-1)*TR; %initialize time array

% Initialize window length
windowLength = diff(dataIndices);
windowLength = windowLength(1)+1;

% Initialize number of segments
numSegments = size(dataIndices,2);

% Initialize number of parameters
num_beta_parameters = 1+1+N_out+2*C_out+2*R_out+4*X; % number of parameters in design matrix

% Initialize Z_hat
Z_hat = zeros(windowLength,num_beta_parameters,numSegments); % optimal design matrix per segment

% Create design matrices
for n = 1:numSegments
    [Z, ~] = HRAN_createDesignMatrix(T(dataIndices(1,n):dataIndices(2,n)),predictedBoldResponse(dataIndices(1,n):dataIndices(2,n),:),w_hr_hat(n),w_rr_hat(n),R_out,C_out,N_out,X);
    Z_hat(:,:,n) = Z;
end

%% Respecify W_freq appropriately

% If W_freq requested...
if nargout > 1
    
    % Iterate through Z_hat and estimate frequencies
    W_freq = zeros(num_beta_parameters,numSegments); % optimal design matrix per segment
    for n = 1:numSegments
        W_f = HRAN_estimateAliasedFreq(TR,Z,w_hr_hat(n),w_rr_hat(n),R_out,C_out,N_out);
        W_freq(:,n) = W_f;
    end
    
    varargout{1} = W_freq;
end