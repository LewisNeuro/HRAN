function W_freq = HRAN_estimateAliasedFreq(TR,Z,w_hr,w_rr,R,C,N)
% HRAN_ESTIMATEALIASEDFREQ estimates frequencies in Z (even if aliased)
% 
%% Input
%   TR: TR of fMRI data
%   Z: design matrix (time x #regressors)
%   w_hr: fundamental cardiac frequency (Hz)
%   w_rr: fundamental respiratory frequency (Hz)
%   R: respiratory model order
%   C: cardiac model order
%   N: number of neural regressors
%
%% Output
%   W_freq: frequencies of each regressor to speed up estimation of prior
%   covariance (# regressors x 1)


%% 0) Initialize W_freq
W_freq = zeros(size(Z,2),1); % initialize vector to store values
fs = 1/TR; % sampling frequency

%% 1) Iterate through cardiac frequencies and determine aliased frequency 
for z = 1:C
    
    % If less than Nyquist, then keep as is
    if z*w_hr < (fs/2)
        fAlias = z*w_hr;
        
    % Otherwise, aliased frequency computed as follows
    else
        fs_int = round(z*w_hr/fs);
        fAlias = abs((fs*fs_int) - z*w_hr);
    end
    
    % Store in W_freq
    W_freq(1+1+N+(z-1)*2+1) = fAlias;
    W_freq(1+1+N+(z-1)*2+2) = fAlias;
end

%% 2) Iterate through respiratory frequencies and determine aliased frequency 
for z = 1:R
    
    % If less than Nyquist, then keep as is
    if z*w_rr < (fs/2)
        fAlias = z*w_rr;
        
    % Otherwise, aliased frequency computed as follows
    else
        fs_int = round(z*w_rr/fs);
        fAlias = abs((fs*fs_int) - z*w_rr);
    end
    
    % Store in W_freq
    W_freq(1+1+N+2*C+(z-1)*2+1) = fAlias;
    W_freq(1+1+N+2*C+(z-1)*2+2) = fAlias;
end

%% 3) If any other terms specified (eg interaction terms), then iterate through
% design matrix and determine frequency with multitaper spectrum
for z = 1+1+N+2*C+2*R+1:size(Z,2)
    
    % Initialize parameters
    params.Fs = fs;
    params.tapers = [3 5];
    
    % Determine frequency
    Zvec = detrend(Z(:,z));
    [P_z, f]=mtspectrumc(Zvec, params);
    [~,zFreqIndex] = max(P_z);
    w_z = f(zFreqIndex);
    
    W_freq(z) = w_z;
     
end

end

