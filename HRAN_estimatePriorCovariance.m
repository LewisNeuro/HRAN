function W = HRAN_estimatePriorCovariance(TR,y,Z,W_f,N)
% HRAN_ESTIMATEPRIORCOVARIANCE estimates prior covariance of design matrix
% regressors
%
%% Input
%   TR: TR of fMRI data (s)
%   y: data
%   Z: design matrix
%   W_f: frequencies corresponding to design matrix
%   N: number of neural regressors
% 
%% Output 
%   W: prior covariance of each regressor
% 

%% 1) If there are neural regressors, then use glm to approximate signal without them
params.Fs = 1/TR;
params.tapers = [2 3];

if N
    [~,~,stats] = glmfit(Z(:,1+1+1:1+1+N),y-mean(y));
    [P_y, f]=mtspectrumc(stats.resid, params);
else
    [P_y, f]=mtspectrumc(y-mean(y), params);
end

%% 2) Compute moving mean
mv_k= find(f>.25,1); %.25
P_movmean = movmean(P_y,mv_k,'Endpoints','shrink');
[P_y, f]=mtspectrumc(y-mean(y), params);

%% 3) Determine prior covariance at each frequency

% Initialize W
W = zeros(size(Z,2));

% Baseline term (assume high prior covariance)
W(1,1) = 10^3;
W(2,2) = 10^3;

% Neural terms (assume high prior covariance)
for h = 1:N
    W(1+1+h,1+1+h) = 10^3;
end

% Physiological terms (determine appropriate weighting)
for z = 1+1+N+1:length(W_f)
    
    % Find frequency band around the one in design matrix
    w_y = [find(f<W_f(z),1,'last') find(f>W_f(z),1,'first')];
    
    % Find max difference between the two
    W_z = max(max(P_y(w_y)-P_movmean(w_y)));
    
    % If negative, then initialize as very low
    if W_z < 0
        W_z = 10^-3;
    end
    
    % Store in W
    W(z,z) = W_z;
end

end

