function [Z,Z_illConditioned] = HRAN_createDesignMatrix(t,neuralZ,w_hr,w_rr,R,C,N,X)
% HRAN_CREATEDESIGNMATRIX creates design matrices with physio noise
%
%% Input
%   t: vector of time points (e.g. [0:TR:length(data)-TR])
%   neuralZ: neural regressors (time x # regressors)
%   w_hr: cardiac frequency (Hz)
%   w_rr: respiratory frequency (Hz)
%   R: order of respiratory term
%   C: order of cardiac term
%   N: number of neural regressors (e.g. size(neuralZ,2))
%   X: order of interaction term
%
%% Output
%   Z: design matrix (time x # regressors)
%   Z_illConditioned: flag to check if matrix rank deficient (1 if
%   ill-conditioned)

%% 0) Initialize design matrix
% baseline + linear trend + neural regressors + cardiac terms + respiratory terms
% + interaction terms
num_beta_parameters = 1 + 1 + N + 2*C + 2*R + 4*X;
T = length(t); % Number of time points
Z = zeros(T,num_beta_parameters);

%% 1) Baseline trend (mu + linear)
Z(:,1) = ones(T,1); % mu_V
Z(:,2) = linspace(0,1,T); % linear term

%% 2) Neural regressors (demeaned, normalized)
for z_i = 1:N
    Z(:,z_i+1+1) = neuralZ(:,z_i) - mean(neuralZ(:,z_i)); % baseline
    Z(:,z_i+1+1) = Z(:,z_i+1+1)./max(Z(:,z_i+1+1)); % normalize
end

%% 3) Cardiac Terms
for z_i = 1:C
    z_col = (1 + 1 + N) + z_i*2 - 1; %assign harmonic to appropriate column
    Z(:,z_col) = cos(z_i*w_hr*2*pi*t'); %all the cos harmonics
    Z(:,z_col+1) = sin(z_i*w_hr*2*pi*t');%all the sin harmonics
end

%% 4) Respiratory Terms
for z_i = 1:R
    z_col = (1 + 1 + N + C * 2) + z_i*2 - 1; %assign harmonic to appropriate column
    Z(:,z_col) = cos(z_i*w_rr*2*pi*t'); %all the cos harmonics
    Z(:,z_col+1) = sin(z_i*w_rr*2*pi*t');%all the sin harmonics
end

%% 5) Interaction terms
for z_i = 1:X
    z_col = (1 + 1+ N + C * 2 + R * 2) + 1; %assign harmonic to appropriate column
    Z(:,z_col) = sin(z_i*w_hr*2*pi*t' + z_i*w_rr*2*pi*t');
    Z(:,z_col+1) = cos(z_i*w_hr*2*pi*t' + z_i*w_rr*2*pi*t');
    Z(:,z_col+2) = sin(z_i*w_hr*2*pi*t' - z_i*w_rr*2*pi*t');
    Z(:,z_col+3) = cos(z_i*w_hr*2*pi*t' - z_i*w_rr*2*pi*t');
end

%% 6) Output boolean indicating whether there is a problem with Z
Z_illConditioned = 0;

% If degenerate (eg some of the matrix columns are the same)
if rank(round(Z,3)) < num_beta_parameters
    Z_illConditioned = 1;
end

% Because of numerical errors, rank doesn't always work. So, check if the
% peak frequency of any columns is the same (and if so, ill conditioned).
% Find out what the neural frequency would be estimated as with given freq resolution
params.Fs = 1/(t(2)-t(1));
params.tapers = [3 5];
if size(Z,1) > 10
    [P,f] = mtspectrumc(Z(:,2:end),params);
    [~,maxInd] = max(P,[],1);
    if length(unique(maxInd(1:2:end))) < R+C+(N-1)+(X-2) || length(unique(maxInd(2:2:end))) < R+C+(N-1)+(X-2)
        Z_illConditioned = 1;
    end
else
    Z_illConditioned = 0;
end

end