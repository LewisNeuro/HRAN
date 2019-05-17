function [Beta_hat,alpha_hat,sigmaSquared_hat,concentratedLikelihood,Q_inv] = HRAN_CyclicDescentAlgorithm(Z, P, y, Q_inv, W)
% HRAN_CYCLICDESCENTALGORITHM separates data into estimated signal + AR
% components
%
%% Input
%   Z: design matrix (time x #regressors)
%   P: AR order (e.g. 1-4)
%   y: data
%   Q_inv: AR covariance matrix (time x time)
%   W: estimated prior covariance (#regressors x #regressors)
%
%% Output
%   Beta_hat: optimal beta estimates
%   alpha_hat: estimated AR coefficients
%   sigmaSquared_hat: estimated AR sigmaSquared
%   concentratedLikelihood: likelihood of model fit (for given w_hr, w_rr)
%   Q_inv: updated AR covariance matrix (time x time)

N = size(Z,1); %number of time points
sigma_squared_w = []; %specific sigma squared for this w_hr and w_rr

%% 1) First iteration through algorithm

% First iteration through algorithm -- use Q_inv = I
Beta_hat = ((Z'*Q_inv*Z + W^-1)^-1)*Z'*Q_inv*y; % generalized least squares estimation for beta
%St = (y - Z*Beta_hat)'*Q_inv*(y - Z*Beta_hat); % mean squared error
V_hat = y - Z*Beta_hat; % AR series estimate

% Burg Algorithm to solve for alpha, residual variances and Cholesky
% decomposition (Q = U*E*U')
[alpha_hat,Eout,U_inv,Sigma_inv] = HRAN_arburg(V_hat, P);
sigma_squared_w = [sigma_squared_w Eout]; % store AR(P) residual variance to determine convergence

% Compute Q_inv efficiently (as in Krishnaswamy et al)
F = U_inv(1:P+1,1:P+1);
G = U_inv(P+1+1:end,1:P+1);
Hq = U_inv(P+1+1:end,P+1+1:end);
D = Sigma_inv(1:P+1,1:P+1);
c = 1/Eout;
Q_inv = [F'*D*F + c*(G'*G), c*G'*Hq; c*Hq'*G, c*(Hq'*Hq)];

%% 2) Second iteration through algorithm
Beta_hat = ((Z'*Q_inv*Z + W^-1)^-1)*Z'*Q_inv*y; % generalized least squares estimation for beta
St = (y - Z*Beta_hat)'*Q_inv*(y - Z*Beta_hat); % mean squared error
V_hat = y - Z*Beta_hat; % AR series estimate

% Burg Algorithm to solve for alpha, residual variances and Cholesky
% decomposition (Q = U*E*U')
[alpha_hat,Eout,U_inv,Sigma_inv] = HRAN_arburg(V_hat, P);
sigma_squared_w = [sigma_squared_w Eout]; % store AR(P) residual variance to determine convergence

% Compute Q_inv efficiently (as in Krishnaswamy et al)
F = U_inv(1:P+1,1:P+1);
G = U_inv(P+1+1:end,1:P+1);
Hq = U_inv(P+1+1:end,P+1+1:end);
D = Sigma_inv(1:P+1,1:P+1);
c = 1/Eout;
Q_inv = [F'*D*F + c*(G'*G), c*G'*Hq; c*Hq'*G, c*(Hq'*Hq)];

%% 3) Remainder of iterations in while loop (until convergence)
while abs(sigma_squared_w(end-1) - sigma_squared_w(end))/sigma_squared_w(end-1) > .0001
    Beta_hat = ((Z'*Q_inv*Z + W^-1)^-1)*Z'*Q_inv*y; % generalized least squares estimation for beta
    St = (y - Z*Beta_hat)'*Q_inv*(y - Z*Beta_hat); % mean squared error
    V_hat = y - Z*Beta_hat; % AR series estimate
    
    % Burg Algorithm to solve for alpha, residual variances and Cholesky
    % decomposition (Q = U*E*U')
    [alpha_hat,Eout,U_inv,Sigma_inv] = HRAN_arburg(V_hat, P);
    sigma_squared_w = [sigma_squared_w Eout]; % store AR(P) residual variance to determine convergence
    
    % Compute Q_inv efficiently (as in Krishnaswamy et al)
    F = U_inv(1:P+1,1:P+1);
    G = U_inv(P+1+1:end,1:P+1);
    Hq = U_inv(P+1+1:end,P+1+1:end);
    D = Sigma_inv(1:P+1,1:P+1);
    c = 1/Eout;
    Q_inv = [F'*D*F + c*(G'*G), c*G'*Hq; c*Hq'*G, c*(Hq'*Hq)];
    
    % If for some reason convergence not achieved, skip
    if length(sigma_squared_w) > 200
        Beta_hat = nan;
        sigmaSquared_hat = nan;
        concentratedLikelihood = nan;
        alpha_hat = nan;
        break
    end
end

% Compute concentrated likelihood
sigmaSquared_hat = sigma_squared_w(end);
[L,p] = chol(Q_inv);

% Check to make sure Q_inv positive definite
if p
    concentratedLikelihood = nan;
else
    logDetQ_inv = 2*sum(log(diag(L)));
    concentratedLikelihood = N * log(sigmaSquared_hat) - logDetQ_inv + St/sigmaSquared_hat;
end
end