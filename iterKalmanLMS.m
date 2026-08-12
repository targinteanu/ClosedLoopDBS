function [Y,P,w] = iterKalmanLMS(Y,artInd,a,P,q,w,stepsize,nLMS)
% 
% Use a combination Least Mean Squares (LMS) adaptive filter + Kalman
% filter technique to remove artifact from a signal. The LMS removes runs
% first, removing some artifact and providing an estimate of observation
% noise for the Kalman filter. Currently only supports single-channel
% output where state space is defined by a provided AR model. 
% 
% Inputs: 
%   Y: data with time-samples as columns; 
%      the first <max(N,M)> samples are assumed to already be artifact free
%   artInd: index of Y where artifact starts 
%   a: AR model coefficients, order M, ROW vector, s.t. a*Y ~ Ynext
%   P: latest estimate of process (state) covariance 
%   q: process noise variance (scalar)
%   w: latest weights of LMS filter, order N, COLUMN vector
%   stepsize: LMS filter learning rate 
%   nLMS: flag to normalize LMS filter weight update 
% 
% Outputs: 
%   Y: artifact-removed data 
%   P: updated estimate of process (state) covariance 
%   w: updated weights of LMS filter, order N, COLUMN vector
% 

L = length(Y); % signal duration
M = length(a); % AR model order / state x dimension
N = length(w); % LMS filter length (# taps)

% LMS setup 
noiseRef = zeros(1,L); % noise ref signal = timing of stim
noiseRef(artInd(1)) = 1;

for t = (max(N,M)+1):L

    %{
    if noiseRef(t+N-1)
        n = 1;
    else
        n = n+1;
    end
    %}
    y = Y(t,:)';
    xest = Y(t-(1:M));

    % LMS est observer noise 
    g = noiseRef(t-(1:N));
    noiseLMS = g*w; 
    e = y' - noiseLMS; % LMS "error" signal
    dw = g'*e; % LMS weight change
    if nLMS
        dw = dw./(g*g' + eps);
    end
    w = w + stepsize*dw;
    gprev = g; eprev = e;
    y = e'; % subtract LMS estimate of noise -- consider removing 
    R = diag(noiseLMS.^2); % Kalman observer noise

    % Kalman predict 
    Ppri = [[P(2:end,2:end), P(2:end,:)*a']; [a*P(:,2:end), a*P*a'+q]];
    xpri = [xest(2:end,:); a*xest];

    % Kalman update 
    S = Ppri(end,end) + R; 
    K = Ppri(:,end)*S^-1;
    xpos = xpri + K*(y-xpri(end));
    Ppos = eye(M); Ppos(:,end) = Ppos(:,end)-K; Ppos = Ppos*Ppri;

    % Kalman forward project and store 
    P = Ppos; 
    xest = xpos; 
    yest = xest(end);
    Y(t,:) = yest';

end

end