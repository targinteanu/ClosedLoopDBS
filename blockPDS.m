function [t2phi, i2phi, phi_inst, f_inst] = ...
    blockPDS(pastData, futureData, fs, phi, tmin, fmin, fmax)
% Determine the time (s) and # samples to next desired phase phi from a
% block of data sampled at a constant rate fs (Hz). Also return the current
% inst. phase phi_inst (rad) and frequency f_inst (Hz).
% Block data should include some length of pastData and
% (forecasted/predicted) futureData to minimize edge effects at the present
% timepoint, which is the last element of pastData. Data should be input as
% columns. 
% phi [desired] is in radians, i.e. phi=0 for peak, phi=pi for trough
% frequency will be clipped within range [fmin, fmax] (Hz) 

N = size(pastData,1); M = size(futureData,1);
blockData = [pastData; futureData];

[phi_block, f_block] = instPhaseFreq(blockData, fs);
phi_inst = phi_block(N,:);
f_block = max(f_block, fmin); 
f_block = min(f_block, fmax);
fwinlen = floor(.03*N); fwinlen = min(fwinlen, M);
fwin = N + ((-fwinlen):fwinlen);
f_inst = mean(f_block(fwin,:));
T=1/f_inst;

% time to next [desired] phi 
t2phi = zeros(size(phi)); i2phi = t2phi;
for p = 1:length(phi)
    phi_ = phi(p);
    t = (mod(phi_+2*pi-phi_inst,2*pi)./f_inst)/(2*pi); 

    % account for minimum delay time tmin 
    nT = (tmin-t)/T; % how many periods needed to add 
    t = t + ceil(nT)*T; 

    t2phi(p) = t;
    i2phi(p) = floor(fs*t2phi(p));
end