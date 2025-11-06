function [MI, P, bcent, fig1] = calcPAC(xPhi, xAmp, nbin, hax_PhaseAmpPlot)
%
% Get the Phase-Amplitude Coupling (PAC) of a signal x. The signal should
% already be filtered in the phase (xPhi) and amplitude (xAmp) ranges, and
% these should be aligned in time. Input nbin = number of phase bins;
% default 18 if empty or omitted.
% Optional fourth input enables output of the Phase-Amplitude plot. If an
% axis handle is passed, the plot will be put on that axis. If true is
% passed, a new figure will be created. Plot will be polar unless a
% cartesian axis handle is input. No plotting if input is empty, false, or
% omitted. 
% 
% Outputs: 
%   MI: modulation index 
%   P: normalized amplitude distribution 
%   bcent: phase bin centers (radians) corresponding to P
%   fig1: handle to figure if created; otherwise empty
% 
% based on https://doi.org/10.1152/jn.00106.2010
% 

if nargin < 4
    hax_PhaseAmpPlot = [];
    if nargin < 3
        nbin = [];
    end
end
if isempty(nbin)
    nbin = 18;
end
fig1 = [];
L = length(xPhi);
if L ~= length(xAmp)
    error('Signals need to be same length and aligned in time.')
end
if nbin < 1
    error('Wrong number of bins requested.')
end

% hilbert 
phi = angle(hilbert(xPhi));
amp = abs(hilbert(xAmp));

% focus on the central 80% of data to avoid edge effects
t1 = floor(.1*L); t2 = ceil(.9*L);
phi = phi(t1:t2); amp = amp(t1:t2);

[MI, P, bcent, fig1] = calcPAChelper(phi, amp, nbin, hax_PhaseAmpPlot);

end