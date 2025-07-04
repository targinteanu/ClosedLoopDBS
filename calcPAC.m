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

% bin based on phase 
bedge = linspace(-pi, pi, nbin+1);
bcent = bedge(2:end) + bedge(1:end-1); bcent = .5*bcent;
bind = discretize(phi, bedge); % could also use histcounts

% normalized mean amplitude
binamp = arrayfun(@(b) mean(amp(bind==b)), 1:nbin);
P = binamp/sum(binamp); % Normalized amplitude distribution
% Phase-Amplitude Plot (fig. 1e)
if ~isempty(hax_PhaseAmpPlot)
    cl = class(hax_PhaseAmpPlot);
    cl = lower(cl);
    if contains(cl, 'axis') || contains(cl, 'axes')
        % put the plot on the provided axis handle
        if contains(cl, 'polar')
            polarplot(hax_PhaseAmpPlot, [bcent, bcent(1)], [P, P(1)], 'LineWidth',1);
        else
            bar(hax_PhaseAmpPlot, bcent, P); 
            xlabel(hax_PhaseAmpPlot, 'Phase (radians)');
            ylabel(hax_PhaseAmpPlot, 'Normalized amplitude distribution');
        end
    elseif islogical(hax_PhaseAmpPlot) || isnumeric(hax_PhaseAmpPlot)
        if hax_PhaseAmpPlot
        % if 'true' was passed in, put the plot in a new figure
        fig1 = figure; polarplot([bcent, bcent(1)], [P, P(1)], 'LineWidth',1);
        title('Phase-Amplitude Plot');
        subtitle('Normalized amplitude distribution');
        end
    end
end

% modulation index 
H = -sum(P.*log(P)); % Shannon entropy 
DKL = log(nbin) - H; % KL distance compared to uniform distribution
MI = DKL / log(nbin); 

end