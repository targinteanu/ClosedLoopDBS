function [cpod_device, cerestim_stimulator] = ...
    stimShutdown_cerestim_cpod(cpod_device, cerestim_stimulator)
% shut down cerestim AND cpod for cerestim in trigger mode

if nargin < 2
    cerestim_stimulator = [];
    if nargin < 1
        cpod_device = [];
    end
end

if ~isempty(cpod_device)
    cpod_device = stimShutdown_cpod(cpod_device);
end

if ~isempty(cerestim_stimulator)
    cerestim_stimulator = stimShutdown_cerestim(cerestim_stimulator);
end

end