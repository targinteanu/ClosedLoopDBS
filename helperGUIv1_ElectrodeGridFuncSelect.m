function [showElecGrid, elecGridFunc] = helperGUIv1_ElectrodeGridFuncSelect(...
    selstr, FilterSetUp, locutoff, hicutoff)
% to be called when selection is made 

if nargin < 3
    locutoff = []; hicutoff = [];
    if nargin < 2
        FilterSetUp = [];
    end
end
if isempty(FilterSetUp)
    FilterSetUp = false;
end
if isempty(locutoff)
    locutoff = nan;
end
if isempty(hicutoff)
    hicutoff = nan;
end

showElecGrid = ~strcmp(selstr, 'None');
if showElecGrid
    if strcmp(selstr, 'PAC')
        % ***** TO DO: phase amplitude coupling: handles.elecGridFunc = ...
    else
        % __ band power
        if strcmp(selstr, 'Selected Band Power')
            if ~FilterSetUp
                error('Filter must be set for this selection.')
            end
            fbnd = [locutoff, hicutoff];
        elseif strcmp(selstr, 'Beta Power')
            fbnd = [13, 30]; % Hz
        elseif strcmp(selstr, 'Theta Power')
            fbnd = [4, 9]; % Hz
        elseif strcmp(selstr, 'Gamma Power')
            fbnd = [50, 200]; % Hz
        end
        elecGridFunc = @(data, fs) bandpower(data, fs, fbnd);
    end
else
    elecGridFunc = @(~,~) nan;
end

end