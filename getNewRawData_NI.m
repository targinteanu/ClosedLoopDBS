function [newContData, chname, chnum] = getNewRawData_NI(selChanInfo, DAQobj)

chsel = cellfun(@(ifo) ifo.IDnumber, selChanInfo);

[continuousData, time] = read(DAQobj, "all", OutputFormat="Matrix");

if isempty(continuousData)
    error('No continuous data; ensure that data acquisition has been enabled.')
end

chnum = 1:size(DAQobj.Channels, 2); % there are no numeric IDs, only indexes
%fs = DAQobj.Rate; fs = repmat(fs, size(chnum));
chname = {DAQobj.Channels.Name};
    % {DAQobj.Channels.Name} = Dev#_ai#
    % {DAQobj.Channels.ID} omits Dev#_

if isempty(chsel)
    % select all channels
    chsel = chnum;
end

newContData = cell(1,length(chsel));

if width(chsel) == width(chnum)
    % assume order is the same for timing 
    for ch = 1:length(chsel)
        newContData{ch} = ...
            [time, continuousData(:,ch)];
    end
else

% if not same size, can't be same order
for ch = 1:length(chsel)
    chInd = find(chnum == chsel(ch)); 
    if ~isempty(chInd)
        if length(chInd) > 1
            error(['Raw data ID ',num2str(chsel(ch)),' is not unique.'])
        end
        newContData{ch} = ...
            [time, continuousData(:,chInd)];
    end
end

end

end