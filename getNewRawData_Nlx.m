function [newContData, chname, chnum] = getNewRawData_Nlx(selChanInfo)

W = length(selChanInfo);
newContData = cell(1,W);
for ch = 1:W
    ifo = selChanInfo{ch};
    try
        newContData{ch} = getNewContinuousData_Nlx(ifo.Name, ifo.Type);
    catch ME
        warning(['Could not update channel ',ifo.Name,': ',ME.message]);
    end
end

chname = cellfun(@(ifo) ifo.Name, selChanInfo, 'UniformOutput',false);
chnum = cellfun(@(ifo) ifo.IDnumber, selChanInfo);

end