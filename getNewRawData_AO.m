function [newContData, chname, chnum] = getNewRawData_AO(chsel)
% chsel (IDs, NOT indexes) must not be empty; it will set the channels in
% order 

W = length(chsel);
[Results,continuousData,DataCapture,time] = AO_GetAlignedData(chsel);
if Results
    msg = ['Failed to acquire data with error code ',num2str(Results)];
    error(msg); % consider trying again until success 
end
continuousData = continuousData(1:DataCapture); 
L = DataCapture/W;
continuousData = reshape(continuousData, ...
    L, W); % columns = channels 

newContData = cell(1,W);
for ch = 1:W
    contData_ch = continuousData(:,ch);
    newContData{ch} = [nan(size(contData_ch)), contData_ch];
    newContData{ch}(1,1) = time;
end

chname = {}; chnum = chsel;
% if function is used correctly, these outputs should not matter.
% Otherwise, the following code will need to be implemented, plus a search
% for which chnum IDs are in chsel, which might cause delays 
%{
[Results, channelsData] = AO_GetAllChannels();
if Results
    error(['Failed to obtain channel info with error code ',num2str(Results)])
end
chnum = [channelsData.channelID];
chname = {channelsData.channelName};
%}

end