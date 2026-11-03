function [emptyData, contData, buffData, chanInfo, startTic] = ...
    initRawData_NI(chsel, bufferSize, DAQobj)
% Initialize the multichannel raw data structure using NI DAQ
% interface. 
%
% Inputs: chsel is selected channel "ID" (index) numbers from NI / DAQ toolbox.
% bufferSize is the size(s) of each corresponding channel buffer (samples). 
% 
% Output data structure is a cell array with columns to each channel. Each
% channel is represented by the matrix [time column, data column]. Time
% column contains time stamps in machine time seconds and nan elsewhere. 
% 
% chanInfo has fields: SampleRate, Name, Unit, IDnumber

%% access DAQ

startTic = tic; 
pause(1);

[continuousData, time] = read(DAQobj, "all", OutputFormat="Matrix");

if isempty(continuousData)
    error('No continuous data; ensure that data acquisition has been enabled.')
end

%% setup 

chnum = 1:size(DAQobj.Channels, 2); % there are no numeric IDs, only indexes
fs = DAQobj.Rate; fs = repmat(fs, size(chnum));
chname = {DAQobj.Channels.Name};
    % {DAQobj.Channels.Name} = Dev#_ai#
    % {DAQobj.Channels.ID} omits Dev#_

if isempty(chsel)
    % select all channels
    chsel = chnum;
end

emptyData = cell(1,length(chsel)); 
contData = emptyData; 
buffData = emptyData;
chanInfo = emptyData;

%% handle/check inputs 
if length(bufferSize) < length(chsel)
    if length(bufferSize) == 1
        % assume the one input applies to all channels.
        bufferSize = repmat(bufferSize, size(chsel));
    else
        error('Incompatible input dimensions.')
    end
end

%% assign data to the structure 

% Does this loop work with removing channels? See differences in
% initRawData_AO
for ch = 1:length(chsel)
    chInd = find(chnum == chsel(ch)); 

    if ~isempty(chInd)

        % Create raw data buffer of zeros of the correct length
        L = size(continuousData,1); 
        emptyData{ch} = [nan(bufferSize(ch),1), zeros(bufferSize(ch),1)];
        if ~isempty(emptyData{ch})
            %emptyData{ch}(1,1) = time(end) - (bufferSize)/fs(chInd);
            emptyData{ch}(end,1) = time(end) - 1/fs(chInd); % ??
        end
        contData{ch} = [time, continuousData(:,chInd)];

        % check units, scale
        unitname = DAQobj.Channels(chInd).Range.Units;
        rangeMin = DAQobj.Channels(chInd).Range.Min;
        rangeMax = DAQobj.Channels(chInd).Range.Max;

        % Channel Info
        ud.SampleRate = fs(chInd);
        ud.Name = chname{chInd};
        ud.Unit = unitname;
        ud.IDnumber = chnum(chInd);
        ud.MinDigital = rangeMin; 
        ud.MinAnalog = rangeMin;
        ud.Resolution = 0.5;
        chanInfo{ch} = ud;

        buffData{ch} = bufferData(emptyData{ch}, contData{ch});

    end

end