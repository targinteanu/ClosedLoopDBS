function [emptyData, contData, buffData, chanInfo] = ...
    initRawData_cbmex(chsel, bufferSize)
% Initialize the multichannel raw data structure using BlackRock cbmex
% interface. 
%
% Inputs: chsel is selected channel ID numbers output from BlackRock.
% bufferSize is the size(s) of each corresponding channel buffer (samples). 
% 
% Output data structure is a cell array with columns to each channel. Each
% channel is represented by the matrix [time column, data column]. Time
% column contains time stamps in machine time seconds and nan elsewhere. 
% 
% chanInfo has fields: SampleRate, Name, Unit, IDnumber

%% handle/check inputs 
if length(bufferSize) < length(chsel)
    if length(bufferSize) == 1
        % assume the one input applies to all channels.
        bufferSize = repmat(bufferSize, size(chsel));
    else
        error('Incompatible input dimensions.')
    end
end

%% access cbmex 

[spikeEvents, time, continuousData] = cbmex('trialdata',1);

if isempty(continuousData)
    error('No continuous data; ensure that data acquisition has been enabled.')
end

%% setup 

chnum = [continuousData{:,1}]';
fs = [continuousData{:,2}]';
chname = spikeEvents(:,1);
chname = chname(chnum);

if isempty(chsel)
    % select all channels
    chsel = chnum;
end

emptyData = cell(1,length(chsel)); 
contData = emptyData; 
buffData = emptyData;
chanInfo = emptyData;

%% assign data to the structure 

for ch = 1:length(chsel)
    chInd = find(chnum == chsel(ch)); 

    if ~isempty(chInd)

        % Create raw data buffer of zeros of the correct length
        L = length(continuousData{chInd,3}); 
        emptyData{ch} = [nan(bufferSize(ch),1), zeros(bufferSize(ch),1)];
        %emptyData{ch}(1,1) = time - (bufferSize)/fs(chInd);
        emptyData{ch}(end,1) = time - 1/fs(chInd); % ??
        contData{ch} = [nan(L,1), continuousData{chInd,3}];
        contData{ch}(1,1) = time;

        % check units 
        config = cbmex('config', chInd);
        unitname_is = lower(config{11,1});
        if contains(unitname_is, 'unit')
            unitname = config{11,2};
        else
            unitname_is = contains(lower(config(:,1)), 'unit');
            unitname_is = find(unitname_is);
            unitname_is = unitname_is(1); 
            unitname = config{unitname_is,2};
        end

        % Channel Info
        ud.SampleRate = fs(chInd);
        ud.Name = chname(chInd);
        ud.Unit = unitname;
        ud.IDnumber = chnum(chInd);
        chanInfo{ch} = ud;

        buffData{ch} = bufferData(emptyData{ch}, contData{ch});

    end
    
end 

end