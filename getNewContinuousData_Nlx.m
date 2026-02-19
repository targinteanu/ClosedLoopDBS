function [tX, chnum, fs] = getNewContinuousData_Nlx(chname, chtype)

if strcmp(chtype, 'CscAcqEnt')
    % Continuously Sampled Channel
    [succeeded, dataArray, timeStampArray, channelNumberArray, samplingFreqArray, numValidSamplesArray, numRecordsReturned, numRecordsDropped ] = NlxGetNewCSCData(chname);
    t = nan(size(dataArray)); 
    irec = cumsum([1, numValidSamplesArray]); % ?? starting index of data corresponding to each entry of time
    t(irec) = timeStampArray;
    fs = samplingFreqArray(end); % should something be done with the other entries ??

elseif strcmp(chtype, 'SEScAcqEnt')
    % Single Electrode
    numSubChannels = 1;
    [succeeded, dataArray, timeStampArray, channelNumberArray, cellNumberArray, featureArray, numRecordsReturned, numRecordsDropped ] = NlxGetNewSEData(chname);
    [dataArray, t] = unInterleave(dataArray, timeStampArray, numSubChannels);
    fs = nan; % fix this!!

elseif strcmp(chtype, 'STScAcqEnt')
    % Stereotrode
    numSubChannels = 2; %number of stereotrode channels
    [succeeded, dataArray, timeStampArray, channelNumberArray, cellNumberArray, featureArray, numRecordsReturned, numRecordsDropped ] = NlxGetNewSTData(chname);
    [dataArray, t] = unInterleave(dataArray, timeStampArray, numSubChannels);
    fs = nan; % fix this!!

elseif strcmp(chtype, 'TTScAcqEnt')
    % Tetrode
    numSubChannels = 4; %number of subchannels in a tetrode
    [succeeded, dataArray, timeStampArray, channelNumberArray, cellNumberArray, featureArray, numRecordsReturned, numRecordsDropped ] = NlxGetNewTTData(chname);
    [dataArray, t] = unInterleave(dataArray, timeStampArray, numSubChannels);
    fs = nan; % fix this!!

else
    error(['Continuous data from channel type ',chtype,' is not supported.'])
end

if ~succeeded
    error(['Nlx Continuous DAQ failed for object ',chname,'.'])
end
if numRecordsDropped > 0
    warning(['Nlx dropped ',num2str(numRecordsDropped),' records since last DAQ.']);
end

chnum = channelNumberArray(end); % should something be done with the other entries ??
tX = [t', dataArray];

    function [X,T] = unInterleave(x, t, n)
        L = length(t); % bufferSize
        W = length(x)/(L*n); % spikeSampleWindowSize
        X = reshape(x, [n,L*W]);
        T = nan(1,L*W);
        T(1:W:end) = t; % ??
    end

end