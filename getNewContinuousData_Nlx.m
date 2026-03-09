function [tX, chnum, fs] = getNewContinuousData_Nlx(chname, chtype)
fs = nan; chnum = nan;

if strcmp(chtype, 'CscAcqEnt')
    % Continuously Sampled Channel
    [succeeded, dataArray, timeStampArray, channelNumberArray, samplingFreqArray, numValidSamplesArray, numRecordsReturned, numRecordsDropped ] = NlxGetNewCSCData(chname);
    t = nan(size(dataArray)); 
    irec = cumsum([1, numValidSamplesArray(1:(end-1))]); % ?? starting index of data corresponding to each entry of time
    % TO DO: can check here whether data size matches numValidSamples
    if ~isempty(timeStampArray)
        t(irec) = double(timeStampArray)/1e6; % us -> s
    end
    if ~isempty(samplingFreqArray)
        fs = samplingFreqArray(end); % should something be done with the other entries ??
    end

elseif strcmp(chtype, 'SEScAcqEnt')
    % Single Electrode
    numSubChannels = 1;
    [succeeded, dataArray, timeStampArray, channelNumberArray, cellNumberArray, featureArray, numRecordsReturned, numRecordsDropped ] = NlxGetNewSEData(chname);
    timeStampArray = double(timeStampArray)/1e6; % us -> s
    [dataArray, t] = unInterleave(dataArray, timeStampArray, numSubChannels);
    fs = nan; % fix this!!

elseif strcmp(chtype, 'STScAcqEnt')
    % Stereotrode
    numSubChannels = 2; %number of stereotrode channels
    [succeeded, dataArray, timeStampArray, channelNumberArray, cellNumberArray, featureArray, numRecordsReturned, numRecordsDropped ] = NlxGetNewSTData(chname);
    timeStampArray = double(timeStampArray)/1e6; % us -> s
    [dataArray, t] = unInterleave(dataArray, timeStampArray, numSubChannels);
    fs = nan; % fix this!!

elseif strcmp(chtype, 'TTScAcqEnt')
    % Tetrode
    numSubChannels = 4; %number of subchannels in a tetrode
    [succeeded, dataArray, timeStampArray, channelNumberArray, cellNumberArray, featureArray, numRecordsReturned, numRecordsDropped ] = NlxGetNewTTData(chname);
    timeStampArray = double(timeStampArray)/1e6; % us -> s
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

if ~isempty(channelNumberArray)
    chnum = channelNumberArray(end); % should something be done with the other entries ??
end
tX = [t', double(dataArray)'];

    function [X,T] = unInterleave(x, t, n)
        L = length(t); % bufferSize
        W = length(x)/(L*n); % spikeSampleWindowSize
        X = reshape(x, [n,L*W]);
        T = nan(1,L*W);
        T(1:W:end) = t; % ??
    end

end