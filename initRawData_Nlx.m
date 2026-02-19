function [emptyData, contData, buffData, chanInfo, startTic] = ...
    initRawData_Nlx(chsel, bufferSize)
% Initialize the multichannel raw data structure using Neuralynx NetCom
% interface. 
%
% Inputs: chsel is selected channel ID numbers. 
% bufferSize is the size(s) of each corresponding channel buffer (samples). 
% 
% Output data structure is a cell array with columns to each channel. Each
% channel is represented by the matrix [time column, data column]. Time
% column contains time stamps in machine time seconds and nan elsewhere. 
% 
% chanInfo has fields: SampleRate, Name, Unit, IDnumber, (Type)

%% access NetCom 

% Nlx should have been connected already 
startTic = tic;
pause(1);

%get a list of all objects in the DAS, along with their types.
[succeeded, dasObjects, dasTypes] = NlxGetDASObjectsAndTypes;
if succeeded == 0
    error('FAILED get DAS objects and types');
else
    fprintf('Retrieved %d objects from the DAS\n', length(dasObjects));
end

%open up a stream for all objects that can stream date
objsconnected = false(size(dasObjects));
exclTypes = {'AcqSource', 'VTAcqEnt', 'EventAcqEnt'}; %EXCLUDE these obj types
for index = 1:length(dasObjects)
    %beginning in Cheetah 5.7.0 and Pegasus 2.0.0, the AcqSource data type
    %was included in the DAS object list. AcqSource objects cannot stream
    %data, but can be used to control the DAS
    if ~sum(strcmp(char(dasTypes(index)), exclTypes)) 
        succeeded = NlxOpenStream(dasObjects(index));
        if succeeded == 0
            warning('FAILED to open stream for %s', char(dasObjects(index)));
        else
            objsconnected(index) = true;
        end
    end
end;
if succeeded == 1
    fprintf('Streams opened for DAS objects\n');
end

dasObjects = dasObjects(objsconnected); 
dasTypes = dasTypes(objsconnected);

%% setup 

% init data streaming 
% -GetDASState requires Cheetah v5.7.0 or Pegasus v2.0.0 or newer. 
[succeeded, reply] = NlxSendCommand('-GetDASState');
if succeeded == 0
    fprintf('Failed to get DAS state\n');
else
    if strcmp(reply, 'Idle') == 1
        [succeeded, ~] = NlxSendCommand('-StartAcquisition');
        if succeeded == 0
            error('Failed to start acquisition');
        end
    end
end

% get initial data 
objsconnected = true(size(dasObjects));
chnum = nan(size(dasObjects)); 
Fs = nan(size(dasObjects));
continuousData = cell(size(dasObjects));
for ch = 1:length(dasObjects)
    chname = char(dasObjects{ch});
    chtype = char(dasTypes{ch});
    try
        [continuousData{ch}, chnum(ch), Fs(ch)] = getNewContinuousData_Nlx(chname, chtype);
    catch ME
        warning(['Could not init channel ',chname,': ',ME.message]);
        objsconnected(ch) = false;
    end
end
chnum = chnum(objsconnected);
Fs = Fs(objsconnected);
dasObjects = dasObjects(objsconnected); 
dasTypes = dasTypes(objsconnected);

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

%% assign data to structure 

for ch = 1:length(chnum)
    chnum_ch = chnum(ch);
    chInd = find(chsel == chnum_ch);

    if ~isempty(chInd)
        if length(chInd) > 1
            error('Non-unique channel ID(s).')
        end
        chname_ch = char(dasObjects{ch});
        chtype_ch = char(dasTypes{ch});

        % Create raw data buffer of zeros of the correct length
        time = min(continuousData{chInd(:,1)}); 
        emptyData{ch} = [nan(bufferSize(ch),1), zeros(bufferSize(ch),1)];
        if ~isempty(emptyData{ch})
            %emptyData{ch}(1,1) = time - (bufferSize)/fs(chInd);
            emptyData{ch}(end,1) = time - 1/fs(chInd); % ??
        end
        contData{ch} = continuousData{chInd};

        % Channel Info 
        ud.SampleRate = Fs(ch);
        ud.Name = chname_ch;
        ud.Unit = ''; % ??
        ud.Type = chtype_ch;
        ud.IDnumber = chnum_ch;
        ud.Resolution = 1; % ??
        ud.MinDigital = 0; % to negate offset calculation - ??
        ud.MinAnalog = 0;  % to negate offset calculation - ??
        chanInfo{chInd} = ud; 

        buffData{ch} = bufferData(emptyData{ch}, contData{ch});

    end

end


end