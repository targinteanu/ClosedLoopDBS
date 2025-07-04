function [emptyData, contData, buffData, chanInfo, startTic] = ...
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

%% access cbmex 

connect_cbmex(); startTic = tic;
pause(1);

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

for ch = 1:length(chsel)
    chInd = find(chnum == chsel(ch)); 

    if ~isempty(chInd)

        % Create raw data buffer of zeros of the correct length
        L = length(continuousData{chInd,3}); 
        emptyData{ch} = [nan(bufferSize(ch),1), zeros(bufferSize(ch),1)];
        if ~isempty(emptyData{ch})
            %emptyData{ch}(1,1) = time - (bufferSize)/fs(chInd);
            emptyData{ch}(end,1) = time - 1/fs(chInd); % ??
        end
        contData{ch} = [nan(L,1), continuousData{chInd,3}];
        contData{ch}(1,1) = time;

        % check units 
        config = cbmex('config', chInd);
        unitname_is = lower(config{10,1});
        if contains(unitname_is, 'unit')
            unitname = config{10,2};
        else
            unitname_is = contains(lower(config(:,1)), 'unit');
            unitname_is = find(unitname_is);
            unitname_is = unitname_is(1); 
            unitname = config{unitname_is,2};
        end

        % unit/scale conversion 
        searchRow = [11; 12; 13; 14]; 
        searchTerm = {'max', 'analog'; 
                      'max', 'digi'; 
                      'min', 'analog'; 
                      'min', 'digi'};
        searchResult = nan(size(searchRow)); 
        for s = 1:length(searchResult)
            searchRow(s) = min(searchRow(s), height(config));
            searchres_is = lower(config{searchRow(s),1});
            if contains(searchres_is, searchTerm{s,1}) && contains(searchres_is, searchTerm{s,2})
                searchResult(s) = config{searchRow(s),2};
            else
                searchres_is = contains(lower(config(:,1)), searchTerm{s,1});
                searchres_is = searchres_is & contains(lower(config(:,1)), searchTerm{s,2});
                searchres_is = find(searchres_is); 
                if ~isempty(searchres_is)
                    if length(searchres_is) > 1
                        warning(['Multiple config fields containing ',...
                            searchTerm{s,1},' and ',searchTerm{s,2}])
                    end
                    searchResult(s) = config{searchres_is(1),2};
                else
                    warning(['No config fields containing ',...
                        searchTerm{s,1},' and ',searchTerm{s,2}])
                end
            end
        end

        % Channel Info
        ud.SampleRate = fs(chInd);
        ud.Name = chname{chInd};
        ud.Unit = unitname;
        ud.IDnumber = chnum(chInd);
        ud.MinDigital = searchResult(4);
        ud.MinAnalog = searchResult(3);
        MaxDigital = searchResult(2); MaxAnalog = searchResult(1);
        ud.Resolution = (MaxAnalog - ud.MinAnalog) / (MaxDigital - ud.MinDigital);
        chanInfo{ch} = ud;

        buffData{ch} = bufferData(emptyData{ch}, contData{ch});

    end
    
end 

end