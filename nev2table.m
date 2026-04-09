function T = nev2table(nevOrFilename)
% Extract data from blackrock nev files using the NPMK and output data as 
% evt table T. Input nevOrFilename can be the data structure obtained by
% running the NPMK (i.e. openNEV) or a filepath to a .ev or .mat file; If
% omitted, the user will be prompted to choose a .nev or .mat file
% graphically. Also output a flag packetLoss for whether packets needed to
% be stitched together. 

if nargin < 1
    nevOrFilename = [];
end
cl = class(nevOrFilename);

%% load file

if isstruct(nevOrFilename)
    NEV = nevOrFilename;
else
    if isempty(nevOrFilename)
        [fn,fp] = uigetfile({'*.nev'; '*.mat'});
        [~,fn,fe] = fileparts(fn);
    elseif ischar(nevOrFilename) || isstring(nevOrFilename)
        [fp,fn,fe] = fileparts(nevOrFilename);
    end
    if strcmpi(fe,'.mat')
        load(fullfile(fp,[fn,fe]), 'nev', 'NEV');
        for vtry = {'nev', 'NEV'}
            if exist(vtry{:}, 'var')
                NEV = eval(vtry{:});
            end
        end
    else
        NEV = openNEV(fullfile(fp,[fn,fe]));
    end
    clear vtry
end

%% electrode labels
lbl = {NEV.ElectrodesInfo.ElectrodeLabel};
lbl = upper(lbl);

% limit names to characters that will behave well as plot labels 
for l = 1:length(lbl)
    L = lbl{l}'; 
    L = L((L >= 33)&(L <= 126));
    L(L=='_') = ' ';
    lbl{l} = L;
end

% handle duplicate labels 
for l = 1:length(lbl)
    L = lbl{l};
    dups = strcmp(lbl, L);
    dups = find(dups);
    if length(dups) > 1
        warning(['Channel name ',L,' is duplicated.'])
        n = 1;
        for d = dups
            lbl{d} = [L,'(',num2str(n),')'];
            n = n+1;
        end
    end
end

%% interpret data from loaded file 
% obtain usable data variables and other information from the file. 

t0 = datetime(NEV.MetaTags.DateTime, 'TimeZone','UTC');

SamplingFreq = getFieldFromOpts(NEV.MetaTags, {'SampleRes', 'TimeRes'});
SamplingFreq = double(SamplingFreq);

datanames = fieldnames(NEV.Data);

% serial digital IO
Tsdio = [];
try
    SDIOname = 'SerialDigitalIO';
    SDIO = getfield(NEV.Data, SDIOname); 
    SDIOind = strcmp(datanames, SDIOname);
    if isfield(SDIO, 'TimeStampSec')
        SDIOtime = SDIO.TimeStampSec;
    else
        if SamplingFreq
            SDIOtime = double(SDIO.TimeStamp)/SamplingFreq;
        else
            error('Sampling rate could not be determined.')
        end
    end
    SDIOtime = SDIOtime';
    SDIOval = getFieldFromOpts(SDIO, {'UnparsedData', 'Value'});
    if isequal(size(SDIOval), size(SDIOtime))
        SDIOval = string(SDIOval);
    elseif numel(SDIOval)==1
        SDIOval = repmat("", size(SDIOtime));
    else
        error('Mismatched data and timestamps size.')
    end
    SDIOval = SDIOname+": "+SDIOval;
    SDIOtime = seconds(SDIOtime) + t0;
    if numel(SDIOtime)
        Tsdio = eventtable(SDIOtime, "EventLabels",SDIOval, "EventEnds",SDIOtime);
    end
    datanames = datanames(~SDIOind);
catch ME1
    warning(['Serial Digital IO not extracted due to error: ',ME1.message]);
end
clear SDIO SDIOtime SDIOval SDIOname SDIOind

% spikes 
Tspk = [];
try
    SPKname = 'Spikes';
    SPK = getfield(NEV.Data, SPKname);
    SPKind = strcmp(datanames, SPKname);
    if isfield(SPK, 'TimeStampSec')
        SPKtime = SPK.TimeStampSec;
    else
        if SamplingFreq
            SPKtime = double(SPK.TimeStamp)/SamplingFreq;
        else
            error('Sampling rate could not be determined.')
        end
    end
    SPKtime = SPKtime';
    SPKchan = SPK.Electrode; SPKchan = SPKchan';
    SPKlbl = arrayfun(@(ch) lbl{ch}, SPKchan, 'UniformOutput',false);
    SPKlbl = "Spike: "+string(SPKlbl);
    % currently ignoring waveform and unit fields!
    SPKtime = seconds(SPKtime) + t0;
    if numel(SPKtime)
        Tspk = eventtable(SPKtime, "EventLabels",SPKlbl, "EventEnds",SPKtime);
    end
    datanames = datanames(~SPKind);
catch ME2
    warning(['Spikes not extracted due to error: ',ME2.message]);
end
clear SPK SPKtime SPKchan SPKlbl SPKname SPKind

% comments 
Tcom = [];
try
    COMname = 'Comments';
    COM = getfield(NEV.Data, COMname);
    COMind = strcmp(datanames, COMname);
    if isfield(COM, 'TimeStampSec')
        COMtime2 = COM.TimeStampSec;
    else
        if SamplingFreq
            COMtime2 = double(COM.TimeStamp)/SamplingFreq;
        else
            error('Sampling rate could not be determined.')
        end
    end
    if isfield(COM, 'TimeStampStartedSec')
        COMtime1 = COM.TimeStampStartedSec;
    else
        if SamplingFreq
            COMtime1 = double(COM.TimeStampStarted)/SamplingFreq;
        else
            error('Sampling rate could not be determined.')
        end
    end
    COMval = COM.Text; COMval = string(COMval);
    COMtime1 = COMtime1'; COMtime2 = COMtime2'; 
    COMtime1 = seconds(COMtime1) + t0; 
    COMtime2 = seconds(COMtime2) + t0; 
    if numel(COMtime1)
        Tcom = eventtable(COMtime1, "EventLabels",COMval, "EventEnds",COMtime2);
    end
    datanames = datanames(~COMind);
catch ME3
    warning(['Comments not extracted due to error: ',ME3.message]);
end
clear COM COMtime1 COMtime2 COMval COMname COMind

% check remaining fields 
for dataname = datanames'
    Dempty = true;
    D = getfield(NEV.Data, dataname{:});
    if isstruct(D)
        Dnames = fieldnames(D);
        for Dname = Dnames'
            DD = getfield(D, Dname{:});
            Dempty = Dempty && isempty(DD);
        end
    end
    if ~Dempty
        warning(['field ',dataname{:},' has data that was not extracted!']);
    end
end

T = [Tsdio; Tspk; Tcom];

%% helpers 

    function val = getFieldFromOpts(S, opts)
        optsTried = 0; val = false;
        while optsTried < length(opts)
            optsTried = optsTried + 1;
            if isfield(S, opts{optsTried})
                val = getfield(S, opts{optsTried});
                return
            end
        end
    end

end