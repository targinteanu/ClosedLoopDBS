%% params 

fbnd = [13,30]; % hz (beta)

%% load data 
disp('locate data file')
[fn,fp] = uigetfile('*.*');
[~,fn,fe] = fileparts(fn);
EVtbl = [];
if contains(lower(fe), '.ns')
    NS = openNSx(fullfile(fp,[fn,fe]), 'uV');
    NEV = openNEV(fullfile(fp,[fn,'.nev']));
    gridsize = [21,3]; % 21x3 ECoG grid 
    Tbl = ns2timetable(NS); EVtbl = nev2table(NEV);
elseif contains(lower(fe),'.mat')
    load(fullfile(fp,[fn,fe]), 'tblsBaseline');
    Tbl = tblsBaseline{1}; EVtbl = Tbl.Properties.Events;
    gridsize = [21,3]; % 21x3 ECoG grid 
elseif contains(lower(fe),'.xls') || strcmpi(fe, '.csv')
    opts = detectImportOptions(fullfile(fp,[fn,fe]));
    out = readmatrix(fullfile(fp,[fn,fe]),opts);
    sigin = out(:,2);
    timein = out(:,1);
    Fs = 1000; % Hz
    times_secs = timein/Fs;
    gridsize = [1,1];
    Tbl = array2timetable(sigin, "RowTimes",seconds(times_secs));
else
    error(['Unrecognized file type: ',fe])
end

%%
chlbls = Tbl.Properties.VariableNames;

Fs = Tbl.Properties.SampleRate;
if isnan(Fs)
    Fs = 1/median(seconds(diff(Tbl.Time)));
end

BPF = fir1(1023, fbnd/(Fs/2)); 
X = Tbl.Variables; 
if width(X) > prod(gridsize)
    warning("Limiting to first "+prod(gridsize)+" channels.")
    X = X(:,1:prod(gridsize));
elseif width(X) < prod(gridsize)
    warning("Attempting to fill missing channels.")
    chlblnum = cellfun(@(str) str((str>=48)&(str<=57)), chlbls, 'UniformOutput',false);
    chlblnum = cellfun(@str2num, chlblnum);
    X_ = zeros(height(X), prod(gridsize));
    X_(:, chlblnum) = X;
    X = X_;
    chlbls_ = convertStringsToChars(repmat("x", 1, width(X)));
    chlbls_(chlblnum) = chlbls;
    chlbls = chlbls_;
end
tRel = seconds(Tbl.Time - Tbl.Time(1));
tReg = tRel(1):(1/Fs):tRel(end); tReg = tReg';

% regularize 
X = interp1(tRel,X,tReg, "nearest","extrap");
% plug NaNs 
inan = isnan(X); 
if any(inan)
    tNan = tReg(~inan); xNan = X(~inan);
    X = interp1(tNan,xNan,tReg, "nearest","extrap");
end

Xf = filtfilt(BPF,1,X); Xe = envelope(Xf);
pwr = movmean(Xe, 1000); 

%%
useBehav = false;
if ~isempty(EVtbl)
    NEVlbl = EVtbl.EventLabels;
    NEVlbl = (contains(NEVlbl, "SerialDigitalIO: "));
    EVtbl = EVtbl(NEVlbl,:); NEVlbl = EVtbl.EventLabels;
    if ~isempty(NEVlbl)
        useBehav = questdlg('Define On/Off Using:', 'Behavior/Ephys Selection', 'behavior', 'ephys', 'behavior');
        useBehav = strcmp(useBehav, 'behavior');
        SrlVal = arrayfun(@(str) sscanf(str, 'SerialDigitalIO: %f'), NEVlbl);
        % find times between 255 and subsequent 253
        SrlStart = SrlVal == 255; SrlEnd = SrlVal == 253;
        SrlStart = find(SrlStart); SrlEnd = find(SrlEnd);
        SrlEnd = SrlEnd(2:2:end); %SrlStart = SrlStart(1:(end-1)); % TO DO: start/end matching should be more robust
        tStart = EVtbl.Time(SrlStart); tEnd = EVtbl.Time(SrlEnd);
        tStartEnd = [tStart, tEnd];
        iOn = false(height(Tbl),1);
        for ti = 1:height(tStartEnd)
            iOn = iOn | ( Tbl.Time>=tStartEnd(ti,1) & (Tbl.Time<=tStartEnd(ti,2)) );
        end
        iOn = repmat(iOn,1,width(X));
    end
end

%% view grid of all chans 
if ~useBehav
    iOn = false(size(pwr)); pwrThresh = zeros(1,width(pwr));
    for ch = 1:width(pwr)
        [iOn_, pwrThresh(ch)] = midcross(pwr(:,ch));
        %iOn_ = round(iOn_);
        %iOn(iOn_,ch) = true;
        iOn(:,ch) = pwr(:,ch) >= pwrThresh(ch);
    end
end
pwrOn = maskedAvg(pwr, iOn);
pwrOff = maskedAvg(pwr, ~iOn);
if useBehav
    pwrThresh = 0.5*(pwrOn + pwrOff);
end
pwrOn = reshape(pwrOn, gridsize); 
pwrOff = reshape(pwrOff, gridsize);
pwrDiff = pwrOn - pwrOff; % pwrDiff = pwrDiff./pwrOn;
figure('Units','normalized', 'Position',[.1,.1,.8,.8]); 
subplot(1,3,1); imagesc(pwrOn); colorbar; title('ON Power');
subplot(1,3,2); imagesc(pwrOff); colorbar; title('OFF Power');
subplot(1,3,3);
imagesc(pwrDiff); colorbar; title('Power Difference (On - Off)');
hold on;
[chlblY, chlblX] = meshgrid(1:gridsize(1), 1:gridsize(2));
chlblX = chlblX'; chlblY = chlblY';
chlblY = chlblY(:); chlblX = chlblX(:);
text(chlblX, chlblY, chlbls(1:prod(gridsize)), ...
    "HorizontalAlignment","center", "VerticalAlignment","middle");

%% view channel(s) spectrogram

chsel = listdlg("PromptString","Select Channel(s) to inspect", "SelectionMode","multiple", ...
    "ListString",chlbls);

for ch = chsel
    figure('Units','normalized', 'Position',[.1,.1,.8,.8])
    x = X(:,ch); xf = Xf(:,ch); Th = pwrThresh(ch); 
    xname = chlbls{ch};
    iOn_ = iOn(:,ch);

    % extract start/end indexes 
    iStart = diff(iOn_) > 0; iEnd = diff(iOn_) < 0;
    iStart = find(iStart); iEnd = find(iEnd);

    % plot unadjusted spectrogram
    [S,fS,tS] = spectrogram(x,1*Fs,[],[],Fs,"yaxis","power");
    ax(1) = subplot(3,1,1); 
    img = imagesc(tS, fS(2:end), 20*log10(abs(S)));
    img.Parent.YDir = 'normal';
    title([xname,' spectrogram (dB)']);
    ylabel('Frequency (Hz)'); xlabel('time (s)');

    % correct pink noise
    [~,k1,c2] = pinkcorrect(mean(abs(S),2),fS);
    Anoise = k1*fS.^c2; Anoise(1)=eps;
    SS = abs(S)./Anoise;
    
    % saturate out outliers for better display
    %SSall = log(SS(:)+eps);
    SSall = SS(:);
    [~,~,OLthresh] = isoutlier(SSall, 'median', 'ThresholdFactor',10);
    OLthresh = max(SSall(SSall<OLthresh));
    %OLthresh = exp(OLthresh);    

    % plot adjusted spectrogram 
    ax(2) = subplot(3,1,2);
    img = imagesc(tS, fS(2:end), (SS(2:end,:)), [0,OLthresh]); %colorbar
    img.Parent.YDir = 'normal';
    title([xname,' adjusted spectrogram']);
    ylabel('Frequency (Hz)'); xlabel('time (s)');

    % plot signal 
    ax(3) = subplot(3,1,3);
    yyaxis("left"); plot(tReg, x); 
    yyaxis("right"); plot(tReg, xf, 'b'); hold on;
    plot([tReg(1), tReg(end)], Th*ones(1,2), ':k');
    plot([tReg(1), tReg(end)], -Th*ones(1,2), ':k');
    yl = ylim(); yl = .75*yl;
    stem(tReg(iStart), yl(1)*ones(size(iStart)), '.g');
    stem(tReg(iStart), yl(2)*ones(size(iStart)), '.g');
    stem(tReg(iEnd), yl(1)*ones(size(iEnd)), '.r');
    stem(tReg(iEnd), yl(2)*ones(size(iEnd)), '.r');

    linkaxes(ax, 'x');

end

%% helper(s)

function avg = maskedAvg(data, mask)
data(~mask) = nan;
avg = median(data,1,'omitnan');
end

function [A, k1, c2] = pinkcorrect(A,f)
% correct for noise that obeys Anoise = k1*f^c2
% i.e. ln(Anoise) = c2*ln(f) + c2*ln(k1)
if f(1) < 2*eps
    f0 = 0; f = f(2:end);
    A0 = A(1,:); A = A(2:end,:);
else
    f0 = zeros(0,width(f)); A0 = zeros(0,width(A));
end
lnA = log(A); lnf = log(f); F = [ones(size(lnf)), lnf];
c = F\lnA; 
% c1 = c2*ln(k1), i.e. k1 = exp(c1/c2)
c2 = c(2); k1 = exp(c(1)/c(2));
lnAnoise = F*c;
lnA = lnA - lnAnoise; A = exp(lnA);
A = [A0; A];
end