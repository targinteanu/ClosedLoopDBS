%% params 

fbnd = [13,30]; % hz (beta)
gridsize = [21,3]; % 21x3 ECoG grid 

%% load data 
disp('locate .ns* file')
[fn,fp] = uigetfile('*.ns*');
NS = openNSx(fullfile(fp,fn), 'uV');
[~,fn,fe] = fileparts(fn);
NEV = openNEV(fullfile(fp,[fn,'.nev']));

%%
NStbl = ns2timetable(NS); NEVtbl = nev2table(NEV);
chlbls = NStbl.Properties.VariableNames;

Fs = NStbl.Properties.SampleRate;
if isnan(Fs)
    Fs = 1/median(seconds(diff(NStbl.Time)));
end

BPF = fir1(1023, fbnd/(Fs/2)); 
X = NStbl{:,1:prod(gridsize)}; 
tRel = seconds(NStbl.Time - NStbl.Time(1));
tReg = tRel(1):(1/Fs):tRel(end);

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
if ~isempty(NEVtbl)
    NEVlbl = NEVtbl.EventLabels;
    NEVlbl = (contains(NEVlbl, "SerialDigitalIO: "));
    NEVtbl = NEVtbl(NEVlbl,:); NEVlbl = NEVtbl.EventLabels;
    if ~isempty(NEVlbl)
        useBehav = true; % TO DO: allow user to select instead
        SrlVal = arrayfun(@(str) sscanf(str, 'SerialDigitalIO: %f'), NEVlbl);
        % find times between 255 and subsequent 253
        SrlStart = SrlVal == 253; SrlEnd = SrlVal == 255;
        SrlStart = find(SrlStart); SrlEnd = find(SrlEnd);
        SrlStart = SrlStart(1:2:end); % TO DO: start/end matching should be more robust
        tStart = NEVtbl.Time(SrlStart); tEnd = NEVtbl.Time(SrlEnd);
        tStartEnd = [tStart, tEnd];
        iOn = false(height(NStbl),1);
        for ti = 1:height(tStartEnd)
            iOn = iOn | ( NStbl.Time>=tStartEnd(ti,1) & (NStbl.Time<=tStartEnd(ti,2)) );
        end
        iOn = repmat(iOn,1,width(X));
    end
end

%%
if ~useBehav
    [iOn, pwrThresh] = midcross(pwr);
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
imagesc(pwrDiff); colorbar; title('Power Difference');
hold on;
[chlblY, chlblX] = meshgrid(1:gridsize(1), 1:gridsize(2));
chlblX = chlblX'; chlblY = chlblY';
chlblY = chlblY(:); chlblX = chlblX(:);
text(chlblX, chlblY, chlbls(1:prod(gridsize)), ...
    "HorizontalAlignment","center", "VerticalAlignment","middle");

%% helper(s)

function avg = maskedAvg(data, mask)
data(mask) = nan;
avg = median(data,1,'omitnan');
end