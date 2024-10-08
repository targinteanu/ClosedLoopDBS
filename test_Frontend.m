%% setup

minfac         = 1;    % this many (lo)cutoff-freq cycles in filter
min_filtorder  = 15;   % minimum filter length
srate = 1000; 
loco = 13; hico = 30; % Hz 
filtorder = minfac*fix(srate/loco);
TimeShiftFIR = filtorder/(2*srate); % seconds
IndShiftFIR = ceil(filtorder/2); % samples

dT = .001; % s between data requests 
PDSwin = 1000; % # samples ahead to forecast
buffSize = 20000; % samples

connect_cbmex(); 
t0 = datetime - seconds(cbmex('time'));
pause(10*dT);

try
    [rawH, rawT, rawB, rawN] = initRawData_cbmex([], buffSize);
catch ME
    warning(['Error on first attempt: ',ME.message]);
    pause(1);
    [rawH, rawT, rawB, rawN] = initRawData_cbmex([], buffSize);
end
fltH = initFilteredData(rawH, repmat(IndShiftFIR, size(rawH))); 
[forH, forT, forB] = initForecastData(fltH, repmat(PDSwin, size(fltH)));

fltT = cellfun(@(D) [1,0].*D, rawT, 'UniformOutput',false);
fltB = cellfun(@(D) [1,0].*D, rawB, 'UniformOutput',false);

timeBuffs = cell(size(rawH));
for ch = 1:length(rawH)
    t_ch = rawH{ch}(:,1);
    dt_ch = 1/rawN{ch}.SampleRate;
    for it = 2:length(t_ch)
        if isnan(t_ch(it))
            t_ch(it) = t_ch(it-1) + dt_ch;
        end
    end
    timeBuffs{ch} = t_ch;% + t0;
end
%forBuffs = cellfun(@(X) t0+seconds(nan(size(X,1),2)), timeBuffs, 'UniformOutput',false);
forBuffs = cellfun(@(X) (nan(size(X,1),2)), timeBuffs, 'UniformOutput',false);

chInd = 33;

forBuffs = forBuffs(chInd);

fig = figure; 

ax(2) = subplot(2,1,2); 
tPlt = timeBuffs{chInd}; tPltDisp = tPlt; td0 = tPltDisp(end); tDisp1 = tic;
pltTime = stem(t0 + seconds(tPlt), [nan; diff(tPlt)], '--');
grid on; hold on;
pltTimeDisp = stem(t0 + seconds(tPltDisp), [nan; diff(tPltDisp)]);
ylabel('Dur (s)');

% add plotting code below to be output from data2timetable? 
ax(1) = subplot(2,1,1);
rawPlt = data2timetable(rawB(chInd),rawN(chInd),t0); rawPlt = rawPlt{1};
fltPlt = data2timetable(fltB(chInd),rawN(chInd),t0); fltPlt = fltPlt{1};
forPlt = data2timetable(forB(chInd),rawN(chInd),t0); forPlt = forPlt{1};
pltRaw = plot(rawPlt.Time, rawPlt.Variables); 
hold on; grid on; 
xlabel('time'); ylabel(rawPlt.Properties.VariableNames{1});
pltFlt = plot(fltPlt.Time, fltPlt.Variables);
pltFor = plot(forPlt.Time, forPlt.Variables);

tPk = forBuffs{1}(:,1); tTr = forBuffs{1}(:,2);
pltPk = plot(t0 + seconds(tPk), zeros(size(tPk)), '^m');
pltTr = plot(t0 + seconds(tTr), zeros(size(tTr)), 'vm');

%linkaxes(ax, 'x')

%% start background parallel pool 
disconnect_cbmex();

% Start parallel pool if not already started
pool = gcp('nocreate');
if isempty(pool)
    pool = parpool;
end

% Set up a data queue to communicate between the background task and the main thread
dataQueue = parallel.pool.PollableDataQueue;

% Execute the data acquisition function asynchronously
f = parfeval(pool, @test_Backend, 0, dataQueue);  % 0 indicates no output

%% loop 
cont = isvalid(fig);
while cont 
    pause(dT)
    tDisp2 = td0 + toc(tDisp1); %tDisp1 = tic;

    [sentData, ok] = poll(dataQueue, .5);
    if ok
        if strcmpi(class(sentData), 'MException')
            getReport(sentData)
            cont = false;
        else
        tPltDisp = bufferData(tPltDisp, tDisp2);

        rawD1 = sentData(1,1); rawD4 = sentData(2,1);
        fltD1 = sentData(1,2); fltD4 = sentData(2,2);
        forD1 = sentData(1,3); forD4 = sentData(2,3);
        timeBuff = sentData{3,1}; forBuff = sentData{3,3};

        rawPlt = data2timetable(rawD4,rawD1,t0); rawPlt = rawPlt{1};
        fltPlt = data2timetable(fltD4,fltD1,t0); fltPlt = fltPlt{1};
        forPlt = data2timetable(forD4,forD1,t0); forPlt = forPlt{1};
        tPlt = timeBuff;
        tPk = forBuff(:,1); tTr = forBuff(:,2);
        tPltRng = [rawPlt.Time; fltPlt.Time; forPlt.Time];
        tPltRng = [min(tPltRng), max(tPltRng)];
        tPltRng = tPltRng + [-1,1]*.1*diff(tPltRng);

        cont = cont && isvalid(fig);
        if cont
            figure(fig);
            pltRaw.YData = rawPlt.Variables; pltRaw.XData = rawPlt.Time;
            pltFlt.YData = fltPlt.Variables; pltFlt.XData = fltPlt.Time;
            pltFor.YData = forPlt.Variables; pltFor.XData = forPlt.Time;
            % subplot(2,1,1); tPltRng = xlim();
            pltPk.XData = t0 + seconds(tPk); pltPk.YData = zeros(size(tPk));
            pltTr.XData = t0 + seconds(tTr); pltTr.YData = zeros(size(tTr));
            pltTime.YData = [nan; diff(tPlt)]; pltTime.XData = t0 + seconds(tPlt);
            pltTimeDisp.YData = [nan; diff(tPltDisp)]; pltTimeDisp.XData = t0 + seconds(tPltDisp);
            if sum(~isnat(tPltRng))
                subplot(2,1,1); xlim(tPltRng);
                subplot(2,1,2); xlim(tPltRng);
            end
        end
        end
    end
    cont = cont && isvalid(fig);
end

%% close 
cancel(f); % stop background parallel pool
disconnect_cbmex();
svname = [mfilename,'_output_',datestr(t0,'yyyymmdd_HHMMSS'),'.mat'];
PeakTrough = forBuffs{1};
save(svname,"PeakTrough");