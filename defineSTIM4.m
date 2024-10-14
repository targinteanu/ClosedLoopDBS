function stimulator = defineSTIM4(channel1, channel2, amp1, amp2, ...
        width1, width2, interphase, frequency, pulses)
% channel1 = Cathode; choose the correct channel combination on one of the headboxes 
% channel2 = Anode

%% connect 

stimulator = cerestim96();
%Check for stimulators
DeviceList = stimulator.scanForDevices();

if isempty(DeviceList)
    error('No CereStim found. Try restarting device.')
end

%Select a stimulator
stimulator.selectDevice(DeviceList(1));
%Connect to the stimulator
stimulator.connect;

%% check 

channel1 = channel1(~isnan(channel1)); channel2 = channel2(~isnan(channel2));
if isempty(channel1) || isempty(channel2)
    error('No channel selected.')
end
N1 = length(channel1); N2 = length(channel2);

    if (sum(sum(channel1 == channel2')) > 0) || (...
            (sum(sum(channel1 == channel1')) > N1) || ...
            (sum(sum(channel2 == channel2')) > N2) )
        proc = questdlg({'Selected Cathode and Annode channels are the same!',...
                         'Do you wish to proceed anyway?'}, ...
                        'Warning: same channel', ...
                        'Yes', 'No', 'No');
        if isempty(proc) || strcmp(proc, 'No')
            error('Retry with new stim parameters.')
        end
    end

[ampMin, ampMax] = stimulator.getMinMaxAmplitude();
if (amp1 < ampMin) || (amp2 < ampMin)
    warning('Amplitude is lower than minimum allowed value.')
end
if (amp1 > ampMax) || (amp2 > ampMax)
    warning('Amplitude is greater than maximum allowed value.')
end

if (sum(channel1 < 1) || sum(channel1 > 96)) || ...
        (sum(channel2 < 1) || sum(channel2 > 96))
    warning('Selected channel is not available for this device.')
end

%% define 

% psuedobipolar stimulation setting
%Program our waveform (stim pattern)
stimulator.setStimPattern('waveform',1,...%We can define multiple waveforms and distinguish them by ID
    'polarity',0,...%0=CF, 1=AF
    'pulses',pulses,...%Number of pulses in stim pattern
    'amp1',amp1/N1,...%Amplitude in uA
    'amp2',amp2/N1,...%Amplitude in uA
    'width1',width1,...%Width for first phase in us
    'width2',width2,...%Width for second phase in us
    'interphase',interphase,...%Time between phases in us
    'frequency',frequency);%Frequency determines time between biphasic pulses


stimulator.setStimPattern('waveform',2,...%We can define multiple waveforms and distinguish them by ID
    'polarity',1,...%0=CF, 1=AF %%% we use this
    'pulses',pulses,...%Number of pulses in stim pattern
    'amp1',amp1/N2,...%Amplitude in uA
    'amp2',amp2/N2,...%Amplitude in uA
    'width1',width1,...%Width for first phase in us
    'width2',width2,...%Width for second phase in us
    'interphase',interphase,...%Time between phases in us
    'frequency',frequency);%Frequency determines time between biphasic pulses
%Create a program sequence using any previously defined waveforms (we only
%have one)

stimulator.beginSequence();%Begin program definitio
stimulator.beginGroup();% What is until endGroup is payed simultaneously
for ch = channel1
    stimulator.autoStim(ch,1);%autoStim(Channel,waveformID)
end
for ch = channel2
    stimulator.autoStim(ch,2);
end
stimulator.endGroup();
% stimulator.wait(65) % ms;
stimulator.endSequence();%End program definition
% end of defintions 

%{
stimpara.pulses=pulses;
stimpara.amp1= amp1 ;  % 1500 - 1.5 ma 
stimpara.amp2 =amp2;
stimpara.width1 =width1;
stimpara.width2 =width2;
stimpara.interphase= interphase ;
stimpara.frequency= frequency ;

h.stimparama=stimpara;
%}
end