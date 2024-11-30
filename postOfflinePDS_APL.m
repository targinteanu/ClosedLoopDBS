%% output CSV file 
[fn,fp] = uigetfile('*.csv'); 
neuromodulation_output_visualization(fullfile(fp,fn)); % APL internal eval 
tbl = readtable(fullfile(fp,fn));

%% interpret APL data
SamplingFreqAPL = 1000; % Hz
dataAPL = (tbl.data)';
tAPL = (tbl.dataTimestamp)'/SamplingFreqAPL; % s
StimIndAPL = (tbl.stimOut > 0);
phaseAPL = (tbl.projPhase)';

%% compare predicted phase, frequency with actual - APL data 
offline_PhaseDetect(dataAPL,[],SamplingFreqAPL,tAPL,'APL');