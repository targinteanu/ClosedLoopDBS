function [siginf,asiginf] = neuromodulation_output_visualization(filename, Hmag_Threshold)
    % This function parses and plots the data and stimulation output
    % from the `output/neuromod_output_YYYY-MM-DD_HH_MM_SS.csv` output
    % log produced by the neuromod_app.py application
    %
    % It runs a zero-phase moving average filter over the data to subract
    % a running DC-component, then it uses the same beta-band FIR filter
    % as in the FPGA to filter the centered-data and computes the hilbert
    % transform to get a baseline of expected 0-phase times.
    %
    % The stimulation timestamps are pulled and compared to the hilbert
    % phases to produce a polar histogram of the stimulation pulses versus
    % phase
    %
    % Returns:
    %   siginf - zero-mean and beta-band filtered version of the ADC input
    %   asiginf - the phase angle of siginf: angle(hilbert(siginf))

    Fs = 1000;

    if(nargin < 2)
        Hmag_Threshold = 100;
    end

    opts = detectImportOptions(filename);
    out = readmatrix(filename,opts);

    %% Signal Data Processing
    sigin = out(:,2);
    timein = out(:,1);
    sigin_zeromean = sigin - mean(sigin);

    % The beta-band filter we are using:
    h = [-0.00235822679806431,-0.00222172563171287,-0.00204188612994641,-0.00183038495484977,-0.00160034390923488,-0.00136581269659799,-0.00114119804128303,-0.000940660566584552,-0.000777502797457442,-0.000663572746103515,-0.00060870768315786,-0.000620241860603209,-0.000702600139262313,-0.000856996727262349,-0.00108125463821507,-0.00136975714725822,-0.00171353761034857,-0.00210050869551155,-0.00251582655346347,-0.00294237994222154,-0.00336138903478806,-0.00375309379670696,-0.00409750762578007,-0.00437520858421974,-0.00456813818108741,-0.00466037640254251,-0.00463886162140641,-0.00449402518340973,-0.00422031285552714,-0.00381656887412917,-0.00328626294222335,-0.00263754604641679,-0.00188312720632949,-0.00103997000959327,-0.000128814776282121,0.000826460826970075,0.00179962321217186,0.00276298068623895,0.0036883666800297,0.00454817311860771,0.00531639560056796,0.00596965038000671,0.0064881229614292,0.00685640947577641,0.00706421489728421,0.00710687651962473,0.00698568681140144,0.0067079966366596,0.00628708762755459,0.00574181096389299,0.00509599864561013,0.00437766221518625,0.00361800246321179,0.00285026159970675,0.00210845637912515,0.00142603643810278,0.000834516390501735,0.000362132821790785,3.25780897138549e-05,-0.000136138309429324,-0.000132653120598098,4.70052982061992e-05,0.000398925699844695,0.000911104654665073,0.00156353765180404,0.00232861857064927,0.00317184445611504,0.00405281127112782,0.00492647522528875,0.00574464379723933,0.00645765107504029,0.00701616390251436,0.00737305886804188,0.00748530568041385,0.00731579015430977,0.00683501001270845,0.00602257906226249,0.00486847998420728,0.00337401289754879,0.00155239580123291,-0.000571016280403731,-0.00295891447561934,-0.00556262940660626,-0.0083231287473506,-0.0111723991318278,-0.0140351789751382,-0.0168309960244307,-0.0194764519931052,-0.0218876868709724,-0.0239829478127192,-0.0256851821892519,-0.0269245716658481,-0.0276409241857489,-0.027785843533624,-0.0273246016808124,-0.0262376472292206,-0.0245216937382922,-0.022190344220844,-0.0192742222365224,-0.0158205953431237,-0.0118924926872715,-0.0075673346993105,-0.0029351086641472,0.00190386116090241,0.00684148672355834,0.0117647588792982,0.0165587765132936,0.0211098722437602,0.0253087409838873,0.0290534767491559,0.032252425249128,0.0348267649435504,0.0367127372169914,0.0378634568836062,0.0382502470358433,0.0378634568836062,0.0367127372169914,0.0348267649435504,0.032252425249128,0.0290534767491559,0.0253087409838873,0.0211098722437602,0.0165587765132936,0.0117647588792982,0.00684148672355834,0.00190386116090241,-0.0029351086641472,-0.0075673346993105,-0.0118924926872715,-0.0158205953431237,-0.0192742222365224,-0.022190344220844,-0.0245216937382922,-0.0262376472292206,-0.0273246016808124,-0.027785843533624,-0.0276409241857489,-0.0269245716658481,-0.0256851821892519,-0.0239829478127192,-0.0218876868709724,-0.0194764519931052,-0.0168309960244307,-0.0140351789751382,-0.0111723991318278,-0.0083231287473506,-0.00556262940660626,-0.00295891447561934,-0.000571016280403731,0.00155239580123291,0.00337401289754879,0.00486847998420728,0.00602257906226249,0.00683501001270845,0.00731579015430977,0.00748530568041385,0.00737305886804188,0.00701616390251436,0.00645765107504029,0.00574464379723933,0.00492647522528875,0.00405281127112782,0.00317184445611504,0.00232861857064927,0.00156353765180404,0.000911104654665073,0.000398925699844695,4.70052982061992e-05,-0.000132653120598098,-0.000136138309429324,3.25780897138549e-05,0.000362132821790785,0.000834516390501735,0.00142603643810278,0.00210845637912515,0.00285026159970675,0.00361800246321179,0.00437766221518625,0.00509599864561013,0.00574181096389299,0.00628708762755459,0.0067079966366596,0.00698568681140144,0.00710687651962473,0.00706421489728421,0.00685640947577641,0.0064881229614292,0.00596965038000671,0.00531639560056796,0.00454817311860771,0.0036883666800297,0.00276298068623895,0.00179962321217186,0.000826460826970075,-0.000128814776282121,-0.00103997000959327,-0.00188312720632949,-0.00263754604641679,-0.00328626294222335,-0.00381656887412917,-0.00422031285552714,-0.00449402518340973,-0.00463886162140641,-0.00466037640254251,-0.00456813818108741,-0.00437520858421974,-0.00409750762578007,-0.00375309379670696,-0.00336138903478806,-0.00294237994222154,-0.00251582655346347,-0.00210050869551155,-0.00171353761034857,-0.00136975714725822,-0.00108125463821507,-0.000856996727262349,-0.000702600139262313,-0.000620241860603209,-0.00060870768315786,-0.000663572746103515,-0.000777502797457442,-0.000940660566584552,-0.00114119804128303,-0.00136581269659799,-0.00160034390923488,-0.00183038495484977,-0.00204188612994641,-0.00222172563171287,-0.00235822679806431];

    siginf = filtfilt(h,1,(sigin - movmean(sigin,1024)));
    siginfh = hilbert(siginf);
    phase_angle_siginf = angle(siginfh);
    mag_siginf = abs(siginfh);
    mavg_mag_siginf = filtfilt(ones(64,1)/64,1,mag_siginf);

    %% Stimulation Time Handling
    stim_times = out(:,13);
    u_stim_times = unique(stim_times);
    % remove 0 entries
    u_stim_times = u_stim_times(u_stim_times ~= out(1,13));
    stims = zeros(length(siginf),1);
    stims(u_stim_times - timein(1) + 1) = 1;
    stims = stims(1:length(siginf)); % truncate incase we have predictions off the end of plot

    less_than_Hamp_Threshold_idx = find(mavg_mag_siginf < Hmag_Threshold);
    stim_high_power_only = stims; stim_high_power_only(less_than_Hamp_Threshold_idx) = 0;

    asiginf = angle(hilbert(siginf));

    %%
    times_secs = timein/Fs;
    figure;
    hold on;
    yyaxis("right");
    plot(times_secs, asiginf,'-', 'DisplayName','θ Beta-filtered Signal','Color',[1 .8 .6]);
    ylabel('Instantaneous Phase (radians)')
    ylim_right=get(gca,'YLim');
    ylim([-1*max(abs(ylim_right)) max(abs(ylim_right))]);
    yyaxis("left");
    plot(times_secs, sigin_zeromean, '-','DisplayName',' Input Signal (0-mean)','Color',[.75 .75 .75]);
    plot(times_secs,siginf,'b-','LineWidth',2,'DisplayName', 'Beta-filtered Signal');
    % xline(times(logical(stims)),'r-',[],'LineWidth',2);
    stem(times_secs, stims*(max(siginf)*1.1),'r-','LineWidth',2,'Marker','None','DisplayName','Calculated Stimulation');
    ylabel('Signal Value')
    xlabel('Time (s)');
    ylim_left=get(gca,'YLim');
    ylim([-1*max(abs(ylim_left)) max(abs(ylim_left))]);
    % plot(times_secs, unwrap(angle(hilbert(siginf))));
    % plot(times_secs, diff(unwrap(angle(hilbert(siginf)))));
    hold off;
    zoom xon;

    %%
    polar_edges = deg2rad(-5:+10:355);
    figure;polarhistogram(phase_angle_siginf(logical(stims)), polar_edges);
    %%%% OR Select just the part we xlim the plot to %%%%%
    % min_sample_idx = 512;
    % max_sample_idx = 1900;
    % xlim([min_sample_idx max_sample_idx]/Fs);
    % phase_angle_siginf_selected = phase_angle_siginf(min_sample_idx+1024:max_sample_idx+1024);
    % figure;polarhistogram(phase_angle_siginf_selected(logical(stims(min_sample_idx+1024:max_sample_idx+1024))), polar_edges);

    title('Polar Histogram of Phases at Calculated Stimulation Times (10° bins)');


    % [xc,lags] = xcorr(siginf.*siginf.^2,stims);
    % figure;plot(lags,xc);
end