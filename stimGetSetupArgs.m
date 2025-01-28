function StimSetupArgs = stimGetSetupArgs(UserArgs)
% 
% Input UserArgs can be e.g. GUIDE handles.
% The popup channel lists in UserArgs (i.e. UserArgs.pop_channel1, etc.)
% MUST have list entries in the format: '<Channel ID [int]>: <Channel Name>'
% 

    amp1 = eval(UserArgs.txt_amp1.String); 
    amp2 = eval(UserArgs.txt_amp2.String); 
    width1 = eval(UserArgs.txt_width1.String); 
    width2 = eval(UserArgs.txt_width2.String); 
    interphase = eval(UserArgs.txt_interphase.String); 
    frequency = 100; 
    pulses = 1;

    channel = nan(1,5); p = 1;
    for pop_ = [UserArgs.pop_channel1, ...
                UserArgs.pop_channel2, ...
                UserArgs.pop_channel3, ...
                UserArgs.pop_channel4, ...
                UserArgs.pop_channel5]
        chanind = pop_.Value;
        if chanind <= size(pop_.String,1)
            chPname = pop_.String{chanind};
            chPid = sscanf(chPname, '%d: '); 
            if numel(chPid)
                channel(p) = chPid(1);
            end
        end
        p = p+1;
    end
    channel1 = channel(1);
    channel2 = channel(2:end);

    StimSetupArgs = struct(...
        'channel1', channel1, ...
        'channel2', channel2, ...
        'amp1', amp1, ...
        'amp2', amp2, ...
        'width1', width1, ...
        'width2', width2, ...
        'interphase', interphase, ...
        'frequency', frequency, ...
        'pulses', pulses);

end