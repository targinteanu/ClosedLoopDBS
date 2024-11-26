function stimulator = setupStimulator_cerestim(UserArgs)

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
            channel(p) = str2double(pop_.String(chanind,:));
        end
        p = p+1;
    end
    channel1 = channel(1);
    channel2 = channel(2:end);

    stimulator = defineSTIM4(channel1, channel2, amp1, amp2, ...
        width1, width2, interphase, frequency, pulses);

end