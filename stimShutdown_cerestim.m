function stimulator = stimShutdown_cerestim(stimulator, ~)

        stimulator.stop();
        stimulator.disconnect;
        delete(stimulator);
        % unclear if this always works 
        % try searching for any connected cerestims and disconnecting them

end