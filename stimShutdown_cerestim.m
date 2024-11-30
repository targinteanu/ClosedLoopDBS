function stimulator = stimShutdown_cerestim(stimulator, ~)

        stimulator.stop();
        stimulator.disconnect;

end