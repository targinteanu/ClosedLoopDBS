function stimulator = shutdownStimulator_cerestim(stimulator, ~)

        stimulator.stop();
        stimulator.disconnect;

end