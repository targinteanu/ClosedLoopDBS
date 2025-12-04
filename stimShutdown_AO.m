function StimArgs = stimShutdown_AO(StimArgs, ~)

Results = AO_StopStimulation(-1);

if Results
    error(['Failed to shut down stimulator with error code ',num2str(Results)])
end

end