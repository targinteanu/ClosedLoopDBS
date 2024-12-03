function [svN, forBuff, stimBuff] = ...
    savebg_PhaseDetect_v1(ForStimSv, forBuff, stimBuff, svname, svN)
% 
% Save buffer data coming from the background 
% 

        % This relies on the forecast buffers being 2xN and the stim buffer
        % being 1xN; should be made more robust. 
        if isempty(ForStimSv)
            stimBuffSv = ForStimSv; forBuffSv = ForStimSv; 
        elseif width(ForStimSv) < 2
            stimBuffSv = ForStimSv; forBuffSv = [];
        elseif width(ForStimSv) > 2
            stimBuffSv = ForStimSv(:,3:end); forBuffSv = ForStimSv(:,1:2);
        else
            stimBuffSv = []; forBuffSv = ForStimSv;
        end

        if ~isempty(forBuffSv)
            PeakTrough = forBuffSv;
            save([svname,num2str(svN),'.mat'], 'PeakTrough');
            svN = svN+1;
            forBuff = [forBuffSv; forBuff];
        end
        if ~isempty(stimBuffSv)
            Stim = stimBuffSv;
            save([svname,num2str(svN),'.mat'], 'Stim');
            svN = svN+1;
            stimBuff = [stimBuffSv; stimBuff];
        end

end