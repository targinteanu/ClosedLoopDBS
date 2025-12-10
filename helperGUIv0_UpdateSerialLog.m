function [handles, newsrl] = helperGUIv0_UpdateSerialLog(handles)
    % update serial log 
    % TO DO: there should be a better way to do this; serial callback
    % should trigger an event or listener that logs the info 
    ReceivedData = handles.srl.UserData.ReceivedData; 
    newsrl = ~strcmp(ReceivedData, handles.srlLastMsg);
    if newsrl
        ud = handles.srl.UserData; 
        ud.TimeStamp = handles.lastSampleProcTime;
        if handles.srlP1 <= length(handles.srlStorage1)
            handles.srlStorage1(handles.srlP1) = ud;
            handles.srlP1 = handles.srlP1+1;
        else
            ud = handles.udBlank;
            % storage full; save
            SerialLog = handles.srlStorage1;
            svfn = [handles.SaveFileLoc,filesep,'SaveFile',num2str(handles.SaveFileN),'.mat'];
            disp(['Saving Serial to ',svfn])
            save(svfn,'SerialLog');
            handles.SaveFileN = handles.SaveFileN + 1;
            handles.srlP1 = 1;
            handles.srlStorage1 = repmat(ud, size(handles.srlStorage1));
        end
    end
    handles.srlLastMsg = ReceivedData;
end