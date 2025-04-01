function handles = helperGUIv1_plotRaw(handles, tNow, rawPlt, timeBuff, common_xlim)
    set(handles.h_rawDataTrace,'YData',rawPlt.Variables);
    set(handles.h_timingTrace,'YData',[nan; diff(timeBuff)]);
    set(handles.h_timeDispTrace,'YData',[nan; diff(handles.timeDispBuff)]);
    if handles.check_polar.Value
        set(handles.h_rawDataTrace,'XData',rawPlt.Time - tNow);
        set(handles.h_timingTrace,'XData', ...
            handles.time0 + seconds(timeBuff) - tNow );
        set(handles.h_timeDispTrace,'XData', ...
            handles.time0 + seconds(handles.timeDispBuff) - tNow );
        if ~sum(isnan(common_xlim)) % why is it sometimes nan??
            set(handles.ax_raw, 'XLim', common_xlim);
            set(handles.ax_timing, 'XLim', common_xlim);
        end
    end
end