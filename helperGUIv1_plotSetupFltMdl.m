function handles = helperGUIv1_plotSetupFltMdl(handles, tNow, ...
    fltPlt, forPlt, forBuff, tSt, common_xlim, unitname)

        % initiate filtered data plot
        common_xdiff = diff(common_xlim); 
        ext_xdiff = common_xdiff * handles.ax_filt.InnerPosition(3) / ...
            handles.ax_raw.InnerPosition(3); 
        ext_xlim = [0, ext_xdiff] + common_xlim(1); % align left 
        axes(handles.ax_filt); hold off; 
        if isempty(fltPlt)
            error('Filter was not actually set up. Something is wrong in the code.')
        end
        handles.h_filtDataTrace = plot(fltPlt.Time - tNow, fltPlt.Variables); 
        grid on; hold on; 
        title('Filtered & Predicted Data'); xlabel('time (s)'); ylabel(unitname);
        xlim(ext_xlim);

        % initiate prediction & peak/trough indicators overlayed on
        % filtered plot
        if handles.MdlSetUp
            if isempty(forPlt)
                error('Forecast was not actually set up. Something is wrong in the code.')
            end
            tPk = forBuff(:,1); tTr = forBuff(:,2);
            handles.h_peakTrace = plot(handles.time0 + seconds(tPk) - tNow, zeros(size(tPk)), ...
                '^', 'Color',"#EDB120"); 
            handles.h_trouTrace = plot(handles.time0 + seconds(tTr) - tNow, zeros(size(tTr)), ...
                'v', 'Color',"#EDB120"); 
            handles.h_stimTrace = plot(handles.time0 + seconds(tSt) - tNow, zeros(size(tSt)), ...
                '*', 'Color','r'); 
            handles.h_predTrace = plot(forPlt.Time - tNow, forPlt.Variables, ':');
            %handles.h_sineTrace = plot(xValues3,handles.sineDataBuffer,'--');

            % initiate polar histogram 
            bedge = (-1:2:35)*pi/18; 
            axes(handles.ax_polar); hold off;
            handles.h_peakPhase = polarhistogram(nan, bedge);
            hold on; 
            handles.h_trouPhase = polarhistogram(nan, bedge);
            title('Actual Phase')
        end

end