function handles = helperGUIv1_plotFlt(handles, tNow, fltPlt, ext_xlim)

        % update filtered data plot
        if numel(fltPlt)
            set(handles.h_filtDataTrace,'YData',fltPlt.Variables);
            set(handles.h_filtDataTrace,'XData',fltPlt.Time - tNow);
        end
        if ~sum(isnan(ext_xlim))
            set(handles.ax_filt, 'XLim', ext_xlim);
        end

end