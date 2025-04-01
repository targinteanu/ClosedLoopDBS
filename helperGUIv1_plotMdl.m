function handles = helperGUIv1_plotMdl(handles, tNow, fltPlt, forPlt, forBuff, tSt, artPlt)

            % update model-forecasted data plot
            if handles.check_polar.Value
                if numel(forPlt)
                set(handles.h_predTrace,'YData',forPlt.Variables);
                set(handles.h_predTrace,'XData',forPlt.Time - tNow);
                end
            end
            tPk = forBuff(:,1); tTr = forBuff(:,2); % time to peak, trough (s)
            set(handles.h_peakTrace,'YData',zeros(size(tPk)));
            set(handles.h_trouTrace,'YData',zeros(size(tTr)));
            set(handles.h_stimTrace,'YData',zeros(size(tSt)));
            if handles.check_polar.Value
            set(handles.h_peakTrace,'XData',handles.time0 + seconds(tPk) - tNow);
            set(handles.h_trouTrace,'XData',handles.time0 + seconds(tTr) - tNow);
            set(handles.h_stimTrace,'XData',handles.time0 + seconds(tSt) - tNow);
            end

            % update artifact-removed plot
            if handles.check_artifact.Value
                if ~isempty(artPlt)
                    set(handles.h_artDataTrace,'YData',artPlt.Variables);
                    if handles.check_polar.Value
                    set(handles.h_artDataTrace,'XData',artPlt.Time - tNow);
                    end
                end
            end

            % plot sine wave 
            if handles.check_polar.Value
                % set(handles.h_sineTrace,'YData', ...
                % set(handles.h_sineTrace,'XData', ...
            end

            % evaluate accuracy of above --> polar histogram
            if handles.check_polar.Value

                % row indexes of peak events 
                rowPk = nan(size(tPk)); 
                for r = 1:height(tPk)
                    % find time of current event relative to time now
                    tPk_r = handles.time0 + seconds(tPk(r)) ;
                    if (tPk_r >= min(fltPlt.Time)) && (tPk_r <= max(fltPlt.Time))
                        % current event is in time range shown on screen,
                        % so let row index be the nearest 
                        [~,rowPk(r)] = min(abs( tPk_r - fltPlt.Time ));
                    end
                end
                rowPk = rowPk(~isnan(rowPk));

                % row indexes of trough events 
                rowTr = nan(size(tTr)); 
                for r = 1:height(tTr)
                    % find time of current event relative to time now
                    tTr_r = handles.time0 + seconds(tTr(r)) ;
                    if (tTr_r >= min(fltPlt.Time)) && (tTr_r <= max(fltPlt.Time))
                        % current event is in time range shown on screen,
                        % so let row index be the nearest 
                        [~,rowTr(r)] = min(abs( tTr_r - fltPlt.Time ));
                    end
                end
                rowTr = rowTr(~isnan(rowTr));

                % calc phase and histogram
                phi = instPhaseFreq(fltPlt.Variables, handles.fSample);
                phiPk = phi(rowPk); phiTr = phi(rowTr);
                set(handles.h_peakPhase,'Data',phiPk);
                set(handles.h_trouPhase,'Data',phiTr);
            end

end