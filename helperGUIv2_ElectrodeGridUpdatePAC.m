function elecGridCData = helperGUIv2_ElectrodeGridUpdatePAC(...
    elecGridImageObj, elecGridFunc, ...
    channelIDlist, rawIDs, fltD4, bufferSizeGrid, fSamples)
% to be called at each iteration of main loop, whenever there is new data
% used for PAC or similar function that requires 2 concurrent data streams
% for each channel

            elecGridCData = elecGridImageObj.CData; 
            for ch = 1:min(numel(elecGridCData), length(channelIDlist))
                chID = channelIDlist(ch); 
                xInd = find(rawIDs == chID); 
                if numel(xInd)
                    x = fltD4{1,xInd}(:,2); y = fltD4{2,xInd}(:,2);
                    l = min(length(x), length(y));
                    x = x(1:l); y = y(1:l); % these should now be same length and time-aligned
                    if length(bufferSizeGrid) > 1
                        L = bufferSizeGrid(xInd);
                    else
                        L = bufferSizeGrid;
                    end
                    if height(x) > L
                        x = x((end-L+1):end, :);
                        y = y((end-L+1):end, :);
                    end
                    if height(x) < L
                        warning(['Channel ',num2str(ch),...
                            ' / ID ',num2str(chID),...
                            ' Electrode Grid buffer is not full length!'])
                    end
                    fSample_ch = fSamples(xInd);
                    if isempty(x)
                        elecGridCData(ch) = nan; % or just don't update (do nothing)?
                    else
                        elecGridCData(ch) = elecGridFunc(x, y, fSample_ch);
                    end
                else
                    elecGridCData(ch) = nan;
                end
            end

end