function elecGridCData = helperGUIv1_ElectrodeGridUpdate(...
    elecGridImageObj, elecGridFunc, ...
    channelIDlist, rawIDs, rawD4, bufferSizeGrid, fSamples)
% to be called at each iteration of main loop, whenever there is new data

            elecGridCData = elecGridImageObj.CData; 
            for ch = 1:min(numel(elecGridCData), length(channelIDlist))
                chID = channelIDlist(ch); 
                xInd = find(rawIDs == chID); 
                if numel(xInd)
                    x = rawD4{xInd}(:,2); 
                    if length(bufferSizeGrid) > 1
                        L = bufferSizeGrid(xInd);
                    else
                        L = bufferSizeGrid;
                    end
                    if height(x) > L
                        x = x((end-L+1):end, :);
                    end
                    if height(x) < L
                        warning(['Channel ',num2str(ch),...
                            ' / ID ',num2str(chID),...
                            ' Electrode Grid buffer is not full length!'])
                    end
                    fSample_ch = fSamples(xInd);
                    elecGridCData(ch) = elecGridFunc(x, fSample_ch);
                else
                    elecGridCData(ch) = nan;
                end
            end

end