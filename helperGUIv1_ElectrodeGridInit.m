function [elecGridImageObj, gridChannelList, gridChannelX, gridChannelY] = ...
    helperGUIv1_ElectrodeGridInit(ax_elecgrid, channelList, ...
    nrow, ncol, gridminval, gridmaxval)
% should be called on connection with hardware 

if nargin < 5
    gridminval = []; gridmaxval = [];
    if nargin < 3
        nrow = []; ncol = [];
    end
end
if isempty(nrow)
    % default to 21-by-3 grid 
    nrow = 21;
end
if isempty(ncol)
    % default to 21-by-3 grid
    ncol = 3;
end
if isempty(gridminval)
    gridminval = 0;
end
if isempty(gridmaxval)
    gridmaxval = 1e9;
end

img = nan(nrow, ncol);
[X,Y] = meshgrid(1:ncol, 1:nrow);
elecGridImageObj = imagesc(ax_elecgrid, img, [gridminval, gridmaxval]);
colormap(ax_elecgrid, 'parula'); colorbar(ax_elecgrid);
hold(ax_elecgrid,"on");
chL_ = channelList(1:min((ncol*nrow), length(channelList)));
X = X(:); X = X(1:length(chL_));
Y = Y(:); Y = Y(1:length(chL_));
text(ax_elecgrid, X,Y, chL_, ...
    'HorizontalAlignment','center', ...
    'VerticalAlignment','middle', ...
    'FontWeight', 'bold', ...
    'Color',[.8 0 0]);
gridChannelList = chL_; gridChannelX = X; gridChannelY = Y;

end