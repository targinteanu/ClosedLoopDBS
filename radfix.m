function radvals = radfix(radvals)
radvals = mod(radvals, 2*pi); % range [0, 2*pi]
radvals(radvals > pi) = radvals(radvals > pi)-2*pi; % range [-pi, pi]
end