function times = extract_elapsed_times(filename)
%EXTRACT_ELAPSED_TIMES Extracts numeric elapsed times from file.
%   times = EXTRACT_ELAPSED_TIMES(filename) reads the text file specified
%   by filename and returns a column vector of numeric values found in
%   lines matching "Elapsed time is _ seconds." The search is case
%   sensitive for that exact phrase structure but allows arbitrary numeric
%   formats (integer, decimal, scientific) for the underscore.
%
%   Example:
%     t = extract_elapsed_times('full hardware timing GUI_PD_v0.txt');

fid = fopen(filename, 'r');
if fid == -1
    error('Could not open file: %s', filename);
end
cleanup = onCleanup(@() fclose(fid));

% Read entire file as text
txt = fread(fid, '*char')';
% Regular expression to capture numbers in lines like:
% "Elapsed time is 12.34 seconds." Allow optional spaces and number formats.
expr = 'Elapsed time is\s*([+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?)\s*seconds\.';
tokens = regexp(txt, expr, 'tokens');

% Flatten tokens and convert to numbers
if isempty(tokens)
    times = zeros(0,1);
else
    nums = cellfun(@(c) str2double(c{1}), tokens);
    times = nums(:);
end
end