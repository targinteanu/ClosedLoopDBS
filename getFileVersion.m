function commitHash = getFileVersion(filepath)
% Return a unique identifier string of the input file that can be used to
% identify its version. 
% Returns the latest commit hash
% by reading .git files (no system(), no Java)
% OR gives a unique date-time string identifier 

if nargin < 1
    filepath = '';
end

[dirpath, filename, filext] = fileparts(filepath);

try
    dirpath = fullfile(dirpath,'.git');

    % Step 1: Read .git/HEAD to get current branch ref
    headPath = fullfile(dirpath, 'HEAD');
    if ~isfile(headPath)
        error('This folder does not appear to be a Git repository.');
    end

    fid = fopen(headPath, 'r');
    headContent = strtrim(fread(fid, '*char')');
    fclose(fid);

    % HEAD can contain a direct ref or a detached commit hash
    if startsWith(headContent, 'ref:')
        refPath = strtrim(extractAfter(headContent, 'ref:'));
        refFile = fullfile(dirpath, refPath);
        if ~isfile(refFile)
            error('Ref file %s not found.', refFile);
        end

        fid = fopen(refFile, 'r');
        commitHash = strtrim(fread(fid, '*char')');
        fclose(fid);
    else
        % Detached HEAD — the hash is directly in HEAD
        commitHash = headContent;
    end

    % shorten for readability 
    commitHash = commitHash(1:14);

catch ME
    warning(['Unable to use git hash due to error: ',ME.message,...
        ' - returning unique date-time string instead.']);
    commitHash = char(datetime('now','TimeZone','UTC'), 'yyyyMMddHHmmss');
    pause(2); % ensures the above is unique on one machine unless clock is wrong
end
end