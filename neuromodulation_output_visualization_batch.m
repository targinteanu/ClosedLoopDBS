% User selects folder; MATLAB runs all csv files 
filepath = uigetdir; 
csvFiles = dir(fullfile(filepath,'*.csv'));
for f = csvFiles'
    fpfn = fullfile(f.folder, f.name);
    [~,~,~,figTime,figPolr] = neuromodulation_output_visualization(fpfn);
    [fp,fn,fe] = fileparts(fpfn);
    fpfn = fullfile(fp,fn); % extension removed 
    usrans = questdlg('Save these figures?');
    if strcmpi(usrans, 'Yes')
        saveas(figTime, [fpfn,'_time'], 'fig');
        saveas(figTime, [fpfn,'_time'], 'png');
        saveas(figPolr, [fpfn,'_polar'], 'fig');
        saveas(figPolr, [fpfn,'_polar'], 'png');
    end
end