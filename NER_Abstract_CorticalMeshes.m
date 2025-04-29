
%PD23N003 - 'pds002.ns2'
%PD22N009 - 'datafile001.ns2'
%PD22N008 - 'PD22N008.ns2'
%PD23N009 - 'PDS001.ns2'

PD22N008_data = 'PD22N008.ns2';

PD22N009_data = 'datafile001.ns2'; %ns2 file name

%PD22N009%%%%%%%%%%%%%%%%%%%%%%

PD22N009_Phase_Data = obtainPhase(PD22N009_data, 'beta', 1, 100000); % (Baseline data 1-10,000) Get the phase data from every channel (
PD22N009_Channel_Names = getChannelNames(PD22N009_data);

PD22N009_PLV_VectorMatrix = compute_PLV_VectorMatrix(PD22N009_Phase_Data); %Get the PLV vector matrix for every channel (64x64)
PD22N009_PLV_Average = getAveragePLV(PD22N009_PLV_VectorMatrix); %Get the Average PLV for every channel (64x1)
ref_channel = 34;
PD22N009_PLV_Relative_ref = getPLV_wrt_Ref(PD22N009_PLV_VectorMatrix, ref_channel); %Get PLV relative to ref channel (64x1)

plot_PLV_Relative(PD22N009_Phase_Data,PD22N009_PLV_Relative_ref, ref_channel, PD22N009_Channel_Names);  %Plot a 21x3 grid - color corresponds to PLV wrt ref channel, # corresponds to avg. phase difference wrt ref channel

load_ACPC_FR_mesh('SubjectPD22N009_elec_acpc_fr.mat', 'PD22N009\PD22N009\freesurfer\surf\lh.pial', PD22N009_PLV_Relative_ref);

load_FSAVG_FRS_mesh('SubjectPD22N009_elec_fsavg_frs.mat', 'PD22N009\PD22N009\freesurfer\surf\lh.pial', PD22N009_PLV_Relative_ref);

min_val = min(PD22N009_PLV_Average);
max_val = max(PD22N009_PLV_Average);
PD22N009_PLV_Average_Normalized = (PD22N009_PLV_Average - min_val) / (max_val - min_val);
PD22N009_PLV_Average_Normalized(PD22N009_PLV_Average_Normalized >= 1) = 0.999;

load_ACPC_FR_mesh('SubjectPD22N009_elec_acpc_fr.mat', 'PD22N009\PD22N009\freesurfer\surf\lh.pial', PD22N009_PLV_Average_Normalized);

central_channel = 34;

best_neighbor = getBestNeighborByPLV(PD22N009_PLV_VectorMatrix, central_channel);

comparison_channels = [best_neighbor, 12,14,33,35,54,55,56];

compareRosePlots(PD22N009_Phase_Data, central_channel, comparison_channels);



%PD22N008%%%%%%%%%%%%%%%%%%%%%%

% PD22N008_Phase_Data = obtainPhase(PD22N008_data, 'beta', 1, 10000); %Get the phase data from every channel
% PD22N008_Channel_Names = getChannelNames(PD22N008_data);
% 
% PD22N008_PLV_VectorMatrix = compute_PLV_VectorMatrix(PD22N008_Phase_Data); %Get the PLV vector matrix for every channel (64x64)
% PD22N008_PLV_Average = getAveragePLV(PD22N008_PLV_VectorMatrix); %Get the Average PLV for every channel (64x1)
% ref_channel = 34;
% PD22N008_PLV_Relative_ref = getPLV_wrt_Ref(PD22N008_PLV_VectorMatrix, ref_channel); %Get PLV relative to ref channel (64x1)
% 
% plot_PLV_Relative(PD22N008_Phase_Data,PD22N008_PLV_Relative_ref, ref_channel, PD22N008_Channel_Names);  %Plot a 21x3 grid - color corresponds to PLV wrt ref channel, # corresponds to avg. phase difference wrt ref channel
% 
% load_ACPC_FR_mesh('SubjectPD22N009_elec_acpc_fr.mat', 'PD22N009\PD22N008\freesurfer\surf\lh.pial', PD22N009_PLV_Relative_ref);
% 
% load_FSAVG_FRS_mesh('SubjectPD22N009_elec_fsavg_frs.mat', 'PD22N009\PD22N008\freesurfer\surf\lh.pial', PD22N009_PLV_Relative_ref);
% 
% min_val = min(PD22N008_PLV_Average);
% max_val = max(PD22N008_PLV_Average);
% PD22N008_PLV_Average_Normalized = (PD22N008_PLV_Average - min_val) / (max_val - min_val);
% PD22N008_PLV_Average_Normalized(PD22N008_PLV_Average_Normalized >= 1) = 0.999;
% 
% load_ACPC_FR_mesh('SubjectPD22N008_elec_acpc_fr.mat', 'lh.pial.T1', PD22N008_PLV_Average_Normalized);



function load_ACPC_FR_mesh(acpc_path, pial_lh_path, data)

    % === Step 1: Load Electrode Positions and Mesh ===
    elec_acpc_fr = load(acpc_path).elec_acpc_fr;
    pial_lh = ft_read_headshape(pial_lh_path);
    pial_lh.coordsys = 'acpc';

    pial_lh = ft_convert_units(pial_lh, 'mm');
    elec_acpc_fr = ft_convert_units(elec_acpc_fr, 'mm');

    % === Step 2: Identify reference channels ===
    ref_idx = isnan(data);          % Logical array: true if it's reference
    data_clean = data;              % Copy data
    data_clean(ref_idx) = 0;         % Set NaNs to 0 for interpolation

    % === Step 3: Build source structure ===
    source = struct();
    source.pos = elec_acpc_fr.chanpos;
    source.pow = data_clean;
    source.dimord = 'pos';
    source.coordsys = 'acpc';

    % === Step 4: Interpolate electrode data ===
    cfg = [];
    cfg.parameter = 'pow';
    cfg.interpmethod = 'sphere_weighteddistance';
    cfg.sphereradius = 3;
    interp_source = ft_sourceinterpolate(cfg, source, pial_lh);

    % === Step 5: Plot the Interpolated Data ===
    cfg = [];
    cfg.funparameter = 'pow';
    cfg.funcolorlim = [0 1];            % normal 0-1 scaling
    cfg.method = 'surface';
    cfg.interpmethod = 'nearest';
    cfg.sphereradius = 3;
    cfg.camlight = 'no';

    ft_sourceplot(cfg, interp_source, pial_lh);

    view([-55 10]);
    material dull;
    lighting gouraud;
    camlight;

    % === Step 6: Normal color map (pure parula) ===
    colormap(parula);   % NO special red added

    % === Step 7: Plot electrodes ===
    % Plot normal electrodes
    elec_plot_colors = data_clean;   % Use interpolated data (0–1)
    ft_plot_sens(elec_acpc_fr, 'elecsize', 35, 'facecolor', elec_plot_colors, 'edgecolor', 'black');

    % === Step 8: Overlay reference electrodes manually ===
    if any(ref_idx)
        ref_positions = elec_acpc_fr.chanpos(ref_idx, :);
        hold on;
        scatter3(ref_positions(:,1), ref_positions(:,2), ref_positions(:,3), ...
            100, 'r', 'filled', 'MarkerEdgeColor', 'k', 'LineWidth', 1.5);
        hold off;
    end

end



function load_FSAVG_FRS_mesh(fsavg_path, pial_lh_path, data)

    % === Step 1: Load Electrode Positions and Mesh ===
    elec_fsavg = load(fsavg_path).elec_fsavg_frs;
    pial_lh = ft_read_headshape('lh.pial');
    pial_lh.coordsys = 'fsaverage';

    pial_lh = ft_convert_units(pial_lh, 'mm');
    elec_fsavg = ft_convert_units(elec_fsavg, 'mm');

    % === Step 2: Prepare Electrode Data ===
    ref_idx = isnan(data);           % Find reference channels
    elec_values = data;

    % Replace NaNs with a very large dummy value
    dummy_value = 2;  % well above PLV 1
    elec_values(ref_idx) = dummy_value;  % mark reference channels

    % === Step 3: Build source structure ===
    source = struct();
    source.pos = elec_fsavg.chanpos;
    source.inside = true(size(source.pos,1),1);
    source.pow = elec_values;
    source.dimord = 'pos';
    source.coordsys = 'fsaverage';

    % === Step 4: Interpolate electrode data ===
    cfg = [];
    cfg.parameter = 'pow';
    cfg.interpmethod = 'sphere_weighteddistance';
    cfg.sphereradius = 8;
    interp_source = ft_sourceinterpolate(cfg, source, pial_lh);

    % === Step 5: After interpolation, restore reference channel to a clean value ===
    % Any interpolated value >1.5 is from dummy (reference channels)
    ref_mask = interp_source.pow > 1.5;
    interp_source.pow(ref_mask) = 1.001;  % set a standard "red" value

    % === Step 6: Plot brain surface ===
    cfg = [];
    cfg.funparameter = 'pow';
    cfg.funcolorlim = [0 1];  % PLVs go between 0 and 1
    cfg.method = 'surface';
    cfg.interpmethod = 'nearest';
    cfg.sphereradius = 3;
    cfg.camlight = 'no';

    ft_sourceplot(cfg, interp_source, pial_lh);

    view([-55 10]);
    material dull;
    lighting gouraud;
    camlight;

    % === Step 7: Build custom colormap ===
    base_cmap = parula(256);
    red_color = [1 0 0];  % pure red
    cmap_with_red = [base_cmap; red_color];
    colormap(gca, cmap_with_red);

    % === Step 8: Also plot electrodes ===
    elec_plot_colors = data;
    elec_plot_colors(ref_idx) = 1.001;  % manually set reference channels to red
    ft_plot_sens(elec_fsavg, 'elecsize', 35, 'facecolor', elec_plot_colors, 'edgecolor', 'black');

end


function chNames = getChannelNames(raw_data)
    [dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
        channelName, channelIndex, channelIndexStim, channelNames] = ...
        getRecordedData_NS(raw_data,1);

    chNames = channelNames;
end


function phaseData = obtainPhase(raw_data, freq_range, start_idx, stop_idx)
   
    [dataOneChannel, StimTrainRec, dataAllChannels, SamplingFreq, t, tRel, ...
        channelName, channelIndex, channelIndexStim, channelNames] = ...
        getRecordedData_NS(raw_data,1);

    numChannels = size(dataAllChannels, 1) - 1;

    phaseData = cell(1, numChannels);
    for idx = 1:numChannels
        ch_raw = dataAllChannels(idx, start_idx:stop_idx);
        if strcmp(freq_range, 'theta')
            ch_filtered = Myeegfilt(ch_raw, SamplingFreq, 4, 9, 0, 1024);
        elseif strcmp(freq_range, 'beta')
            ch_filtered = Myeegfilt(ch_raw, SamplingFreq, 13, 30, 0, 1024);
        else
            error('Unknown frequency range: %s', freq_range);
        end
        [ch_phase, ~] = instPhaseFreq(ch_filtered, SamplingFreq);
        phaseData{idx} = ch_phase;
    end

end


function PLV_Relative_Vector = getPLV_wrt_Ref(PLV_VectorMatrix, ref_channel)
    % Extracts PLV of all channels with respect to the reference channel
    % Inputs:
    %   - PLV_VectorMatrix: NxN PLV matrix
    %   - ref_channel: reference channel index (scalar)
    % Output:
    %   - PLV_Relative_Vector: Nx1 vector of PLV values relative to ref_channel

    PLV_Relative_Vector = abs(PLV_VectorMatrix(:, ref_channel)); % take absolute value in case of complex numbers
    PLV_Relative_Vector(ref_channel) = NaN; % ignore self-PLV (you could also set to 0 if you want)
end



function PLV_VectorMatrix = compute_PLV_VectorMatrix(phase_data)
    
    numChannels = length(phase_data);
    PLV_vector_matrix = zeros(numChannels, numChannels);
    for ch1 = 1:numChannels
        for ch2 = ch1+1:numChannels
            phase_diff = phase_data{ch1} - phase_data{ch2};
            PLV_vector_matrix(ch1, ch2) = mean(exp(1j * phase_diff));
            PLV_vector_matrix(ch2, ch1) = conj(PLV_vector_matrix(ch1, ch2));
        end
    end

    PLV_VectorMatrix = PLV_vector_matrix;
end


function PLV_Average = getAveragePLV(vector_matrix)
    PLV_Average = mean(abs(vector_matrix), 2);
end


function plot_PLV_Relative(phase_data, PLV_relative_vector, ref_channel, channelNames)
    numChannels = length(phase_data);
    avg_phase_diff = zeros(numChannels, 1);

    for ch = 1:numChannels
        phase_diff = angle(exp(1j * (phase_data{ch} - phase_data{ref_channel}))); % circular diff
        avg_phase_diff(ch) = mean(phase_diff); % Mean circular phase difference
    end
    
    % Best neighbor determination (from PLV vector)
    PLV_copy = PLV_relative_vector;
    PLV_copy(ref_channel) = -1; % exclude self-comparison
    [max_plv, best_neighbor_ch] = max(PLV_copy);
    best_neighbor_phase_diff = avg_phase_diff(best_neighbor_ch);
    
    % Display results
    fprintf('Best neighbor of channel %d is channel %d\n', ref_channel, best_neighbor_ch);
    fprintf('Highest PLV: %.3f\n', max_plv);
    fprintf('Phase difference with channel %d: %.4f radians\n', ref_channel, best_neighbor_phase_diff);
    
    % --- Plot ---
    figure; hold on;
    cmap = parula(256);

    % Normalize PLV for color mapping
    plv_norm = PLV_relative_vector;
    plv_norm(isnan(plv_norm)) = 0; % Replace NaN with 0 for plotting
    plv_norm = plv_norm / max(plv_norm); % Normalize to [0,1]

    % Layout for plotting
    grid_width = 3;
    grid_height = 21;
    num_channels = grid_width * grid_height;
    electrode_positions = cell(num_channels, 1);
    channel_num = 1;
    for x = 0:grid_width-1
        for y = 0:grid_height-1
            electrode_positions{channel_num} = [x, y];
            channel_num = channel_num + 1;
        end
    end
        
       % === Plot filled boxes and text FIRST ===
    for i = 1:numChannels
        x_i = electrode_positions{i}(1);
        y_i = electrode_positions{i}(2);
    
        if i == ref_channel
            rectangle('Position', [x_i, y_i, 1, 1], ...
                      'FaceColor', 'r', 'EdgeColor', 'k'); % green fill for ref
        else
            color_idx = round(plv_norm(i) * 255) + 1;
            color_idx = min(max(color_idx, 1), 256);
    
            rectangle('Position', [x_i, y_i, 1, 1], ...
                      'FaceColor', cmap(color_idx, :), 'EdgeColor', 'k');
        end
    
        % Add channel name and phase value
        chan_name = string(channelNames{i});
        phase_val = avg_phase_diff(i);
    
        text(x_i + 0.5, y_i + 0.65, chan_name, ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 7, 'FontWeight', 'bold', 'Color', 'k');
    
        text(x_i + 0.5, y_i + 0.35, sprintf('%.3f', phase_val), ...
            'HorizontalAlignment', 'center', ...
            'VerticalAlignment', 'middle', ...
            'FontSize', 6, 'FontWeight', 'normal', 'Color', 'k');
    end

    % === Now separately highlight reference and best neighbor ON TOP ===
    % Draw bright green border for reference channel
    x_ref = electrode_positions{ref_channel}(1);
    y_ref = electrode_positions{ref_channel}(2);
    rectangle('Position', [x_ref, y_ref, 1, 1], ...
              'EdgeColor', 'r', 'LineWidth', 1, 'LineStyle', '-');
    
    % Draw red border for best neighbor
    x_best = electrode_positions{best_neighbor_ch}(1);
    y_best = electrode_positions{best_neighbor_ch}(2);
    rectangle('Position', [x_best, y_best, 1, 1], ...
              'EdgeColor', [0 0.8 0], 'LineWidth', 2, 'LineStyle', '-');

    % --- Plot final settings ---
    axis equal; xlim([0, grid_width]); ylim([0, grid_height]);
    set(gca, 'XTick', [], 'YTick', []);
    xlabel(''); ylabel('');
    title(sprintf('Avg Phase Diff (vs Ch %d) + PLV Heatmap', ref_channel));
    set(gca, 'YDir', 'normal');

    colormap(cmap);
    cbar = colorbar;
    ylabel(cbar, sprintf('PLV with Channel %d', ref_channel));

    hold off;
end


function best_neighbor = getBestNeighborByPLV(PLV_VectorMatrix, ref_channel)
% Returns the best neighbor (highest PLV) for a given reference channel
%
% Inputs:
%   PLV_VectorMatrix - NxN PLV matrix (complex or real)
%   ref_channel - scalar index of the reference channel
%
% Output:
%   best_neighbor - index of the best neighbor channel (highest PLV with ref_channel)

    plv_vector = abs(PLV_VectorMatrix(:, ref_channel)); % Take abs to handle complex numbers
    plv_vector(ref_channel) = -Inf; % Exclude self (don't pick self-PLV)
    [~, best_neighbor] = max(plv_vector); % Find max
end




function compareRosePlots(phaseData, ref_channel, comparison_channels)
% Compare reference channel to others and plot rose histograms with 95% CI and ±1 STD shading

    figure;
    numComparisons = numel(comparison_channels);

    % --- Decide layout ---
    if numComparisons <= 4
        nRows = 1;
        nCols = numComparisons;
    else
        nRows = 2;
        nCols = ceil(numComparisons / 2);
    end

    for idx = 1:numComparisons
        comp_ch = comparison_channels(idx);

        ref_phase = angle(exp(1j * phaseData{ref_channel}));
        comp_phase = angle(exp(1j * phaseData{comp_ch}));

        phase_diff = angle(exp(1j * (comp_phase - ref_phase)));
    
        [mu, ul, ll] = circ_mean(phase_diff, [], 2);   % Mean, upper, lower 95% CI
        [~, sigma] = circ_std(phase_diff, [], [], 2);  % Circular standard deviation

        % --- Plot ---
        subplot(nRows, nCols, idx);

        % Rose histogram
        h = polarhistogram(phase_diff, 75, 'Normalization', 'probability');
        max_r = max(h.Values); 
        rlim([0 max_r * 1.1]);
        hold on;

        % --- Shaded ±1 STD Cone ---
        theta_fill = linspace(mu - sigma, mu + sigma, 100);
        r_fill = ones(size(theta_fill)) * max_r * 0.8;
        polarplot(theta_fill, r_fill, 'Color', [0.6 0.8 1], 'LineWidth', 6);

        % --- Mean Direction ---
        polarplot([mu mu], [0 max_r], 'r-', 'LineWidth', 2);

        % --- 95% Confidence Interval Lines ---
        polarplot([ul ul], [0 max_r * 0.8], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
        polarplot([ll ll], [0 max_r * 0.8], '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);

        % --- Title ---
        title({
            sprintf('Reference %d vs Neighbor %d', ref_channel, comp_ch), ...
            sprintf('\\mu = %.3f rad', mu), ...
            sprintf('\\sigma = %.3f rad', sigma)
            }, 'FontSize', 10);

    end
end
