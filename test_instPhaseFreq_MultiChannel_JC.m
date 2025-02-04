function [phi_block, f_block] = test_instPhaseFreq_MultiChannel_JC(bdMatrix, fs)
% Provided a matrix of data sampled at constant rate fs, calculate the
% instantaneous phase phi and frequency f using the Hilbert transform.
% Each row of bdMatrix is a separate signal.
    % Ensure bdMatrix is a numeric array
    if ~isnumeric(bdMatrix)
        error('Input bdMatrix must be a numeric matrix.');
    end

    % Apply Hilbert transform to each row (each signal)
    H_all = hilbert(bdMatrix.').';  % Compute the analytic signal for each channel
    
    % Compute the average analytic signal across all channels
    H_avg = mean(H_all, 1);  
    
    % Extract instantaneous phase from the averaged analytic signal
    phi_block = angle(H_avg);  
    
    % Compute instantaneous frequency
    f_block = gradient(unwrap(phi_block)) * fs / (2 * pi);  

end