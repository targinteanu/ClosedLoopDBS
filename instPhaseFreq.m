function [phi_block, f_block, A_block] = instPhaseFreq(blockData, fs)
% Provided a block of data sampled at constant rate fs, calculate the
% instantaneous phase phi, frequency f, and amplitude A using the Hilbert 
% transform. 

% use hilbert transform to determine inst. phase
H_block = hilbert(blockData); 
A_block = abs(H_block);
phi_block = angle(H_block); 

% calculate inst. freq. 
f_block = gradient(unwrap(phi_block)')' *fs/(2*pi);
%f_block = abs(f_block); % ?? why does it go negative??

end