function [carrier, tvec, sig_mod, sig_super] = simulate_LIA_modulation(sig, tsig, lia_freqs__hz, out_fs__hz)
arguments
  sig           (:,:)
  tsig          (1,:)
  lia_freqs__hz (1,:) {mustBePositive}
  out_fs__hz    (1,1) {mustBePositive}
end
%% 
%   [carrier, tvec, sig_mod, sig_super] = simulate_LIA_modulation(sig, tsig, lia_freqs__hz, out_fs__hz)
%
% Simulates LIA modulation of "sig" with lia_freqs__hz carrier frequencies.
% 
% Input:
%   sig           - nsig x time of signals to modulate
%   tsig          - corresponding 1 x time vector 
%   lia_freqs__hz - 1 x nsig carrier frequencies
%   out_fs__hz    - sampling rate of the modulation (and output)
%
% Output:
%   carrier       - 1 x time vector of the combined LIA modulated output
%   tvec          - corresponding timestamp vector
%   sig_mod       - individually LIA modulated signal
%   sig_super     - original sig up-sampled to out_fs__hz


    %% Input checking
    assert(size(sig,1) == numel(lia_freqs__hz), 'Number of carrier frequencies does not match the number of signals.');

    %% LIA modulate
    
    % up sample
    [sig_super, tvec] = resample(sig', tsig, out_fs__hz);
    sig_super = sig_super';

    % carrier sine waves
    % lia_sines = (sin(2*pi*lia_freqs__hz'.*tvec) + 1.0) / 2;  % in phase
    % lia_sines = (sin(2*pi*lia_freqs__hz'.*tvec + pi/2) + 1.0) / 2;  % in phase

    % random
    lia_sines = (sin(2*pi*lia_freqs__hz'.*tvec + rand(size(lia_freqs__hz'))*2*pi) + 1.0) / 2;
    
    % modulated signals
    sig_mod = sig_super .* lia_sines;
    carrier = sum(sig_mod,1);

end