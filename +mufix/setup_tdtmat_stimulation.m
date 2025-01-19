function tdtmat = setup_tdtmat_stimulation(tdtmat, nbursts, init_delay_s, burst_interval_s, burst_len_s, pulse_hz, pulse_on_s, fs_mod)
%% 
%   tdtmat = setup_tdtmat_stimulation(tdtmat, nbursts, init_delay_s, burst_interval_s, burst_len_s, pulse_hz, pulse_on_s, fs_mod)
%   tdtmat = setup_tdtmat_stimulation(tdtmat, pulse_hz, pulse_on_s, fs_mod)
%
% Simulates optogenetic pulse-train stimulation delivery.
% 
% Input:
%   tdtmat           - tdt recording data structure
%   nbursts          - number of stimulation bursts
%   init_delay_s     - time in seconds to the first burst
%   burst_interval_s - time in seconds between start of each burst
%   burst_len_s      - duration of each burst in seconds
%   pulse_hz         - frequency of of pulse train in each burst
%   pulse_on_s       - duration of each pulse in seconds
%   fs_mod           - sampling frequency of the stimulation delivery.
%                      Exact timing of pulses are aligned to sampling frequency.

% arguments
%   tdtmat           (1,1) struct
%   nbursts          (1,1) {mustBePositive, mustBeInteger}
%   init_delay_s     (1,1) {mustBePositive}
%   burst_interval_s (1,1) {mustBePositive}
%   burst_len_s      (1,1) {mustBePositive}
%   pulse_hz         (1,1) {mustBePositive}
%   pulse_on_s       (1,1) {mustBePositive}
%   fs_mod           (1,1) {mustBePositive}
% end


  % Set up stim bursts (epochs)
  if ~isempty(nbursts)
    tdtmat.epocs.Brst.onset  = fsalign((0:nbursts-1)' * burst_interval_s + init_delay_s);
    tdtmat.epocs.Brst.offset = tdtmat.epocs.Brst.onset + fsalign(burst_len_s);
    tdtmat.epocs.Brst.data   = (1:numel(tdtmat.epocs.Brst.onset))';
  else
    burst_len_s = min(tdtmat.epocs.Brst.offset - tdtmat.epocs.Brst.onset);
  end

  % Set up stim pulses
  inter_pulse = fsalign(1 / pulse_hz);
  nPulse = burst_len_s * pulse_hz;
  tdtmat.epocs.Pls_.onset  = (0:nPulse-1)' * inter_pulse;
  tdtmat.epocs.Pls_.onset  = repmat(tdtmat.epocs.Pls_.onset, 1, length(tdtmat.epocs.Brst.onset));
  tdtmat.epocs.Pls_.onset  = tdtmat.epocs.Pls_.onset + tdtmat.epocs.Brst.onset';
  tdtmat.epocs.Pls_.onset  = fsalign(tdtmat.epocs.Pls_.onset(:));
  tdtmat.epocs.Pls_.offset = tdtmat.epocs.Pls_.onset + fsalign(pulse_on_s);
  tdtmat.epocs.Pls_.data   = repmat((1:nPulse)', length(tdtmat.epocs.Brst.onset), 1);

  function t = fsalign(t)
    % turn this off to better simulate misalignment in TDT timing
    % t = round(t * fs_mod) / fs_mod;
  end


end