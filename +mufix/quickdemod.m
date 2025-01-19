function [ydemod, cfg, figs] = quickdemod(y,fs,frqs,cfg_)
% Demodulate sinusoidal response in signal
%
%   ydemod = quickdemod(y,fs,frqs,cfg)
%
% Extract sinusoidal response at specified frqs from signal y.
%
% Inputs:
%   y    - data vector
%   fs   - sampling rate for y (Hz)
%   frqs - vector of sinusoidal frequencies to extract (Hz)
%   cfg  - optional settings struct (see below)
%
% Output:
%   ydemod - demodulated signals (nY x nfrqs)
%   cfg    - configuration struct with intermediate outputs
%
% cfg configuration struct fields:
%   lowpass_cutoff__Hz - Butterworth filter lowpass cut-off frequency (Hz)
%                        (default: 3 Hz)
%   lowpass_order      - Butterworth filter order for lowpass filtering (default: 8)
%   lowpass_zerophase  - Performs filtfilt zero-phase filtering (default: true)
%   downsample         - downsampling factor (default: round(fs/1000) aims for around 1kHz)
%   plot               - Plot output (0: off, 1: plot ydemod)
%
% Sun Lab: spencer.chen@rutgers.edu

VERSION = 'v1.5, 2024.6.25';

%% Configuration

figs = [];

cfg = struct();
cfg.lowpass_cutoff__Hz  = 3;
cfg.lowpass_order       = 8;  % 8 to match TDT-Fi2r demod, and 12 to match Fi1r (OP421)
cfg.lowpass_zerophase   = true;
cfg.downsample          = round(fs/1000);
cfg.plot                = 0;
cfg.frqs__Hz            = frqs;

if exist('cfg_','var')
  FF = fieldnames(cfg_);
  for ff = 1:numel(FF)
    cfg.(FF{ff}) = cfg_.(FF{ff});
  end
end

%% Demodulation

t = (0:numel(y)-1)/fs;
F = exp(1j*2*pi*frqs(:)*t).';
D = double(y(:)) .* F;


%% Filter and decimate demodulated data

% this approach has a singularity problem
% [filtB, filtA] = butter(cfg.lowpass_order,cfg.lowpass_cutoff/fs*2);

% this filter similates TDT filtering better
[cfg.filtZ, cfg.filtP, cfg.filtK] = butter(cfg.lowpass_order,cfg.lowpass_cutoff__Hz/fs*2);
[cfg.filtSOS, cfg.sosgain] = zp2sos(cfg.filtZ, cfg.filtP, cfg.filtK);

if cfg.lowpass_zerophase
  filtfunc = @(x) filtfilt(cfg.filtSOS, cfg.sosgain, x);
else
  filtfunc = @(x) sosfilt(cfg.filtSOS, x) * cfg.sosgain;
end

for ii = 1:size(D,2)
  D(:,ii) = filtfunc(D(:,ii));
end
ydemod = 4 * abs(downsample(D,cfg.downsample));
cfg.demodfs__Hz = fs/cfg.downsample;

if cfg.plot
  tr = t(1:cfg.downsample:end);
  figs(1) = figure;
  hold on
  for ii = 1:size(ydemod,2)
    plot(tr,ydemod(:,ii))
  end
  legstr = cellstr(num2str(frqs(:)));
  legend(legstr{:});
  xlabel('time (s)');
  ylabel('signal (V)');
  title('quickdemod');
end