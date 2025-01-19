function [yhat, cfg, figs] = fillmissingInterp(y,fs,method,tseg,cfg_)
%% Interpolate data segments using sinusoids
%
%   yhat = fillmissingInterp(y,fs,frqs,tseg,cfg)
%
% Algorithm fills data in the time segments specified by tseg with
% combination of sinusoids specified by frqs.
%
% This is used to remove saturated segments in the photometry sensor data
% with underlying oscillatory excitation drives.
%
% Inputs:
%   y      - data vector
%   fs     - sampling frequency (Hz)
%   method - fillmissing method
%   tseg   - [start end] column vector of time segments to interpolate
%            (seconds relative to start of y)
%   cfg    - additional configuration struct (see below)
%
% Output:
%   yhat   - interpolated data vector
%   cfg    - full configuration struct with intermediate processing info
%
% cfg struct fields:
%   plot      - enable plot processing (0: off (default), 1: final output, 2: segment fit)
%   label     - signal name label

VERSION = 'v1.5, 2024.6.25';

%% Process config

figs = [];

% Defaults
cfg = struct();
cfg.plot      = 0;
cfg.crop_pre__sec  = 0.01;
cfg.crop_post__sec = 0.01;
cfg.label     = '';
cfg.waitbar   = false;

% Overwrite with user settings
if exist('cfg_','var')
  FF = fieldnames(cfg_);
  for ff = 1:numel(FF)
    cfg.(FF{ff}) = cfg_.(FF{ff});
  end
end

yhat = y;

if isempty(tseg), return; end

idx_tseg = tseg * fs;
idx_tseg(:,1) = floor(idx_tseg(:,1));
idx_tseg(:,2) = ceil (idx_tseg(:,2));

% beyond the end of the data
idx_tseg = min(idx_tseg,numel(y));
idx_tseg = idx_tseg(idx_tseg(:,1) ~= idx_tseg(:,2),:);

idx_plot = round([-2*cfg.crop_pre__sec 2*cfg.crop_post__sec] * fs);
idx_msd  = idx_plot + [0 max(diff(idx_tseg,[],2))];

nTseg = size(idx_tseg,1);

%% Fitting loop
t = (0:numel(y)-1)/fs;

% Debug plot
if cfg.plot >= 2
  figs(2) = figure;
  tf = tiledlayout(2,1);
  ax1 = nexttile;
  ax2 = nexttile;
  hold(ax1,'on');
  hold(ax2,'on');
  xlabel(ax2,'index');
  ylabel(ax1,'original signal (V)');
  ylabel(ax2,'interp signal (V)');
elseif cfg.waitbar
  wf = waitbar(0, "Removing " + string(cfg.label) + " cross-talk ...");
end

for ii = 1:nTseg

  tframe  = idx_tseg(ii,1):idx_tseg(ii,2);
  
  % beyond the end of the signal
  tframe  = tframe (tframe  <= numel(t));
  yhat(tframe) = nan;
  cfg.tframe{ii}  = tframe;  

end
yhat = fillmissing(yhat,method);

mseg0 = zeros(1,diff(idx_msd)+1);
sseg0 = zeros(1,diff(idx_msd)+1);
mseg  = zeros(1,diff(idx_msd)+1);
sseg  = zeros(1,diff(idx_msd)+1);
nseg  = 0;
for ii = 1:nTseg
  if cfg.plot >= 1
    midx = (idx_tseg(ii,1)+idx_msd(1)) : (idx_tseg(ii,1)+idx_msd(2));
    if midx(end) <= numel(y)  % ignore segments beyond the end of the signal
      mseg0 = mseg0 + y   (midx);
      sseg0 = sseg0 + y   (midx).^2;
      mseg  = mseg  + yhat(midx);
      sseg  = sseg  + yhat(midx).^2;     
      nseg  = nseg  + 1;
    end
  end
  if cfg.plot >= 2
    plotidx = (idx_tseg(ii,1)+idx_plot(1)) : (idx_tseg(ii,2)+idx_plot(2));
    plotidx = plotidx(plotidx  <= numel(y));    
    xx_     = plotidx - idx_tseg(ii,1);
    plot(ax1, xx_, y   (plotidx));
    plot(ax2, xx_, yhat(plotidx));
    title(tf, sprintf('Interpolating Segment %d/%d',ii,nTseg));
    drawnow
  elseif exist('wf','var')
    waitbar(ii/nTseg,wf);
  end

end

if exist('wf','var')
  close(wf)
end

if cfg.plot >= 1
  figs(1) = figure;
  tiledlayout(2,1);
  
  nexttile
  hold on; 
  plot(t,y); 
  plot(t,yhat);
  legend('original','interp');
  xlabel('time (s)');
  ylabel('signal (V)');
  title(['fillmissing: ', method]);
  
  nexttile
  mseg0 = mseg0 / nseg;
  mseg  = mseg  / nseg;
  sseg0 = sseg0 / nseg - mseg0.^2;
  sseg  = sseg  / nseg - mseg.^2;
  sseg0 = sqrt(sseg0);
  sseg  = sqrt(sseg);
  xx    = (idx_msd(1):idx_msd(2)) / fs * 1000;
  pxx   = [xx, fliplr(xx)]; 
  cmap  = lines(2);
  hold on;
  patch(pxx, [mseg0, fliplr(mseg0)] + [sseg0, -fliplr(sseg0)], cmap(1,:), 'LineStyle', 'none', 'FaceAlpha', 0.5);
  patch(pxx, [mseg,  fliplr(mseg) ] + [sseg,  -fliplr(sseg)],  cmap(2,:), 'LineStyle', 'none', 'FaceAlpha', 0.5);
  h1 = plot(xx, mseg0, 'Color', cmap(1,:), 'LineWidth', 2);
  h2 = plot(xx, mseg,  'Color', cmap(2,:), 'LineWidth', 2);
  legend([h1, h2], 'original','interp');
  xlabel('pulse time (ms)');
  ylabel('average pulse signal (V)');  
  
end

