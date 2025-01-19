function [tdtmat, cfgout, dmod, figs] = rmCrosstalk_tdtmat(tdtmat, ph_ch_map, cfg_)
% Remove crosstalk from recording
%
%   cfg = rmCrosstalk(matfile, ph_ch_map, cfg)
%   cfg = rmCrosstalk(tdtmat, ph_ch_map, cfg)
%
% Inputs:
%   matfile   - path to the extracted TDT mat file
%   tdtmat    - struct in the format extracted by TDTbin2mat.m
%   ph_ch_map - table mapping photometry stream store names to sensor source channels
%               to be processed. Original data is copied to orgdata field in output.
%               e.g.:
%               StoreName    SourceStore    Channel    SeqID
%               ---------    -----------    -------    -----
%               "RCa2"       "Fi2r"         2          1
%               ...
%
% Output:
%   cfg    - configuration struct with intermediate outputs
%
% cfg configuration struct fields:
%   plot   - Plot output
%            0: off (default)
%            1: grand summary
%            2-n: enable n-1 level plots at calling functions
%   poststim_grace__sec - time after end of each pulse to crop out (default: 0.004)
%   dmod_shift__sec     - time shift adjustment to align with original TDT demodulated data
%                         (default: 0)
%   dmod_scale          - scaling adjustment to align with original TDT demodulated data
%                         (default: 1.1413)
%   xtalk_method        - "mufix"  = multi-frequency sinusoidal fit (sinePiecewiseInterp)
%                         "linear" = linear fit (linPiecewiseInterp)
%                         "sine1"  = single sinudoidal fit (sigPiecewiseInterp)
%                         "fill_linear" = calls fillmissing with "linear" interpolation
%                         "fill_spline" = calls fillmissing with "spline" interpolation
%                         "fill_pchip"  = calls fillmissing with "pchip" interpolation
%                         "fill_makima" = calls fillmissing with "makima" interpolation
%                         "none"        = calls fillmissing with "makima" interpolation
%
%   For xtalk removal:
%   max_amp        - vector of maximum amplitude for fitted sinusoids (default: Inf)
%   crop_pre__sec  - Amount of data to crop before tseg for fitting (seconds)
%                    (default: 2.5 cycles of the lowest frqs)
%   crop_post__sec - Amount of data to crop after tseg for fitting (seconds)
%                    (default: 2.5 cycles of the lowest frqs)
%
%   From quickdemod():
%   lowpass_cutoff - Butterworth filter lowpass cut-off frequency (Hz)
%                    (default: 3 Hz)
%   lowpass_order  - Butterworth filter order for lowpass filtering (default: 8)
%   lowpass_zerophase - Performs filtfilt zero-phase filtering (default: true)

VERSION = 'v1.5, 2024.6.25';
PHCFG_FREQ_IDX = [2 6 10];

import mufix.*

%% Configuration

% Default returns
figs = [];
dmod = {};
cfgout = struct();

cfg = struct();
cfg.xtalk_method = "mufix";
cfg.plot         = 0;
cfg.verbose      = 1;
cfg.dmod_scale   = 1.413/4;
cfg.dmod_shift__sec     = 0;
cfg.poststim_grace__sec = 0.004;  % response delay (~1ms) + saturation rebound (~2ms)

if exist('cfg_','var')
  FF = fieldnames(cfg_);
  for ff = 1:numel(FF)
    cfg.(FF{ff}) = cfg_.(FF{ff});
  end
end

switch cfg.verbose
  case 1
    vprintf = @(varargin) fprintf(varargin{:});
  otherwise
    vprintf = @(varargin) [];
end


%% Load TDT mat
if ischar(tdtmat) || isstring(tdtmat)
  mfile = tdtmat;
  
  % Check file exists
  assert(exist(mfile,'file'), 'MAT file not found at: %s', mfile);
  
  tdtmat = load(mfile);
end

% Include checks on required tdtmat fields here


%% Remove crosstalk and demodulate signal

if ~exist("ph_ch_map",'var') 
  ph_ch_map = guessPhMap(tdtmat);
end
if isempty(ph_ch_map), return; end

[uSrc,~,ridx] = unique(ph_ch_map(:,{'SourceStore','Channel'}));
if isfield(tdtmat,'epocs') &&  isfield(tdtmat.epocs,'Pls_')
  tseg = [tdtmat.epocs.Pls_.onset(:), tdtmat.epocs.Pls_.offset(:) + cfg.poststim_grace__sec];
else
  tseg = [];
end

i = 1;
while(i <= length(uSrc.SourceStore))
    if ~any(strcmp(fieldnames(tdtmat.streams), uSrc.SourceStore(i)))
        uSrc(i,:) = [];
    else
        i = i + 1;
    end
end

if cfg.plot
  figs(1) = figure;
  tf = tiledlayout(2,1);
  ax1 = nexttile;
  ax2 = nexttile;
  hold(ax1,'on');
  hold(ax2,'on');
  title(tf, 'rmCrosstalk');
  xlabel(ax2,'time (s)');  
  ylabel(ax1,'signal (mV)');
  ylabel(ax2,'\Delta signal (mV)');
  cmap = lines(height(ph_ch_map));
  plotorder = [];
end

dmod = cell(height(uSrc),1);

for ss = 1:height(uSrc)
  
  vprintf('Processing "%s" Channel %d ...\n', uSrc.SourceStore(ss), uSrc.Channel(ss));
  chstr = sprintf('%s-%d', uSrc.SourceStore(ss), uSrc.Channel(ss));
  
  % The corresponding config store has the 'i' prefix: e.g. Fi2r --> Fi2i
  phcfg_store = char(uSrc.SourceStore(ss));
  phcfg_store(4) = 'i';
  
  phsrc  = tdtmat.streams.(uSrc.SourceStore(ss));
  phcfg  = tdtmat.scalars.(phcfg_store);  
  srcdat = phsrc.data(uSrc.Channel(ss),:);
  fs     = phsrc.fs;
  
  nDrv   = size(phcfg.data,1) / 4;
  frqs   = phcfg.data(PHCFG_FREQ_IDX(1:nDrv),1);

  dmod{ss} = struct('SourceStore',uSrc.SourceStore(ss),'Channel',uSrc.Channel(ss));
  
  % remove crosstalk
  if cfg.xtalk_method ~= "none"
    interpcfg = cfg;
    interpcfg.plot = max(0,cfg.plot-1);
    interpcfg.label = sprintf("%s:ch%d", uSrc.SourceStore(ss), uSrc.Channel(ss));

    switch cfg.xtalk_method
      case "linear"
        [dmod{ss}.interpdat, dmod{ss}.interpcfg, figs] = linPiecewiseInterp(srcdat, fs, tseg, interpcfg);  
      case "sine1"
        [dmod{ss}.interpdat, dmod{ss}.interpcfg, figs] = sigPiecewiseInterp(srcdat, fs, frqs, tseg, interpcfg);  
      case "fill_linear"
        [dmod{ss}.interpdat, dmod{ss}.interpcfg, figs] = fillmissingInterp(srcdat, fs, 'linear', tseg, interpcfg);  
      case "fill_spline"
        [dmod{ss}.interpdat, dmod{ss}.interpcfg, figs] = fillmissingInterp(srcdat, fs, 'spline', tseg, interpcfg);  
      case "fill_pchip"
        [dmod{ss}.interpdat, dmod{ss}.interpcfg, figs] = fillmissingInterp(srcdat, fs, 'pchip',  tseg, interpcfg);  
      case "fill_makima"
        [dmod{ss}.interpdat, dmod{ss}.interpcfg, figs] = fillmissingInterp(srcdat, fs, 'makima', tseg, interpcfg);  
      case "mufix1"
        [dmod{ss}.interpdat, dmod{ss}.interpcfg, figs] = sinePiecewiseInterp(srcdat, fs, frqs(1), tseg, interpcfg);  
      otherwise  % "mufix"
        [dmod{ss}.interpdat, dmod{ss}.interpcfg, figs] = sinePiecewiseInterp(srcdat, fs, frqs, tseg, interpcfg);  
    end
    dmod{ss}.interpcfg.tseg = tseg;
    dmod{ss}.interpcfg.fs   = fs;

    if numel(figs) >= 1
      ax = findobj(get(figs(1),'Children'),'-class','matlab.graphics.axis.Axes');
      ax(end).Title.String = [chstr,': ',ax(end).Title.String];
    end
    if numel(figs) >= 2
      ax = findobj(get(figs(2),'Children'),'-class','matlab.graphics.axis.Axes');
      ax(end).Title.String = [chstr,': ',ax(end).Title.String];
    end
  end
  
  % demodulate
  phrows = find(ridx == ss);
  drv    = ph_ch_map.SeqID(phrows);  
  frqs2  = frqs(drv);
  
  dmodcfg = cfg;
  dmodcfg.plot = max(0,cfg.plot-1);
  dmodcfg.downsample = 1; % do not downsample here
  if cfg.xtalk_method == "none"
    [dmod{ss}.dmoddat, dmod{ss}.dmodcfg, figs] = quickdemod(srcdat, fs, frqs2, dmodcfg);
  else
    [dmod{ss}.dmoddat, dmod{ss}.dmodcfg, figs] = quickdemod(dmod{ss}.interpdat, fs, frqs2, dmodcfg);
  end
  if numel(figs) >= 1
    ax = findobj(get(figs(1),'Children'),'-class','matlab.graphics.axis.Axes');
    ax.Title.String = [chstr,': ',ax.Title.String];
  end
  
  % fudge it like TDT
  dmod{ss}.dmoddat = 1000 * cfg.dmod_scale * dmod{ss}.dmoddat;
  dmod_shift = cfg.dmod_shift__sec * fs;
  if dmod_shift > 0
    dmod{ss}.dmoddat = [zeros(dmod_shift,size(dmod{ss}.dmoddat,2)); dmod{ss}.dmoddat];
  elseif dmod_shift < 0
    dmod{ss}.dmoddat = dmod{ss}.dmoddat((-dmod_shift+1):end,:);
  end
  
  % selected config output
  cfgout = cfg;
  cfgout.VERSION   = VERSION;
  cfgout.dmodcfg   = dmod{ss}.dmodcfg;
  % cfgout.dmodcfg   = subset_struct(dmod{ss}.dmodcfg, ...
  %   {'lowpass_cutoff__Hz','lowpass_order','lowpass_zerophase','downsample'});
  if cfg.xtalk_method ~= "none"
    cfgout.interpcfg = subset_struct(dmod{ss}.interpcfg, ...
      {'max_amp','crop_pre__sec','crop_post__sec','fitparams','frqs__Hz'});
  end
  
  % assign the output to tdtmat
  for rr = 1:numel(phrows)
    phname = ph_ch_map.StoreName(phrows(rr));
    vprintf('  Extracting "%s" ...\n', phname);
    
    phdat  = tdtmat.streams.(phname);
    if ~isfield(phdat,'orgdata')
      phdat.orgdata = phdat.data;
    end
    phdat.data    = dmod{ss}.dmoddat(:,rr)';
    phdat.dmodcfg = cfgout;
    
    % resample to match existing data 
    if ~isfield(phdat,'fs')
      % If no sampling frequency exist, aim for 1kHz
      ds = round(fs/1000);
      phdat.fs = fs / ds;
    end
    phdat.dmodcfg.downsample = round(fs/phdat.fs);  % this should end up to be an integer anyway
    phdat.data = downsample(phdat.data, phdat.dmodcfg.downsample);
    tdtmat.streams.(phname)  = phdat;
    
    % plot result
    if cfg.plot
      nOrg = numel(phdat.orgdata);
      nNew = numel(phdat.data);
      tt1 = (0:(nOrg-1)) * phdat.fs;
      tt2 = (0:(nNew-1)) * phdat.fs;
      plot(ax1, tt1, phdat.orgdata, 'Color', cmap(phrows(rr),:).^(0.3));
      plot(ax1, tt2, phdat.data,    'Color', cmap(phrows(rr),:));
      n = min(nOrg,nNew);
      plot(ax2, tt1(1:n), phdat.orgdata(1:n)-phdat.data(1:n), 'Color', cmap(phrows(rr),:));
      plotorder(end+1) = phrows(rr);
    end
  end
  
  % legends
  if cfg.plot
    legstr2 = ph_ch_map.StoreName(plotorder);
    legstr1 = [legstr2+" (org)",  legstr2+" (interp)"];
    legstr1 = legstr1';    
    legend(ax1,legstr1(:));
    legend(ax2,legstr2);
  end
end

tdtmat.info.phmap = ph_ch_map;

end

function S = subset_struct(S,keep_fields)
  FF = fieldnames(S);
  S  = rmfield(S, setdiff(FF,keep_fields));
end