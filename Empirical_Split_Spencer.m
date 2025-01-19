% Modular code should allow to adjust the following variables:
% Lenght of after-pulse artifact, amount of before and after data to model
% from, frequencies chosen, number of frequencies, 

SZDB    = 'seizure_list.xls';
SZ_LIST = readtable(SZDB,"TextType","string");

RECORDROOT = 'Recordings';
RECORDINGS = [ ...
              "OP2718-211103-111406", "OP2718-211104-095749", "OP2718-211115-134907", "OP2718-211117-130838"...
              "OP275-211103-115508", "OP275-211104-110601", "OP275-211115-125031", "OP275-211117-140342"
             ];
RECFILES  = fullfile(RECORDROOT, RECORDINGS);

% Frequencies
FREQ.f1 = 210; % Frequency of signal A in Hz
FREQ.f2 = 330; % Frequency of signal B in Hz
FREQ.f3 = 530; % Frequency of signal C in Hz

% % Correlation
CORR.extraEp_s = [-2 4];


%% Base simulation

STIM = struct();
STIM.pulse_hz = 20;
STIM.pulse_on_s = 5/1000;
STIM.burst_len_s = 30;

[T2,  SErr, SErr_Time] = mufix.sim_empirical(RECFILES,STIM,FREQ,CORR,["none","mufix"],SZ_LIST);

% Save output
currtime = datetime();
T2.TimeStamp(:) = string(currtime,'yyyyMMdd_HHmmss');
filename        = sprintf('%s_%dHz_%dms_empirtrials', ...
  char(currtime,'yyyyMMdd'),round(STIM.pulse_hz), round(1000*STIM.pulse_on_s));
xlfile = [filename,'.xls'];

% filename        = [char(currtime,'yyyyMMdd'),'_empirtrials'];
writetable(T2, xlfile, 'Sheet','Corr','WriteMode', 'append');

% Configuration
config = [[fieldnames(STIM), struct2cell(STIM)];
          [fieldnames(CORR), struct2cell(CORR)]];
config = [cell(height(config),2) config];

writecell({'TimeStamp' T2.TimeStamp(1)}, xlfile, 'Sheet','Config','WriteMode', 'append');
writecell(config, xlfile, 'Sheet','Config','WriteMode', 'append');

save(filename, 'T2', 'STIM', 'CORR', 'FREQ', 'SZ_LIST');
save([filename,'_SErr'], '-v7.3', 'SErr', 'SErr_Time', 'STIM', 'CORR', 'FREQ', 'SZ_LIST');


%% 20Hz, multiple pulse widths
PULSE_WIDTHS = [15 25 33 34 35];

STIM = struct();
STIM.pulse_hz = 20;

for pp = 1:numel(PULSE_WIDTHS)
  fprintf('Running 20Hz Empirical simulation with pulse width %d ms ...\n', PULSE_WIDTHS(pp));

  STIM.pulse_on_s = PULSE_WIDTHS(pp)/1000;    
  T2 = mufix.sim_empirical(RECFILES,STIM,FREQ,CORR,["none","mufix"],SZ_LIST);

  % Save output
  currtime = datetime();
  T2.TimeStamp(:) = string(currtime,'yyyyMMdd_HHmmss');
  filename        = sprintf('%s_%dHz_%dms_empirtrials', ...
    char(currtime,'yyyyMMdd'),round(STIM.pulse_hz), round(1000*STIM.pulse_on_s));
  
  save(filename, 'T2', 'STIM', 'CORR', 'FREQ', 'SZ_LIST');
end


%% 1Hz, multiple pulse widths
PULSE_WIDTHS = [50 100 150 200 400 900];

STIM = struct();
STIM.pulse_hz = 1;

for pp = 1:numel(PULSE_WIDTHS)
  fprintf('Running 1Hz Empirical simulation with pulse width %d ms ...\n', PULSE_WIDTHS(pp));

  STIM.pulse_on_s = PULSE_WIDTHS(pp)/1000;    
  T2 = mufix.sim_empirical(RECFILES,STIM,FREQ,CORR,"mufix",SZ_LIST);

  % Save output
  currtime = datetime();
  T2.TimeStamp(:) = string(currtime,'yyyyMMdd_HHmmss');
  filename        = sprintf('%s_%dHz_%dms_empirtrials', ...
    char(currtime,'yyyyMMdd'),round(STIM.pulse_hz), round(1000*STIM.pulse_on_s));
  
  save(filename, 'T2', 'STIM', 'CORR', 'FREQ', 'SZ_LIST');
end




%% Simulate 1 trial for plotting

XTALK_METHOD = ["none","mufix","linear"]; %,"sine1","fill_linear","fill_spline"];

% Crosstalk simulation
STIM1 = STIM;
STIM1.pulse_hz   = 20;
STIM1.sim_bursts = 1;

PRE_POST_CYCLES  = 1/8;

[T2,  SErr, SErr_Time, tdtextr, carrier, dmodout] = mufix.sim_empirical(RECFILES(1), STIM1, FREQ, CORR, XTALK_METHOD, SZ_LIST, PRE_POST_CYCLES);

% Plot dmod output
figure;
T = tiledlayout(3, numel(XTALK_METHOD),'tileindexing','columnmajor');
title(T,"Demodulation Output Illustration");
tt = (0:size(dmodout{1}.original,2)-1) / tdtextr.streams.RCa2.fs;
tt = tt - tdtextr.epocs.Brst.onset(1);
ax = [];
for ii = 1:numel(XTALK_METHOD)

  % Demod against original
  ax(1,ii) = nexttile;
  plot(tt, dmodout{1}.(XTALK_METHOD(ii))(2,:));
  hold on
  plot(tt, dmodout{1}.original(2,:),'k--');
  title(XTALK_METHOD(ii))
  xlabel('s')
  ylabel('mV')

  % Demod Error
  ax(2,ii) = nexttile;
  plot(tt, dmodout{1}.(XTALK_METHOD(ii))(2,:) - dmodout{1}.original(2,:));
  xlabel('s')
  ylabel('\DeltamV')
  title(sprintf('r_{cr} = %.2f, r_{dm} = %.2f, mse_{dm} = %f', ...
    T2.("Carr_"  +XTALK_METHOD(ii))(1), ...
    T2.("Corr_B_"+XTALK_METHOD(ii))(1), ...
    T2.("MSE_B_" +XTALK_METHOD(ii))(1)))
  
  % Interp error
  tframe = carrier{1}.mufix_tframe;
  % nidx   = cellfun(@numel, tframe);
  for jj = 1:numel(tframe)
    t_   = (tframe{jj}(1)-20):(tframe{jj}(end)+20);
    yorg = carrier{1}.original(t_);
    if XTALK_METHOD(ii) == "none"
      yhat = carrier{1}.xtalk(t_);
    else
      yhat = carrier{1}.(XTALK_METHOD(ii))(t_);
    end
    ydiff(jj,1:numel(t_)) = (yhat - yorg).^2;
  end
  mydiff = mean(ydiff,1);
  ax(3,ii) = nexttile;
  tts = 1000*((1:numel(mydiff)) - 20) / tdtextr.streams.Fi2r.fs;
  plot(tts, 1e6*mydiff);
  xlabel('ms')
  ylabel('mV^2')
  title(sprintf('R2_{fit} = %.2f [%.2f,%.2f], mse_{fit} = %f', ...
    mean(T2.("R2_Pulse_"  +XTALK_METHOD(ii))(1,:)), ...
    min (T2.("R2_Pulse_"  +XTALK_METHOD(ii))(1,:)), ...
    max (T2.("R2_Pulse_"  +XTALK_METHOD(ii))(1,:)), ...
    1e6*mean(T2.("MSE_Pulse_" +XTALK_METHOD(ii))(1,:))))

end
linkaxes(ax(1:2,:),'x')
xlim(ax(1,1),[-10 40])