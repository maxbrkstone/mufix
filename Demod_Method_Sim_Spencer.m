%import recording to use the PhMap
TDTMAT_TEMPLATE = "tdtmat_sim3ph_template";

% Crosstalk simulation
STIM.burst_len_s          = 30;
STIM.pulse_hz             = 20;
STIM.pulse_on_s           = 5/1000;
STIM.burst_interval_s     = 120;
STIM.init_delay_s         = 60;
STIM.sim_bursts           = 10;

% Correlation
CORR.extraEp_s            = [-2 4];


%% Simulate all interp methods
XTALK_METHOD = ["none","mufix", "linear","sine1","fill_linear","fill_spline"];

% Ground truth configurations
GNDTRUTH = struct();
GNDTRUTH.padding_s        = 120;
GNDTRUTH.recording_s      = 1800;
GNDTRUTH.nodes_per_second = 3;
GNDTRUTH.sig_range        = [15 30; 15 20; 5 10];


starttime = datetime;
[T2, SErr, SErr_Time, tdtextr] = mufix.sim_artificial(TDTMAT_TEMPLATE, 1, STIM, GNDTRUTH, CORR, XTALK_METHOD);
endtime   = datetime;

% Output file
currtime = datetime();
T2.TimeStamp(:) = string(currtime,'yyyyMMdd_HHmmss');
filename        = [char(currtime,'yyyyMMdd'),'_method_trials'];
xlfile          = [filename,'.xlsx'];
writetable(T2, xlfile, 'Sheet','Corr','WriteMode', 'append');

% Configuration
config = [[fieldnames(GNDTRUTH), struct2cell(GNDTRUTH)];
          [fieldnames(STIM),     struct2cell(STIM)];
          [fieldnames(CORR),     struct2cell(CORR)]];
config = [cell(height(config),2) config];

writecell({'TimeStamp' T2.TimeStamp(1)}, xlfile, 'Sheet','Config','WriteMode', 'append');
writecell(config, xlfile, 'Sheet','Config','WriteMode', 'append');

save(filename, 'T2', 'tdtextr', 'GNDTRUTH', 'STIM', 'CORR');
save([filename,'_SErr'], '-v7.3', 'SErr', 'SErr_Time', 'tdtextr', 'GNDTRUTH', 'STIM', 'CORR');


%% Simulate 1 trial for plotting

XTALK_METHOD = ["none","mufix","linear"]; %,"sine1","fill_linear","fill_spline"];

% Ground truth configurations
GNDTRUTH = struct();
GNDTRUTH.padding_s        = 0;
GNDTRUTH.recording_s      = 200;
GNDTRUTH.nodes_per_second = 3;
GNDTRUTH.sig_range        = [15 30; 15 30; 15 30];

% Crosstalk simulation
STIM1 = STIM;
STIM1.pulse_hz   = 21;
STIM1.sim_bursts = 1;

PRE_POST_CYCLES  = 1/2;

[T2, SErr, SErr_Time, tdtextr, carrier, dmodout] = mufix.sim_artificial(TDTMAT_TEMPLATE, 1, STIM1, GNDTRUTH, CORR, XTALK_METHOD, PRE_POST_CYCLES);

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


%% Simulate different signal complexity
NODES_PER_S = [1 40 60 80 100 120 150];
XTALK_METHOD = ["none", "mufix"];
            
% Ground truth configurations
GNDTRUTH = struct();
GNDTRUTH.padding_s        = 120;
GNDTRUTH.recording_s      = 1800;
GNDTRUTH.sig_range        = [15 30; 15 20; 5 10];

for pp = 1:numel(NODES_PER_S)
  fprintf('Running Articial ground truth simulation with dynamics %d nodes per second ...\n', NODES_PER_S(pp));

  GNDTRUTH.nodes_per_second = NODES_PER_S(pp);    
  [T2, ~, ~, tdtextr] = mufix.sim_artificial(TDTMAT_TEMPLATE, 20, STIM, GNDTRUTH, CORR, XTALK_METHOD);

  % Save output
  currtime = datetime();
  T2.TimeStamp(:) = string(currtime,'yyyyMMdd_HHmmss');
  filename        = sprintf('%s_%dnps_sim_dynamics', ...
    char(currtime,'yyyyMMdd'),GNDTRUTH.nodes_per_second);
  
  save(filename, 'T2', 'tdtextr', 'GNDTRUTH', 'STIM', 'CORR');
end

%% Simulate different fitting window
NODES_PER_S = 3;
XTALK_METHOD = ["none", "mufix"];
PRE_POST_WIN = [(1/8)];
            
% Ground truth configurations
GNDTRUTH = struct();
GNDTRUTH.padding_s        = 120;
GNDTRUTH.recording_s      = 1800;
GNDTRUTH.sig_range        = [15 30; 15 20; 5 10];

for pp = 1:numel(PRE_POST_WIN)
  fprintf('Running Articial ground truth simulation with %d fitting cycles ...\n', PRE_POST_WIN(pp));

  GNDTRUTH.nodes_per_second = NODES_PER_S;    
  [T2, ~, ~, tdtextr] = mufix.sim_artificial(TDTMAT_TEMPLATE, 20, STIM, GNDTRUTH, CORR, XTALK_METHOD, PRE_POST_WIN(pp));

  % Save output
  currtime = datetime();
  T2.TimeStamp(:) = string(currtime,'yyyyMMdd_HHmmss');
  filename        = sprintf('%s_%gcycle_sim_dynamics', ...
    char(currtime,'yyyyMMdd'),125);
  
  save(filename, 'T2', 'tdtextr', 'GNDTRUTH', 'STIM', 'CORR');
end


