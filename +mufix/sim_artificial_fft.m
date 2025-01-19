function [T2, SErr, SErr_Time, tdtextr, carrier_out, dmoddat_out] = sim_artificial_fft(tdt_template, recording_trials, STIM, GNDTRUTH, CORR, XTALK_METHOD)
  
  REVIEW_PLOT  = 0;

  %% Load template and configure stimulus

  tdtextr = load(tdt_template);
  fs      = tdtextr.streams.GCa2.fs; % Sampling frequency of photometry output
  fs_mod  = tdtextr.streams.Fi2r.fs; % Sampling frequency of modulated signal
  tdtmat = mufix.setup_tdtmat_stimulation(tdtextr, STIM.sim_bursts, STIM.init_delay_s + GNDTRUTH.padding_s/2, STIM.burst_interval_s, STIM.burst_len_s, STIM.pulse_hz, STIM.pulse_on_s, fs_mod);

  % Photometry signal map
  PHMAP = mufix.guessPhMap(tdtextr);

  %% Main Simulation
          
  T2 = {};
  SErr = {};
  SErr_Time = [];
  
  cfg = struct();
  cfg.plot = 0;
  cfg.dmod_scale = 1;
  cfg.lowpass_zerophase = 1;    
  
  if REVIEW_PLOT
    revfig = figure;
  end
  
  ffts = [];

  %%
  for recording = 1:recording_trials
  
      fprintf('Simulation Run #%d ...\n', recording);

      % Generate ground truth
      % [t, ~, ch1_pre_pre, ch2_pre_pre, ch3_pre_pre] = simulate_recording(fs,fs_mod,f1,f2,f3,arange,brange,crange,numSamples, p);
      [gndtruth, t] = mufix.generate_groundtruth(3,GNDTRUTH.recording_s + GNDTRUTH.padding_s,fs,GNDTRUTH.nodes_per_second,GNDTRUTH.sig_range);
  
  
      %% Calculate correlations
      for burst = 1:STIM.sim_bursts

        compare_start = round(tdtmat.epocs.Brst.onset(burst)*fs + CORR.extraEp_s(1)*fs);
        compare_end   = round(tdtmat.epocs.Brst.offset(burst)*fs + CORR.extraEp_s(2)*fs);            
        gt_short_rca      = gndtruth(1,compare_start:compare_end);
        gt_short_gca      = gndtruth(2,compare_start:compare_end);
        gt_short_isb      = gndtruth(3,compare_start:compare_end);

        [ffttemp, f] = pwelch((gt_short_rca - mean(gt_short_rca)), fs, [], [], fs);
        ffttemp = [ffttemp' zeros([1 (2^18 - length(ffttemp))])];
        f = [f' zeros([1 (2^18 - length(f))])];
        ffts = [ffts; ffttemp];

        [ffttemp, f] = pwelch((gt_short_gca - mean(gt_short_gca)), fs, [], [], fs);
        ffttemp = [ffttemp' zeros([1 (2^18 - length(ffttemp))])];
        f = [f' zeros([1 (2^18 - length(f))])];
        ffts = [ffts; ffttemp];

        [ffttemp, f] = pwelch((gt_short_isb - mean(gt_short_isb)), fs, [], [], fs);
        ffttemp = [ffttemp' zeros([1 (2^18 - length(ffttemp))])];
        f = [f' zeros([1 (2^18 - length(f))])];
        ffts = [ffts; ffttemp];
      end
      fprintf("recording %g done\n", recording)
  end
end