function [T2, SErr, SErr_Time, tdtextr, carrier_out, dmoddat_out] = sim_artificial(tdt_template, recording_trials, STIM, GNDTRUTH, CORR, XTALK_METHOD, pre_post_win)
  
  REVIEW_PLOT  = 0;

  if ~exist('pre_post_win', 'var')
    pre_post_win = 2.5;
  end

  %% Load template and configure stimulus

  tdtextr = load(tdt_template);
  fs      = tdtextr.streams.GCa2.fs; % Sampling frequency of photometry output
  fs_mod  = tdtextr.streams.Fi2r.fs; % Sampling frequency of modulated signal
  tdtextr = mufix.setup_tdtmat_stimulation(tdtextr, STIM.sim_bursts, STIM.init_delay_s + GNDTRUTH.padding_s/2, STIM.burst_interval_s, STIM.burst_len_s, STIM.pulse_hz, STIM.pulse_on_s, fs_mod);

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
  

  %%
  for recording = 1:recording_trials
  
      fprintf('Simulation Run #%d ...\n', recording);

      % Generate ground truth
      % [t, ~, ch1_pre_pre, ch2_pre_pre, ch3_pre_pre] = simulate_recording(fs,fs_mod,f1,f2,f3,arange,brange,crange,numSamples, p);
      [gndtruth, t] = mufix.generate_groundtruth(3,GNDTRUTH.recording_s + GNDTRUTH.padding_s,fs,GNDTRUTH.nodes_per_second,GNDTRUTH.sig_range);
  
      % Rotate ground truth between 3 channels
      carrier_out = {};
      dmoddat_out = {};
      for rotation = 1:3
  
          % reset
          dmodtime  = struct();
          dmoddat   = struct();
          interpdat = struct();
  
          % Analysis performed 3 times swapping carrier frequencies
          if rotation == 1
              rotidx = [1 2 3];
          elseif rotation == 2
              rotidx = [3 1 2];
          elseif rotation == 3
              rotidx = [2 3 1];
          end
  
          % Simulate LIA modulation
          [carrier, t_super, ch_mod, ch_super] = mufix.simulate_LIA_modulation(gndtruth(rotidx,:)/1000, t, PHMAP.Freq__Hz', fs_mod);
  
          % Demodulate original
          tdtmat = tdtextr;
          tdtmat.streams.Fi2r.data = carrier;
          cfg.xtalk_method = "none";
          cfg.lowpass_cutoff__Hz = max(3,GNDTRUTH.nodes_per_second);

          tic;
          tdtmat = mufix.rmCrosstalk_tdtmat(tdtmat, PHMAP, cfg);
          dmodtime.original = toc;        
          dmoddat.original  = cat(1,tdtmat.streams.RCa2.data, tdtmat.streams.GCa2.data, tdtmat.streams.IGa2.data);
  
  
          %% Demodulate with different approaches
  
          % SIMULATE CROSSTALK
          [carrierX, xtalk_cfg] = mufix.simulate_crosstalk(tdtmat, carrier, fs_mod);
          tdtmat.streams.Fi2r.data = carrierX;
          interpdat.xtalk_tframe = arrayfun(@(x,y) x:y, xtalk_cfg.xtalk_on, xtalk_cfg.xtalk_off, 'UniformOutput', false);
  
          for method = 1:numel(XTALK_METHOD)
              cfg.xtalk_method   = XTALK_METHOD(method);
              cfg.crop_pre__sec  = pre_post_win / min(PHMAP.Freq__Hz);
              cfg.crop_post__sec = pre_post_win / min(PHMAP.Freq__Hz);
              tic;
              [tdtmat,~,tmpdat] = mufix.rmCrosstalk_tdtmat(tdtmat, PHMAP, cfg);
              dmodtime.(cfg.xtalk_method)  = toc;        
              dmoddat.(cfg.xtalk_method)   = cat(1,tdtmat.streams.RCa2.data, tdtmat.streams.GCa2.data, tdtmat.streams.IGa2.data);
              if isfield(tmpdat{1},'interpdat')
                interpdat.(cfg.xtalk_method)           = tmpdat{1}.interpdat;
                interpdat.(cfg.xtalk_method+"_tframe") = tmpdat{1}.interpcfg.tframe;
             end
          end
          carrier_out{rotation} = interpdat;
          carrier_out{rotation}.original = carrier;
          carrier_out{rotation}.xtalk    = carrierX;
          dmoddat_out{rotation} = dmoddat;
          
          %% Plot dmod output
          if REVIEW_PLOT
              clf(revfig);
              T = tiledlayout(revfig, size(gndtruth,1), 1);
              title(T,sprintf('Recording %d, Rotation %d', recording, rotation));
              revax = [];          
              for ii = 1:size(gndtruth,1)
                  revax(ii) = nexttile(T);
                  hold(revax,'on');
                  METHOD_ = fieldnames(dmoddat);
                  for method = 1:numel(METHOD_)
                    plot(revax(ii), t(1:end-1), dmoddat.(METHOD_{method})(ii,:));
                  end
                  plot(revax(ii), t, gndtruth(rotidx(ii),:),'k')
                  title("Ch"+char('A'+ii-1))
              end
              legend(revax(ii),[METHOD_;'Ground Truth'],'Location','NorthEastOutside')
              linkaxes(revax,'x')
              pause
          end
  
  
          %% Calculate correlations
          for burst = 1:STIM.sim_bursts
              
              %% Organize output as a struct
              trial_dat = struct();
  
              % Save trial identification data
              trial_dat.ID      = strcat("SIM_REC_", num2str(recording), "_TRIAL_", num2str(burst));
              trial_dat.rotID   = rotation;
              trial_dat.epoch   = burst;
              trial_dat.sigA    = PHMAP.StoreName(rotidx(1));
              trial_dat.sigB    = PHMAP.StoreName(rotidx(2));
              trial_dat.sigC    = PHMAP.StoreName(rotidx(3));
              trial_dat.frqA    = PHMAP.Freq__Hz(1);
              trial_dat.frqB    = PHMAP.Freq__Hz(2);
              trial_dat.frqC    = PHMAP.Freq__Hz(3);
  
              trial_dat.burst_interval_s = STIM.burst_interval_s;
              trial_dat.burst_len_s      = STIM.burst_len_s;
              trial_dat.pulse_hz         = STIM.pulse_hz;
              trial_dat.pulse_on_s       = STIM.pulse_on_s;
              trial_dat.nodes_per_s      = GNDTRUTH.nodes_per_second;
              trial_dat.pre_post_win     = pre_post_win;
  
              err_dat = trial_dat();
  
              % Timing data
              METHOD_ = fieldnames(dmoddat);
              for method = 1:numel(METHOD_)
                  mm = METHOD_{method};
                  trial_dat.("Time_"+mm) = dmodtime.(mm);
              end     
              
              % Carrier correlation values
              compare_start = round(tdtmat.epocs.Brst.onset (burst)*fs_mod + CORR.extraEp_s(1)*fs_mod);
              compare_end   = round(tdtmat.epocs.Brst.offset(burst)*fs_mod + CORR.extraEp_s(2)*fs_mod);            
              carr_gt       = carrier(compare_start:compare_end);
              for method = 1:numel(XTALK_METHOD)
                  mm = XTALK_METHOD(method);
                  if mm == "none"
                    carrX = carrierX(compare_start:compare_end);
                  else
                    carrX = interpdat.(mm)(compare_start:compare_end);
                  end
                  trial_dat.("Carr_"+mm) = corr(carr_gt',carrX');
  
                  err2 = (carr_gt-carrX).^2;
                  fitmask = err2 ~= 0;
                  trial_dat.("R2_"+mm) = 1 - mean(err2(fitmask)) / var(carr_gt);
              end     
  
              % Demod Correlation values
              compare_start = round(tdtmat.epocs.Brst.onset (burst)*fs + CORR.extraEp_s(1)*fs);
              compare_end   = round(tdtmat.epocs.Brst.offset(burst)*fs + CORR.extraEp_s(2)*fs);            
              gt_short      = gndtruth(rotidx,compare_start:compare_end);
              for method = 1:numel(METHOD_)
                  mm = METHOD_{method};
                  for ch = 1:size(gt_short,1)
                    trial_dat.("Corr_"+char(ch+'A'-1)+"_"+mm) = corr(gt_short(ch,:)', dmoddat.(mm)(ch,compare_start:compare_end)');
  
                    varGT = var(gt_short(ch,:));
                    err_dat.("Err2_"+char(ch+'A'-1)+"_"+mm)   = (gt_short(ch,:) - dmoddat.(mm)(ch,compare_start:compare_end)).^2 / varGT;
                    trial_dat.("MSE_"+char(ch+'A'-1)+"_"+mm)  = mean(err_dat.("Err2_"+char(ch+'A'-1)+"_"+mm));
                  end
              end
              
              % Pulse by pulse fits
              for method = 1:numel(XTALK_METHOD)

                  mm = XTALK_METHOD(method);
                  if mm == "none"
                    tframe = interpdat.xtalk_tframe;
                  else
                    tframe = interpdat.(cfg.xtalk_method+"_tframe");
                  end

                  for ff = 1:numel(tframe)
                      if mm == "none"
                        carrX = carrierX(tframe{ff});
                      else
                        carrX = interpdat.(mm)(tframe{ff});
                      end
                      carr_gt = carrier(tframe{ff});
      
                      err2    = (carr_gt-carrX).^2;
                      fitmask = err2 ~= 0;
                      mse     = mean(err2(fitmask));
                      trial_dat.("MSE_Pulse_"+mm)(ff) = mse;
                      trial_dat.("R2_Pulse_"+mm)(ff)  = 1 - mse / var(carr_gt);
                  end
              end                           

              T2{end+1} = trial_dat;
              SErr{end+1} = err_dat;
              SErr_Time = t(compare_start:compare_end) - tdtmat.epocs.Brst.onset(burst);
              
          end          
       
      end
      fprintf("recording %g done\n", recording)
  end
  

  %% Convert structs to a table
  T2 = cat(1,T2{:});
  T2 = struct2table(T2);
  SErr = cat(1,SErr{:});
  SErr = struct2table(SErr);
  
  
  %% Calculate Fisher Z transform
  T2_Fields  = T2.Properties.VariableNames;
  corrFields = T2_Fields(startsWith(T2_Fields,'Corr'));
  zFields    = strrep(corrFields,'Corr','Z');
  T2(:,zFields) = atanh(T2(:,corrFields));
  
  corrFields = T2_Fields(startsWith(T2_Fields,'Carr'));
  zFields    = strrep(corrFields,'Carr','CarZ');
  T2(:,zFields) = atanh(T2(:,corrFields));
  
  
end