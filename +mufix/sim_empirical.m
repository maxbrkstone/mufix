function [T2, SErr, SErr_Time, tdtmat_, carrier_out_, dmoddat_out_] = sim_empirical(RECORDINGS, STIM, FREQ, CORR, XTALK_METHOD, SZ_LIST, pre_post_win)

if ~exist('pre_post_win', 'var')
  pre_post_win = 2.5;
end
  
%%
T2 = cell(length(RECORDINGS),1);
SErr = cell(length(RECORDINGS),1);
SErr_Time = cell(length(RECORDINGS),1);

% Set up parallel processing
% - resample() does not work with thread pools, have to set up process pools
% - 4 recordings: parfor=517s, single=1407s.
% - all recordings (43), parfor=2866s
if length(RECORDINGS) >= 2
  pargrp = gcp('nocreate');
  if isempty(pargrp) || isa(pargrp,'parallel.ThreadPool')
    delete(pargrp);
    parpool('Processes',8);
  end
end

tdtmat_      = cell(length(RECORDINGS),1);
carrier_out_ = cell(length(RECORDINGS),1);
dmoddat_out_ = cell(length(RECORDINGS),1);
extra_out    = length(RECORDINGS) == 1;

tic
parfor recording_index = 1:length(RECORDINGS)

    % inside each thread
    tT2 = {};
    tSErr = {};
    tSErr_Time = [];

    filename = strcat(RECORDINGS(recording_index), '.mat');
    fprintf('Processing "%s" ... \n', filename)
    tdtextr  = load(filename);

    map = mufix.guessPhMap(tdtextr);
    map = sortrows(map,'SeqID');
    if any(map.SourceStore == "Fi2r")
      map = map(["RCa2","GCa2","IGa2"],:);
    else
      map = map(["RCa1","GCa1","Isb1"],:);
    end
    phscalar = char(map.SourceStore(1));
    phscalar = [phscalar(1:3),'i'];

    cfg = struct();
    cfg.plot = 0;
    cfg.lowpass_zerophase = 1;
    
    % Extract the demodulated as the ground truth
    cfg.xtalk_method = "mufix";
    [tdtextr, ~] = mufix.rmCrosstalk_tdtmat(tdtextr, map, cfg);
    gndtruth = arrayfun(@(x) tdtextr.streams.(x).data, map.StoreName, 'UniformOutput', false);
    gndtruth = cat(1, gndtruth{:});
    fs       = tdtextr.streams.(map.StoreName(1)).fs; % Sampling frequency in Hz
    fs_mod   = tdtextr.streams.(map.SourceStore(1)).fs; % Sampling frequency in Hz
    

    % Modify for simulation
    tdtmat = tdtextr;
    tdtmat = mufix.setup_tdtmat_stimulation(tdtmat, [], [], [], [], STIM.pulse_hz, STIM.pulse_on_s, fs_mod);
    tdtmat.streams.(map.SourceStore(1)).data = tdtmat.streams.(map.SourceStore(1)).data(1,:);
    tdtmat.scalars.(phscalar).data(2,:)  = FREQ.f1;
    tdtmat.scalars.(phscalar).data(6,:)  = FREQ.f2;
    tdtmat.scalars.(phscalar).data(10,:) = FREQ.f3;
    map_sim = mufix.guessPhMap(tdtmat);
    if any(map_sim.SourceStore == "Fi2r")
      map_sim = map_sim(["RCa2","GCa2","IGa2"],:);
    else
      map_sim = map_sim(["RCa1","GCa1","Isb1"],:);
    end

    % Save tdtmat for debugging, but remove data to save memory
    if extra_out
      tdtmat_{recording_index} = tdtmat;
      fstreams = fieldnames(tdtmat_{recording_index}.streams);
      for ff = 1:numel(fstreams)
        tdtmat_{recording_index}.streams.(fstreams{ff}).data = [];
      end
    end

    %% Gather experiment info

    SZINFO    = strsplit(extractAfter(RECORDINGS(recording_index),filesep),"-");
    if endsWith(SZINFO{1},'A')
      SZINFO{1} = extractBefore(SZINFO{1},'A');
    end
    SZPREFIX  = strcat("SZ_", SZINFO{1}, "_D", SZINFO{2}, "-", SZINFO{3}(1:4), "_T");
    sz_list   = SZ_LIST(startsWith(SZ_LIST.SeizureID, SZPREFIX),:);
    sz_list.Properties.RowNames = string(sz_list.Epoch);    

    
    %%
    rotidx = [];
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
        t = (0:(size(gndtruth,2)-1)) / fs; % Time vector    
        [carrier, t_super, ch_mod, ch_super] = mufix.simulate_LIA_modulation(gndtruth(rotidx,:)/1000, t, map_sim.Freq__Hz', fs_mod);

        % Demodulate original
        tdtmat.streams.(map.SourceStore(1)).data = carrier;
        cfg.dmod_scale = 1;
        cfg.xtalk_method = "none";
        %   crop_pre__sec  - Amount of data to crop before tseg for fitting (seconds)
        %                    (default: 2.5 cycles of the lowest frqs)
        %   crop_post__sec - Amount of data to crop after tseg for fitting (seconds)
        %                    (default: 2.5 cycles of the lowest frqs)
        cfg.crop_pre__sec  = pre_post_win / min(map.Freq__Hz);
        cfg.crop_post__sec = pre_post_win / min(map.Freq__Hz);

        tdtmat = mufix.rmCrosstalk_tdtmat(tdtmat, map_sim, cfg);
        dmoddat.original = arrayfun(@(x) tdtmat.streams.(x).data, map.StoreName, 'UniformOutput', false);
        dmoddat.original = cat(1,dmoddat.original{:});


        %% Demodulate with different approaches

        % SIMULATE CROSSTALK
        if STIM.burst_len_s ~= 30
            for burst = 1:tdtmat.epocs.Brst.size    
                % insert_on = tdtmat.epocs.Pls_.onset(pre_end) + 0.05 : 0.05 : tdtmat.epocs.Pls_.onset(pre_end) + 0.05 * (598 + 599*ninetysec_stim);
                % insert_off = tdtmat.epocs.Pls_.offset(pre_end) + 0.05 : 0.05 : tdtmat.epocs.Pls_.offset(pre_end) + 0.05 * (598 + 599*ninetysec_stim);
                % 
                % % if burst ~= tdtmat.epocs.Brst.size
                % %     post_start = find(tdtmat.epocs.Pls_.onset == tdtmat.epocs.Brst.onset(burst+1));
                % % 
                % %     tdtmat.epocs.Pls_.onset = [tdtmat.epocs.Pls_.onset(1:pre_end); insert_on'; tdtmat.epocs.Pls_.onset(post_start:end);];
                % %     tdtmat.epocs.Pls_.offset = [tdtmat.epocs.Pls_.offset(1:pre_end); insert_off'; tdtmat.epocs.Pls_.offset(post_start:end);];
                % % else
                % %     tdtmat.epocs.Pls_.onset = [tdtmat.epocs.Pls_.onset(1:pre_end); insert_on';];
                % %     tdtmat.epocs.Pls_.offset = [tdtmat.epocs.Pls_.offset(1:pre_end); insert_off';];
                % % end
                tdtmat.epocs.Pls_.onset = zeros(1,(STIM.pulse_hz*STIM.burst_len_s));
                tdtmat.epocs.Pls_.offset = zeros(1,(STIM.pulse_hz*STIM.burst_len_s));
                
                for pulse = 0:(length(tdtmat.epocs.Pls_.onset) - 1)
                    tdtmat.epocs.Pls_.onset(pulse + STIM.pulse_on_s*STIM.pulse_hz*(burst-1)) = tdtmat.epocs.Brst.onset(burst) + pulse*STIM.pulse_on_s + pulse*(1/STIM.pulse_hz - STIM.pulse_on_s);
                    tdtmat.epocs.Pls_.offset(pulse + STIM.pulse_on_s*STIM.pulse_hz*((burst-1)) = tdtmat.epocs.Brst.onset(burst) + (pulse + 1)*STIM.pulse_on_s;
                end

                tdtmat.epocs.Brst.offset(burst) =  tdt(end);
            end
            tdtextr.epocs.Brst.onset = tdtmat.epocs.Brst.onset;
            tdtextr.epocs.Brst.offset = tdtmat.epocs.Brst.offset;
            tdtextr.epocs.Pls_.onset = tdtmat.epocs.Pls_.onset;
            tdtextr.epocs.Pls_.offset = tdtmat.epocs.Pls_.offset;
        end

        [carrierX, xtalk_cfg]  = mufix.simulate_crosstalk(tdtmat, carrier, fs_mod, STIM);
        tdtmat.streams.(map.SourceStore(1)).data = carrierX;
        interpdat.xtalk_tframe = arrayfun(@(x,y) x:y, xtalk_cfg.xtalk_on, xtalk_cfg.xtalk_off, 'UniformOutput', false);

        for method = 1:numel(XTALK_METHOD)
            cfg.xtalk_method = XTALK_METHOD(method);
            [tdtmat,~,tmpdat] = mufix.rmCrosstalk_tdtmat(tdtmat, map_sim, cfg);
            dmoddat.(cfg.xtalk_method)  = arrayfun(@(x) tdtmat.streams.(x).data, map.StoreName, 'UniformOutput', false);
            dmoddat.(cfg.xtalk_method)  = cat(1,dmoddat.((cfg.xtalk_method) ){:});
            if isfield(tmpdat{1},'interpdat')
                interpdat.(cfg.xtalk_method) = tmpdat{1}.interpdat;
                interpdat.(cfg.xtalk_method+"_tframe") = tmpdat{1}.interpcfg.tframe;
            end
        end
        carrier_out{rotation} = interpdat;
        carrier_out{rotation}.original = carrier;
        carrier_out{rotation}.xtalk    = carrierX;
        dmoddat_out{rotation} = dmoddat;

                
        [~,stim] = mufix.getStimMode(tdtextr);

        % Calculate demodulated correlations
        for burst = 1:length(tdtmat.epocs.Brst.onset)

            trial_dat = struct();

            % Save trial identification data
            trial_dat.ID      = strcat(extractAfter(RECORDINGS(recording_index),filesep), "_EPOCH_", num2str(burst));
            trial_dat.Mouse   = sz_list.Mouse(1);
            trial_dat.epoch   = burst;
            trial_dat.rotID   = rotation;
            trial_dat.sigA    = map_sim.StoreName(rotidx(1));
            trial_dat.sigB    = map_sim.StoreName(rotidx(2));
            trial_dat.sigC    = map_sim.StoreName(rotidx(3));
            trial_dat.frqA    = map_sim.Freq__Hz(1);
            trial_dat.frqB    = map_sim.Freq__Hz(2);
            trial_dat.frqC    = map_sim.Freq__Hz(3);

            trial_dat.original_frqA    = map.Freq__Hz(1);
            trial_dat.original_frqB    = map.Freq__Hz(2);
            trial_dat.original_frqC    = map.Freq__Hz(3);

            trial_dat.burst_interval_s = stim.burstperiod;
            trial_dat.burst_len_s      = stim.burstdur;
            trial_dat.pulse_hz         = STIM.pulse_hz;
            trial_dat.pulse_on_s       = STIM.pulse_on_s;
            trial_dat.original_pulse_hz   = stim.pulsefreq;
            trial_dat.original_pulse_on_s = stim.pulsedur/1000;

            % Seizure information
            % - Get seizure type / epoch type
            if ismember(burst,sz_list.Epoch)
              szid = string(burst);
              trial_dat.SeizureID   = sz_list.SeizureID(szid);
              trial_dat.paroxtime_s = sz_list.ParoxysmalStart__sec(szid) - sz_list.EpochStart__sec(szid);
              trial_dat.hasSeizure  = true;
            else
              trial_dat.SeizureID   = "";
              trial_dat.paroxtime_s = nan;
              trial_dat.hasSeizure  = false;
            end            

            err_dat = trial_dat();

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

            %get original nps value
            %trial_dat.mean_freq = meanfreq(gt_short(1,:), fs);

            trial_dat.meanA  = mean(gt_short(1,:));
            trial_dat.meanB  = mean(gt_short(2,:));
            trial_dat.meanC  = mean(gt_short(3,:));
            trial_dat.sdA    = std(gt_short(1,:));
            trial_dat.sdB    = std(gt_short(2,:));
            trial_dat.sdC    = std(gt_short(3,:));

            METHODS_ = fieldnames(dmoddat);
            for method = 1:numel(METHODS_)
                mm = METHODS_{method};
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
            
            tT2{end+1} = trial_dat;
            tSErr{end+1} = err_dat;
            tSErr_Time = t(compare_start:compare_end) - tdtmat.epocs.Brst.onset(burst);
        end
    end

    % thread output assignments
    T2{recording_index} = cat(1,tT2{:});
    SErr{recording_index} = cat(1,tSErr{:});
    SErr_Time{recording_index} = tSErr_Time;

    if extra_out
      carrier_out_{recording_index} = carrier_out;
      dmoddat_out_{recording_index} = dmoddat_out;
    end
    
end

%Tidy up
T2 = cat(1, T2{:});
T2 = struct2table(T2);
SErr = cat(1,SErr{:});
SErr = struct2table(SErr);
SErr_Time = SErr_Time{1};

toc

if extra_out
  tdtmat_      = tdtmat_{1};
  carrier_out_ = carrier_out_{1};
  dmoddat_out_ = dmoddat_out_{1};
end


%% Calculate Fisher Z transform
T2_Fields  = T2.Properties.VariableNames;
corrFields = T2_Fields(startsWith(T2_Fields,'Corr'));
zFields    = strrep(corrFields,'Corr','Z');
T2(:,zFields) = atanh(T2(:,corrFields));

corrFields = T2_Fields(startsWith(T2_Fields,'Carr'));
zFields    = strrep(corrFields,'Carr','CarZ');
T2(:,zFields) = atanh(T2(:,corrFields));


end