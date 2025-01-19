function [T2, SErr, SErr_Time] = sim_empirical_fft(RECORDINGS, STIM, FREQ, CORR, XTALK_METHOD, SZ_LIST)

    ffts_rca = [];
    ffts_gca = [];
    ffts_isb = [];
    
    tic
    for recording_index = 1:length(RECORDINGS)
    
        filename = strcat(RECORDINGS(recording_index), '.mat');
        fprintf('Processing "%s" ... \n', filename)
        tdtmat  = load(filename);
    
        map = mufix.guessPhMap(tdtmat);
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
        [tdtmat, ~] = mufix.rmCrosstalk_tdtmat(tdtmat, map, cfg);
        gndtruth = arrayfun(@(x) tdtmat.streams.(x).data, map.StoreName, 'UniformOutput', false);
        gndtruth = cat(1, gndtruth{:});
        fs       = tdtmat.streams.(map.StoreName(1)).fs; % Sampling frequency in Hz
    
        %% Gather experiment info
    
        SZINFO    = strsplit(extractAfter(RECORDINGS(recording_index),filesep),"-");
        if endsWith(SZINFO{1},'A')  
          SZINFO{1} = extractBefore(SZINFO{1},'A');
        end
        SZPREFIX  = strcat("SZ_", SZINFO{1}, "_D", SZINFO{2}, "-", SZINFO{3}(1:4), "_T");
        sz_list   = SZ_LIST(startsWith(SZ_LIST.SeizureID, SZPREFIX),:);
        sz_list.Properties.RowNames = string(sz_list.Epoch);    
    
        for burst = 1:length(tdtmat.epocs.Brst.onset)

            index = sz_list.Epoch == burst;
            if sz_list(index, :).RacineScore <= 1
                continue;
            end

            compare_start = round(tdtmat.epocs.Brst.onset(burst)*fs + CORR.extraEp_s(1)*fs);
            compare_end   = round(tdtmat.epocs.Brst.offset(burst)*fs + CORR.extraEp_s(2)*fs);            
            gt_short_rca      = gndtruth(1,compare_start:compare_end);
            gt_short_gca      = gndtruth(2,compare_start:compare_end);
            gt_short_isb      = gndtruth(3,compare_start:compare_end);

            [ffttemp, f] = pwelch((gt_short_rca - mean(gt_short_rca)), fs, [], [], fs);
            ffttemp = [ffttemp' zeros([1 (2^18 - length(ffttemp))])];
            f = [f' zeros([1 (2^18 - length(f))])];
            ffts_rca = [ffts_rca; ffttemp];

            [ffttemp, f] = pwelch((gt_short_gca - mean(gt_short_gca)), fs, [], [], fs);
            ffttemp = [ffttemp' zeros([1 (2^18 - length(ffttemp))])];
            f = [f' zeros([1 (2^18 - length(f))])];
            ffts_gca = [ffts_gca; ffttemp];

            [ffttemp, f] = pwelch((gt_short_isb - mean(gt_short_isb)), fs, [], [], fs);
            ffttemp = [ffttemp' zeros([1 (2^18 - length(ffttemp))])];
            f = [f' zeros([1 (2^18 - length(f))])];
            ffts_isb = [ffts_isb; ffttemp];
        end
    end
    % ffts_rca = ffts_rca/length(RECORDINGS);
    % figure;
    % hold on;
    % plot(f,10*log10(ffts_rca));
    % xlim([0 10]);
end