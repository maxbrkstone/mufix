function carrier = simulate_crosstalk(tdtmat, carrier, fs_mod)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

    %ctcenter = rand*(7.5-4.5) + 4.5;
    %ctmin = ctcenter - 2.5;
    %ctmax = ctcenter + 2.5;
    %ctfreq = 230;

    %carrier(:, stim_on:stim_off) = ...
            %    (rand(1, length(stim_on:stim_off))*(ctmax-ctmin) + ctmin).*(sin(2*pi*ctfreq*t_super(stim_on:stim_off)) + 1.0)/2;

    ctcenter  = 10;
    rebound_s = 2/1000;

    % % stim_len = round((max(tdtmat.epocs.Pls_.offset - tdtmat.epocs.Pls_.onset) + rebound_s) * fs_mod);
    % % stim_on  = round(tdtmat.epocs.Pls_.onset * fs_mod);
    % % stim_off = stim_on + stim_len - 1;
    % % rebx     = get_rebound();
    % % xtalk    = ctcenter * ones(1,stim_len);
    % % xtalk(1:numel(rebx)) = fliplr(rebx);
    % % xtalk    = fliplr(xtalk);

    stim_len = round(max(tdtmat.epocs.Pls_.offset - tdtmat.epocs.Pls_.onset + 0.0005) * fs_mod);
    rebx     = get_rebound();
    xtalk    = ctcenter * tukeywin(stim_len, (fs_mod/1000)/stim_len);
    ndelay   = round(fs_mod/2000);
    xtalk    = [zeros(1,ndelay) xtalk']; %, rebx];

    stim_on  = round(tdtmat.epocs.Pls_.onset * fs_mod);
    stim_off = stim_on + numel(xtalk) - 1;

    for stim_num = 1:length(stim_on)           
        carrier(:, stim_on(stim_num):stim_off(stim_num)) = xtalk;
    end
    % carrier = min(carrier,10);

  function reboundx = get_rebound(stim_len)
      reboundx(1) = -1.0;
      reboundx(2) = -1.5;
      reboundx(3) = -0.5;
      reboundx(4) = ctcenter/8;
      reboundx(5) = ctcenter/4;
      reboundx(6) = ctcenter/3;
      n = round(rebound_s * fs_mod);
      reboundx = interp1(linspace(0,1,numel(reboundx)), reboundx, linspace(0,1,n));
  end

end