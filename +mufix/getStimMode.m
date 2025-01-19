function [modestr, S] = getStimMode(mat)

  S = struct();
  modestr = '';
  if ~isfield(mat,'epocs') || ~isfield(mat.epocs,'Brst') || ~isfield(mat.epocs,'Pls_')
    return;
  end

  S.burststart  = round(mat.epocs.Brst.onset(1));
  S.nbursts     = numel(mat.epocs.Brst.offset);  % completed bursts
  S.burstdur    = mat.epocs.Brst.offset(1:numel(mat.epocs.Brst.onset)) - mat.epocs.Brst.onset;
  S.burstdur    = median(round(S.burstdur));
  S.burstperiod = median(round(diff(mat.epocs.Brst.onset)));
  S.pulsedur    = mat.epocs.Pls_.offset(1:numel(mat.epocs.Pls_.onset)) - mat.epocs.Pls_.onset;
  S.pulsedur    = median(round(1000*S.pulsedur));
  S.pulsefreq   = median(round(1./diff(mat.epocs.Pls_.onset)));

  modestr = sprintf('B%d_P%dF%d_I%d', ...
    S.burstdur, S.pulsedur, S.pulsefreq, S.burstperiod);

end