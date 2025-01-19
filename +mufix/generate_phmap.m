function phmap = generate_phmap(signames, lia_freqs__hz)
%%
%   phmap = generate_phmap(template_name)
%   phmap = generate_phmap(signames, lia_freqs__hz)
%
% Generate photometry configuration table as produced by sunlab.tdt.guessPhMap
% for the purpose of simulating LIA demodulation.


%% Constants
PHMAP_TITLE    = {'StoreName','SourceStore','Channel','SeqID','Freq__Hz','PowerOn__sec','PowerOff__sec','ACPower__mA','DCPower__mA'};
POWER_ON__sec  = {5};
POWER_OFF__sec = {0};
ACPOWER__mA    = {100};
DCPOWER__mA    = {10};

% Defaults 
PowerOn__sec  = [];
PowerOff__sec = [];
ACPower__mA   = [];
DCPower__mA   = [];


%% Check template
template = [];
if isstring(signames) && numel(signames) == 1
  switch signames
    case "3ph"
      template  = "3ph";
      signames  = ["RCa2", "GCa2", "IGa2"];
      sigsource = repmat("Fi2r",1,3);      
      channel   = [1 1 1];
      seqid     = [1 2 3];
      lia_freqs__hz = [210 330 530];
  end
end

%% Custom setup
if isempty(template)
  if numel(signames) ~= numel(lia_freqs__hz)
    error('Number of carrier frequencys does not match the number of signals.');
  end
  sigsource = repmat("Fi2r",1,numel(signames));      
  channel   = ones(1,numel(signames));
  seqid     = 1:numel(signames);  
end

%% Fill in the rest
nsig = numel(signames);
if isempty(PowerOn__sec)
  PowerOn__sec = repmat(POWER_ON__sec,1,nsig);
end
if isempty(PowerOff__sec)
  PowerOff__sec = repmat(POWER_OFF__sec,1,nsig);
end
if isempty(ACPower__mA)
  ACPower__mA = repmat(ACPOWER__mA,1,nsig);
end
if isempty(DCPower__mA)
  DCPower__mA = repmat(DCPOWER__mA,1,nsig);
end


%% Create table
phmap = table(signames(:), sigsource(:), channel(:), seqid(:), lia_freqs__hz(:), PowerOn__sec(:), PowerOff__sec(:), ACPower__mA(:), DCPower__mA(:), ...
  'VariableNames', PHMAP_TITLE, 'RowNames', signames);



