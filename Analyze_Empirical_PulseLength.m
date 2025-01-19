%% Analyze pulse duration

SRCDAT.H20{1} = load('20240701_20Hz_5ms_empirtrials.mat');
SRCDAT.H20{2} = load('20240705_20Hz_15ms_empirtrials.mat');
SRCDAT.H20{3} = load('20240705_20Hz_25ms_empirtrials.mat');
SRCDAT.H20{4} = load('20240705_20Hz_33ms_empirtrials.mat');
SRCDAT.H20{5} = load('20240705_20Hz_34ms_empirtrials.mat');
SRCDAT.H20{6} = load('20240705_20Hz_35ms_empirtrials.mat');

SRCDAT.H1{1} = load('20240705_1Hz_50ms_empirtrials.mat');
SRCDAT.H1{2} = load('20240705_1Hz_100ms_empirtrials.mat');
SRCDAT.H1{3} = load('20240705_1Hz_150ms_empirtrials.mat');
SRCDAT.H1{4} = load('20240705_1Hz_200ms_empirtrials.mat');
SRCDAT.H1{5} = load('20240705_1Hz_300ms_empirtrials.mat');
SRCDAT.H1{6} = load('20240705_1Hz_400ms_empirtrials.mat');
SRCDAT.H1{7} = load('20240705_1Hz_500ms_empirtrials.mat');
SRCDAT.H1{8} = load('20240705_1Hz_750ms_empirtrials.mat');
SRCDAT.H1{9} = load('20240705_1Hz_900ms_empirtrials.mat');


%% Combine tables

T2.H20    = cellfun(@(x) x.T2, SRCDAT.H20, 'UniformOutput', false);
T2.H1     = cellfun(@(x) x.T2, SRCDAT.H1,  'UniformOutput', false);
T2_Fields = cellfun(@(x) x.Properties.VariableNames, [T2.H20, T2.H1], 'UniformOutput', false);

% find minimal set
for ii = 2:numel(T2_Fields)
  T2_Fields{1} = intersect(T2_Fields{1},T2_Fields{ii},'stable');
end
T2_Fields = T2_Fields{1};

T2.H20 = cellfun(@(x) x(:,T2_Fields), T2.H20, 'UniformOutput', false);
T2.H20 = cat(1, T2.H20{:});
T2.H1  = cellfun(@(x) x(:,T2_Fields), T2.H1,  'UniformOutput', false);
T2.H1  = cat(1, T2.H1{:});

% Extract and Stack A B C channels
T2_Corr.H20   = {};
T2_Corr.H20{1} = T2.H20(:,{'ID','rotID','Mouse','sigA','epoch','hasSeizure','frqA','meanA','sdA','paroxtime_s','pulse_on_s','Z_A_mufix'});
T2_Corr.H20{2} = T2.H20(:,{'ID','rotID','Mouse','sigB','epoch','hasSeizure','frqB','meanB','sdB','paroxtime_s','pulse_on_s','Z_B_mufix'});
T2_Corr.H20{3} = T2.H20(:,{'ID','rotID','Mouse','sigC','epoch','hasSeizure','frqC','meanC','sdC','paroxtime_s','pulse_on_s','Z_C_mufix'});

T2_Corr.H20{2}.Properties.VariableNames = T2_Corr.H20{1}.Properties.VariableNames;
T2_Corr.H20{3}.Properties.VariableNames = T2_Corr.H20{1}.Properties.VariableNames;
T2_Corr.H20 = cat(1, T2_Corr.H20{:});
T2_Corr.H20.sigMouse = T2_Corr.H20.Mouse + "-" + T2_Corr.H20.sigA;

% Extract and Stack A B C channels
T2_Corr.H1   = {};
T2_Corr.H1{1} = T2.H1(:,{'ID','rotID','Mouse','sigA','epoch','hasSeizure','frqA','meanA','sdA','paroxtime_s','pulse_on_s','Z_A_mufix'});
T2_Corr.H1{2} = T2.H1(:,{'ID','rotID','Mouse','sigB','epoch','hasSeizure','frqB','meanB','sdB','paroxtime_s','pulse_on_s','Z_B_mufix'});
T2_Corr.H1{3} = T2.H1(:,{'ID','rotID','Mouse','sigC','epoch','hasSeizure','frqC','meanC','sdC','paroxtime_s','pulse_on_s','Z_C_mufix'});

T2_Corr.H1{2}.Properties.VariableNames = T2_Corr.H1{1}.Properties.VariableNames;
T2_Corr.H1{3}.Properties.VariableNames = T2_Corr.H1{1}.Properties.VariableNames;
T2_Corr.H1 = cat(1, T2_Corr.H1{:});
T2_Corr.H1.sigMouse = T2_Corr.H1.Mouse + "-" + T2_Corr.H1.sigA;


%% 20Hz Anova
sdmask = T2_Corr.H20.sdA > 0.3;
T2_Corr_ = T2_Corr.H20(sdmask,:);

[CorrZ_P.H20,~,CorrZ_Stats.H20] = anovan(T2_Corr_.Z_A_mufix, ...
  {T2_Corr_.frqA, T2_Corr_.pulse_on_s, log(T2_Corr_.sdA)}, ...
  'continuous', 3, ...
  'varnames', {'freq', 'pdur', 'sd'}, ...
  'model', 'interaction' ...
  );  

figure
[CorrZ_comps.H20,CorrZ_means.H20,~,CorrZ_mean_names.H20] = multcompare(CorrZ_Stats.H20,'Dimension',[1 2]);

figure
[CorrZ_comps2.H20,CorrZ_means2.H20,~,CorrZ_mean_names2.H20] = multcompare(CorrZ_Stats.H20,'Dimension',2);

figure
[CorrZ_comps1.H20,CorrZ_means1.H20,~,CorrZ_mean_names1.H20] = multcompare(CorrZ_Stats.H20,'Dimension',1, 'Alpha', 0.01);

gsum = groupsummary(T2_Corr_,'pulse_on_s',{'mean','std','median','min', 'max', @(x) quantile(x,[0.25 0.75]), @(x) quantile(x,[0.025 0.975])},'Z_A_mufix');
m    = tanh(gsum(:,{'median_Z_A_mufix','min_Z_A_mufix','max_Z_A_mufix','fun1_Z_A_mufix','fun2_Z_A_mufix'})) 

figure;
T = tiledlayout(2,1,'tileindexing','columnmajor');
title(T,'20Hz pulse width')

T2_Corr_.frqPulse = "F" + T2_Corr_.frqA + "_P" + num2str(1000*T2_Corr_.pulse_on_s,'%03d');
nexttile
swarmchart(categorical(T2_Corr_.frqPulse), tanh(T2_Corr_.Z_A_mufix), 10, T2_Corr_.frqA, 'filled');
hold on
ylabel('Correlation (\rho)');
set(gca,'TickLabelInterpreter','none')

nexttile
swarmchart(categorical(T2_Corr_.frqPulse), T2_Corr_.Z_A_mufix, 10, T2_Corr_.frqA, 'filled');
hold on
ylabel('Correlation (Z)');
set(gca,'TickLabelInterpreter','none')


%% 1Hz Anova
sdmask = T2_Corr.H1.sdA > 0; %0.3;
T2_Corr_ = T2_Corr.H1(sdmask,:);

[CorrZ_P.H1,~,CorrZ_Stats.H1] = anovan(T2_Corr_.Z_A_mufix, ...
  {T2_Corr_.frqA, T2_Corr_.pulse_on_s, log(T2_Corr_.sdA)}, ...
  'continuous', 3, ...
  'varnames', {'freq', 'pdur', 'sd'}, ...
  'model', 'interaction' ...
  );  

figure
[CorrZ_comps.H1,CorrZ_means.H1,~,CorrZ_mean_names.H1] = multcompare(CorrZ_Stats.H1,'Dimension',[1 2]);

figure
[CorrZ_comps2.H1,CorrZ_means2.H1,~,CorrZ_mean_names2.H1] = multcompare(CorrZ_Stats.H1,'Dimension',2);

figure
[CorrZ_comps1.H1,CorrZ_means1.H1,~,CorrZ_mean_names1.H1] = multcompare(CorrZ_Stats.H1,'Dimension',1, 'Alpha', 0.01);


gsum = groupsummary(T2_Corr_,'pulse_on_s',{'mean','std','median','min', 'max', @(x) quantile(x,[0.25 0.75]), @(x) quantile(x,[0.025 0.975])},'Z_A_mufix');
m    = tanh(gsum(:,{'median_Z_A_mufix','min_Z_A_mufix','max_Z_A_mufix','fun1_Z_A_mufix','fun2_Z_A_mufix'})) 

figure;
T = tiledlayout(2,1,'tileindexing','columnmajor');
title(T,'1Hz pulse width')

T2_Corr_.frqPulse = "F" + T2_Corr_.frqA + "_P" + num2str(1000*T2_Corr_.pulse_on_s,'%03d');
nexttile
swarmchart(categorical(T2_Corr_.frqPulse), tanh(T2_Corr_.Z_A_mufix), 10, T2_Corr_.frqA, 'filled');
hold on
ylabel('Correlation (\rho)');
set(gca,'TickLabelInterpreter','none')

nexttile
swarmchart(categorical(T2_Corr_.frqPulse), T2_Corr_.Z_A_mufix, 10, T2_Corr_.frqA, 'filled');
hold on
ylabel('Correlation (Z)');
set(gca,'TickLabelInterpreter','none')
