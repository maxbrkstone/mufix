%% Analyze pulse duration

SRCDAT{1} = load('20240707_1nps_sim_dynamics.mat');
SRCDAT{2} = load('20240707_method_trials.mat');
SRCDAT{3} = load('20240706_6nps_sim_dynamics.mat');
SRCDAT{4} = load('20240706_10nps_sim_dynamics.mat');
SRCDAT{5} = load('20240706_20nps_sim_dynamics.mat');
SRCDAT{6} = load('20240707_30nps_sim_dynamics.mat');
SRCDAT{7} = load('20240707_40nps_sim_dynamics.mat');
SRCDAT{8} = load('20240706_50nps_sim_dynamics.mat');
SRCDAT{9} = load('20240707_60nps_sim_dynamics.mat');
SRCDAT{10} = load('20240707_70nps_sim_dynamics.mat');
SRCDAT{11} = load('20240707_80nps_sim_dynamics.mat');
SRCDAT{12} = load('20240707_90nps_sim_dynamics.mat');
SRCDAT{13} = load('20240707_100nps_sim_dynamics.mat');
SRCDAT{14} = load('20240707_120nps_sim_dynamics.mat');


%% Combine tables

T2 = cellfun(@(x) x.T2, SRCDAT, 'UniformOutput', false);


% find minimal set
T2_Fields = cellfun(@(x) x.Properties.VariableNames, T2, 'UniformOutput', false);
for ii = 1:numel(T2_Fields)
  T2_Fields{1} = intersect(T2_Fields{1},T2_Fields{ii},'stable');
end
T2_Fields = T2_Fields{1};

T2 = cellfun(@(x) x(:,T2_Fields), T2, 'UniformOutput', false);
T2 = cat(1, T2{:});

% Extract and Stack A B C channels
T2_Corr   = {};
T2_Corr{1} = T2(:,{'ID','rotID','sigA','epoch','nodes_per_s','frqA','Z_A_original','Z_A_mufix'});
T2_Corr{2} = T2(:,{'ID','rotID','sigB','epoch','nodes_per_s','frqB','Z_B_original','Z_B_mufix'});
T2_Corr{3} = T2(:,{'ID','rotID','sigC','epoch','nodes_per_s','frqC','Z_C_original','Z_C_mufix'});

T2_Corr{2}.Properties.VariableNames = T2_Corr{1}.Properties.VariableNames;
T2_Corr{3}.Properties.VariableNames = T2_Corr{1}.Properties.VariableNames;
T2_Corr = cat(1, T2_Corr{:});

% Stack by method
T2_Corr = stack(T2_Corr,T2_Fields(startsWith(T2_Fields,'Z_A_')), ...
  'IndexVariableName', 'Method', 'NewDataVariableName', 'CorrZ');


%% Anova

[CorrZ_P,~,CorrZ_Stats] = anovan(T2_Corr.CorrZ, ...
  {T2_Corr.frqA, T2_Corr.sigA, T2_Corr.Method, T2_Corr.nodes_per_s}, ...
  'varnames', {'freq', 'sig', 'Method', 'dyn'}, ...
  'model', 'interaction' ...
  );  

figure
[CorrZ_comps,CorrZ_means,~,CorrZ_mean_names] = multcompare(CorrZ_Stats,'Dimension',[1 3 4]);

figure
multcompare(CorrZ_Stats,'Dimension',[3 4]);

gsum = groupsummary(T2_Corr,{'nodes_per_s','Method'},{'mean','std','median','min', 'max', @(x) quantile(x,[0.25 0.75]), @(x) quantile(x,[0.025 0.975])},'CorrZ');
m    = tanh(gsum(:,{'median_CorrZ','min_CorrZ','max_CorrZ','fun1_CorrZ','fun2_CorrZ'})) ;
[gsum(:,1:2) m]

gsum_frq = groupsummary(T2_Corr,{'nodes_per_s','Method','frqA'},{'mean','std','median','min', 'max', @(x) quantile(x,[0.25 0.75]), @(x) quantile(x,[0.025 0.975])},'CorrZ');
m_frq    = tanh(gsum_frq(:,{'median_CorrZ','min_CorrZ','max_CorrZ','fun1_CorrZ','fun2_CorrZ'})) ;
[gsum_frq(:,1:3) m_frq]

figure;
T = tiledlayout(2,1,'tileindexing','columnmajor');
title(T,'Signal Dynamics')

T2_Corr.dynfrq = "N" + num2str(T2_Corr.nodes_per_s,'%02d') + "_F" + T2_Corr.frqA + "_" + char(T2_Corr.Method);
nexttile
swarmchart(categorical(T2_Corr.dynfrq), tanh(T2_Corr.CorrZ), 10, T2_Corr.frqA, 'filled');
hold on
ylabel('Correlation (\rho)');
set(gca,'TickLabelInterpreter','none')

nexttile
swarmchart(categorical(T2_Corr.dynfrq), T2_Corr.CorrZ, 10, T2_Corr.frqA, 'filled');
hold on
ylabel('Correlation (Z)');
set(gca,'TickLabelInterpreter','none')

