%% Analyze artificial simulation results

load('20240707_method_trials.mat');


%% Timing analysis

% Extract timing data
T2_Fields = T2.Properties.VariableNames;
T2_Timing = T2(:,['ID','rotID',T2_Fields(startsWith(T2_Fields,'Time_'))]);
T2_Timing = stack(T2_Timing, T2_Fields(startsWith(T2_Fields,'Time_')), ...
  'IndexVariableName', 'Method', 'NewDataVariableName', 'Time');

% Anova
[Timing_P,~,Timing_Stats] = anovan(T2_Timing.Time, ...
  {categorical(T2_Timing.ID), T2_Timing.rotID, T2_Timing.Method}, ...
  'varnames', {'TrialID', 'rotID', 'Method'});
figure
[Timing_comps,Timing_means,~,Timing_mean_names] = multcompare(Timing_Stats,'Dimension',3);

gsum = groupsummary(T2_Timing,'Method',{'mean','std','median','min', 'max', @(x) quantile(x,[0.25 0.75]), @(x) quantile(x,[0.025 0.975])},'Time');
m  = 1000*gsum.mean_Time / (STIM.pulse_hz*STIM.burst_len_s)
sd = 1000*gsum.std_Time  / (STIM.pulse_hz*STIM.burst_len_s)
dm = gsum.mean_Time - gsum.mean_Time(2)

% Plot
figure;
swarmchart(T2_Timing.Method, T2_Timing.Time)
hold on
errorbar(Timing_means(:,1),Timing_means(:,2),'.','MarkerSize',30)
title('Timing Analysis')
ylabel('Total Execution Time (s)');
set(gca,'TickLabelInterpreter','none')


%% Correlation analysis

% Extract and Stack A B C channels
T2_Fields = T2.Properties.VariableNames;
T2_Corr   = {};
T2_Corr{1} = T2(:,['ID','rotID','sigA','frqA',T2_Fields(startsWith(T2_Fields,'Z_A_'))]);
T2_Corr{2} = T2(:,['ID','rotID','sigB','frqB',T2_Fields(startsWith(T2_Fields,'Z_B_'))]);
T2_Corr{3} = T2(:,['ID','rotID','sigC','frqC',T2_Fields(startsWith(T2_Fields,'Z_C_'))]);

T2_Corr{2}.Properties.VariableNames = T2_Corr{1}.Properties.VariableNames;
T2_Corr{3}.Properties.VariableNames = T2_Corr{1}.Properties.VariableNames;
T2_Corr = cat(1, T2_Corr{:});

% Stack by method
T2_Corr = stack(T2_Corr,T2_Fields(startsWith(T2_Fields,'Z_A_')), ...
  'IndexVariableName', 'Method', 'NewDataVariableName', 'CorrZ');

% Anova
T2_Corr.SigID = T2_Corr.ID + "_" + T2_Corr.sigA;
[CorrZ_P,~,CorrZ_Stats] = anovan(T2_Corr.CorrZ, ...
  {categorical(T2_Corr.ID), T2_Corr.frqA, T2_Corr.sigA, T2_Corr.Method}, ...
  'varnames', {'TrialID', 'freq', 'sig', 'Method'}, ...
  'model', [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 0 1 1 0; 0 1 0 1; 0 0 1 1; 0 1 1 1]);
figure
[CorrZ_comps,CorrZ_means,~,CorrZ_mean_names] = multcompare(CorrZ_Stats,'Dimension',[2 3 4]);

figure
multcompare(CorrZ_Stats,'Dimension',[2 4]);

figure
multcompare(CorrZ_Stats,'Dimension',4);

figure
multcompare(CorrZ_Stats,'Dimension',2,'Alpha',0.0001);

gsum_method = groupsummary(T2_Corr,'Method',{'mean','std','median','min', 'max', @(x) quantile(x,[0.25 0.75]), @(x) quantile(x,[0.025 0.975])},'CorrZ');
m_method    = tanh(gsum_method(:,{'median_CorrZ','min_CorrZ','max_CorrZ','fun1_CorrZ','fun2_CorrZ'})) ;
[gsum_method(:,1:2) m_method]

gsum_frqmethod = groupsummary(T2_Corr,{'frqA','Method'},{'mean','std','median','min', 'max', @(x) quantile(x,[0.25 0.75]), @(x) quantile(x,[0.025 0.975])},'CorrZ');
m_frqmethod    = tanh(gsum_frqmethod(:,{'median_CorrZ','min_CorrZ','max_CorrZ','fun1_CorrZ','fun2_CorrZ'}));
[gsum_frqmethod(:,1:2) m_frqmethod]

gsum_sigmethod = groupsummary(T2_Corr,{'sigA','Method'},{'mean','std','median','min', 'max', @(x) quantile(x,[0.25 0.75]), @(x) quantile(x,[0.025 0.975])},'CorrZ');
m_sigmethod    = tanh(gsum_sigmethod(:,{'median_CorrZ','min_CorrZ','max_CorrZ','fun1_CorrZ','fun2_CorrZ'}));
[gsum_sigmethod(:,1:2) m_sigmethod]

T2_Corr.MethodG = T2_Corr.Method;
T2_Corr.MethodG(ismember(T2_Corr.MethodG,categorical(["Z_A_linear","Z_A_fill_spline"]))) = "Z_A_Alternate";
gsum_frqmethod2 = groupsummary(T2_Corr,{'frqA','MethodG'},{'mean','std','median','min', 'max', @(x) quantile(x,[0.25 0.75]), @(x) quantile(x,[0.025 0.975])},'CorrZ');
m_frqmethod2    = tanh(gsum_frqmethod2(:,{'median_CorrZ','min_CorrZ','max_CorrZ','fun1_CorrZ','fun2_CorrZ'})); 
[gsum_frqmethod2(:,1:2) m_frqmethod2]

% Plot
T2_Corr.frqMethod = string(T2_Corr.Method) + "_" + T2_Corr.sigA + "_" + T2_Corr.frqA;
figure;
tiledlayout(2,1)
nexttile
swarmchart(categorical(T2_Corr.frqMethod), tanh(T2_Corr.CorrZ))
hold on
title('Correlation Analysis')
ylabel('Correlation (\rho)');
set(gca,'TickLabelInterpreter','none')

nexttile
swarmchart(categorical(T2_Corr.frqMethod), T2_Corr.CorrZ)
hold on
title('Correlation Analysis')
ylabel('Correlation (Z)');
set(gca,'TickLabelInterpreter','none')


%% Mean Squared error analysis

% Extract and Stack A B C channels
T2_Fields = T2.Properties.VariableNames;
T2_MSE   = {};
T2_MSE{1} = T2(:,['ID','rotID','sigA','frqA',T2_Fields(startsWith(T2_Fields,'MSE_A_'))]);
T2_MSE{2} = T2(:,['ID','rotID','sigB','frqB',T2_Fields(startsWith(T2_Fields,'MSE_B_'))]);
T2_MSE{3} = T2(:,['ID','rotID','sigC','frqC',T2_Fields(startsWith(T2_Fields,'MSE_C_'))]);

T2_MSE{2}.Properties.VariableNames = T2_MSE{1}.Properties.VariableNames;
T2_MSE{3}.Properties.VariableNames = T2_MSE{1}.Properties.VariableNames;
T2_MSE = cat(1, T2_MSE{:});

% Stack by method
T2_MSE = stack(T2_MSE,T2_Fields(startsWith(T2_Fields,'MSE_A_')), ...
  'IndexVariableName', 'Method', 'NewDataVariableName', 'MSE');

% Anova
T2_MSE.SigID = T2_MSE.ID + "_" + T2_MSE.sigA;
[CorrZ_P,~,CorrZ_Stats] = anovan(T2_MSE.MSE, ...
  {categorical(T2_MSE.ID), T2_MSE.frqA, T2_MSE.sigA, T2_MSE.Method}, ...
  'varnames', {'TrialID', 'freq', 'sig', 'Method'}, ...
  'model', [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 0 1 1 0; 0 1 0 1; 0 0 1 1; 0 1 1 1]);
figure
[CorrZ_comps,CorrZ_means,~,CorrZ_mean_names] = multcompare(CorrZ_Stats,'Dimension',[2 3 4]);

% Plot
T2_MSE.frqMethod = string(T2_MSE.Method) + "_" + T2_MSE.sigA + "_" + T2_MSE.frqA;
figure;
swarmchart(categorical(T2_MSE.frqMethod), T2_MSE.MSE)
title('MSE Analysis')
ylabel('MSE (mV^2)');
set(gca,'TickLabelInterpreter','none','YScale','log')



%% Least Squared error time course

load('20240707_method_trials_SErr.mat');

% Extract and Stack A B C channels
LS_Fields = SErr.Properties.VariableNames;
LS_tc   = {};
LS_tc{1} = SErr(:,['ID','rotID','sigA','frqA',LS_Fields(startsWith(LS_Fields,'Err2_A_'))]);
LS_tc{2} = SErr(:,['ID','rotID','sigB','frqB',LS_Fields(startsWith(LS_Fields,'Err2_B_'))]);
LS_tc{3} = SErr(:,['ID','rotID','sigC','frqC',LS_Fields(startsWith(LS_Fields,'Err2_C_'))]);

LS_tc{2}.Properties.VariableNames = LS_tc{1}.Properties.VariableNames;
LS_tc{3}.Properties.VariableNames = LS_tc{1}.Properties.VariableNames;
LS_tc = cat(1, LS_tc{:});

Err2_fields = LS_Fields(startsWith(LS_Fields,'Err2_A_'));
for ii = 1:numel(Err2_fields)
  if iscell(LS_tc.(Err2_fields{ii}))
    dat = LS_tc.(Err2_fields{ii});
    n = min(cellfun(@numel, dat));
    dat = cellfun(@(x) x(1:n), dat, 'UniformOutput', false);
    LS_tc.(Err2_fields{ii}) = cat(1, dat{:});
  end
end
LS_tc_sum = groupsummary(LS_tc, {'sigA','frqA'}, 'mean', LS_Fields(startsWith(LS_Fields,'Err2_A_')));

% Plot
uSig = unique(LS_tc_sum.sigA);
uFrq = unique(LS_tc_sum.frqA);
Meth = LS_Fields(startsWith(LS_Fields,'Err2_A_'));
figure
T = tiledlayout(numel(uSig),numel(uFrq));
title(T,'Demod Error Time Course')
ax = [];
for ss = 1:numel(uSig)
  for ff = 1:numel(uFrq)
    ax(ss,ff) = nexttile;
    tc_ = LS_tc_sum(LS_tc_sum.sigA == uSig(ss) & LS_tc_sum.frqA == uFrq(ff),:);
    hold on;
    for mm = 1:numel(Meth)
      plot(SErr_Time,tc_.("mean_"+Meth(mm)));
    end
    set(ax(ss,ff),'YScale','log')
    title(sprintf('%s/%d',uSig(ss),uFrq(ff)))
  end
end
legend(ax(1), Meth,'Location','NorthWestOutside','interpreter','none')
xlabel('time (s)')
ylabel('Mean Squred Error (mV^2/variance)')
linkaxes(ax)

%% Mufix Only
figure
T = tiledlayout(numel(uSig),1);
title(T,'MuFIX Error Time Course')
ax = [];
for ss = 1:numel(uSig)
  ax(ss) = nexttile;
  tc_ = LS_tc_sum(LS_tc_sum.sigA == uSig(ss),:);
  plot(SErr_Time,tc_.mean_Err2_A_mufix);
  set(ax(ss),'YScale','log')
  title(uSig(ss))
end
legend(ax(1), string(uFrq), 'Location','NorthWestOutside')
xlabel('time (s)')
ylabel('Mean Squared Error (mV^2/variance)')
linkaxes(ax)