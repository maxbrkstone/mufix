%% Analysis of the Empirical data

load('20240701_20Hz_5ms_empirtrials.mat');


%% Correlation analysis

% Extract and Stack A B C channels
T2_Fields = T2.Properties.VariableNames;
T2_Corr   = {};
T2_Corr{1} = T2(:,['ID','rotID','Mouse','sigA','epoch','hasSeizure','frqA','meanA','sdA','paroxtime_s',T2_Fields(startsWith(T2_Fields,'Z_A_'))]);
T2_Corr{2} = T2(:,['ID','rotID','Mouse','sigB','epoch','hasSeizure','frqB','meanB','sdB','paroxtime_s',T2_Fields(startsWith(T2_Fields,'Z_B_'))]);
T2_Corr{3} = T2(:,['ID','rotID','Mouse','sigC','epoch','hasSeizure','frqC','meanC','sdC','paroxtime_s',T2_Fields(startsWith(T2_Fields,'Z_C_'))]);

T2_Corr{2}.Properties.VariableNames = T2_Corr{1}.Properties.VariableNames;
T2_Corr{3}.Properties.VariableNames = T2_Corr{1}.Properties.VariableNames;
T2_Corr = cat(1, T2_Corr{:});
T2_Corr.sigMouse = T2_Corr.Mouse + "-" + T2_Corr.sigA;

% Stack by method
T2_Corr = stack(T2_Corr,T2_Fields(startsWith(T2_Fields,'Z_A_')), ...
  'IndexVariableName', 'Method', 'NewDataVariableName', 'CorrZ');

% Anova
T2_Corr.SigID = T2_Corr.Mouse + "_" + T2_Corr.sigA;
[CorrZ_P,~,CorrZ_Stats] = anovan(T2_Corr.CorrZ, ...
  {T2_Corr.frqA, T2_Corr.sigMouse, T2_Corr.Method, log(T2_Corr.sdA)}, ...
  'continuous', 4, ...
  'varnames', {'freq', 'sig', 'Method', 'sd'}, ...
  'model', [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 1 1 0 0; 1 0 1 0; 1 0 0 1; 0 1 1 0; 0 0 1 1] ...
  );  

figure
[CorrZ_comps1,CorrZ_means1,~,CorrZ_mean_names1] = multcompare(CorrZ_Stats,'Dimension',[1 2 3]);

figure
[CorrZ_comps2,CorrZ_means2,~,CorrZ_mean_names2] = multcompare(CorrZ_Stats,'Dimension',[2 3]);

gsum_method = groupsummary(T2_Corr,'Method',{'mean','std','median','min', 'max', @(x) quantile(x,[0.25 0.75]), @(x) quantile(x,[0.025 0.975])},'CorrZ');
m_method    = tanh(gsum_method(:,{'median_CorrZ','min_CorrZ','max_CorrZ','fun1_CorrZ','fun2_CorrZ'})) 

gsum_frqmethod = groupsummary(T2_Corr,{'frqA','Method'},{'mean','std','median','min', 'max', @(x) quantile(x,[0.25 0.75]), @(x) quantile(x,[0.025 0.975])},'CorrZ');
m_frqmethod    = tanh(gsum_frqmethod(:,{'median_CorrZ','min_CorrZ','max_CorrZ','fun1_CorrZ','fun2_CorrZ'})) 

gsum_sigmethod = groupsummary(T2_Corr,{'sigA','Method'},{'mean','std','median','min', 'max', @(x) quantile(x,[0.25 0.75]), @(x) quantile(x,[0.025 0.975])},'CorrZ');
m_sigmethod    = tanh(gsum_sigmethod(:,{'median_CorrZ','min_CorrZ','max_CorrZ','fun1_CorrZ','fun2_CorrZ'})) 


% Plot
cmap = lines(2);

T2_Corr.frqMethod = string(T2_Corr.Method) + "_" + T2_Corr.sigA + "_" + T2_Corr.frqA;
figure;
T = tiledlayout(2,2,'tileindexing','columnmajor');
title(T,'Correlation Analysis')

nexttile
swarmchart(categorical(T2_Corr.frqMethod), tanh(T2_Corr.CorrZ), 10, T2_Corr.hasSeizure, 'filled');
hold on
title('Freq x Signal Source x Method')
ylabel('Correlation (\rho)');
set(gca,'TickLabelInterpreter','none')

nexttile
swarmchart(categorical(T2_Corr.frqMethod), T2_Corr.CorrZ, 10, T2_Corr.hasSeizure, 'filled');
hold on
ylabel('Correlation (Z)');
set(gca,'TickLabelInterpreter','none')

T2_Corr.sigMethod = string(T2_Corr.Method) + "_" + T2_Corr.sigMouse;
nexttile
swarmchart(categorical(T2_Corr.sigMethod), tanh(T2_Corr.CorrZ), 10,  T2_Corr.hasSeizure, 'filled');
hold on
h(1) = plot(nan(1),nan(1),'.','Color',cmap(1,:));
h(2) = plot(nan(1),nan(1),'.','Color',cmap(2,:));
legend(h,'Non-Seizure','Seizure','Location','NorthEastOutside')
title('Signal Source x Method')
ylabel('Correlation (\rho)');
set(gca,'TickLabelInterpreter','none')

nexttile
swarmchart(categorical(T2_Corr.sigMethod), T2_Corr.CorrZ, 10,  T2_Corr.hasSeizure, 'filled');
hold on
ylabel('Correlation (Z)');
set(gca,'TickLabelInterpreter','none')

colormap(cmap)



%% Mean Squared error analysis

% Extract and Stack A B C channels
T2_Fields = T2.Properties.VariableNames;
T2_MSE   = {};
T2_MSE{1} = T2(:,['ID','rotID','Mouse','sigA','frqA','hasSeizure','frqA','meanA','sdA','paroxtime_s',T2_Fields(startsWith(T2_Fields,'MSE_A_'))]);
T2_MSE{2} = T2(:,['ID','rotID','Mouse','sigB','frqB','hasSeizure','frqB','meanB','sdB','paroxtime_s',T2_Fields(startsWith(T2_Fields,'MSE_B_'))]);
T2_MSE{3} = T2(:,['ID','rotID','Mouse','sigC','frqC','hasSeizure','frqC','meanC','sdC','paroxtime_s',T2_Fields(startsWith(T2_Fields,'MSE_C_'))]);

T2_MSE{2}.Properties.VariableNames = T2_MSE{1}.Properties.VariableNames;
T2_MSE{3}.Properties.VariableNames = T2_MSE{1}.Properties.VariableNames;
T2_MSE = cat(1, T2_MSE{:});
T2_MSE.sigMouse = T2_MSE.Mouse + "-" + T2_MSE.sigA;

% Stack by method
T2_MSE = stack(T2_MSE,T2_Fields(startsWith(T2_Fields,'MSE_A_')), ...
  'IndexVariableName', 'Method', 'NewDataVariableName', 'MSE');

% Anova
T2_MSE.SigID = T2_MSE.ID + "_" + T2_MSE.sigA;
[CorrZ_P,~,CorrZ_Stats] = anovan(T2_MSE.MSE, ...
  {T2_MSE.frqA, T2_MSE.sigA, T2_MSE.Method, T2_MSE.hasSeizure}, ...
  'varnames', {'freq', 'sig', 'Method', 'seizure'}, ...
  'model', [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1; 0 1 1 0; 0 1 0 1; 0 0 1 1; 0 1 1 1] ...
  );
figure
[CorrZ_comps,CorrZ_means,~,CorrZ_mean_names] = multcompare(CorrZ_Stats,'Dimension',[2 3 4]);


% Plot
cmap = lines(2);
T2_MSE.frqMethod = string(T2_MSE.Method) + "_" + T2_MSE.sigA + "_" + T2_MSE.frqA;
figure;
swarmchart(categorical(T2_MSE.frqMethod), T2_MSE.MSE, 10, T2_MSE.hasSeizure, 'filled')
hold on
h(1) = plot(nan(1),nan(1),'.','Color',cmap(1,:));
h(2) = plot(nan(1),nan(1),'.','Color',cmap(2,:));
legend(h,'Non-Seizure','Seizure','Location','NorthEastOutside')
title('MSE Analysis')
ylabel('MSE (mV^2)');
set(gca,'TickLabelInterpreter','none','YScale','log')
colormap(cmap)


%% Relationship between signal SD
uMethod = struct();
uMethod.Corr = unique(T2_Corr.Method);
uMethod.MSE  = unique(T2_MSE.Method);
cmap    = lines(3);

figure
T = tiledlayout(3,numel(uMethod.Corr),'tileindexing','columnmajor');
title(T,'Relationship between SD and MSE')
for mm = 1:numel(uMethod.Corr)

  T2_Corr_ = T2_Corr(T2_Corr.Method == uMethod.Corr(mm) & T2_Corr.frqA == FREQ.f1,:);

  nexttile
  scatter(T2_Corr_.sdA,tanh(T2_Corr_.CorrZ),10,categorical(T2_Corr_.sigA))
  xlabel('Signal SD (mV)')
  ylabel('Correlation (\rho)')
  title(uMethod.Corr(mm),'interpreter','none')
  set(gca,'XScale','log')

  nexttile
  scatter(T2_Corr_.sdA,T2_Corr_.CorrZ,10,categorical(T2_Corr_.sigA))
  xlabel('Signal SD (mV)')
  ylabel('Correlation (Z)')
  set(gca,'XScale','log')


  T2_MSE_ = T2_MSE(T2_MSE.Method == uMethod.MSE(mm) & T2_MSE.frqA == FREQ.f1,:);
  nexttile
  scatter(T2_MSE_.sdA,T2_MSE_.MSE,10,categorical(T2_MSE_.sigA))
  xlabel('Signal SD (mV)')
  ylabel('MSE (mV^2)')
  set(gca,'XScale','log','YScale','log')
  
end

colormap(cmap)


%% Least Squared error time course

% Extract and Stack A B C channels
LS_Fields = SErr.Properties.VariableNames;
LS_tc   = {};
LS_tc{1} = SErr(:,['ID','rotID','sigA','frqA',LS_Fields(startsWith(LS_Fields,'Err2_A_'))]);
LS_tc{2} = SErr(:,['ID','rotID','sigB','frqB',LS_Fields(startsWith(LS_Fields,'Err2_B_'))]);
LS_tc{3} = SErr(:,['ID','rotID','sigC','frqC',LS_Fields(startsWith(LS_Fields,'Err2_C_'))]);

LS_tc{2}.Properties.VariableNames = LS_tc{1}.Properties.VariableNames;
LS_tc{3}.Properties.VariableNames = LS_tc{1}.Properties.VariableNames;
LS_tc = cat(1, LS_tc{:});

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


