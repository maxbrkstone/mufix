% Modular code should allow to adjust the following variables:
% Lenght of after-pulse artifact, amount of before and after data to model
% from, frequencies chosen, number of frequencies, 

SZDB    = 'seizure_list.xls';
SZ_LIST = readtable(SZDB,"TextType","string");

RECORDROOT = 'Recordings';
RECORDING  = "OP275-211027-143224";

              % "OP2718-211103-111406", ...
              % "OP2718-211104-095749", "OP2718-211115-134907", "OP2718-211117-130838"...
              % "OP2718-211118-124325", "OP2718-211119-112816", "OP2718-211122-141834", "OP2718-211123-115331"...
              % "OP275-211103-115508", "OP275-211104-110601", "OP275-211115-125031", "OP275-211117-140342" ...
              % "OP275-211118-132745", "OP275-211119-122324", "OP275-211122-131852", "OP275-211123-125257"...
              % "OP275-211124-125329", "OP275-211129-154529", "OP275-211201-113714", "OP275-211203-104522" ... 
              % "OP275-211130-110840", ...
              % "OP193-210721-130721", "OP193-210722-123926", "OP193-210723-112957", "OP193-210726-110721"...
              % "OP193-210727-115314", "OP193-210728-113710", "OP193-210729-102931", "OP193-210802-143423", ...
              % "OP193-210803-132527", "OP193-210804-120004", "OP193-210805-112840", ... 
              % "OP193A-210811-120318", "OP193A-210812-120126", ...              
              % "OP191-210721-112640", "OP191-210722-105613", "OP191-210723-094254", ...
              % "OP191-210726-091257", "OP191-210727-100409", "OP191-210728-095734", ...
              % "OP191-210802-125502", "OP191-210803-114236", "OP191-210804-101237", ...

% Correlation
CORR.extraEp_s = [-2 4];

%% Perform MuFIX

filename = fullfile(RECORDROOT, strcat(RECORDING, '.mat'));
tdtextr  = load(filename);

map = mufix.guessPhMap(tdtextr);
map = sortrows(map,'SeqID');
phscalar = char(map.SourceStore(1));
phscalar = [phscalar(1:3),'i'];

cfg = struct();
cfg.plot = 0;
cfg.lowpass_zerophase = 1;
    
% Demodulate without MuFIX
cfg.xtalk_method = "none";
tdtextr = mufix.rmCrosstalk_tdtmat(tdtextr, map, cfg);

% Demodulate with MuFIX
cfg.xtalk_method = "mufix";
[tdtmufix,~,dmodat] = mufix.rmCrosstalk_tdtmat(tdtextr, map, cfg);


%% Characterize Crosstalk

fs = tdtextr.streams.(dmodat{1}.SourceStore).fs;
pstart = round((tdtextr.epocs.Pls_.onset - 0.01) * fs);
pdur   = round(0.05 * fs);
pidx   = (0:pdur);
[xx,yy] = meshgrid(pstart,pidx);
xy = xx + yy;
phdat = tdtextr.streams.(dmodat{1}.SourceStore).data(2,xy);
phdat = reshape(phdat,numel(pidx),[]);
mphdat = mean(phdat,2);

figure;
plot(1000*pidx/fs,phdat);
xlabel('ms');
ylabel('V')

% save crosstalk template


%% Calculate R2 between original and interpdat
R2 = {};
chname = {};
for ii = 1:numel(dmodat)
  chname{ii} = dmodat{ii}.SourceStore + "_Ch" + dmodat{ii}.Channel;
  orgdat     = tdtextr.streams.(dmodat{ii}.SourceStore).data(dmodat{ii}.Channel,:);
  interpdat  = dmodat{ii}.interpdat;
  tidx       = dmodat{ii}.interpcfg.tseg * dmodat{ii}.interpcfg.fs;
  tidx(:,1)  = floor(tidx(:,1));
  tidx(:,2)  = ceil(tidx(:,2));
  R2{ii}     = nan(size(tidx,1),1);
  for jj = 1:size(tidx,1)
    inds = tidx(jj,1):tidx(jj,2);
    R2{ii}(jj) = 1 - mean(interpdat(inds) - orgdat(inds)).^2 / var(orgdat(inds));
  end
end

chname
nR2    = cellfun(@numel, R2)
midR2  = cellfun(@median, R2)
minR2  = cellfun(@min, R2)
maxR2  = cellfun(@max, R2)


%% Calculate correlation between no-mufix and mufix demodulated output
TIME_CORR = [-2 4];

try
  fs    = tdtextr.streams.GCa1.fs;
catch
  fs    = tdtextr.streams.GCa2.fs;
end
epstart = floor((tdtextr.epocs.Brst.onset  + TIME_CORR(1)) * fs);
epend   = floor((tdtextr.epocs.Brst.offset + TIME_CORR(2)) * fs);
for ii = 1:height(map)
  orgdat   = tdtextr .streams.(map.StoreName(ii)).data;
  mufixdat = tdtmufix.streams.(map.StoreName(ii)).data;
  cr = nan(1,numel(epstart));
  for jj = 1:numel(epstart)
    inds   = epstart(jj):epend(jj);
    cr(jj) = corr(orgdat(inds)',mufixdat(inds)');    
  end
  map.Corr{ii} = cr;
end

midCR  = cellfun(@(x) tanh(median(atanh(x))), map.Corr)
minCR  = cellfun(@min, map.Corr)
maxCR  = cellfun(@max, map.Corr)
nCR  = cellfun(@numel, map.Corr)

