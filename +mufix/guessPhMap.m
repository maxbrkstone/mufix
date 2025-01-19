function phmap = guessPhMap(tdtmat)
% Guesstimate the best photometry channel source mapping for re-demodulation
%   phmap = guessPhMap(tdtmat)
%
% tdtmat - structure returned by TDT Matlab API TDTBin2Mat.m
% Specifically required the followinf fields:
%     info    - recording meta data: mouse, date, Experiment (protocol)
%     streams - substruct with photometry signals: GCa1, GCa2, RCa2, etc.
%               and photometry sources: Fi1r and Fi2r.
%     scalars - substruct with photometry config: Fi1i and Fi2i.
%
% phmap  - Photometry signal-source mapping table for rmCrosstalk.m
%           e.g.:
%           StoreName    SourceStore    Channel    SeqID
%           ---------    -----------    -------    -----
%           "RCa2"       "Fi2r"         2          1
%           ...

VERSION = 'v1.5, 2024.6.25';

phmap = table();

if isempty(tdtmat) || ~all(isfield(tdtmat,{'streams','scalars'})) || isempty(tdtmat.streams) || isempty(tdtmat.scalars)
  return;
end

% Extract photometry info
FF = fieldnames(tdtmat.streams);
phstores = FF(startsWith(FF, ["Isb","IG","IR","IC","GC","RC"]));
phsrc    = FF(startsWith(FF, ["Fi1r","Fi2r"]));
phinfo   = cellfun(@(x) [x(1:3),'i'], phsrc, 'UniformOutput', false); 
phsrc    = phsrc(ismember(phinfo,fieldnames(tdtmat.scalars)));
hasFi1   = any(strcmp(phsrc,"Fi1r"));
hasFi2   = any(strcmp(phsrc,"Fi2r"));

% Number of light sources in Fi1i/Fi2i
nlightsrc = cellfun(@(x) size(tdtmat.scalars.([x(1:3),'i']).data,1)/4, phsrc); 
nchan     = cellfun(@(x) size(tdtmat.streams.(x).data,1), phsrc); 

if isempty(phstores) || isempty(phsrc), return; end

phmap.StoreName      = string(phstores);
phmap.SourceStore(:) = string(phsrc(1));
phmap.Channel(:)     = 1;
phmap.SeqID(:)       = 1;
phmap.Properties.RowNames = phmap.StoreName;

if numel(phsrc) > 1
  bstores = endsWith(phstores,["B","2"]);
  phmap.SourceStore(bstores) = {'Fi2r'};
end


%% Guess Channel and SeqID values based on general heuristic 
rcamp = startsWith(phstores,["RC","IR"]);
phmap.Channel(rcamp) = arrayfun(@(x) size(tdtmat.streams.(x).data,1), phmap.SourceStore(rcamp));
    % account for the actual number of channels in the source signal

gcamp2 = startsWith(phstores,"GC");
phmap.SeqID(gcamp2)  = 2;


%% Specific mapping based on mouse/protocol/date

% mouse    = string(tdtmat.info.mouse);
% protocol = string(tdtmat.info.Experiment); 
date     = datetime(tdtmat.info.date);

%flip actually happened on 9/14/21, recording sheet updated on 9/30/21:
%before, s1 iso s2 gca s3 rca 
%after, s1 rca s2 gca s3 iso
if(date >= datetime(2021, 09, 14))
    rcamp = startsWith(phstores,"RC");
    phmap.SeqID(rcamp) = 1;
    
    isb = startsWith(phstores,"I");
    phmap.SeqID(isb)    = 3;
else
    rcamp = startsWith(phstores,"RC");
    phmap.SeqID(rcamp) = 3;

    isb = startsWith(phstores,"I");
    phmap.SeqID(isb)    = 1;
end

fi1r = strcmp(phmap.SourceStore,'Fi1r');
fi2r = strcmp(phmap.SourceStore,'Fi2r');
frqs1 = sum(strcmp(phmap.SourceStore,'Fi1r'), 'all');
frqs2 = sum(strcmp(phmap.SourceStore,'Fi2r'), 'all');
if(frqs1 == 2)
    phmap.SeqID(fi1r)    = arrayfun(@(x) max(1, x-1), phmap.SeqID(fi1r));
end
if(frqs2 == 2)
    phmap.SeqID(fi2r)    = arrayfun(@(x) max(1, x-1), phmap.SeqID(fi2r));
end
%has _duet suffix if two mice; add this check in protocol for a double ipsi recording setup

%% Read Fi1i and Fi2i for information

phmap.Freq__Hz(:)      = 0;
phmap.PowerOn__sec(:)  = {{}};
phmap.PowerOff__sec(:) = {{}};
phmap.ACPower__mA(:)   = {{}};
phmap.DCPower__mA(:)   = {{}};

sname = strrep(phmap.SourceStore,'r','i');
for ii = 1:height(phmap)
  fsrc = tdtmat.scalars.(sname(ii));
  fdat = fsrc.data((1:4) + 4*(phmap.SeqID(ii)-1),:);
  onmask = fdat(1,:) ~= 0;

  phmap.Freq__Hz(ii)      = fdat(2,1);
  phmap.PowerOff__sec{ii} = fsrc.ts(~onmask);
  phmap.PowerOn__sec{ii}  = fsrc.ts(onmask);
  phmap.ACPower__mA{ii}   = fdat(3,onmask);
  phmap.DCPower__mA{ii}   = fdat(4,onmask);
end


%% Final output
phmap = sortrows(phmap, {'SourceStore','Channel','SeqID'});

