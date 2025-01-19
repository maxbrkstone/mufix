%% Save a copy of the recording with necessary info

RECORDROOT = 'S:\RUData\CA1_ChR2_RU';
RECORDINGS = [...
              "OP2718-211103-111406", "OP2718-211104-095749", "OP2718-211115-134907", "OP2718-211117-130838"...
              "OP2718-211118-124325", "OP2718-211119-112816", "OP2718-211122-141834", "OP2718-211123-115331"...
              "OP275-211103-115508", "OP275-211104-110601", "OP275-211115-125031", "OP275-211117-140342" ...
              "OP275-211118-132745", "OP275-211119-122324", "OP275-211122-131852", "OP275-211123-125257"...
              "OP275-211124-125329", "OP275-211129-154529", "OP275-211201-113714", "OP275-211203-104522" ... 
              "OP275-211130-110840", ...
              "OP193-210721-130721", "OP193-210722-123926", "OP193-210723-112957", "OP193-210726-110721"...
              "OP193-210727-115314", "OP193-210728-113710", "OP193-210729-102931", "OP193-210802-143423", ...
              "OP193-210803-132527", "OP193-210804-120004", "OP193-210805-112840", ... 
              "OP191-210721-112640", "OP191-210722-105613", "OP191-210723-094254", ...
              "OP191-210726-091257", "OP191-210727-100409", "OP191-210728-095734", ...
              "OP191-210802-125502", "OP191-210803-114236", "OP191-210804-101237", ...
              "OP193A-OP166B-210811-120318\OP193A-210811-120318", "OP193A-OP166B-210812-120126\OP193A-210812-120126", ...
              "OP275-211027-143224", "OP2718-211027-132637"
             ];  

REC2 = ["OP275-211027-143224", "OP2718-211027-132637"];


OUTPUT_FOLDER = 'Recordings';


%% Loop through and copy the recordings to local folder

for recording_index = 44:45 %1:length(RECORDINGS)
    
  % Load data
  recname = split(RECORDINGS(recording_index),'\');
  filename = fullfile(RECORDROOT, recname(1), strcat(recname(end), '.mat'));
  fprintf('Loading "%s" ... ', filename)
  tdtextr = load(filename);
  
  % Remove previously demodulated data
  sfields = fieldnames(tdtextr.streams);
  for ii = 1:numel(sfields)
    if isfield(tdtextr.streams.(sfields{ii}),'orgdata')
     tdtextr.streams.(sfields{ii}).data = tdtextr.streams.(sfields{ii}).orgdata;
     tdtextr.streams.(sfields{ii}) = rmfield(tdtextr.streams.(sfields{ii}),{'orgdata','dmodcfg'});
    end
  end

  % Keep only channel 1 of EEG
  filteeg = sfields(startsWith(sfields,'EEG'));
  for ii = 1:numel(filteeg)
    tdtextr.streams.(filteeg{ii}).data = tdtextr.streams.(filteeg{ii}).data(1,:);
  end

  % Remove Raw
  raweeg  = sfields(startsWith(sfields,'Raw'));
  tdtextr.streams = rmfield(tdtextr.streams,raweeg);

  % Remove other non-necessary data
  tdtextr = rmfield(tdtextr,'snips');
  efields = fieldnames(tdtextr.epocs);
  cam     = efields(startsWith(efields,'Cam'));
  tdtextr.epocs = rmfield(tdtextr.epocs,cam);

  % Special treatment to align recordings in REC2
  % - Assigning photometry store names to match with other OP275 and OP2718
  %  recordings
  if ismember(RECORDINGS(recording_index),REC2)
    tdtextr.streams.IRa2 = tdtextr.streams.IRa1;
    tdtextr.streams.RCa2 = tdtextr.streams.RCa1;
    tdtextr.streams.IGa2 = tdtextr.streams.IGa1;
    tdtextr.streams.GCa2 = tdtextr.streams.GCa1;
    tdtextr.streams.Fi2r = tdtextr.streams.Fi1r;
    tdtextr.scalars.Fi2i = tdtextr.scalars.Fi1i;
    tdtextr.streams = rmfield(tdtextr.streams, {'IRa1','RCa1','IGa1','GCa1','Fi1r'});
    tdtextr.scalars = rmfield(tdtextr.scalars, 'Fi1i');
  end

  % Save output
  fileout = fullfile(OUTPUT_FOLDER, strcat(recname(end), '.mat'));
  fprintf(' saving ...\n')
  save(fileout,'-struct','tdtextr');

end
