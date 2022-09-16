%% ABRs analysis script - Estelle Herv� - 2022 - %80PRIME Project

%% Variables to enter manually before running the code
clear all;

addpath(genpath('/Users/anne-sophiedubarry/Documents/4_Software/eeglab2020_0'));

%Set filepath (must contain .bdf and .txt files from recording)
INDIR = '/Users/anne-sophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data';

% Reads all folders that are in INDIR 
d = dir(INDIR); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

%Set variables for filtering
hp = 80; %value for high-pass filter (Hz)
lp = 3000; %value for low-pass filter (Hz)

%Rejection treshold for bad epochs
rej_low = -45;
rej_high = 45;

%Epoch window
epoch_timew =  [-0.04, 0.2] ; 
baseline_timew = [-39 0]  ; 

% Loop through subjects
% for jj=1:1%length(subjects) 
  
jj=find(ismember(subjects,'pilote_salome')) ; 

fname= dir(fullfile(INDIR,subjects{jj},'*.bdf'));

bdf_filename = fname.name ; 
txt_filename = strrep(bdf_filename,'.bdf','.txt');
filename = strrep(bdf_filename,'.bdf','');


%% From loading file to get ABR + trigg channels

% Open eeglab
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

% Select bdf file in the folder
EEG = pop_biosig(fullfile(INDIR,subjects{jj},bdf_filename));

% Rename file
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 0,'setname',filename,'gui','off');
EEG = eeg_checkset( EEG );

% Select channels to keep
EEG = eeg_checkset( EEG );
EEG = pop_select( EEG, 'channel',{'Erg1','Left','Right'});
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'gui','off'); 

eeglab redraw;
%%
% Compute ABR channel from Left / Right channels
% forumla : {(Left)+(Right)}/-2 = Ref - {(LA+RA)/2}
EEG.nbchan = EEG.nbchan+1;
if ~isempty(EEG.chanlocs)
	EEG.chanlocs(end+1).labels = 'ABR';
    EEG.chanlocs(end).type = 'EEG';
end

for m = 1:size(EEG.data,2)
    EEG.data(4,m) = ((EEG.data(1,m))+(EEG.data(2,m)))/2; 
end
%%
% Keep only ABR channel and triggers (select channel)
EEG = eeg_checkset( EEG );
EEG = pop_select( EEG, 'channel',{'ABR', 'Erg1'});
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'gui','off'); 

% Extract event from trigger channel (Erg1)
EEG = pop_chanevent(EEG, 1,'oper','X>20000','edge','leading','edgelen',1);
[ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
EEG = eeg_checkset( EEG );
pop_eegplot( EEG, 1, 1, 1);

eeglab redraw;

%%
% Reject bad events --> pour l'instant pas ok pour tous les cas de figure, ex joseph_T0
%pas besoin pour DVL_003_T3
% bad_events = [];
% ind = 1;
% for ii = 1:size(EEG.event,2)-1
%     if EEG.event(ii).latency-EEG.event(ii+1).latency > -1000 
%         bad_events(ind) = EEG.event(ii).urevent;
%         ind = ind+1;
%     elseif EEG.event(ii).latency-EEG.event(ii+1).latency < -150000 
%         bad_events(ind) = EEG.event(ii).urevent;
%         ind = ind+1;
%         end
%     ii=ii+1;
% end
% 
% new_latency = [];
% new_type = [];
% new_urevent = [];
% indx = 1;
% 
% for pp = 1:size(EEG.event,2)
%    if ismember(EEG.event(pp).urevent, bad_events(:))== 0
%        new_latency(indx) = EEG.event(pp).latency;
%        new_type(indx) = EEG.event(pp).type;
%        new_urevent(indx) = EEG.event(pp).urevent;
%        indx = indx+1;
%    end
% end
% 
% new_latency = num2cell(new_latency);
% new_type = num2cell(new_type);
% new_urevent = num2cell(new_urevent);
% 
% EEG.event = struct('latency', (new_latency(:))', ...
%                 'type', (new_type(:))',...
%                 'urevent', (new_urevent(:))');
            

%% Replace event type by infos from .txt
% Read .txt
my_events = readtable(fullfile(INDIR,subjects{jj}, txt_filename), 'ReadVariableNames', 0);

%Insert info from .txt into EEG.event
my_events = table2array(my_events);

temp = struct('latency', {EEG.event(:).latency}, ...
                'type', (my_events(:))',...
                'urevent', {EEG.event(:).urevent});

EEG.event = temp;          
eeglab redraw;

%% Filtering
% Filter data
% EEG = pop_eegfiltnew(EEG, 'locutoff',80,'hicutoff',3000,'channels',{'ABR'});%values chosen from Skoe and Kraus, 2010 and stimuli features (Praat) --> cf. Word document 'filtering_choices_ABRs.docx'
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'gui','off'); 
% eeglab redraw;

% Plot data
% EEG = eeg_checkset( EEG );
% pop_eegplot( EEG, 1, 1, 1);

%% Extract epochs for HF

% EEG = pop_epoch( EEG, {  'HF'  }, [-0.01         0.1], 'newname', 'clement_pilote_rej resampled epochs', 'epochinfo', 'yes');
% [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 6,'setname','clement_pilote_rej resampled epochs HF ABR','gui','off'); 
% EEG = eeg_checkset( EEG );
% eeglab redraw;

EEG = pop_epoch( EEG, {  'HF'  }, epoch_timew, 'newname', 'epochs', 'epochinfo', 'yes');
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'gui','off') 
EEG = eeg_checkset( EEG );

%Remove baseline
EEG = pop_rmbase( EEG, baseline_timew,[]);
[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 6,'gui','off'); 

%Epoch rejection
EEG = eeg_checkset( EEG );
%EEG = pop_eegthresh(EEG,1,1,rej_low,rej_high,-0.039978,0.059937,0,1);

binf = EEG.times(1)/1000 ; 
bsup = EEG.times(end)/1000;

EEG = pop_eegthresh(EEG,1,1,rej_low,rej_high,binf,bsup,0,1);

[ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 7,'setname','epochs_rej','gui','off'); 

eeglab redraw;

%% Filtering

% Filter data
EEG  = pop_basicfilter( EEG,  1 , 'Cutoff', [ 80 3000], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2 ); % GUI: 11-Apr-2022 12:47:48
% EEG  = pop_basicfilter( EEG,  elec , 'Boundary', 'boundary', 'Cutoff', [hp lp], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2, 'RemoveDC', 'on' );
[ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 6, 'savenew', filename_filter,'gui','off');
 
%  
% %Filtering
% [ALLEEG EEG CURRENTSET] = eeg_store(ALLEEG, EEG, CURRENTSET);
% EEG = pop_eegfiltnew(EEG, 'locutoff',hp,'hicutoff',lp);
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 7,'gui','off'); 

%Extract mean activity (erp) and replace data
abr = mean(EEG.data(1,:,:),3);
EEG.data = abr;

eeglab redraw;

% Add tube delay (27 cm x 340 m/s ) 
nsample_delay = fix(EEG.srate * (0.27 / 340) ) ; 

abr_shifted = circshift(abr,nsample_delay) ;



%% Export ABR data into .txt file
fname_out = fullfile(filepath,strrep(bdf_filename,'.bdf','abr_shifted_data.txt')) ;
fid = fopen(fname_out,'w');
fprintf(fid,'%c\n',abr_shifted);
fclose(fid);

bt_txt2avg(fname_out, EEG.srate, epoch_timew(1)*1000, epoch_timew(2)*1000);
