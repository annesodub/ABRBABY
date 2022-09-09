%% ERPs analysis script - Estelle Herv� - 2022 - %80PRIME Project

%% Variables to enter manually before running the code

% Load EEGLAB 
% addpath(genpath('/Users/anne-sophiedubarry/Documents/4_Software/eeglab2020_0'));
tmp = pwd ; 
cd '/Users/anne-sophiedubarry/Documents/4_Software/eeglab2020_0' ; 
% Open eeglab
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
run('/Users/anne-sophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/dev/signal_processing/biosig4octmat-3.8.0/biosig_installer.m') ; 

cd(tmp) ; 

% Set filepath (must contain .bdf and .txt files from recording)
% INDIR = '/Users/anne-sophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data';
INDIR = '/Users/anne-sophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data/DEVLANG_data' ;

% Reads all folders that are in INDIR 
d = dir(INDIR); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Set variables for filtering
% hp = 0.1; %value for high-pass filter (Hz) (APICE)
% lp = 40; %value for low-pass filter (Hz) (APICE) 
hp = 1; %value for high-pass filter (Hz) (APICE)
lp = 30; %value for low-pass filter (Hz) (APICE) 


% Rejection treshold for bad epochs
rej_low = -150; %150 infants; 120 adults
rej_high = 150; %150 infants; 120 adults

% List of channel labels to reref with 
mastos = {'Lmon','Rmon','MASTOG','MASTOD'};
trig = {'Erg1'};

baseline = [-99, 0] ; 
win_of_interest = [-0.1, 0.5] ; 
conditions = {'STD','DEV1','DEV2'} ; 
elec = 1:16 ; 

% FOR SANITY CHECK
for jj=find(ismember(subjects,'DVL_007_T8'))

% Loop through subjects
% for jj=1:length(subjects) 
        
%     fprintf(strcat(subjects{jj}, '...\n'));
%     jj=find(ismember(subjects,'DVL_007_T8')) ; 
    
    %% IMPORT
    % Get BDF file
    fname= dir(fullfile(INDIR,subjects{jj},'*.bdf'));
 
    % Select bdf file in the folder
    EEG = pop_biosig(fullfile(INDIR, subjects{jj}, fname.name));

    % Find REF electrodes indices by labels 
    ref_elec = find(ismember({EEG.chanlocs.labels},mastos)); 
    
    % Save a first dataset in EEGLAB 
    [~,filename,~] = fileparts(fname.name);    
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',filename,'gui','off');
 
    %% RE-REF (excluding trig channel)
    % Find TRIG electrodes indices by labels 
    trigg_elec = find(ismember({EEG.chanlocs.labels},trig)); 

    % Re-reference data and rename new file
    EEG = pop_reref(EEG, ref_elec, 'exclude',trigg_elec, 'keepref','on');
    
    %% EVENTS 
    % Extract event from trigger channel (Erg1)
    EEG = pop_chanevent(EEG, trigg_elec,'oper','X>20000','edge','leading','edgelen',1);
    
    % Removes events outliers (e.g. boundaries) or too close events 
    EEG.event(find(diff([EEG.event.latency])<1000)+1) = [] ;
    EEG.event(find(diff([EEG.event.latency])>100000))= [] ; % threshold to adjust

    % Relabels events with condition name (defined in txt file <SUBJECT>.txt)
    EEG.event = read_custom_events(strrep(fullfile(fname.folder,fname.name),'.bdf','.txt'),EEG.event) ;
    
    %% FILTERS the data with ERPLab
%     EEG  = pop_basicfilter(EEG,  elec , 'Boundary', 'boundary', 'Cutoff', [hp lp], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2, 'RemoveDC', 'on' );
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );

    %% SAVE DATASET BEFORE EPOCHING
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', strcat(filename,'_filtered'),'savenew', fullfile(INDIR,subjects{jj}, strcat(filename,'_filtered')),'gui','off');
    CURR_FILTERED = CURRENTSET ; 
    
    %% EPOCHING 
    for cc=1:length(conditions)
        
        [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'retrieve',CURR_FILTERED,'study',0); 
        EEG = eeg_checkset( EEG );
    
        % Build condition name
        cond_filename = strrep(fname.name,'.bdf',conditions{cc});
        
        % Extract epochs 
        EEG = pop_epoch(EEG, {conditions{cc}}, win_of_interest, 'newname', strcat(filename,'_',conditions{cc}), 'epochinfo', 'yes');
        
        % Remove baseline
        EEG = pop_rmbase( EEG, baseline,[] );
    
        % Reject data epochs 
%         EEG = pop_eegthresh(EEG,1,elec ,rej_low, rej_high,-0.099609,0.49805,0,1);
        EEG = pop_eegthresh(EEG,1,elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    
        % Add channels information
        EEG=pop_chanedit(EEG, 'lookup','/Users/anne-sophiedubarry/Documents/4_Software/eeglab2020_0/plugins/dipfit/standard_BEM/elec/standard_1005.elc');
 
         %  Plot average
        figure; pop_plottopo(EEG, elec , conditions{cc}, 0, 'ydir',-1);

        % Save dataset in EEGLAB 
%         [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET,'savenew',strcat(filename,'_',conditions{cc}),'gui','off');
        [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname',strcat(filename,'_',conditions{cc}),'savenew', fullfile(INDIR,subjects{jj},strcat(filename,'_',conditions{cc})),'gui','off');
  
    end
    
    %% SANITY CHECK DISPLAY F3-FZ-F4 , C3-CZ-C4 ou une topo complete? (16 elec)
%     % ii= 5 corresponds to electrode F3 
% grd_STD_F3 = squeeze(mean(STD_subj(:,5,:),1)) ; 
% grd_DEV1_F3 = squeeze(mean(DEV1_subj(:,5,:),1)) ; 
% grd_DEV2_F3 = squeeze(mean(DEV2_subj(:,5,:),1)) ; 
% 
% figure ; 
% plot(EEG_dev2.times,grd_STD_F3,'k','Linewidth',1.5); hold on ;set(gca,'YDir','reverse') ;
% plot(EEG_dev2.times,grd_DEV1_F3,'r','Linewidth',1.5);  hold on; set(gca,'YDir','reverse') ;
% plot(EEG_dev2.times, grd_DEV2_F3,'b','Linewidth',1.5); set(gca,'YDir','reverse') ;
% grid on ; 
% 
% legend('STD','DEV1','DEV2');
% xlabel('Times (ms)'); ylabel('uV'); title ('Grand average F3');
    
    %%%% TO FINISH (following APICE?)
    
end

%--------------------------------------------------------------
% FUNCTION that reads events from text file and output 
% an EEGLAB events structure 
%--------------------------------------------------------------
function out_event = read_custom_events(fname, in_event) 

% Read .txt 
my_events = readtable(fname, 'ReadVariableNames', 0);

% Insert info from .txt into EEG.event
my_events = table2array(my_events);

out_event = struct('latency', {in_event(:).latency}, ...
                'type', (my_events(:))',...
                'urevent', {in_event(:).urevent});

end