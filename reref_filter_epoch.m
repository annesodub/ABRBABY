function [EEG] = reref_filter_epoch(indir, hp,lp, mastos, trig,baseline,win_of_interest)
% ERPs sanity check script - 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

% Reads all folders that are in indir 
d = dir(indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% %% Variables to enter manually before running the code
% 
% % Set variables for filtering
% % hp = 0.1; %value for high-pass filter (Hz) (APICE)
% % lp = 40; %value for low-pass filter (Hz) (APICE) 
% hp = filter.hp; %value for high-pass filter (Hz) (APICE)
% lp = 30; %value for low-pass filter (Hz) (APICE) 
% 
% elec_to_disp_labels = {'F3','Fz','F4';'C3','Cz','C4'};
% 
% % Rejection treshold for bad epochs
% rej_low = -150; %150 infants; 120 adults
% rej_high = 150; %150 infants; 120 adults
% 
% % List of channel labels to reref with 
% mastos = {'Lmon','Rmon','MASTOG','MASTOD'};
% trig = {'Erg1'};
% 
% baseline = [-99, 0] ; 
% win_of_interest = [-0.1, 0.5] ; 
% conditions = {'STD','DEV1','DEV2'} ; 
% cond_sylab = {'BA','GA'} ; 
% 
% elec = 1:16 ; 

%Colors for plots
STD_color = [0.4941 0.1019 0.8863]; %purple
DEV1_color = [1 0.7686 0]; %light orange
DEV2_color = [1 0.4 0]; %dark orange
DEV_colors = {DEV1_color, DEV2_color};
DIFF_color = [0 0 0]; %black

% FOR SANITY CHECK
%for jj=find(ismember(subjects,'DVL_010_T24'))

    %Loop through subjects
    for jj=1:length(subjects) 
        
     fprintf(strcat(subjects{jj}, '...\n'));
     
    %% IMPORT
    % Get BDF file
    %[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
    fname= dir(fullfile(indir,subjects{jj},'*.bdf'));
 
    % Select bdf file in the folder
    EEG = pop_biosig(fullfile(indir, subjects{jj}, fname.name));

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
    
    % Identifies outliers events (e.g. boundaries) or too close events 
    %idx_to_remove = [   find(diff([EEG.event.latency])<0.1*EEG.srate)+1,... % minimum intretrial duration = 220 ms
     %                   find(diff([EEG.event.latency])>2*EEG.srate) ];        
    idx_to_remove = [   find(diff([EEG.event.latency])<0.1*EEG.srate),... % minimum intretrial duration = 220 ms
                        find(diff([EEG.event.latency])>2*EEG.srate) ];  
    % Removes outliers events
    EEG.event(idx_to_remove) = [] ;  EEG.urevent(idx_to_remove) = [] ; 
    
    % Relabels events with condition name (defined in txt file <SUBJECT>.txt)
    EEG.event = read_custom_events(strrep(fullfile(fname.folder,fname.name),'.bdf','.txt'),EEG.event) ;
    EEG.orig_events = EEG.urevent ; EEG.urevent = EEG.event;
    
    %% FILTERS the data with ERPLab
    EEG  = pop_basicfilter(EEG,  1:16 , 'Boundary', 'boundary', 'Cutoff', [hp lp], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2, 'RemoveDC', 'on' );
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );

    %% SAVE DATASET BEFORE EPOCHING
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', strcat(filename,'_filtered'),'savenew', fullfile(indir,subjects{jj}, strcat(filename,'_filtered')),'gui','off');
    CURR_FILTERED = CURRENTSET ; 
    
    % Extract ALL conditions epochs
    EEG = pop_epoch(EEG, conditions, win_of_interest, 'newname', strcat(filename,'_ALL'), 'epochinfo', 'yes');

    % Remove baseline
    EEG = pop_rmbase( EEG, baseline,[] );
    
    % Add channels information
    %EEG=pop_chanedit(EEG, 'lookup','/Users/anne-sophiedubarry/Documents/4_Software/eeglab2020_0/plugins/dipfit/standard_BEM/elec/standard_1005.elc');
    EEG=pop_chanedit(EEG, 'lookup',chan_dir);
    
    end
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