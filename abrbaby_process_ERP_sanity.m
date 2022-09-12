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
INDIR = '/Users/anne-sophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data/DEVLANG_data/' ;

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

% Fontsize for graphics
FONTSZ = 12 ; 

elec_to_disp_labels = {'F3','Fz','F4';'C3','Cz','C4'};

% Rejection treshold for bad epochs
rej_low = -150; %150 infants; 120 adults
rej_high = 150; %150 infants; 120 adults

% List of channel labels to reref with 
mastos = {'Lmon','Rmon','MASTOG','MASTOD'};
trig = {'Erg1'};

baseline = [-99, 0] ; 
win_of_interest = [-0.1, 0.5] ; 
conditions = {'STD','DEV1','DEV2'} ; 
cond_sylab = {'BA','GA'} ; 

elec = 1:16 ; 

% FOR SANITY CHECK
for jj=find(ismember(subjects,'DVL_005_T18'))

% Loop through subjects
% for jj=1:length(subjects) 
        
%     fprintf(strcat(subjects{jj}, '...\n'));
%     jj=find(ismember(subjects,'DVL_007_T8')) ; 
    
   
    %% IMPORT
    % Get BDF file
    fname= dir(fullfile(INDIR,subjects{jj},'*.bdf'));
 
    % Print a message if wrong file name 
    if isempty(fname) ; error('File does not exist, please check path and name'); end ;
            
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
    
    % Identifies outliers events (e.g. boundaries) or too close events 
    idx_to_remove = [   find(diff([EEG.event.latency])<0.1*EEG.srate),... % minimum intretrial duration = 100 ms
                        find(diff([EEG.event.latency])>2*EEG.srate) ];        
    
    % Removes outliers events
    EEG.event(idx_to_remove) = [] ;  EEG.urevent(idx_to_remove) = [] ; 
    
    % Relabels events with condition name (defined in txt file <SUBJECT>.txt)
    EEG.event = read_custom_events(strrep(fullfile(fname.folder,fname.name),'.bdf','.txt'),EEG.event) ;
    EEG.orig_events = EEG.urevent ; EEG.urevent = EEG.event;
    
    %% FILTERS the data with ERPLab
%     EEG  = pop_basicfilter(EEG,  elec , 'Boundary', 'boundary', 'Cutoff', [hp lp], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2, 'RemoveDC', 'on' );
    [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
    EEG = eeg_checkset( EEG );

    %% SAVE DATASET BEFORE EPOCHING
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, CURRENTSET, 'setname', strcat(filename,'_filtered'),'savenew', fullfile(INDIR,subjects{jj}, strcat(filename,'_filtered')),'gui','off');
    CURR_FILTERED = CURRENTSET ; 
    
    % Extract ALL conditions epochs
    EEG = pop_epoch(EEG, conditions, win_of_interest, 'newname', strcat(filename,'_ALL'), 'epochinfo', 'yes');

    % Remove baseline
    EEG = pop_rmbase( EEG, baseline,[] );
    
    % Add channels information
    EEG=pop_chanedit(EEG, 'lookup','/Users/anne-sophiedubarry/Documents/4_Software/eeglab2020_0/plugins/dipfit/standard_BEM/elec/standard_1005.elc');

    % Select DEV
    [EEG_DEV1,target_indices1] = pop_selectevent(EEG,'type','DEV1'); % list of DEV1 relative to 900 trials 
    [EEG_DEV2,target_indices2] = pop_selectevent(EEG,'type','DEV2'); % list of DEV2 relative to 900 trials 
    [EEG_STD,target_indices_std] = pop_selectevent(EEG,'type','STD'); % list of STD relative to 900 trials 
    
    idx_std1 = target_indices_std(ismember(target_indices_std,target_indices1-1)); % list of STD preceeding DEV1 relative to 900 trials 
    idx_std2 = target_indices_std(ismember(target_indices_std,target_indices2-1)); % list of STD preceeding DEV2 relative to 900 trials 
    
    % Select all STD events in data
    [EEG_STD1,target_indices_std1] = pop_selectevent(EEG,'event',idx_std1); 
    [EEG_STD2,target_indices_std2] = pop_selectevent(EEG,'event',idx_std2);
   
    % Remove artifactual trials (STD, DEV)
    [EEG_STD1_thresh,idx_std1_rej] = pop_eegthresh(EEG_STD1,1,elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    [EEG_STD2_thresh,idx_std2_rej] = pop_eegthresh(EEG_STD2,1,elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    [EEG_DEV1_thresh,idx_dev1_rej] = pop_eegthresh(EEG_DEV1,1,elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
    [EEG_DEV2_thresh,idx_dev2_rej] = pop_eegthresh(EEG_DEV2,1,elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
   
    % Remove artifactual trials from all STD
    [EEG_STD_thresh, idx_removed] = pop_eegthresh(EEG_STD,1,elec ,rej_low, rej_high, win_of_interest(1), win_of_interest(2),0,1);
   
    std_good = setdiff(1:900,target_indices_std(idx_removed)); 
    begining_of_block = repelem((1:30:900)-1,3)+repmat(1:3,1,30); 
    
     % If nubmber of STD1 < number of DEV1 : randomly select other STD
    if length(EEG_DEV1_thresh.event)~=length(EEG_STD1_thresh.event)
       % Find number of trial to add in STD
       ntrials =  length(EEG_DEV1_thresh.event)-length(EEG_STD1_thresh.event) ;
       %Apply function to balance number of STD trial to reach the number of DEV trials 
       [EEG_STD1_thresh] = balance_number_of_STD(EEG,ntrials,std_good,target_indices1,begining_of_block,target_indices_std1,idx_std1_rej,idx_std1, idx_dev1_rej) ;
       
    end
    
    % If nubmber of STD2 < number of DEV2 : randomly select other STD
    if length(EEG_DEV2_thresh.event)~=length(EEG_STD2_thresh.event)
       % Find number of trial to add in STD
       ntrials =  length(EEG_DEV2_thresh.event)-length(EEG_STD2_thresh.event) ;
       %Apply function to balance number of STD trial to reach the number of DEV trials 
       [EEG_STD2_thresh] = balance_number_of_STD(EEG,ntrials,std_good,target_indices2,begining_of_block,target_indices_std2,idx_std2_rej,idx_std2, idx_dev2_rej) ; 

    end
     
    for cc=1:2 
        
        figure('Name',strcat('Subject :',subjects{jj},'Condition :',conditions{cc+1}),'Units','normalized','Position',[0,0,1,1]);
        
        EEG_DEV = eval(sprintf('EEG_DEV%d_thresh',cc)) ;
        EEG_STD = eval(sprintf('EEG_STD%d_thresh',cc)) ;
        nfig = 1 ; 

        for elec_letter=1:size(elec_to_disp_labels,1)
   
            for elec_numb=1:size(elec_to_disp_labels,2)

                hAxes = subplot(size(elec_to_disp_labels,1),size(elec_to_disp_labels,2),nfig) ; 
                nfig = nfig +1 ;

                idx_elec = find(ismember({EEG_DEV.chanlocs.labels},elec_to_disp_labels(elec_letter,elec_numb))) ; 

                % Compute grand average over one electrode
                grd_STD = squeeze(mean(EEG_STD.data(idx_elec,:,:),3)) ; 
                grd_DEV = squeeze(mean(EEG_DEV.data(idx_elec,:,:),3)) ; 
                grd_DIFF = grd_DEV - grd_STD ; 

                % Plot timeseries
                plot(EEG_STD.times,grd_STD,'g','Linewidth',1.5); hold on ;set(gca,'YDir','reverse') ; 
                plot(EEG_STD.times,grd_DEV,'r','Linewidth',1.5);  hold on; set(gca,'YDir','reverse') ;
                plot(EEG_STD.times,grd_DIFF,'k','Linewidth',1.5);  hold on; set(gca,'YDir','reverse') ;
                
                % Plot transparetn halo (+-mad)
                plotHaloPatchMAD(hAxes, EEG_STD.times, squeeze(EEG_STD.data(idx_elec,:,:)), [0,255,0]) ; 
                plotHaloPatchMAD(hAxes, EEG_DEV.times, squeeze(EEG_DEV.data(idx_elec,:,:)), [255,0,0]) ; 
           
                % Adjust graphics
                xlim([EEG_STD.xmin, EEG_STD.xmax]*1000); grid on ; 
                legend('STD (/DA/)',sprintf('DEV (/%s/)',cond_sylab{cc}),sprintf('DEV-STD (/%s/)',cond_sylab{cc}));
                title(elec_to_disp_labels(elec_letter,elec_numb));
                xlabel('Times (ms)'); ylabel('uV');
                set(hAxes,'Fontsize',FONTSZ);
                % To save data in vectoriel
%                 print('-dsvg',fullfile(INDIR,strcat(subjects{jj},'_',conditions{cc+1},'.svg'))) ; 
                
            end
        
        end
        % Export data
        [ALLEEG, EEG_DEV, CURRENTSET] = pop_newset(ALLEEG, EEG_DEV, CURRENTSET, 'setname',strcat(filename,'_','DEV',num2str(cc)),'savenew', fullfile(INDIR,subjects{jj},strcat(filename,'_','DEV',num2str(cc))),'gui','off');
        [ALLEEG, EEG_STD, CURRENTSET] = pop_newset(ALLEEG, EEG_STD, CURRENTSET, 'setname',strcat(filename,'_','STD',num2str(cc)),'savenew', fullfile(INDIR,subjects{jj},strcat(filename,'_','STD',num2str(cc))),'gui','off');
    end

end
    
%--------------------------------------------------------------
% FUNCTION that select from EEG_STD ntrial with no repetition with exisitng
% STD 
%--------------------------------------------------------------
function [EEG_STD_ALL] = balance_number_of_STD(EEG,ntrials,std_good,target_indices,begining_of_block,target_indices_std,idx_std_rej, idx_std,idx_dev_rej)
   
% We need to add some STD 
if ntrials > 0 
        % Pool of STD without those rejected by threshold detection
        pool_std = setdiff(std_good,target_indices);   
        
        % Pool of STD without beginners in block (3 first trials) 
        pool_std_w_no_beginners = setdiff(pool_std,begining_of_block);
        
        % Trials which were already selected 
        idx_std_already_included = setdiff(target_indices_std, target_indices_std(idx_std_rej)) ; 
     
        % Trials to add to balance the number of trial to the same number
        % as DEV
        idx_to_add = pool_std_w_no_beginners(randperm(length(pool_std_w_no_beginners),ntrials));
        
        % Select trial : 1) randomly a number = ntrial and 2) those which
        % were already selected and 'good'
        [EEG_STD_ALL,~] = pop_selectevent(EEG,'event',[idx_std_already_included idx_to_add]);

% We need to remove some STD
else
        %Find indices of rejected DEV1 in the 900 trials referential
        list_of_rej_dev = target_indices(idx_dev_rej);
        
        available_std =  setdiff(list_of_rej_dev-1,target_indices_std(idx_std_rej)) ; 
        
        %Select random indices amont the rejected DEV
        random_std_to_remove = available_std(randperm(length(available_std),abs(ntrials)));
       
        %Remove STD corresponding to rejected DEV
        idx_std_already_included = setdiff(target_indices_std, target_indices_std(idx_std_rej)) ; %problem : if there are rejected STD among the selected STD to remove!!
        new_list_of_std = setdiff(idx_std_already_included,random_std_to_remove);
        
        %Modify EEG struct
        [EEG_STD_ALL,~] = pop_selectevent(EEG,'event',[new_list_of_std]);
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