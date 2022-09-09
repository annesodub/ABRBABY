%% ERPs analysis script - Estelle Herv� - 2022 - %80PRIME Project

%% Variables to enter manually before running the code

% Load EEGLAB 
% addpath(genpath('/Users/anne-sophiedubarry/Documents/4_Software/eeglab2020_0'));
%tmp = pwd ; 
%cd '/Users/anne-sophiedubarry/Documents/4_Software/eeglab2020_0' ;
%cd '\\Filer\home\Invites\herv�\Mes documents\MATLAB' ;

% % Open eeglab
% [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

%cd(tmp) ; 

% Set filepath (must contain .bdf and .txt files from recording)
INDIR = '/Users/anne-sophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data';
% INDIR = '\\Filer\home\Invites\herv�\Mes documents\These\EEG\Data\DEVLANG_data';

% Reads all folders that are in INDIR 
d = dir(INDIR); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Set variables for filtering
hp = 0.1; %value for high-pass filter (Hz) (APICE)
lp = 40; %value for low-pass filter (Hz) (APICE) 

% Rejection treshold for bad epochs
rej_low = -150; %150 infants; 120 adults
rej_high = 150; %150 infants; 120 adults

% List of channel labels to reref with 
mastos = {'Lmon','Rmon'};
trig = {'Erg1'};

baseline = [-99, 0] ; 
elec = [1:16];

% Loop through subjects
%for jj=11:11%length(subjects) 
    
    jj=find(ismember(subjects,'DVL_016_T18')) ; 

    %% IMPORT
    % Open eeglab
    [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
    
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
    %EEG.events = read_custom_events(strrep(fullfile(fname.folder,fname.name),'.bdf','.txt'),EEG.event) ;
    
    txt_filename = [filename, '.txt'];
   
    % Read .txt 
    my_events = readtable(fullfile(INDIR,subjects{jj},txt_filename), 'ReadVariableNames', 0);
   
    % Insert info from .txt into EEG.event
    my_events = table2array(my_events);

    temp = struct('latency', {EEG.event(:).latency}, ...
                'type', (my_events(:))',...
                'urevent', {EEG.event(:).urevent});
    
    EEG.event = temp;   
    
    %% FILTERS
    % Filters the data with ERPLab
    bdf_filename = fname.name ; 
    filename_filter = strrep(bdf_filename,'.bdf','_filtered');
    EEG  = pop_basicfilter(EEG,  elec , 'Boundary', 'boundary', 'Cutoff', [hp lp], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2, 'RemoveDC', 'on' );
    [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 2, 'savenew', filename_filter,'gui','off');
    
    %filename = strrep(bdf_filename,'.bdf','');
    std_filename = strrep(bdf_filename,'.bdf','_STD');
    dev1_filename = strrep(bdf_filename,'.bdf','_DEV1');
    dev2_filename = strrep(bdf_filename,'.bdf','_DEV2');
    
    %% EPOCHING
    %% ERPs for STD

    % Extract epochs for STD
    EEG = pop_epoch( EEG, {  'STD'  }, [-0.1         0.5], 'newname', std_filename, 'epochinfo', 'yes');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 7,'setname',std_filename,'gui','off'); 
    EEG = eeg_checkset( EEG );

    % Baseline removal
    EEG = pop_rmbase( EEG, [-99 0],[] );
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 8,'gui','off');  

    % Reject data epochs (extreme values)
    EEG = eeg_checkset( EEG );
    EEG = pop_eegthresh(EEG,1,elec ,rej_low, rej_high,-0.099609,0.49805,0,1);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 9,'gui','off'); 

    % Indicate channel location
    EEG = eeg_checkset( EEG );
    EEG=pop_chanedit(EEG, 'lookup','\\Filer\\home\\Invites\\herv�\\Mes documents\\MATLAB\\eeglab2021.1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

    % Plot ERPs
    EEG = eeg_checkset( EEG );
    figure; pop_plottopo(EEG, elec , 'STD', 0, 'ydir',-1);
    eeglab redraw;

    %Save dataset
    EEG = pop_saveset( EEG, 'filename',std_filename,'filepath', fname.folder);
    %EEG = pop_saveset( EEG, 'filename','DVL_004_T8_STDs.set','filepath','\\\\Filer\\home\\Invites\\herv�\\Mes documents\\These\\EEG\\Data\\DEVLANG_data\\DVL_004_T8\\');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

    %% ERPs for DEV1

    % Extract epochs for DEV1
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 6,'retrieve',2, 'study',0); 
    EEG = eeg_checkset( EEG );
    EEG = pop_epoch( EEG, {  'DEV1'  }, [-0.1         0.5], 'newname', dev1_filename, 'epochinfo', 'yes');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 10,'setname',dev1_filename,'gui','off'); 
    EEG = eeg_checkset( EEG );

    % Baseline removal
    EEG = pop_rmbase( EEG, [-99 0],[] );
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 11,'gui','off'); 

    % Reject data epochs (extreme values)
    EEG = eeg_checkset( EEG );
    EEG = pop_eegthresh(EEG,1, elec , rej_low, rej_high,-0.099609,0.49805,0,1);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 12,'gui','off');

    % Channel location
    EEG = eeg_checkset( EEG );
    EEG=pop_chanedit(EEG, 'lookup','\\Filer\\home\\Invites\\herv�\\Mes documents\\MATLAB\\eeglab2021.1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

    %Plot ERPs
    EEG = eeg_checkset( EEG );
    figure; pop_plottopo(EEG, elec , 'DEV1', 0, 'ydir',-1)

    %Save dataset
    EEG = pop_saveset( EEG, 'filename',dev1_filename,'filepath',fname.folder);
    %EEG = pop_saveset( EEG, 'filename','DVL_004_T8_STDs.set','filepath','\\\\Filer\\home\\Invites\\herv�\\Mes documents\\These\\EEG\\Data\\DEVLANG_data\\DVL_004_T8\\');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

    %% ERPs for DEV2

    % Extract epochs for DEV2
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 13,'retrieve',2,'study',0); 
    EEG = eeg_checkset( EEG );
    EEG = pop_epoch( EEG, {  'DEV2'  }, [-0.1         0.5], 'newname', dev2_filename, 'epochinfo', 'yes');
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 13,'setname', dev2_filename,'gui','off'); 

    % Baseline removal
    EEG = eeg_checkset( EEG );
    EEG = pop_rmbase( EEG, [-99 0],[] );
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 14,'gui','off'); 

    % Reject data epochs (extreme values)
    EEG = eeg_checkset( EEG );
    EEG = pop_eegthresh(EEG,1,elec, rej_low, rej_high,-0.099609,0.49805,0,1);
    [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 15,'gui','off'); 

    % Indicate channel location
    EEG = eeg_checkset( EEG );
    EEG=pop_chanedit(EEG, 'lookup','\\Filer\\home\\Invites\\herv�\\Mes documents\\MATLAB\\eeglab2021.1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

    %Plot ERPs
    EEG = eeg_checkset( EEG );
    figure; pop_plottopo(EEG, elec, 'DEV2', 0, 'ydir',-1);

    %Save dataset
    EEG = pop_saveset( EEG, 'filename',dev2_filename,'filepath',fname.folder);
    %EEG = pop_saveset( EEG, 'filename','DVL_004_T8_STDs.set','filepath','\\\\Filer\\home\\Invites\\herv�\\Mes documents\\These\\EEG\\Data\\DEVLANG_data\\DVL_004_T8\\');
    [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);[ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);

    eeglab redraw;

     
     %%OLD CODE
%     % Extract epochs for STD
%     EEG = pop_epoch(EEG, {'STD'}, [-0.1, 0.5], 'newname', strcat(filename,'_STD'), 'epochinfo', 'yes');
%     
%     % Save dataset in EEGLAB 
%     [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname',strcat(filename,'_STD'),'gui','off');
%      
%     %% BASELINE
%     % Baseline removal
%     EEG = pop_rmbase( EEG, [-99 0],[] );
%     
%     % Reject data epochs (extreme values)
%     EEG = pop_eegthresh(EEG,1,[1:16] ,rej_low, rej_high,-0.099609,0.49805,0,1);
%     
%     % Indicate channel location
%     EEG=pop_chanedit(EEG, 'lookup','\\Filer\\home\\Invites\\herv�\\Mes documents\\MATLAB\\eeglab2021.1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc');
%     
%     % Plot ERPs
%     figure; pop_plottopo(EEG, [1:16] , 'STD', 0, 'ydir',-1);
%     
%     eeglab redraw;

    %%%%% TO FINISH (following APICE?)
    
%end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% function out_event = read_custom_events(fname, in_event) 
% 
% % Read .txt 
% my_events = readtable(fname, 'ReadVariableNames', 0);
% 
% % Insert info from .txt into EEG.event
% my_events = table2array(my_events);
% 
% out_event = struct('latency', {in_event(:).latency}, ...
%                 'type', (my_events(:))',...
%                 'urevent', {in_event(:).urevent});
% 
% end

%  
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% STARTING FROM HERE : old code
% 
% %Set tresholds for bad event detection:
% %min_thr= XX; %double detection
% %max_thr=XX; %outliers
% 
% %% Set variables based on previous section
% 
% % %Set filnames
% 
% % filename = ['DVL_', subj_n, '_',subj_age];
% % filename_trig = ['DVL_', subj_n, '_', subj_age, '_trig'];
% % filename_filter = ['DVL_', subj_n, '_', subj_age, '_filtered'];
% 
% % filename = 'pilote_salome_170';
% % filename_trig = 'pilote_salome_170_trig';
% % filename_filter = 'pilote_salome_170_filtered';
% % 
% filename = 'BBL_RK_T0';
% filename_trig = 'BBL_RK_T0_trig';
% filename_filter = 'BBL_RK_T0_filtered';
% 
% %txt_filename = ['DVL_', subj_n, '_', subj_age, '.txt'];
% %bdf_filename = ['DVL_', subj_n, '_', subj_age, '.bdf']; 
% 
% % txt_filename = 'pilote_salome_170.txt';
% % bdf_filename = 'pilote_salome_170.bdf'; 
% 
% txt_filename = 'BBL_RK_T0.txt';
% bdf_filename = 'BBL_RK_T0.bdf';
% 
% %% From loading file to extract event from data channel
% 
% % Open eeglab
% [ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;
% 
% % Select bdf file in the folder
% EEG = pop_biosig(fullfile(filepath, bdf_filename));
% 
% % Rename file
% [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 1,'setname',filename,'gui','off');
% EEG = eeg_checkset( EEG );
% 
% % Rereference data and rename new file
% %EEG = pop_reref( EEG, [33 34]); %adapted for pilote_clem
% EEG = pop_reref( EEG, ref_elec);
% [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 2,'setname',filename,'gui','off');
% EEG = eeg_checkset( EEG );
% 
% % Select channels to keep
% EEG = eeg_checkset( EEG );
% %EEG = pop_select( EEG, 'channel',{'Fp1','Fp2','F4','Fz','F3','T7','C3','C4','T8','P4','Pz','P3','O1','Oz','O2','Erg1'}); %here Cz is missing !!
% EEG = pop_select( EEG, 'channel',{'Fp1','Fp2','F4','Fz','F3','T7','C3','Cz','C4','T8','P4','Pz','P3','O1','Oz','O2','Erg1'}); %with Cz
% [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 3, 'setname', filename, 'gui','off');
% 
% % Extract event from trigger channel (Erg1)
% %EEG = pop_chanevent(EEG, trigg_elec,'oper','X>40000','edge','leading','edgelen',1); %better for DVL_005_T18
% EEG = pop_chanevent(EEG, trigg_elec,'oper','X>20000','edge','leading','edgelen',1);
% [ALLEEG, EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% EEG = eeg_checkset( EEG );
% pop_eegplot( EEG, 1, 1, 1);
% [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 4,'setname',filename_trig,'gui','off');
% 
% eeglab redraw;
% %EEG.event(find(diff([EEG.event.latency])<1000)+1) = [] ;
% %EEG.event(find(diff([EEG.event.latency])>A DEFINIR))= [] ;
% 
% %-----STARTING HERE : OLD CODE
% % bad_events = [];
% % ind = 1;
% % for ii = 1:size(EEG.event,2)-1
% %     if EEG.event(ii).latency-EEG.event(ii+1).latency > -3500  % -1000 for 40 ms; -3500 for 170ms
% %         bad_events(ind) = EEG.event(ii).urevent;
% %         ind = ind+1;
% %     elseif EEG.event(ii).latency-EEG.event(ii+1).latency < -150000 
% %         bad_events(ind) = EEG.event(ii).urevent;
% %         ind = ind+1;
% %         end
% %     ii=ii+1;
% % end
% % 
% % new_latency = [];
% % new_type = [];
% % new_urevent = [];
% % indx = 1;
% % 
% % for pp = 1:size(EEG.event,2)
% %    if ismember(EEG.event(pp).urevent, bad_events(:))== 0
% %        new_latency(indx) = EEG.event(pp).latency;
% %        new_type(indx) = EEG.event(pp).type;
% %        new_urevent(indx) = EEG.event(pp).urevent;
% %        indx = indx+1;
% %    end
% % end
% % 
% % new_latency = num2cell(new_latency);
% % new_type = num2cell(new_type);
% % new_urevent = num2cell(new_urevent);
% % 
% % EEG.event = struct('latency', (new_latency(:))', ...
% %                 'type', (new_type(:))',...
% %                 'urevent', (new_urevent(:))');
% % 
% % eeglab redraw;
% %%
% % %% Reject bad events (manually), fonctionne mais � optimiser (josephT0)
% % EEG = eeg_checkset( EEG );
% % EEG = pop_editeventvals(EEG,'changefield',{1,'latency',250},'delete',1658,'delete',1658,'delete',1658,'delete',1658);
% % [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% 

% 
% 
% %% Resampling and filtering
% 
% %Resample data
% EEG = pop_resample( EEG, 512);
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 5,'gui','off'); 
% 
% % Filter data
% EEG  = pop_basicfilter( EEG,  elec , 'Boundary', 'boundary', 'Cutoff', [hp lp], 'Design', 'butter', 'Filter', 'bandpass', 'Order',  2, 'RemoveDC', 'on' );
% [ALLEEG, EEG, CURRENTSET] = pop_newset(ALLEEG, EEG, 6, 'savenew', filename_filter,'gui','off');
% 
% % Plot data
% EEG = eeg_checkset( EEG );
% pop_eegplot( EEG, 1, 1, 1);
% 
% eeglab redraw;
% 
% %% ERPs for STD
% 
% % Extract epochs for STD
% EEG = pop_epoch( EEG, {  'STD'  }, [-0.1         0.5], 'newname', 'clement_pilote_rej resampled epochs', 'epochinfo', 'yes');
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 7,'setname','clement_pilote_rej resampled epochs STD','gui','off'); 
% EEG = eeg_checkset( EEG );
% 
% % Baseline removal
% EEG = pop_rmbase( EEG, [-99 0],[] );
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 8,'gui','off');  
% 
% % Reject data epochs (extreme values)
% EEG = eeg_checkset( EEG );
% EEG = pop_eegthresh(EEG,1,[1:16] ,rej_low, rej_high,-0.099609,0.49805,0,1);
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 9,'gui','off'); 
% 
% % Indicate channel location
% EEG = eeg_checkset( EEG );
% EEG=pop_chanedit(EEG, 'lookup','\\Filer\\home\\Invites\\herv�\\Mes documents\\MATLAB\\eeglab2021.1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc');
% [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% 
% % Plot ERPs
% EEG = eeg_checkset( EEG );
% figure; pop_plottopo(EEG, [1:16] , 'STD', 0, 'ydir',-1);
% eeglab redraw;
% 
% 
% %% ERPs for DEV1
% 
% % Extract epochs for DEV1
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 10,'retrieve',6, 'study',0); 
% EEG = eeg_checkset( EEG );
% EEG = pop_epoch( EEG, {  'DEV1'  }, [-0.1         0.5], 'newname', 'clement_pilote_rej resampled epochs', 'epochinfo', 'yes');
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 10,'setname','clement_pilote_rej resampled epochsDEV1','gui','off'); 
% EEG = eeg_checkset( EEG );
% 
% % Baseline removal
% EEG = pop_rmbase( EEG, [-99 0],[] );
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 11,'gui','off'); 
% 
% % Reject data epochs (extreme values)
% EEG = eeg_checkset( EEG );
% EEG = pop_eegthresh(EEG,1,[1:16] , rej_low, rej_high,-0.099609,0.49805,0,1);
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 12,'gui','off');
% 
% % Channel location
% EEG = eeg_checkset( EEG );
% EEG=pop_chanedit(EEG, 'lookup','\\Filer\\home\\Invites\\herv�\\Mes documents\\MATLAB\\eeglab2021.1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc');
% [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% 
% %Plot ERPs
% EEG = eeg_checkset( EEG );
% figure; pop_plottopo(EEG, [1:16] , 'DEV1', 0, 'ydir',-1)
% 
% %% ERPs for DEV2
% 
% % Extract epochs for DEV2
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 13,'retrieve',6,'study',0); 
% EEG = eeg_checkset( EEG );
% EEG = pop_epoch( EEG, {  'DEV2'  }, [-0.1         0.5], 'newname', 'clement_pilote_rej resampled epochs', 'epochinfo', 'yes');
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 13,'setname','clement_pilote_rej resampled epochsDEV2','gui','off'); 
% 
% % Baseline removal
% EEG = eeg_checkset( EEG );
% EEG = pop_rmbase( EEG, [-99 0],[] );
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 14,'gui','off'); 
% 
% % Reject data epochs (extreme values)
% EEG = eeg_checkset( EEG );
% EEG = pop_eegthresh(EEG,1,[1:16], rej_low, rej_high,-0.099609,0.49805,0,1);
% [ALLEEG EEG CURRENTSET] = pop_newset(ALLEEG, EEG, 15,'gui','off'); 
% 
% % Indicate channel location
% EEG = eeg_checkset( EEG );
% EEG=pop_chanedit(EEG, 'lookup','\\Filer\\home\\Invites\\herv�\\Mes documents\\MATLAB\\eeglab2021.1\\plugins\\dipfit\\standard_BEM\\elec\\standard_1005.elc');
% [ALLEEG EEG] = eeg_store(ALLEEG, EEG, CURRENTSET);
% 
% %Plot ERPs
% EEG = eeg_checkset( EEG );
% figure; pop_plottopo(EEG, [1:16] , 'DEV2', 0, 'ydir',-1);
% 
% eeglab redraw;
% 
% %% Compare DEV1 and STD / DEV2 and STD 
% EEG = eeg_checkset( EEG );
% pop_comperp( ALLEEG, 1, [12 15] ,[9 9] ,'addavg','off','addstd','off','addall','on','subavg','on','diffavg','off','diffstd','off','chans',[1:13 15] ,'tplotopt',{'ydir',-1});
% %pop_comperp( ALLEEG, 1, [7 10] ,[4 4] ,'addavg','off','addstd','off','addall','on','subavg','on','diffavg','off','diffstd','off','chans',[1:13 15] ,'tplotopt',{'ydir',-1});
% 
% %%
% %difference DEV1-STD
% EEG = eeg_checkset( EEG );
% pop_comperp( ALLEEG, 1, 12, 9,'addavg','on','addstd','off','subavg','on','diffavg','on','diffstd','off','tplotopt',{'ydir',-1});
% %difference DEV2-STD
% EEG = eeg_checkset( EEG );
% pop_comperp( ALLEEG, 1, 15, 9,'addavg','on','addstd','off','subavg','on','diffavg','on','diffstd','off','tplotopt',{'ydir',-1});
% 
% %% Supplement
% 
% % erp = mean(EEG.data(1,:,:),3);
% % EEG.data = erp;
% 
% %Export events as txt file:
% fid = fopen('events_KRT0.txt','w');
% fprintf(fid,'%c\n',EEG.event.latency);
% fclose(fid);