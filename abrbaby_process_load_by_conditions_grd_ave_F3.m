% Load EEGLAB 
% addpath(genpath('/Users/anne-sophiedubarry/Documents/4_Software/eeglab2020_0'));
tmp = pwd ; 
cd '/Users/anne-sophiedubarry/Documents/4_Software/eeglab2020_0' ; 
% Open eeglab
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

cd(tmp) ; 

% Modify preferences in order to be able to load multiple datasets 
pop_editoptions( 'option_storedisk', 1);

% Set filepath (must contain .bdf and .txt files from recording)
% INDIR = '/Volumes/ASDUB/DEVLANG_data';
INDIR = '/Users/anne-sophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data/STUDY';

% Reads all folders that are in INDIR 
d = dir(INDIR); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

conditions = {'STD','DEV1','DEV2'} ; 
commands = {}; % initialize STUDY dataset list

% Loop through all of the subjects in the study to create the dataset
for loopnum = 1:length(subjects) %for each subject
    
    STDFile = fullfile(INDIR,subjects{loopnum},strcat(subjects{loopnum},'_STD.set')) ; 
    EEG_std = pop_loadset(STDFile) ; 
    STD_subj(loopnum,:,:)  = mean(EEG_std.data,3) ; 
    
    DEV1File = fullfile(INDIR,subjects{loopnum},strcat(subjects{loopnum},'_DEV1.set')) ; 
    EEG_dev1 = pop_loadset(DEV1File) ; 
    DEV1_subj(loopnum,:,:)  = mean(EEG_dev1.data,3) ; 
    
    DEV2File = fullfile(INDIR,subjects{loopnum},strcat(subjects{loopnum},'_DEV2.set')) ; 
    EEG_dev2 = pop_loadset(DEV2File) ; 
    DEV2_subj(loopnum,:,:)  = mean(EEG_dev2.data,3) ; 

end

% ii= 5 corresponds to electrode F3 
grd_STD_F3 = squeeze(mean(STD_subj(:,5,:),1)) ; 
grd_DEV1_F3 = squeeze(mean(DEV1_subj(:,5,:),1)) ; 
grd_DEV2_F3 = squeeze(mean(DEV2_subj(:,5,:),1)) ; 

figure ; 
plot(EEG_dev2.times,grd_STD_F3,'k','Linewidth',1.5); hold on ;set(gca,'YDir','reverse') ;
plot(EEG_dev2.times,grd_DEV1_F3,'r','Linewidth',1.5);  hold on; set(gca,'YDir','reverse') ;
plot(EEG_dev2.times, grd_DEV2_F3,'b','Linewidth',1.5); set(gca,'YDir','reverse') ;
grid on ; 

legend('STD','DEV1','DEV2');
xlabel('Times (ms)'); ylabel('uV'); title ('Grand average F3');
