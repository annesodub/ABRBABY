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
    DEV1File = fullfile(INDIR,subjects{loopnum},strcat(subjects{loopnum},'_DEV1.set')) ; 
    DEV2File = fullfile(INDIR,subjects{loopnum},strcat(subjects{loopnum},'_DEV2.set')) ; 
    
    commands = {commands{:} ...
    {'index' 3*loopnum-2 'load' STDFile 'subject' subjects{loopnum} 'condition' 'STD'} ...
    {'index' 3*loopnum-1 'load' DEV1File 'subject' subjects{loopnum} 'condition' 'DEV1'} ...
    {'index' 3*loopnum 'load' DEV2File 'subject' subjects{loopnum} 'condition' 'DEV2'}};
end

% Uncomment the line below to select ICA components with less than 15% residual variance
% commands = {commands{:} {'dipselect', 0.15}};
[STUDY, ALLEEG] = std_editset(STUDY, ALLEEG, 'name','ABRBABY_ERP','commands',commands,'updatedat','on');

% Update workspace variables and redraw EEGLAB
CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];
[STUDY, ALLEEG] = std_checkset(STUDY, ALLEEG);
eeglab redraw

% Precompute ERP 
[STUDY, ALLEEG] = std_precomp(STUDY, ALLEEG, {},'savetrials','on','interp','on','recompute','on','erp','on','erpparams',{'rmbase',[-99 0] });

% [STUDY ALLEEG] = std_precomp(STUDY, ALLEEG, {'Fp1','Fp2','F4','Fz','F3','T7','C3','Cz','C4','T8','P4','Pz','P3','O1','Oz','O2','Lmon','Ref','Rmon','Left','Right'}, 'interp', 'on', 'recompute','on','erp','on');
% 
% [STUDY ALLEEG] = std_precomp(STUDY, ALLEEG, 'components', 'erp', 'on', 'rmbase',[-99 0] , 'scalp', 'on');

% Precompute ERP 
tmpchanlocs = ALLEEG(1).chanlocs; 

STUDY = std_erpplot(STUDY, ALLEEG, 'channels', { tmpchanlocs.labels }, 'topoplotopt',{'ydir',-1},'plotconditions', 'together','ydir',-1);

% STUDY = pop_erpparams(STUDY, 'plotconditions','together');
% STUDY = std_erpplot(STUDY,ALLEEG,'channels',{'Fp1','Fp2','F4','Fz','F3','T7','C3','Cz','C4','T8','P4','Pz','P3','O1','Oz','O2','Lmon','Ref','Rmon','Left','Right'}, 'design', 1);
% CURRENTSTUDY = 1; EEG = ALLEEG; CURRENTSET = [1:length(EEG)];


