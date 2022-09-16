% Load EEGLAB 
% addpath(genpath('/Users/anne-sophiedubarry/Documents/4_Software/eeglab2020_0'));
% tmp = pwd ; 
% cd '/Users/anne-sophiedubarry/Documents/4_Software/eeglab2020_0' ; 

%Set channel indice and label
%3=F4; 4=Fz; 5=F3; 7=C3; 8=Cz; 9=C4 
chan_ind = [5];
chan_lab = 'F3';

% Open eeglab
[ALLEEG, EEG, CURRENTSET, ALLCOM] = eeglab;

%cd(tmp) ; 

% Modify preferences in order to be able to load multiple datasets 
pop_editoptions( 'option_storedisk', 1);

% Set filepath (must contain .bdf and .txt files from recording)
% INDIR = '/Volumes/ASDUB/DEVLANG_data';
INDIR = '\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\DEVLANG_data';
%INDIR ='C:\Users\hervé\Documents\data_test';
%INDIR = '\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\DEVLANG_data\DATA_6-10';
%INDIR = '\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\DEVLANG_data\DATA_18-24';
%INDIR = '\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\DEVLANG_data_issues';


% Reads all folders that are in INDIR 
d = dir(INDIR); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

conditions = {'STD1', 'STD2', 'DEV1','DEV2'} ; 
commands = {}; % initialize STUDY dataset list

grpA.suffix = {'_T3','_T6','_T8','_T10'};
grpB.suffix = {'_T18','_T24'};

[grpA.DEV1_avg, grpA.DEV2_avg, grpA.STD1_avg, grpA.STD2_avg] = extract_averages_DEV_STD(INDIR,subjects( contains(subjects,grpA.suffix)));
[grpB.DEV1_avg, grpB.DEV2_avg, grpB.STD1_avg, grpB.STD2_avg] = extract_averages_DEV_STD(INDIR,subjects( contains(subjects,grpB.suffix)));

display_comparison_by_group(grpA,grpB) ;
%TODO