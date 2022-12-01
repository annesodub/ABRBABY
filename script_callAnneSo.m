

%% Set environment 
% Variables to enter manually before running the code
eeglab_path = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/dev/signal_processing/ABRBABY/eeglab2021.1' ; 
erplab_path = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/dev/signal_processing/ABRBABY/erplab8.30';
biosig_installer_path = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/dev/signal_processing/ABRBABY/biosig4octmat-3.8.0/biosig_installer.m' ; 

% Load path and start Matlab : returns ALLEEG (EEGLAB structure)
ALLEEG = prep_and_start_environement(eeglab_path, biosig_installer_path, erplab_path) ;

%% Preprocess 
indir = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data/DEVLANG_data/' ;

% Set variables for preproc
hp = 1; % high-pass (Hz) (APICE)
lp = 30; % low-pass (Hz) (APICE) 
mastos = {'Lmon','Rmon','MASTOG','MASTOD'}; % REF
trig = {'Erg1'};
baseline = [-99, 0] ; 
win_of_interest = [-0.1, 0.5] ; 
conditions = {'STD','DEV1','DEV2'} ; 

EEG = reref_filter_epoch(indir, hp,lp, mastos, trig,baseline,win_of_interest);

plots_dir = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data/png_folder';

STD_number = 1;


abrbaby_process_ERP_sanity_exportdata_allSTD(eeglab_path, biosig_installer_path,indir,plots_dir,chan_dir,STD_number);
% abrbaby_process_ERP_sanity(eeglab_path, biosig_installer_path,indir) ;


FFR_analysis(indir) ;


% HERE I CAN EXECUTE ANY FUNCTION OF THE PACKAGE WITH MY PATHS

% ALLEEG = prep_and_start_environnement() : load path and start Matlab : returns ALLEEG 
% one_subject_process() % sanity + save processed data : en header documenter etape de proc 
% all_subjects_process() : skip subject if data exist 
% display_one_subjects()
% display : creer des fonctions de displays