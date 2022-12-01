

%% Variables to enter manually before running the code
eeglab_path = '/Users/annesophiedubarry/Documents/4_Software/eeglab2021.1' ; 
erplab_path = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/dev/signal_processing/ABRBABY/erplab8.30';
biosig_installer_path = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/dev/signal_processing/ABRBABY/biosig4octmat-3.8.0/biosig_installer.m' ; 
indir = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data/DEVLANG_data/' ;
plots_dir = '/Users/annesophiedubarry/Documents/0_projects/in_progress/ABRBABY_cfrancois/data/png_folder';
chan_dir = fullfile(eeglab_path,'plugins/dipfit/standard_BEM/elec/standard_1005.elc');
STD_number = 1;
abrbaby_process_ERP_sanity_exportdata_allSTD(eeglab_path, biosig_installer_path,indir,plots_dir,chan_dir,STD_number);
% abrbaby_process_ERP_sanity(eeglab_path, biosig_installer_path,indir) ;

FFR_analysis(indir) ;


% HERE I CAN EXECUTE ANY FUNCTION OF THE PACKAGE WITH MY PATHS

% prep_and_start_environnement() : load path and start Matlab : returns ALLEEG 
% one_subject_process() % sanity + save processed data : en header documenter etape de proc 
% all_subjects_process() : skip subject if data exist 
% display_one_subjects()
% display : creer des fonctions de displays