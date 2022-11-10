
% HERE I CAN EXECUTE ANY FUNCTION OF THE PACKAGE WITH MY PATHS

%% Variables to enter manually before running the code
eeglab_path = '\\Filer\home\Invites\hervé\Mes documents\MATLAB\eeglab2021.1';
biosig_installer_path = '\\Filer\home\Invites\hervé\Mes documents\MATLAB\eeglab2021.1\plugins\Biosig3.7.9\biosig_installer.m'; 
indir = '\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\DEVLANG_data';
%indir = '\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\DEVLANG_DATA_NEW';

%% Sanity checks:

%indir = '\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\DEVLANG_DATA_NEW';

abrbaby_process_ERP_sanity_exportdata(eeglab_path, biosig_installer_path,indir) ;

abrbaby_process_ERP_sanity_exportdata_allSTD(eeglab_path, biosig_installer_path,indir);

abrbaby_sanity_check_FFR(eeglab_path, biosig_installer_path,indir);

beep;
disp("Done");

%% Analyses:

%indir = '\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\DEVLANG_data';

%abrbaby_process_visual_comparison_by_group_age(eeglab_path, biosig_installer_path,indir);

%FFR_analysis(indir) ;

%beep;
%disp("Done");


