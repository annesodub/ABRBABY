function [DEV1_subj, DEV2_subj,STD1_subj, STD2_subj, timepoints] = extract_averages_DEV_STD(INDIR,subjects)
%function [DEV1_subj, DEV2_subj,STD1_subj, STD2_subj, grd_STD1, grd_STD2, grd_DEV1, grd_DEV2, diff_DEV1_STD1, diff_DEV2_STD2, mean_STD1_STD2, timepoints] = extract_averages_DEV_STD(INDIR,subjects)

% Loop through all of the subjects in the study to create the dataset
for loopnum = 1:length(subjects) %for each subject

    %dimensions datasets (.set) =  subjects x channels x timepoints

    STD1File = fullfile(INDIR,subjects{loopnum},strcat(subjects{loopnum},'_STD1.set')) ; 
    EEG_std1 = pop_loadset(STD1File) ; 
    STD1_subj(loopnum,:,:)  = mean(EEG_std1.data,3) ;

    STD2File = fullfile(INDIR,subjects{loopnum},strcat(subjects{loopnum},'_STD2.set')) ; 
    EEG_std2 = pop_loadset(STD2File) ; 
    STD2_subj(loopnum,:,:)  = mean(EEG_std2.data,3) ;

    DEV1File = fullfile(INDIR,subjects{loopnum},strcat(subjects{loopnum},'_DEV1.set')) ; 
    EEG_dev1 = pop_loadset(DEV1File) ; 
    DEV1_subj(loopnum,:,:)  = mean(EEG_dev1.data,3) ; 

    DEV2File = fullfile(INDIR,subjects{loopnum},strcat(subjects{loopnum},'_DEV2.set')) ; 
    EEG_dev2 = pop_loadset(DEV2File) ; 
    DEV2_subj(loopnum,:,:)  = mean(EEG_dev2.data,3) ; 

end

% chan_indices = {'Fp1','Fp2','F4','Fz','F3','T7','C3','Cz','C4','T8','P4','Pz','P3','O1','Oz','O2','Lmon','Ref','Rmon','Left','Right'}
% 
% for chan_ind = 1:length(chan_indices)
%     grd_STD1 = squeeze(mean(mean(STD1_subj(:,chan_ind,:),1),2)) ; 
%     grd_STD2 = squeeze(mean(mean(STD2_subj(:,chan_ind,:),1),2)) ;
%     grd_DEV1 = squeeze(mean(mean(DEV1_subj(:,chan_ind,:),1),2)) ; 
%     grd_DEV2 = squeeze(mean(mean(DEV2_subj(:,chan_ind,:),1),2)) ; 
%     
% %Difference DEV1-STD1
% diff_DEV1_STD1 = grd_DEV1 - grd_STD1;
% 
% %Difference DEV2-STD2 for grpA
% diff_DEV2_STD2 = grd_DEV2 - grd_STD2;
% 
% %Mean STD1+STD2 for grpA
% mean_STD1_STD2 = (grd_STD1 + grd_STD2)/2;
% 
% end


timepoints = EEG_std1.times;

end
