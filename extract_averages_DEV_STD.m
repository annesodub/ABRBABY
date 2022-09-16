function [DEV1_subj, DEV2_subj,STD1_subj, STD2_subj] = extract_averages_DEV_STD(INDIR,subjects)

% Loop through all of the subjects in the study to create the dataset
for loopnum = 1:length(subjects) %for each subject

    %dimensions datasets (.set)  subjects x channels x timepoints

    STD1File = fullfile(INDIR,subjects{loopnum},strcat(subjects{loopnum},'_STD1.set')) ; 
    EEG_std1 = pop_loadset(STD1File) ; 
    STD1_subj(loopnum,:,:)  = mean(EEG_std1.data,3) ;

    %apply low pass 10Hz (visual display)

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

% grd_STD = squeeze(mean(STD_subj(:,chan_ind,:),1)) ; 
% grd_DEV1 = squeeze(mean(DEV1_subj(:,chan_ind,:),1)) ; 
% grd_DEV2 = squeeze(mean(DEV2_subj(:,chan_ind,:),1)) ;

%3=F4; 4=Fz; 5=F3; 7=C3; 8=Cz; 9=C4 
chan_ind = [5];
chan_lab = 'F3';

%Mean of F3, Fz & F4: ?
grd_STD1 = squeeze(mean(mean(STD1_subj(:,chan_ind,:),1),2)) ; 
grd_STD2 = squeeze(mean(mean(STD2_subj(:,chan_ind,:),1),2)) ; 
grd_DEV1 = squeeze(mean(mean(DEV1_subj(:,chan_ind,:),1),2)) ; 
grd_DEV2 = squeeze(mean(mean(DEV2_subj(:,chan_ind,:),1),2)) ; 

%Difference DEV-STD
diff_DEV1_STD1 = grd_DEV1-grd_STD1;
diff_DEV2_STD2 = grd_DEV2-grd_STD2;

figure ; 
plot(EEG_dev1.times,grd_STD1,'k','Linewidth',1.5); hold on ;set(gca,'YDir','reverse') ;
plot(EEG_dev1.times, grd_DEV1,'r','Linewidth',1.5); set(gca,'YDir','reverse') ;
plot(EEG_dev1.times, diff_DEV1_STD1,'g','Linewidth',1.5); set(gca,'YDir','reverse') ;
grid on ; 
legend('STD (/DA/)','DEV1 (/BA/)','DEV1-STD (/BA/)');
xlabel('Times (ms)'); ylabel('uV'); %title ([filename, ' Grand average ', chan_lab], 'Interpreter', 'None');


figure ; 
plot(EEG_dev2.times,grd_STD2,'k','Linewidth',1.5); hold on ;set(gca,'YDir','reverse') ;
plot(EEG_dev2.times,grd_DEV2,'b','Linewidth',1.5);  hold on; set(gca,'YDir','reverse') ;
plot(EEG_dev2.times,diff_DEV2_STD2,'c','Linewidth',1.5);  hold on; set(gca,'YDir','reverse') ;
grid on ; 
legend('STD (/DA/)','DEV2 (/GA/)','DEV2-STD (/GA/)');
xlabel('Times (ms)'); ylabel('uV'); %title ([filename, ' Grand average ', chan_lab],'Interpreter', 'None');

% %legend('STD','DEV1-STD','DEV2-STD');
end
