%% FFR analysis script - Estelle HERVE - 80PRIME DEVLANG project

INDIR = '\\Filer\home\Invites\hervé\Mes documents\These\EEG\Data\DEVLANG_data';

% Reads all folders that are in INDIR 
d = dir(INDIR); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

% Get timepoints
FFR_times_file = fullfile(INDIR,'ABR_timepoints.txt') ; 
fileID = fopen(FFR_times_file,'r');
formatSpec = '%f';
timepoints = fscanf(fileID,formatSpec);

% Load .txt files containing FFR for each subject put them in a matrix
% datapoints x subjects
mat = 1;
all_subj = zeros(size(timepoints,1),1);
for loopnum = 1:length(subjects) %for each subject
%for loopnum=find(ismember(subjects,'DVL_004_T10')) ; 
    FFR_file = fullfile(INDIR,subjects{loopnum},strcat(subjects{loopnum},'_abr_shifted_data_HF.txt')) ; 
    fileID = fopen(FFR_file,'r');
    formatSpec = '%f';
    FFR_subj = fscanf(fileID,formatSpec);
    all_subj(:,mat) = FFR_subj(:,1) ;
    mat = mat+1;
end

% % Plot individual FFRs
% for jj = 1:size(subjects,1)
%     figure ; 
%     plot(timepoints,all_subj(:,jj),'b','Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
%     grid on ; 
%     legend('Individual FFR', subjects(jj), 'Interpreter', 'None');
%     xlabel('Times (ms)'); ylabel('uV'); title ([' FFR ', subjects(jj)], 'Interpreter', 'None');
% end

% % Exclude noisy participants from observation
% a = ismember(subjects,noise_list)==1
% for excl = 1:size(subjects,1)
%     if a(excl)==1
%         subjects(excl) = [];
%         all_subj(:,excl) = []; 
%     end
% end

% Exclude noisy participants from max and min values
excluded_subj = {};
ex = 1;
for noise = 1:size(all_subj,2)
    if max(all_subj(:,noise)) > 0.5 || min(all_subj(:,noise)) < -0.5
        excluded_subj(ex) = subjects(noise,1);   
        disp([subjects(noise,1), 'excluded']);
        ex = ex +1;
    end
   
end

% Delete excluded subjects
exclud = ismember(subjects,excluded_subj);
new_all_subj = zeros(size(timepoints,1),1);
new_subjects = {};
n = 1;
for tt = 1:size(exclud,1)
    if exclud(tt)==0
        new_all_subj(:,n) = all_subj(:,tt);
        new_subjects(n) = subjects(tt);
        n = n+1;
    end
end

  
%%
% Compute grand average and plot
grd_FFR = mean(new_all_subj,2);
figure ; 
plot(timepoints,grd_FFR,'r','Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
grid on ; 
legend('Grand average FFR');
xlabel('Times (ms)'); ylabel('uV'); title ('Grand average FFR, 6-24 mo');

% Classify into groups
grpA.suffix = {'_T3','_T6','_T8','_T10'};
grpB.suffix = {'_T18','_T24'};

% Compute grand average by group
grpA.subj = new_subjects(contains(new_subjects,grpA.suffix));
grpB.subj = new_subjects(contains(new_subjects,grpB.suffix));

grpA.data = zeros(size(timepoints,1),1);
grpB.data = zeros(size(timepoints,1),1);

groups = {grpA.subj, grpB.subj};
data_groups = {grpA.data, grpB.data};
for k = 1:length(groups)
    grp_subj = zeros(size(timepoints,1),1);
    indices = find(ismember(new_subjects,groups{k})==1);
    for l = 1:length(indices)
        grp_subj(:,l) = new_all_subj(:,indices(l));
    end
    data_groups{k} = grp_subj;
end

FFR_avg_grpA = mean(data_groups{1,1},2);
FFR_avg_grpB = mean(data_groups{1,2},2);

% Visualization by group

% Temporal
figure ; 
plot(timepoints,FFR_avg_grpA,'b','Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
plot(timepoints,FFR_avg_grpB,'g','Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
grid on ; 
legend('Grand average FFR 6-10 mo', 'Grand average FFR 18-24 mo');
xlabel('Times (ms)'); ylabel('uV'); title ('Grand average FFR group comparison');

%%
% Frequential
%TODO


Fs = 16384;            % Sampling frequency                    
T = 1/Fs;              % Sampling period       
L = 250;               % Length of signal
t = (0:L-1)*T;         % Time vector

Y = fft(grd_FFR);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
plot(f,P1) 
title("Single-Sided Amplitude Spectrum of X(t)")
xlabel("f (Hz)")
ylabel("|P1(f)|")

