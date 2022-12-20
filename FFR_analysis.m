function [] = FFR_analysis(indir) 
% ERPs analysis script - 
% Estelle Herve, A.-Sophie Dubarry - 2022 - %80PRIME Project

% Reads all folders that are in INDIR 
d = dir(indir); 
isub = [d(:).isdir]; % returns logical vector if is folder
subjects = {d(isub).name}';
subjects(ismember(subjects,{'.','..'})) = []; % Removes . and ..

%Color for plot
DIFF_color = [0 0 0]; %black
FFR_color = [0.8902 0 0]; %red
grpA_color = [0.2 0.2 1]; %blue
grpB_color = [0.2 0.7765 0.2]; %green

% Get timepoints
FFR_times_file = fullfile(indir,'ABR_timepoints.txt') ; 
fileID = fopen(FFR_times_file,'r');
formatSpec = '%f';
timepoints = fscanf(fileID,formatSpec);

% Load .txt files containing FFR for each subject put them in a matrix
% datapoints x subjects
mat = 1;
all_subj = zeros(size(timepoints,1),1);
for loopnum = 1:length(subjects) %for each subject
%for loopnum=find(ismember(subjects,'DVL_003_T10')) ; 
    FFR_file = fullfile(indir,subjects{loopnum},strcat(subjects{loopnum},'_abr_shifted_data_HF.txt')) ; 
    if exist(FFR_file,'file')==0 ; error(['File does not exist, please run FFR sanity check on ',subjects{loopnum}]); end
    fileID = fopen(FFR_file,'r');
    formatSpec = '%f';
    FFR_subj = fscanf(fileID,formatSpec);
    all_subj(:,mat) = FFR_subj(:,1) ;
    mat = mat+1;
end

% % Exclude noisy participants from max and min values
% excluded_subj = {};
% ex = 1;
% for noise = 1:size(all_subj,2)
%     if max(all_subj(:,noise)) > 0.5 || min(all_subj(:,noise)) < -0.5
%         excluded_subj(ex) = subjects(noise,1);   
%         disp([subjects(noise,1), 'excluded']);
%         ex = ex +1;
%     end
%    
% end
% 
% % Delete excluded subjects
% exclud = ismember(subjects,excluded_subj);
% new_all_subj = zeros(size(timepoints,1),1);
% new_subjects = {};
% n = 1;
% for tt = 1:size(exclud,1)
%     if exclud(tt)==0
%         new_all_subj(:,n) = all_subj(:,tt);
%         new_subjects(n) = subjects(tt);
%         n = n+1;
%     end
% end

  
%% Time domain : FFR visualization
% Compute grand average and plot
grd_FFR = mean(all_subj,2);
%grd_FFR = mean(new_all_subj,2);
figure ; 
plot(timepoints,grd_FFR,'Color',FFR_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
grid on ; 
legend('Grand average FFR');
xlabel('Times (ms)'); ylabel('uV'); title ('Grand average FFR, 6-24 mo');

% Classify into groups
grpA.suffix = {'_T3','_T6','_T8','_T10'};
grpB.suffix = {'_T18','_T24'};

% Compute grand average by group
%grpA.subj = new_subjects(contains(new_subjects,grpA.suffix));
%sgrpB.subj = new_subjects(contains(new_subjects,grpB.suffix));
grpA.subj = subjects(contains(subjects,grpA.suffix));
grpB.subj = subjects(contains(subjects,grpB.suffix));

grpA.data = zeros(size(timepoints,1),1);
grpB.data = zeros(size(timepoints,1),1);

groups = {grpA.subj, grpB.subj};
data_groups = {grpA.data, grpB.data};
for k = 1:length(groups)
    grp_subj = zeros(size(timepoints,1),1);
    %indices = find(ismember(new_subjects,groups{k})==1);
    indices = find(ismember(subjects,groups{k})==1);
    for l = 1:length(indices)
        %grp_subj(:,l) = new_all_subj(:,indices(l));
        grp_subj(:,l) = all_subj(:,indices(l));
    end
    data_groups{k} = grp_subj;
end

FFR_avg_grpA = mean(data_groups{1,1},2);
FFR_avg_grpB = mean(data_groups{1,2},2);

% Visualization by group

% Temporal
figure ; 
plot(timepoints,FFR_avg_grpA,'Color', grpA_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
plot(timepoints,FFR_avg_grpB,'Color', grpB_color,'Linewidth',0.5); hold on ;set(gca,'YDir','reverse') ;
grid on ; 
legend('Grand average FFR 6-10 mo', 'Grand average FFR 18-24 mo');
xlabel('Times (ms)'); ylabel('uV'); title ('Grand average FFR group comparison');


%% Root mean square

% Adapted from bt_rms.m from bt_gui toolbox (Kraus & Skoe, 2010)
% Note: std(X,1) is a shortcut RMS calc. It is equivalent to RMS
% on a baselined (demeaned to 0) waveform.  This is what we want here.
startPeriod = 0.0610; % (in ms) to obtain 1 in samples
endPeriod = 40; % (in ms)

% response period
RMS_grpA = std(FFR_avg_grpA,1); 
RMS_grpB = std(FFR_avg_grpB,1); 

% Prestim period
Segment_A = FFR_avg_grpA(round(startPeriod*16384/1000):round(endPeriod*16384/1000));
Segment_B = FFR_avg_grpB(round(startPeriod*16384/1000):round(endPeriod*16384/1000));
RMSprestimA = std(Segment_A,1);
RMSprestimB = std(Segment_B,1);

% signal-to-noise
SNR_grpA = RMS_grpA./RMSprestimA;
SNR_grpB = RMS_grpB./RMSprestimB;

% display results
disp(["Group A: RMS total = ",RMS_grpA, "RMS prestim = ", RMSprestimA, "SNR = ", SNR_grpA]);
disp(["Group B: RMS total = ",RMS_grpB, "RMS prestim = ", RMSprestimB, "SNR = ", SNR_grpB]);


%% Stimulus-to-response correlation

%function [time autocorr LAG FFT freqAxis preFFT blocks]= pitchtrack(avgname, block, step, startAnalysis,channel,exportData)

block = 1;
step = 1;
startAnalysis = 1;
%channel = 1;
%exportData = 1;
file.xmin =1; 
file.xmax=240;
file.pnts=2932;

% if block<40
%     display('Block Size is too small. Using default: 40 ms');
%     block = 40;
% end

% %Open average file 
% [file]= openavg(avgname);
%Define time axis 
timeaxis = linspace(file.xmin, file.xmax, file.pnts)';

 % extract signal
 %SIGNAL = file.signal(:,channel);
 SIGNAL = FFR_avg_grpB;
 % get sampling rate
 fs = 16384;

 
 % ---------------------------PRESTIM, ---------------------------
 % extract portion of prestimulus time, and detrend. The size of the prestim portion is dependent on the block size.
 % To do SNR calculations block and prestim must be same the same number of ms.
 % If the block size is 40ms, then only 40ms of the prestim will be
 % extracted. 
 
 %PRESTIM = SIGNAL(1: ms2row(file, block));   %start with the very first point.
 PRESTIM = Segment_B;
 % ramp
 r = hann(size(PRESTIM, 1));  % the entire prestim is ramped.
 % ramp and detrend
 PRESTIM = detrend(PRESTIM.*hann(size(PRESTIM,1)), 'constant');
 % FFT (zero-padded to sampling rate);
 preFFT = abs(fft(PRESTIM, fs));
 %scale preFFT
 preFFT= preFFT*(2./length(PRESTIM));
 preFFT=preFFT(1:1001,1);  %truncate above 1000 Hz
 
   
%  ------------------ FFTS of RESPONSE CHUNKS-----------------------
j = startAnalysis; % each time through loop j increases by step size;

    chunks = 5000;    % an arbitrary maximum number of blocks that the program will create. 
                        
    for k = 1:chunks;   %the program knows to stop once file.xmax is exceeded 

        % variables created: 
        ramptime = (block/1000);   % ramp the entire chunk
        start = j;
        stop = j+block;

        if stop>(file.xmax)  % if stop exceeds the maximum ms time then abort and break out from loop
            k=k-1;
            j=j-step;
            break;
        else
            signal = detrend(SIGNAL(ms2row(file, start):ms2row(file, stop)), 'constant');   % de-mean to zero           
        end
        
        midpoint(k) = mean(ms2row(file, start):ms2row(file,stop));  %calculates the time corresponding to the midpoint of the chunk

        % generate ramp
        ramp = hann(size(signal,1));
        % ramp and de-mean
        signal = detrend(signal.*ramp, 'constant');


        % autocorrelation (see Boersma 1993)
        [c lag]=xcorr(signal, 'coeff');
       
        
        [cwin lagwin]=xcorr(ramp, 'coeff');
        LAG = linspace(-block, block, length(c));
    
        autoc=c./cwin;
        autoc(autoc>1)=1;  %this handles the very rare case that the remainder of the previous step is >1.
        % only plot the first 15 ms; % lowest frequency is ~66 Hz.
        
        
            
            startlag = find(LAG==0);  
            endlag = find(LAG==closestrc(LAG, 15));

        
        % truncate lag and r value matrices to only include values up to first 15 ms.
        autocorr(:,k)=autoc(startlag:endlag)';
        LAG = LAG(startlag:endlag)';
     
        ostartlag = startlag;
        
        ostoplag = endlag;
         
         %%% Now do FFTs;
        % fft, pads to sampling rate
        fftsignal{k} = abs(fft(signal, fs));
        % only go up to 1000 Hz;
        FFT{k} = fftsignal{k}(1:1001,1);
        FFT{k}= FFT{k}*(2./length(signal));
        freqAxis = linspace(0, 1000, 1001);
        
        j = j+step;  % loop through next time chunk
        
       
       
    end
    time = timeaxis(round(midpoint));
    
    blocks = k;
    
    FFT = cell2mat(FFT);
    
    
%     if exportData == 1
%         [fpath fname ext]=fileparts(avgname);
%         FFTfile = [fpath, '\', fname, '-FFTmatrix.xls'];
%         ACfile = [fpath, '\', fname, '-ACmatrix.xls'];
%         
%         xlswrite(FFT, fname, {'FFT matrix'}, FFTfile, 'Sheet1');
%         xlswrite(freqAxis', fname, {'Frequency Axis'}, FFTfile, 'Sheet2');
%         xlswrite(time, fname, {'Time Axis'}, FFTfile,     'Sheet3');
%         
%         xlswrite(autocorr, fname, {'autocorrelation matrix'}, ACfile, 'Sheet1');
%         xlswrite( LAG, fname, {'Lag Axis'} , ACfile, 'Sheet2');
%         xlswrite(time, fname, {'Time Axis'} , ACfile, 'Sheet3');
%        
%      
%     end
%     
%% Response-to-response correlation    

%compare grpA and grpB

%% Response consistency 

%compare 2 subaverages of the same FFR (use all subjects from eacg group)

%% Frequency domain

% Adaptation of Skoe function (bt_fftsc)
 %function [Freq1 Freq2 Freq3 fftFFR HzScale]=bt_fftsc(FILENAME,start,stop,F0_Lo,F0_Hi,F1_Lo,F1_Hi,HF_Lo,HF_Hi, chan)
% bt_fftsc computes frequency-domain amplitudes of three user-defined 
% frequency bins of Brainstem response.  Results are scaled to peak �V.
%
% Usage: [F0 F1 HF] = bt_fftsc('filename.avg',10,40,100,150,300,350,600,800);
%    over the range of 10 to 40 ms, finds average frequency amplitude of
%    100-150 Hz, 300-350 Hz and 600-800 Hz.
% 
% Three variables are returned the workspace:
% (1) Freq1: mean amplitude over F0_Lo-F0_Hi Hz
% (2) Freq2: mean amplitude over F1_Lo-F1_Hi Hz
% (3) Freq3: mean amplitude over HF_Lo to HF_Hi Hz 

% % Dependancies: ms2row, openavg

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Originally developed by E.E. Skoe.  
% Toolbox version by E.E. Skoe & T.G. Nicol
% eeskoe@northwestern.edu tgn@northwestern.edu
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

F0_Lo = 100;
F0_Hi = 110;
F1_Lo = 455;
F1_Hi = 740;
HF_Lo = 1060;
HF_Hi = 2750;

%%%  OPEN UP AVG FILE
%avg = openavg(FILENAME);
%data = grd_FFR;
FS = 16384;
FFR = grd_FFR;
%******** STEP 1. CREATE VARIABLE "FFR" CORRESPONDING TO FFR PERIOD 

%startPoint = ms2row(avg, start);
%endPoint = ms2row(avg, stop);

%FFR = data(startPoint:endPoint);
numPoints = length(FFR);
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                      
%**** STEP 2. FFT OF FFR
%******** STEP 2a. CREATE and APPLY HANNING RAMP 2 msec rise, 2 msec fall
rampMS = 4/1000; % length of ramp (on and off) in seconds
% hanPoints = 26;  %hard coded to be the same as Biologic's settings (December 7, 2005);
hanPoints = rampMS.*FS; % length of ramp in points
hanPoints = 2.*round(hanPoints/2); % force it to be nearest even.
hanHalfPoints = round(hanPoints./2);
numberOfOnes = numPoints - hanPoints;
FFRhan = hann(hanPoints);  
FFRhan = [FFRhan(1:hanHalfPoints); ones(numberOfOnes,1); FFRhan(hanHalfPoints+1:hanPoints)];

% baseline, window, then baseline again
FFR = detrend(FFR, 'constant');
FFR = detrend(FFR.*FFRhan, 'constant');

%******** STEP 2b. Perform FFT
fftFFR = abs(fft(FFR, round(FS)));
fftFFR = fftFFR(1:round(round(FS)/2));
fftFFR = fftFFR.*(2./numPoints); % scale to peak �V
HzScale = [0:1:round(FS/2)]'; % frequency 'axis'
HzScale = HzScale(1:length(fftFFR));

% Plot FFT
figure ;

subplot(2,2,1);
plot(HzScale,fftFFR);
grid on;
title("Single-Sided Amplitude Spectrum of X(t)");
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

subplot(2,2,2); 
plot(HzScale,fftFFR);
xlim([90 110]);
grid on;
title("Single-Sided Amplitude Spectrum of X(t)");
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

subplot(2,2,3); 
plot(HzScale,fftFFR);
xlim([300 500]);
grid on;
title("Single-Sided Amplitude Spectrum of X(t)");
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

subplot(2,2,4); 
plot(HzScale,fftFFR);
xlim([1100 1300]);
grid on;
title("Single-Sided Amplitude Spectrum of X(t)");
xlabel("Frequency (Hz)");
ylabel("Amplitude (µV)");

%% 

