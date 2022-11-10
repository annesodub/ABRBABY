function visualisation = visual_comparison_several_electrodes_3values(condition,elec_to_disp_labels,elec_indices,STD,DEV,MMN,timepoints,fig_name,legend1,legend2,legend3)
   
visualisation = figure('Name', fig_name,'Units','normalized','Position',[0,0,1,1]);
 
span = 1000;

STD_color = [0.4941 0.1019 0.8863]; %purple
DEV1_color = [1 0.7686 0]; %light orange
DEV2_color = [1 0.4 0]; %dark orange
DEV_colors = {DEV1_color, DEV2_color};
DIFF_color = [0 0 0]; %black

nfig = 1 ; 

for elec_letter=1:size(elec_to_disp_labels,1)
    
    for elec_numb=1:size(elec_to_disp_labels,2)
        hAxes = subplot(size(elec_to_disp_labels,1),size(elec_to_disp_labels,2),nfig) ; 
        nfig = nfig +1 ;
        idx_elec = elec_indices(elec_letter,elec_numb);
        
        % Compute grand average over one electrode
        grd_std = squeeze(mean(STD(idx_elec,:),1)) ; 
        grd_dev = squeeze(mean(DEV(idx_elec,:),1)) ; 
        grd_mmn = squeeze(mean(MMN(idx_elec,:),1)) ; 
        
        %Add 10Hz low pass filter for display
        grd_std = smooth(timepoints, grd_std,span); %'span' points moving average (ex from pop_smootherp in ERPLAB)
        grd_dev = smooth(timepoints, grd_dev,span);
        grd_mmn = smooth(timepoints, grd_mmn,span); 
        %grd_std = lowpass(grd_std',10,16384)';
        %grd_dev = lowpass(grd_dev',10,16384)';
        %grd_mmn = lowpass(grd_mmn',10,16384)';
        %grd_std_filt = pop_eegfilt(grd_std, 0, 10, 8);
        %grd_dev_filt = pop_eegfilt(grd_dev, 0, 10, 8);
        %grd_mmn_filt = pop_eegfilt(grd_mmn, 0, 10, 8); 
        
        % Plot timeseries
        plot(timepoints,grd_std_filt,'Color',STD_color,'Linewidth',1.5); hold on ;set(gca,'YDir','reverse') ; 
        plot(timepoints,grd_dev_filt,'Color', DEV_colors{condition},'Linewidth',1.5);  hold on; set(gca,'YDir','reverse') ;
        plot(timepoints,grd_mmn_filt,'Color', DIFF_color,'Linewidth',1.5);  hold on; set(gca,'YDir','reverse') ; hold on ; plot(timepoints,grd_std,'b','Linewidth',1.5)
        
        % % Plot transparetn halo (+-mad)
        %plotHaloPatchMAD(hAxes, timepoints, squeeze(mean((raw_avg1_grpA(:,idx_elec,:)+raw_avg2_grpA(:,idx_elec,:))/2,1)), [0,255,0]) ; 
        %plotHaloPatchMAD(hAxes, timepoints, squeeze(mean((raw_avg1_grpB(:,idx_elec,:)+raw_avg2_grpB(:,idx_elec,:))/2,1)), [255,0,0]) ; 
                 
        % Adjust graphics
         %xlim([grd_STD_grpA.xmin, grd_STD_grpB.xmax]*1000); grid on ; 
         grid on;
         ylim([-4 5]);
         legend(legend1,legend2,legend3);
         title(elec_to_disp_labels(elec_letter,elec_numb));
         xlabel('Times (ms)'); ylabel('uV');
         set(hAxes,'Fontsize',9);
                
         % To save data in vectoriel
         %print('-dsvg',fullfile(INDIR,strcat(subjects{jj},'_',conditions{cc+1},'.svg'))) ; 
    end
    
end

end