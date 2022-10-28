%function visualisation = visual_comparison_several_electrodes_2values(mean_grpA, mean_grpB,raw_avg1_grpA,raw_avg2_grpA, raw_avg1_grpB,raw_avg2_grpB,timepoints,fig_name,legend1,legend2)
function visualisation = visual_comparison_several_electrodes_2values(elec_to_disp_labels,elec_indices,mean_grpA, mean_grpB,timepoints,fig_name,legend1,legend2)

visualisation = figure('Name', fig_name,'Units','normalized','Position',[0,0,1,1]);
  
grpA_color = [0.2 0.2 1]; %blue
grpB_color = [0.2 0.7765 0.2]; %green

nfig = 1 ; 

for elec_letter=1:size(elec_to_disp_labels,1)
    
    for elec_numb=1:size(elec_to_disp_labels,2)
        hAxes = subplot(size(elec_to_disp_labels,1),size(elec_to_disp_labels,2),nfig) ; 
        nfig = nfig +1 ;
        idx_elec = elec_indices(elec_letter,elec_numb);
        
        % Compute grand average over one electrode
        grd_avg_grpA = squeeze(mean(mean_grpA(idx_elec,:),1)) ; 
        grd_avg_grpB = squeeze(mean(mean_grpB(idx_elec,:),1)) ; 
                
        %Add 10 Hz lowpass filter
        grd_avg_grpA = lowpass(grd_avg_grpA, 10, 16384);
        grd_avg_grpB = lowpass(grd_avg_grpB, 10, 16384);
        
        % Plot timeseries
        plot(timepoints,grd_avg_grpA,'Color', grpA_color,'Linewidth',1.5); hold on ;set(gca,'YDir','reverse') ; 
        plot(timepoints,grd_avg_grpB,'Color', grpB_color,'Linewidth',1.5);  hold on; set(gca,'YDir','reverse') ;
        
        % Plot transparetn halo (+-mad)
        %plotHaloPatchMAD(hAxes, timepoints, squeeze(mean((raw_avg1_grpA(:,idx_elec,:)+raw_avg2_grpA(:,idx_elec,:))/2,1)), [0,255,0]) ; 
        %plotHaloPatchMAD(hAxes, timepoints, squeeze(mean((raw_avg1_grpB(:,idx_elec,:)+raw_avg2_grpB(:,idx_elec,:))/2,1)), [255,0,0]) ; 
                 
        % Adjust graphics
         %xlim([grd_STD_grpA.xmin, grd_STD_grpB.xmax]*1000); grid on ; 
         grid on;
         legend(legend1,legend2);
         title(elec_to_disp_labels(elec_letter,elec_numb));
         xlabel('Times (ms)'); ylabel('uV');
         set(hAxes,'Fontsize',12);
                
         % To save data in vectoriel
         %print('-dsvg',fullfile(INDIR,strcat(subjects{jj},'_',conditions{cc+1},'.svg'))) ; 
    end
    
end

end