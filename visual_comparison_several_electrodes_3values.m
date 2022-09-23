function visualisation = visual_comparison_several_electrodes_3values(elec_to_disp_labels,elec_indices,STD,DEV,MMN,timepoints,fig_name,legend1,legend2,legend3)
   
visualisation = figure('Name', fig_name,'Units','normalized','Position',[0,0,1,1]);
     
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
                
        % Plot timeseries
        plot(timepoints,grd_std,'g','Linewidth',1.5); hold on ;set(gca,'YDir','reverse') ; 
        plot(timepoints,grd_dev,'r','Linewidth',1.5);  hold on; set(gca,'YDir','reverse') ;
        plot(timepoints,grd_mmn,'k','Linewidth',1.5);  hold on; set(gca,'YDir','reverse') ;
        
        % % Plot transparetn halo (+-mad)
        %plotHaloPatchMAD(hAxes, timepoints, squeeze(mean((raw_avg1_grpA(:,idx_elec,:)+raw_avg2_grpA(:,idx_elec,:))/2,1)), [0,255,0]) ; 
        %plotHaloPatchMAD(hAxes, timepoints, squeeze(mean((raw_avg1_grpB(:,idx_elec,:)+raw_avg2_grpB(:,idx_elec,:))/2,1)), [255,0,0]) ; 
                 
        % Adjust graphics
         %xlim([grd_STD_grpA.xmin, grd_STD_grpB.xmax]*1000); grid on ; 
         legend(legend1,legend2,legend3);
         title(elec_to_disp_labels(elec_letter,elec_numb));
         xlabel('Times (ms)'); ylabel('uV');
         set(hAxes,'Fontsize',9);
                
         % To save data in vectoriel
         %print('-dsvg',fullfile(INDIR,strcat(subjects{jj},'_',conditions{cc+1},'.svg'))) ; 
    end
    
end

end