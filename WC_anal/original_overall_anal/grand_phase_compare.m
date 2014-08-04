%% grad phase compare

clear all
close all
 

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\phase\phase_data

cell_type(cell_type == 8) = [];
cell_type(cell_type == 7) = [];
cell_type(cell_type == 6) = 4;
cell_type(cell_type == 5) = 4;

lec_cells = find(cell_type==4);
up_dist(lec_cells([2 3 5 6 7 8 9]),:) = [];
cell_type(lec_cells([2 3 5 6 7 8 9])) = [];


phase_range = [-fliplr(phase_range) phase_range];

%% cycle through all cell types and compute the average trig averages
diff_cell_types = unique(cell_type);
cmap = colormap(jet(9));
for c = 1:length(diff_cell_types)
 
    get_cells = find(cell_type == diff_cell_types(c));  
    [ov_mean(c),ul_mean(c),ll_mean(c)] = circ_mean(circ_up_mean(get_cells));
    mp_up_phase(c,:) = nanmean(up_dist(get_cells,:));
    uci(c,:) = mp_up_phase(c,:)+std(up_dist(get_cells,:))/sqrt(length(get_cells));
    lci(c,:) = mp_up_phase(c,:)-std(up_dist(get_cells,:))/sqrt(length(get_cells));
    
    temp = [mp_up_phase(c,:) mp_up_phase(c,:)];
    
    plot(phase_range,temp,'color',cmap(c,:),'linewidth',2)
        hold on
 
end


% load('C:\WC_Germany\overall_calcs\HC_Cort\hc_cor_dir.mat')
% load C:\WC_Germany\overall_calcs\HC_Cort\phase\phase_data
% phase_range = [-fliplr(phase_range) phase_range];
% 
% diff_cell_types = unique(cell_type);
% 
% for c = 5:9
%  
%     get_cells = find(cell_type == diff_cell_types(c-4));  
%     
%     mp_up_phase(c,:) = nanmean(up_dist(get_cells,:));
%     uci(c,:) = mp_up_phase(c,:)+std(up_dist(get_cells,:))/sqrt(length(get_cells));
%     lci(c,:) = mp_up_phase(c,:)-std(up_dist(get_cells,:))/sqrt(length(get_cells));
%     
%     if c ~= 6 & c ~= 8
%     temp = [mp_up_phase(c,:) mp_up_phase(c,:)];
%     
%     plot(phase_range,temp,'color',cmap(c,:),'linewidth',2)
%         hold on
%     end
%     
% end
% 
% legend('MEC Pyr','MEC Ste','MEC 5','LEC','Cort','CA1 Int','DG')
% 
% for c = 5:9
%     if c ~= 6 & c ~= 8
%         
%         temp = [uci(c,:) uci(c,:)];
%         plot(phase_range,temp,'--','color',cmap(c,:))
%         temp = [lci(c,:) lci(c,:)];
%         plot(phase_range,temp,'--','color',cmap(c,:))
%         
%     end
% end


load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\phase\phase_data
phase_range = [-fliplr(phase_range) phase_range];

cell_type(cell_type == 8) = [];
cell_type(cell_type == 7) = [];
cell_type(cell_type == 6) = 4;
cell_type(cell_type == 5) = 4;


lec_cells = find(cell_type==4);
up_dist(lec_cells([2 3 5 6 7 8 9]),:) = [];
cell_type(lec_cells([2 3 5 6 7 8 9])) = [];

diff_cell_types = unique(cell_type);
for c = 1:length(diff_cell_types)
 
    get_cells = find(cell_type == diff_cell_types(c));  
    
    mp_up_phase(c,:) = nanmean(up_dist(get_cells,:));
    uci(c,:) = mp_up_phase(c,:)+std(up_dist(get_cells,:))/sqrt(length(get_cells));
    lci(c,:) = mp_up_phase(c,:)-std(up_dist(get_cells,:))/sqrt(length(get_cells));
        
       temp = [uci(c,:) uci(c,:)];
        plot(phase_range,temp,'--','color',cmap(c,:))
        temp = [lci(c,:) lci(c,:)];
        plot(phase_range,temp,'--','color',cmap(c,:))

end


xlim([-360 360]);grid on
line([180 180],[0 0.08],'Color','k','linewidth',2)
line([-180 -180],[0 0.08],'Color','k','linewidth',2)