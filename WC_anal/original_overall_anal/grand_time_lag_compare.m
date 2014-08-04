%% grad phase compare

clear all
close all
 

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\time_lag\time_lag_data

cell_type(cell_type == 8) = [];
cell_type(cell_type == 7) = [];
cell_type(cell_type == 6) = 4;
cell_type(cell_type == 5) = 4;

%% cycle through all cell types and compute the average trig averages
diff_cell_types = unique(cell_type);
cmap = colormap(jet(9));
for c = 1:length(diff_cell_types)
 
    get_cells = find(cell_type == diff_cell_types(c));  
    
    mp_up_phase(c,:) = nanmean(up_dist(get_cells,:));
    uci(c,:) = mp_up_phase(c,:)+std(up_dist(get_cells,:))/sqrt(length(get_cells));
    lci(c,:) = mp_up_phase(c,:)-std(up_dist(get_cells,:))/sqrt(length(get_cells));
        
    plot(lag_range,mp_up_phase(c,:),'color',cmap(c,:),'linewidth',2)
        hold on
 
end


load('C:\WC_Germany\overall_calcs\HC_Cort\hc_cor_dir.mat')
load C:\WC_Germany\overall_calcs\HC_Cort\time_lag\time_lag_data

diff_cell_types = unique(cell_type);

for c = 5:9
 
    get_cells = find(cell_type == diff_cell_types(c-4));  
    
    mp_up_phase(c,:) = nanmean(up_dist(get_cells,:));
    uci(c,:) = mp_up_phase(c,:)+std(up_dist(get_cells,:))/sqrt(length(get_cells));
    lci(c,:) = mp_up_phase(c,:)-std(up_dist(get_cells,:))/sqrt(length(get_cells));
    
    if c ~= 6 & c ~= 8    
    plot(lag_range,mp_up_phase(c,:),'color',cmap(c,:),'linewidth',2)
        hold on
    end
    
end

legend('MEC Pyr','MEC Ste','MEC 5','LEC','Cort','CA1 Int','DG')

% for c = 5:9
%     if c ~= 6 & c ~= 8
%         
%         plot(lag_range,uci(c,:),'--','color',cmap(c,:))
%         plot(lag_range,lci(c,:),'--','color',cmap(c,:))
%         
%     end
% end
% 
% 
% load('C:\WC_Germany\overall_calcs\overall_dir.mat')
% load C:\WC_Germany\overall_calcs\time_lag\time_lag_data
% 
% cell_type(cell_type == 8) = [];
% cell_type(cell_type == 7) = [];
% cell_type(cell_type == 6) = 4;
% cell_type(cell_type == 5) = 4;
% 
% diff_cell_types = unique(cell_type);
% for c = 1:length(diff_cell_types)
%  
%     get_cells = find(cell_type == diff_cell_types(c));  
%     
%     mp_up_phase(c,:) = nanmean(up_dist(get_cells,:));
%     uci(c,:) = mp_up_phase(c,:)+std(up_dist(get_cells,:))/sqrt(length(get_cells));
%     lci(c,:) = mp_up_phase(c,:)-std(up_dist(get_cells,:))/sqrt(length(get_cells));
%         
%         plot(lag_range,uci(c,:),'--','color',cmap(c,:))
%         plot(lag_range,lci(c,:),'--','color',cmap(c,:))
% 
% end
% 

xlim([-2 2]);grid on
