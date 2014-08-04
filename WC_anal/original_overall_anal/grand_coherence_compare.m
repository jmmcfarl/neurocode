clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\coherence\coherence_data_hipp

%% smooth coherence
% for i = 1:length(over_dir)
%     Cmn3(i,:) = jmm_smooth_1d(Cmn3(i,:),10);
% end

cell_type(cell_type == 8) = [];
cell_type(cell_type == 7) = [];
cell_type(cell_type == 6) = 4;
cell_type(cell_type == 5) = 4;

%% cycle through all cell types and compute the average trig averages
diff_cell_types = unique(cell_type);
cmap = colormap(jet(11));
for c = 1:length(diff_cell_types)
 
    get_cells = find(cell_type == diff_cell_types(c));  
    
    avg_coh(c,:) = nanmean(Cmn3s(get_cells,:));
    uci(c,:) = avg_coh(c,:)+std(Cmn3s(get_cells,:))/sqrt(length(get_cells));
    lci(c,:) = avg_coh(c,:)-std(Cmn3s(get_cells,:))/sqrt(length(get_cells));
        
    plot(f,avg_coh(c,:),'color',cmap(c,:),'linewidth',2)
        hold on
 
end


load('C:\WC_Germany\overall_calcs\HC_Cort\hc_cor_dir.mat')
load C:\WC_Germany\overall_calcs\HC_Cort\coherence\coherence_data_hipp
% for i = 1:length(over_dir)
%     Cmn3s(i,:) = jmm_smooth_1d(Cmn3s(i,:),10);
% end

diff_cell_types = unique(cell_type);

for c = 5:9
 
    get_cells = find(cell_type == diff_cell_types(c-4));  
    
    avg_coh(c,:) = nanmean(Cmn3s(get_cells,:));
    uci(c,:) = avg_coh(c,:)+std(Cmn3s(get_cells,:))/sqrt(length(get_cells));
    lci(c,:) = avg_coh(c,:)-std(Cmn3s(get_cells,:))/sqrt(length(get_cells));
    
    plot(f,avg_coh(c,:),'color',cmap(c,:),'linewidth',2)
        hold on
    
end

legend('MEC Pyr','MEC Ste','MEC 5','LEC','Cort','CA1 pyr','CA1 Int','CA3 pyr','DG')

% for c = 5:9
%         
%         plot(f,uci(c,:),'--','color',cmap(c,:))
%         plot(f,lci(c,:),'--','color',cmap(c,:))
%         
% end
% 
% 
% load('C:\WC_Germany\overall_calcs\overall_dir.mat')
% load C:\WC_Germany\overall_calcs\coherence\coherence_data
% for i = 1:length(over_dir)
%     Cmn3s(i,:) = jmm_smooth_1d(Cmn3s(i,:),5);
% end

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
%     avg_coh(c,:) = nanmean(Cmn3s(get_cells,:));
%     uci(c,:) = avg_coh(c,:)+std(Cmn3s(get_cells,:))/sqrt(length(get_cells));
%     lci(c,:) = avg_coh(c,:)-std(Cmn3s(get_cells,:))/sqrt(length(get_cells));
%         
%         plot(f,uci(c,:),'--','color',cmap(c,:))
%         plot(f,lci(c,:),'--','color',cmap(c,:))
% 
% end
% 

xlim([0 10]);
