clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\time_domain\time_domain_data

cell_type(cell_type == 8) = [];
cell_type(cell_type == 7) = [];
cell_type(cell_type == 6) = 4;
cell_type(cell_type == 5) = 4;

%% cycle through all cell types and compute the average trig averages
diff_cell_types = unique(cell_type);
cmap = colormap(jet(11));
for c = 1:length(diff_cell_types)
 
    get_cells = find(cell_type == diff_cell_types(c));  
    
    mp_acorr(c,:) = nanmean(tot_wcv_acorr(get_cells,:));
    uci(c,:) = mp_acorr(c,:)+std(tot_wcv_acorr(get_cells,:))/sqrt(length(get_cells));
    lci(c,:) = mp_acorr(c,:)-std(tot_wcv_acorr(get_cells,:))/sqrt(length(get_cells));
        
    plot(lags/Fsd,mp_acorr(c,:),'color',cmap(c,:),'linewidth',2)
        hold on
 
end


load('C:\WC_Germany\overall_calcs\HC_Cort\hc_cor_dir.mat')
load C:\WC_Germany\overall_calcs\HC_Cort\time_domain\time_domain_data

diff_cell_types = unique(cell_type);

for c = 5:9
 
    get_cells = find(cell_type == diff_cell_types(c-4));  
    
    mp_acorr(c,:) = nanmean(tot_wcv_acorr(get_cells,:));
    uci(c,:) = mp_acorr(c,:)+std(tot_wcv_acorr(get_cells,:))/sqrt(length(get_cells));
    lci(c,:) = mp_acorr(c,:)-std(tot_wcv_acorr(get_cells,:))/sqrt(length(get_cells));
    
%     if c ~= 6 & c ~= 8    
    plot(lags/Fsd,mp_acorr(c,:),'color',cmap(c,:),'linewidth',2)
        hold on
%     end
    
end

legend('MEC Pyr','MEC Ste','MEC 5','LEC','Cort','Ca1 pyr','CA1 Int','CA3 Pyr','DG')
% 
for c = 5:9
%     if c ~= 6 & c ~= 8
        
        plot(lags/Fsd,uci(c,:),'--','color',cmap(c,:))
        plot(lags/Fsd,lci(c,:),'--','color',cmap(c,:))
        
%     end
end


load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\time_domain\time_domain_data

cell_type(cell_type == 8) = [];
cell_type(cell_type == 7) = [];
cell_type(cell_type == 6) = 4;
cell_type(cell_type == 5) = 4;

diff_cell_types = unique(cell_type);
for c = 1:length(diff_cell_types)
 
    get_cells = find(cell_type == diff_cell_types(c));  
    
    mp_acorr(c,:) = nanmean(tot_wcv_acorr(get_cells,:));
    uci(c,:) = mp_acorr(c,:)+std(tot_wcv_acorr(get_cells,:))/sqrt(length(get_cells));
    lci(c,:) = mp_acorr(c,:)-std(tot_wcv_acorr(get_cells,:))/sqrt(length(get_cells));
        
        plot(lags/Fsd,uci(c,:),'--','color',cmap(c,:))
        plot(lags/Fsd,lci(c,:),'--','color',cmap(c,:))

end
% 

xlim([0 10]);grid on
ylim([-0.4 0.4])
