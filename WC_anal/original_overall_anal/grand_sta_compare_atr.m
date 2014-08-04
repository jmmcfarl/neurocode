clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\sta\sta_data5

% cell_type(cell_type == 8) = [];
% cell_type(cell_type == 7) = [];
% cell_type(cell_type == 6) = 4;
% cell_type(cell_type == 5) = 4;

%% cycle through all cell types and compute the average trig averages
diff_cell_types = unique(cell_type);
cmap = colormap(jet(4));
used_types = [1 2 7 8];
for c = 1:length(used_types)
 
    get_cells = find(cell_type == diff_cell_types(used_types(c)));  
    
    avg_st5(c,:) = nanmean(sta_lf5(get_cells,:));
    uci(c,:) = avg_st5(c,:)+std(sta_lf5(get_cells,:))/sqrt(length(get_cells));
    lci(c,:) = avg_st5(c,:)-std(sta_lf5(get_cells,:))/sqrt(length(get_cells));
        
    plot(lags/Fs,avg_st5(c,:),'color',cmap(c,:),'linewidth',2)
        hold on
 
end


legend('MEC Pyr','MEC Ste','Pyr Atr','Ste Atr')

for c = 1:length(used_types)
        
        plot(lags/Fs,uci(c,:),'--','color',cmap(c,:))
        plot(lags/Fs,lci(c,:),'--','color',cmap(c,:))
        
end

xlim([-0.2 0.2]);grid on

