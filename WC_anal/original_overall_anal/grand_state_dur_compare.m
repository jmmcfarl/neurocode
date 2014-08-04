clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\UDS_synch_state_dur\UDS_synch_state_dur_data

up_range = [0.3 50];
down_range = [0.3 50];
numBins = 100;

for d = 1:length(over_dir)
    
    [up_hist(d,:),up_grid] = log_hist(synch_up_dur{d},up_range,numBins);
    [up_hist8(d,:),up_grid] = log_hist(synch_up_dur8{d},up_range,numBins);
    up_hist(d,:) = jmm_smooth_1d(up_hist(d,:),2);
    up_hist8(d,:) = jmm_smooth_1d(up_hist8(d,:),2);

%     [up_hist(d,:),up_grid] = log_hist(synch_up_dur{d},up_range,numBins);
%     [up_hist8(d,:),up_grid] = log_hist(synch_up_dur8{d},up_range,numBins);
%     up_hist(d,:) = jmm_smooth_1d(up_hist(d,:),2);
%     up_hist8(d,:) = jmm_smooth_1d(up_hist8(d,:),2);
    
end

cell_type(cell_type == 8) = [];
cell_type(cell_type == 7) = [];
cell_type(cell_type == 6) = 4;
cell_type(cell_type == 5) = 4;

%% cycle through all cell types and compute the average trig averages
diff_cell_types = unique(cell_type);
cmap = colormap(jet(9));
for c = 1:length(diff_cell_types)
 
    get_cells = find(cell_type == diff_cell_types(c));  
    
    mp_up_dur(c,:) = nanmean(up_hist(get_cells,:));
    uci(c,:) = mp_up_dur(c,:)+std(up_hist(get_cells,:))/sqrt(length(get_cells));
    lci(c,:) = mp_up_dur(c,:)-std(up_hist(get_cells,:))/sqrt(length(get_cells));
        
    plot(up_grid,cumsum(mp_up_dur(c,:)),'color',cmap(c,:),'linewidth',2)
        hold on
 
end


load('C:\WC_Germany\overall_calcs\HC_Cort\hc_cor_dir.mat')
load C:\WC_Germany\overall_calcs\HC_Cort\UDS_synch_state_dur\UDS_synch_state_dur_data

up_range = [0.3 20];
up_range = [0.3 20];
numBins = 100;

for d = 1:length(over_dir)
    
    [up_hist(d,:),up_grid] = log_hist(synch_up_dur{d},up_range,numBins);
    [up_hist8(d,:),up_grid] = log_hist(synch_up_dur8{d},up_range,numBins);
    up_hist(d,:) = jmm_smooth_1d(up_hist(d,:),2);
    up_hist8(d,:) = jmm_smooth_1d(up_hist8(d,:),2);

%     [up_hist(d,:),up_grid] = log_hist(synch_up_dur{d},up_range,numBins);
%     [up_hist8(d,:),up_grid] = log_hist(synch_up_dur8{d},up_range,numBins);
%     up_hist(d,:) = jmm_smooth_1d(up_hist(d,:),2);
%     up_hist8(d,:) = jmm_smooth_1d(up_hist8(d,:),2);
    
end

diff_cell_types = unique(cell_type);

for c = 5:9
 
    get_cells = find(cell_type == diff_cell_types(c-4));  
    
    mp_up_dur(c,:) = nanmean(up_hist(get_cells,:));
    uci(c,:) = mp_up_dur(c,:)+std(up_hist(get_cells,:))/sqrt(length(get_cells));
    lci(c,:) = mp_up_dur(c,:)-std(up_hist(get_cells,:))/sqrt(length(get_cells));
    
    if c ~= 6 & c ~= 8    
    plot(up_grid,cumsum(mp_up_dur(c,:)),'color',cmap(c,:),'linewidth',2)
        hold on
    end
    
end

legend('MEC Pyr','MEC Ste','MEC 5','LEC','Cort','CA1 Int','DG')

for c = 5:9
    if c ~= 6 & c ~= 8
        
        plot(up_grid,uci(c,:),'--','color',cmap(c,:))
        plot(up_grid,lci(c,:),'--','color',cmap(c,:))
        
    end
end


load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\UDS_synch_state_dur\UDS_synch_state_dur_data

up_range = [0.3 20];
up_range = [0.3 20];
numBins = 100;

for d = 1:length(over_dir)
    
    [up_hist(d,:),up_grid] = log_hist(synch_up_dur{d},up_range,numBins);
    [up_hist8(d,:),up_grid] = log_hist(synch_up_dur8{d},up_range,numBins);
    up_hist(d,:) = jmm_smooth_1d(up_hist(d,:),2);
    up_hist8(d,:) = jmm_smooth_1d(up_hist8(d,:),2);

%     [up_hist(d,:),up_grid] = log_hist(synch_up_dur{d},up_range,numBins);
%     [up_hist8(d,:),up_grid] = log_hist(synch_up_dur8{d},up_range,numBins);
%     up_hist(d,:) = jmm_smooth_1d(up_hist(d,:),2);
%     up_hist8(d,:) = jmm_smooth_1d(up_hist8(d,:),2);
    
end

cell_type(cell_type == 8) = [];
cell_type(cell_type == 7) = [];
cell_type(cell_type == 6) = 4;
cell_type(cell_type == 5) = 4;

diff_cell_types = unique(cell_type);
for c = 1:length(diff_cell_types)
 
    get_cells = find(cell_type == diff_cell_types(c));  
    
    mp_up_dur(c,:) = nanmean(up_hist(get_cells,:));
    uci(c,:) = mp_up_dur(c,:)+std(up_hist(get_cells,:))/sqrt(length(get_cells));
    lci(c,:) = mp_up_dur(c,:)-std(up_hist(get_cells,:))/sqrt(length(get_cells));
        
        plot(up_grid,cumsum(uci(c,:)),'--','color',cmap(c,:))
        plot(up_grid,cumsum(lci(c,:)),'--','color',cmap(c,:))

end
% 

xlim([0 10]);grid on
% set(gca,'xscale','log')
% set(gca,'yscale','log')
