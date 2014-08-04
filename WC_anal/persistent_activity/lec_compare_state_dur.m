clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\UDS_synch_state_dur\UDS_synch_state_dur_data_new_method

maxdur = 50;
up_range = [0.3 maxdur];
down_range = [0.3 maxdur];
% numBins = 300;
% cur_grid = [0.3; maxdur; numBins];
% up_grid = linspace(0.3,maxdur,numBins);
% down_grid = linspace(0.3,maxdur,numBins);
numBins = 50;
up_grid = linspace(0.3,maxdur,numBins);
down_grid = linspace(0.3,maxdur,numBins);

for d = 1:length(over_dir)
    
    avg_up_dur(d) = nanmean(synch_up_dur{d});
    avg_up_dur8(d) = nanmean(synch_up_dur8{d});
    avg_down_dur(d) = nanmean(synch_down_dur{d});
    avg_down_dur8(d) = nanmean(synch_down_dur8{d});
    
    [up_hist(d,:),up_grid] = log_hist(synch_up_dur{d}(synch_up_dur{d} <= maxdur),up_range,numBins);
    [up_hist8(d,:),up_grid] = log_hist(synch_up_dur8{d}(synch_up_dur8{d} <= maxdur),up_range,numBins);
    
%     up_hist(d,:) = jmm_smooth_1d_cor(up_hist(d,:),2);
%     up_hist8(d,:) = jmm_smooth_1d_cor(up_hist8(d,:),2);
    
    [down_hist(d,:),down_grid] = log_hist(synch_down_dur{d}(synch_down_dur{d} <= maxdur),down_range,numBins);
    [down_hist8(d,:),down_grid] = log_hist(synch_down_dur8{d}(synch_down_dur8{d} <= maxdur),down_range,numBins);
    
%     down_hist(d,:) = jmm_smooth_1d_cor(down_hist(d,:),2);
%     down_hist8(d,:) = jmm_smooth_1d_cor(down_hist8(d,:),2);
   
%     if ~find(isempty(synch_up_dur8{d}>5))
%         disp('ERROR!!!')
%     end
%     [up_hist(d,:)] = gpkde(synch_up_dur{d}',-3,cur_grid);
%     [up_hist8(d,:)] = gpkde(synch_up_dur8{d}',-3,cur_grid);
%     
%     [down_hist(d,:)] = gpkde(synch_down_dur{d}',-3,cur_grid);
%     [down_hist8(d,:)] = gpkde(synch_down_dur8{d}',-3,cur_grid);
    
    
end

%%now for EC adjacent areas
load 'C:\WC_Germany\EC_adjacent_areas\UDS_synch_state_dur\UDS_synch_state_dur_data_new_method'
for d = 1:5
    [up_hist_adj(d,:),up_grid] = log_hist(synch_up_dur{d}(synch_up_dur{d} <= maxdur),up_range,numBins);
%     [up_hist8(d,:),up_grid] = log_hist(synch_up_dur8{d}(synch_up_dur8{d} <= maxdur),up_range,numBins);
    [down_hist_adj(d,:),down_grid] = log_hist(synch_down_dur{d}(synch_down_dur{d} <= maxdur),down_range,numBins);
%     [down_hist8(d,:),down_grid] = log_hist(synch_down_dur8{d}(synch_down_dur8{d} <= maxdur),down_range,numBins);  
end


cell_type(cell_type == 8) = [];
cell_type(cell_type == 7) = [];
cell_type(cell_type == 6) = 4;
cell_type(cell_type == 5) = 4;





pyr_up_hist = up_hist(cell_type==1,:);
pyr_up_hist8 = up_hist8(cell_type==1,:);
lec_up_hist = up_hist(cell_type==4,:);
lec_up_hist8 = up_hist8(cell_type==4,:);
pyr_down_hist = down_hist(cell_type==1,:);
pyr_down_hist8 = down_hist8(cell_type==1,:);
lec_down_hist = down_hist(cell_type==4,:);
lec_down_hist8 = down_hist8(cell_type==4,:);

pyr_avg_up_dur = avg_up_dur(cell_type==1);
pyr_avg_up_dur8 = avg_up_dur8(cell_type==1);
pyr_avg_down_dur = avg_down_dur(cell_type==1);
pyr_avg_down_dur8 = avg_down_dur8(cell_type==1);
lec_avg_up_dur = avg_up_dur(cell_type==4);
lec_avg_up_dur8 = avg_up_dur8(cell_type==4);
lec_avg_down_dur = avg_down_dur(cell_type==4);
lec_avg_down_dur8 = avg_down_dur8(cell_type==4);

%for lec layer3 pyramidal only
lec_up_hist([2 3 5 6 7 8 9],:) = [];
lec_down_hist([2 3 4 6 7 8 9],:) = [];
lec_up_hist8([2 3 5 6 7 8 9],:) = [];
lec_down_hist8([2 3 4 6 7 8 9],:) = [];
lec_avg_up_dur([2 3 5 6 7 8 9]) = [];
lec_avg_up_dur8([2 3 5 6 7 8 9]) = [];
lec_avg_down_dur([2 3 5 6 7 8 9]) = [];
lec_avg_down_dur8([2 3 5 6 7 8 9]) = [];


avg_pyr_up = nanmean(pyr_up_hist);
avg_pyr_up8 = nanmean(pyr_up_hist8);
avg_lec_up = nanmean(lec_up_hist);
avg_lec_up8 = nanmean(lec_up_hist8);
avg_adj_up = nanmean(up_hist_adj);
se_pyr_up = nanstd(pyr_up_hist)/sqrt(17);
se_pyr_up8 = nanstd(pyr_up_hist8)/sqrt(17);
se_lec_up = nanstd(lec_up_hist)/sqrt(4);
se_lec_up8 = nanstd(lec_up_hist8)/sqrt(4);
se_adj_up = nanstd(up_hist_adj)/sqrt(5);
avg_pyr_down = nanmean(pyr_down_hist);
avg_pyr_down8 = nanmean(pyr_down_hist8);
avg_lec_down = nanmean(lec_down_hist);
avg_lec_down8 = nanmean(lec_down_hist8);
avg_adj_down = nanmean(down_hist_adj);
se_pyr_down = nanstd(pyr_down_hist)/sqrt(17);
se_pyr_down8 = nanstd(pyr_down_hist8)/sqrt(17);
se_lec_down = nanstd(lec_down_hist)/sqrt(4);
se_lec_down8 = nanstd(lec_down_hist8)/sqrt(4);
se_adj_down = nanstd(down_hist_adj)/sqrt(5);

% figure(1)
% errorbar(up_grid,1-avg_pyr_up,se_pyr_up)
% hold on
% errorbar(up_grid,1-avg_lec_up,se_lec_up,'r')

save C:\WC_Germany\Persistent_activity\lec_state_dur_data_new_method *lec* *pyr* *adj* up_grid down_grid maxdur


