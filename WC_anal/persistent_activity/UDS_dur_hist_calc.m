clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data 

for d = 1:length(dir_array)

%% Calculate distributions of up and down states
up_range = [0.3 10];
down_range = [0.3 10];
numBins = 100;

synch_up_dur = up_state_dur{d};
synch_down_dur = down_state_dur{d};
synch_up_dur8 = up_state_dur8{d};
synch_down_dur8 = down_state_dur8{d};

if min(synch_up_dur < up_range(1))
    disp('error up')
end
if min(synch_up_dur8 < up_range(1))
    disp('error up8')
end
if min(synch_down_dur < down_range(1))
    disp('error down')
end
if min(synch_down_dur8 < down_range(1))
    disp('error down8')
end

[up_hist(d,:),up_grid] = log_hist(synch_up_dur,up_range,numBins);
[up_hist8(d,:),up_grid] = log_hist(synch_up_dur8,up_range,numBins);

num_long_ups(d) = length(find(synch_up_dur > 10));

[down_hist(d,:),down_grid] = log_hist(synch_down_dur,down_range,numBins);
[down_hist8(d,:),down_grid] = log_hist(synch_down_dur8,down_range,numBins);

% stairs(up_grid,up_hist(d,:),'linewidth',2)
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% hold on
% stairs(up_grid,up_hist8(d,:),'r','linewidth',2)
% grid
% title('Up State Duration','FontSize',14)
% legend('MP','LFP')
% ylim([3e-3 0.4])
% xlim([0.2 100])
%     tname = ['C:\WC_Germany\Atropine\UDS_dist\up_state_' f_names{d}];
%     print('-dpng',tname);
%     close
% 
% stairs(down_grid,down_hist(d,:),'linewidth',2)
% set(gca,'xscale','log')
% set(gca,'yscale','log')
% hold on
% stairs(down_grid,down_hist8(d,:),'r','linewidth',2)
% grid
% title('Down State Duration','FontSize',14)
% legend('MP','LFP')
% ylim([3e-3 0.4])
% xlim([0.2 100])
%     tname = ['C:\WC_Germany\Atropine\UDS_dist\down_state_' f_names{d}];
%     print('-dpng',tname);
%     close

    d

end


save C:\WC_Germany\Atropine\UDS_dist\UDS_dist_data up_hist* up_grid down_hist* down_grid