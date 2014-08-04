clear all

load('C:\WC_Germany\Persistent_activity\dir_tree_update.mat')
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data_new_method
load C:\WC_Germany\Persistent_activity\desynch_extract\desynch_points_new_method

% load('C:\WC_Germany\Persistent_activity\dir_tree_update.mat')
% load C:\WC_Germany\Persistent_activity\UDS_dur_raw_no_cut\UDS_raw_data_no_cut
% load C:\WC_Germany\Persistent_activity\desynch_extract_no_cut\desynch_points_no_cut


for d = 1:length(dir_array)

    %% Calculate distributions of up and down states
    up_range = [0.3 50];
    down_range = [0.3 150];
    numBins = 100;

    synch_ups{d} = 1:length(up_trans{d});
    synch_downs{d} = 1:length(down_trans{d});
    synch_ups8{d} = 1:length(up_trans8{d});
    synch_downs8{d} = 1:length(down_trans8{d});
    
    synch_downs8{d}(end) = [];
    synch_downs{d}(end) = [];
    
    % if there are any MP up states during desynch times get rid of them
    if ~isempty(desynch_start{d})

        bad_trans = [];
        bad_trans8 = [];
        for i = 1:length(desynch_start{d})
            bad_trans = [bad_trans ...
                find(up_trans{d} > desynch_start{d}(i) & up_trans{d} < desynch_stop{d}(i))];
            bad_trans = [bad_trans ...
                find(down_trans{d} > desynch_start{d}(i) & down_trans{d} < desynch_stop{d}(i))];

            bad_trans8 = [bad_trans8 ...
                find(up_trans8{d} > desynch_start{d}(i) & up_trans8{d} < desynch_stop{d}(i))];
            bad_trans8 = [bad_trans8 ...
                find(down_trans8{d} > desynch_start{d}(i) & down_trans8{d} < desynch_stop{d}(i))];

        end
        
        bad_trans = unique(bad_trans);
        bad_trans8 = unique(bad_trans8);
        synch_ups{d}(bad_trans) = [];
        synch_downs{d}(bad_trans) = [];
        synch_ups8{d}(bad_trans8) = [];
        synch_downs8{d}(bad_trans8) = [];
        
    end

    synch_up_dur{d} = up_state_dur{d}(synch_ups{d});
    synch_down_dur{d} = down_state_dur{d}(synch_downs{d});
    synch_up_dur8{d} = up_state_dur8{d}(synch_ups8{d});
    synch_down_dur8{d} = down_state_dur8{d}(synch_downs8{d});
    
    tot_short_up(d) = length(find(up_state_dur{d}(synch_ups{d}) < 0.3))/length(synch_ups{d});
    tot_short_down(d) = length(find(down_state_dur{d}(synch_downs{d}) < 0.3))/length(synch_downs{d});
    tot_short_up8(d) = length(find(up_state_dur8{d}(synch_ups8{d}) < 0.3))/length(synch_ups8{d});
    tot_short_down8(d) = length(find(down_state_dur8{d}(synch_downs8{d}) < 0.3))/length(synch_downs8{d});
    
    
    if ~isempty(find(synch_up_dur8{d} > 5))
        disp('Error!')
    end
    
    [up_hist(d,:),up_grid] = log_hist(synch_up_dur{d},up_range,numBins);
    [up_hist8(d,:),up_grid] = log_hist(synch_up_dur8{d},up_range,numBins);


    [down_hist(d,:),down_grid] = log_hist(synch_down_dur{d},down_range,numBins);
    [down_hist8(d,:),down_grid] = log_hist(synch_down_dur8{d},down_range,numBins);

    stairs(up_grid,up_hist(d,:),'linewidth',2)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    hold on
    stairs(up_grid,up_hist8(d,:),'r','linewidth',2)
    grid
    title('Up State Duration','FontSize',14)
    legend('MP','LFP')
    ylim([3e-3 0.4])
    xlim([0 20])
    t_names = ['C:\WC_Germany\Persistent_activity\UDS_synch_state_dur\up_state_new_method_' f_names{d}];
        print('-dpng',t_names);
        close
    
    stairs(down_grid,down_hist(d,:),'linewidth',2)
    set(gca,'xscale','log')
    set(gca,'yscale','log')
    hold on
    stairs(down_grid,down_hist8(d,:),'r','linewidth',2)
    grid
    title('Down State Duration','FontSize',14)
    legend('MP','LFP')
    ylim([3e-3 0.4])
    xlim([0 20])
    t_names = ['C:\WC_Germany\Persistent_activity\UDS_synch_state_dur\down_state_new_method_' f_names{d}];
        print('-dpng',t_names);
        close

    d

end

% figure
% plot(up_grid,nanmean(up_hist))
% hold on
% plot(up_grid,nanmean(up_hist8),'r')
% 
% figure
% plot(down_grid,nanmean(down_hist))
% hold on
% plot(down_grid,nanmean(down_hist8),'r')

save C:\WC_Germany\Persistent_activity\UDS_synch_state_dur\UDS_synch_state_dur_data_new_method synch* *grid *hist*
% save C:\WC_Germany\Persistent_activity\UDS_synch_state_dur_no_cut\UDS_synch_state_dur_data_no_cut synch* *grid *hist*