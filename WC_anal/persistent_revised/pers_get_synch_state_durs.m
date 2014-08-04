clear all

load C:\WC_Germany\persistent_revised\pers_revised_dir
load C:\WC_Germany\persistent_revised\desynch_detect\desynch_times
load C:\WC_Germany\persistent_revised\UDS_dur_raw\UDS_raw_data

Fsd = 2016/8;

% Calculate distributions of up and down states
up_range = [0.3 10];
down_range = [0.3 10];
numBins = 50;


for d = 1:28
   
    cd(dir_array{d})
    load used_data lf8
   
    datalen = length(lf8)/8;
    
    synch_ups{d} = 1:length(up_trans{d});
    synch_downs{d} = 1:length(down_trans{d});
    synch_ups8{d} = 1:length(up_trans8{d});
    synch_downs8{d} = 1:length(down_trans8{d});
    
    desynch_log = zeros(1,datalen);
    
    for w = 1:length(desynch_start_times{d})
       cur_up_id = desynch_start_times{d}(w)*Fsd;
       cur_down_id = desynch_stop_times{d}(w)*Fsd;
       desynch_log(cur_up_id:cur_down_id) = 1;
    end
    
    desynch_log = logical(desynch_log);
    
    bad_trans = [];
    bad_trans8 = [];
    for i = 1:length(up_trans{d})
       if desynch_log(up_trans{d}(i)) | desynch_log(down_trans{d}(i))
          bad_trans = [bad_trans i]; 
       end
    end
    for i = 1:length(up_trans8{d})
       if desynch_log(up_trans8{d}(i)) | desynch_log(down_trans8{d}(i))
          bad_trans8 = [bad_trans8 i]; 
       end
    end

    synch_ups{d}(bad_trans) = [];
    synch_downs{d}(bad_trans) = [];
    synch_ups8{d}(bad_trans8) = [];
    synch_downs8{d}(bad_trans8) = [];
    
    synch_up_dur{d} = (down_trans{d}(synch_downs{d})-up_trans{d}(synch_ups{d}))/Fsd;  
    synch_down_dur{d} = (up_trans{d}(synch_ups{d}(2:end))-down_trans{d}(synch_downs{d}(1:end-1)))/Fsd;  
    synch_up_dur8{d} = (down_trans8{d}(synch_downs8{d})-up_trans8{d}(synch_ups8{d}))/Fsd;  
    synch_down_dur8{d} = (up_trans8{d}(synch_ups8{d}(2:end))-down_trans8{d}(synch_downs8{d}(1:end-1)))/Fsd;  
    
    
    [up_hist(d,:),up_grid] = log_hist(synch_up_dur{d},up_range,numBins);
    [up_hist8(d,:),up_grid] = log_hist(synch_up_dur8{d},up_range,numBins);


    [down_hist(d,:),down_grid] = log_hist(synch_down_dur{d},down_range,numBins);
    [down_hist8(d,:),down_grid] = log_hist(synch_down_dur8{d},down_range,numBins);

%     stairs(up_grid,up_hist(d,:),'linewidth',2)
%     set(gca,'yscale','log')
%     hold on
%     stairs(up_grid,up_hist8(d,:),'r','linewidth',2)
%     grid
%     title('Up State Duration','FontSize',14)
%     legend('MP','LFP')
%     ylim([1e-3 0.2])
%     xlim([0.2 10])
%     t_names = ['C:\WC_Germany\persistent_revised\UDS_synch_state_dur\up_state_' f_names{d}];
%         print('-dpng',t_names);
%         close
%     
%     stairs(down_grid,down_hist(d,:),'linewidth',2)
%     set(gca,'yscale','log')
%     hold on
%     stairs(down_grid,down_hist8(d,:),'r','linewidth',2)
%     grid
%     title('Down State Duration','FontSize',14)
%     legend('MP','LFP')
%     ylim([1e-3 0.2])
%     xlim([0.2 10])
%     t_names = ['C:\WC_Germany\persistent_revised\UDS_synch_state_dur\down_state_' f_names{d}];
%         print('-dpng',t_names);
%         close

net_up_time(d) = sum(synch_up_dur{d});
net_down_time(d) = sum(synch_down_dur{d});
net_up_time8(d) = sum(synch_up_dur8{d});
net_down_time8(d) = sum(synch_down_dur8{d});

mean_up_dur(d) = mean(synch_up_dur{d});
mean_down_dur(d) = mean(synch_down_dur{d});
mean_up_dur8(d) = mean(synch_up_dur8{d});
mean_down_dur8(d) = mean(synch_down_dur8{d});
% duty_cycle8{d} = synch_up_dur8{d}./(synch_up_dur8{d}+synch_down_dur8{d});
% mean_duty_cycle8(d) = mean(duty_cycle8{d});
num_long_pers(d) = length(find(synch_up_dur{d} > 10));
total_number_up_states(d) = length(synch_up_dur{d});
total_number_up_states8(d) = length(synch_up_dur8{d});
total_number_down_states(d) = length(synch_down_dur{d});
total_number_down_states8(d) = length(synch_down_dur8{d});

fraction_chopped_up(d) = tot_chopped_up(d)/total_number_up_states(d);
fraction_chopped_up8(d) = tot_chopped_up8(d)/total_number_up_states8(d);
fraction_chopped_down(d) = tot_chopped_down(d)/total_number_down_states(d);
fraction_chopped_down8(d) = tot_chopped_down8(d)/total_number_down_states8(d);

end

mec_cells = 1:21;
lec_cells = 22:28;
lec_cells(4) = [];
mec_cells(21) = [];

fraction_up_time = net_up_time./(net_up_time+net_down_time);
fraction_up_time8 = net_up_time8./(net_up_time8+net_down_time8);

mec_up_time = fraction_up_time(mec_cells);
lec_up_time = fraction_up_time(lec_cells);
mec_up_time8 = fraction_up_time8(mec_cells);
lec_up_time8 = fraction_up_time8(lec_cells);

mec_up_dur = mean_up_dur(mec_cells);
lec_up_dur = mean_up_dur(lec_cells);
mec_up_dur8 = mean_up_dur8(mec_cells);
lec_up_dur8 = mean_up_dur8(lec_cells);

mec_down_dur = mean_down_dur(mec_cells);
lec_down_dur = mean_down_dur(lec_cells);
mec_down_dur8 = mean_down_dur8(mec_cells);
lec_down_dur8 = mean_down_dur8(lec_cells);


% save C:\WC_Germany\persistent_revised\UDS_synch_state_dur\UDS_synch_state_dur_data synch* *hist* *grid net*