clear all
close all
cd C:\WC_Germany\Cortical_analysis

load C:\WC_Germany\Cortical_analysis\cortical_dir
load C:\WC_Germany\Cortical_analysis\desynch_detect\desynch_times
load C:\WC_Germany\Cortical_analysis\UDS_dur_raw\UDS_raw_data

Fsd = 2016/8;

% Calculate distributions of up and down states
up_range = [0.3 10];
down_range = [0.3 10];
numBins = 50;


for d = 1:length(sess_data)
   
    cd(sess_data(d).directory)
    load used_data lf8 
   
    datalen = length(lf8)/8;
    
    synch_ups{d} = 1:length(up_trans{d});
    synch_downs{d} = 1:length(down_trans{d});
    synch_ups8{d} = 1:length(up_trans8{d});
    synch_downs8{d} = 1:length(down_trans8{d});
%     synch_ups5{d} = 1:length(up_trans5{d});
%     synch_downs5{d} = 1:length(down_trans5{d});

    desynch_log = zeros(1,datalen);
    
    for w = 1:length(desynch_start_times{d})
       cur_up_id = desynch_start_times{d}(w)*Fsd;
       cur_down_id = desynch_stop_times{d}(w)*Fsd;
       desynch_log(cur_up_id:cur_down_id) = 1;
    end
    
    desynch_log = logical(desynch_log);
    
    bad_trans = [];
    bad_trans8 = [];
%     bad_trans5 = [];

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
%     for i = 1:length(up_trans5{d})
%        if desynch_log(up_trans5{d}(i)) | desynch_log(down_trans5{d}(i))
%           bad_trans5 = [bad_trans5 i]; 
%        end
%     end
    synch_ups{d}(bad_trans) = [];
    synch_downs{d}(bad_trans) = [];
    synch_ups8{d}(bad_trans8) = [];
    synch_downs8{d}(bad_trans8) = [];
%     synch_ups5{d}(bad_trans5) = [];
%     synch_downs5{d}(bad_trans5) = [];
    
    synch_up_dur{d} = (down_trans{d}(synch_downs{d})-up_trans{d}(synch_ups{d}))/Fsd;  
    synch_down_dur{d} = (up_trans{d}(synch_ups{d}(2:end))-down_trans{d}(synch_downs{d}(1:end-1)))/Fsd;  
    synch_up_dur8{d} = (down_trans8{d}(synch_downs8{d})-up_trans8{d}(synch_ups8{d}))/Fsd;  
    synch_down_dur8{d} = (up_trans8{d}(synch_ups8{d}(2:end))-down_trans8{d}(synch_downs8{d}(1:end-1)))/Fsd;  
%     synch_up_dur5{d} = (down_trans5{d}(synch_downs5{d})-up_trans5{d}(synch_ups5{d}))/Fsd;  
%     synch_down_dur5{d} = (up_trans5{d}(synch_ups5{d}(2:end))-down_trans5{d}(synch_downs5{d}(1:end-1)))/Fsd;  
    
    [up_hist(d,:),up_grid] = log_hist(synch_up_dur{d},up_range,numBins);
    [up_hist8(d,:),up_grid] = log_hist(synch_up_dur8{d},up_range,numBins);
%     [up_hist5(d,:),up_grid] = log_hist(synch_up_dur5{d},up_range,numBins);


    [down_hist(d,:),down_grid] = log_hist(synch_down_dur{d},down_range,numBins);
    [down_hist8(d,:),down_grid] = log_hist(synch_down_dur8{d},down_range,numBins);
%     [down_hist5(d,:),down_grid] = log_hist(synch_down_dur5{d},down_range,numBins);

    stairs(up_grid,up_hist(d,:),'linewidth',2)
    set(gca,'yscale','log')
    hold on
    stairs(up_grid,up_hist8(d,:),'r','linewidth',2)
%     stairs(up_grid,up_hist5(d,:),'k','linewidth',2)
    grid
    title('Up State Duration','FontSize',14)
    legend('MP','LF8')
    ylim([1e-3 0.2])
    xlim([0.2 10])
    cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
    t_names = ['C:\WC_Germany\Cortical_analysis\UDS_synch_state_dur\up_state_' cell_name];
        print('-dpng',t_names);
        close
    
    stairs(down_grid,down_hist(d,:),'linewidth',2)
    set(gca,'yscale','log')
    hold on
    stairs(down_grid,down_hist8(d,:),'r','linewidth',2)
%     stairs(down_grid,down_hist5(d,:),'k','linewidth',2)
    grid
    title('Down State Duration','FontSize',14)
    legend('MP','LF8')
    ylim([1e-3 0.2])
    xlim([0.2 10])
    cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
    t_names = ['C:\WC_Germany\Cortical_analysis\UDS_synch_state_dur\down_state_' cell_name];
        print('-dpng',t_names);
        close

end

save C:\WC_Germany\Cortical_analysis\UDS_synch_state_dur\UDS_synch_state_dur_data synch* *hist* *grid