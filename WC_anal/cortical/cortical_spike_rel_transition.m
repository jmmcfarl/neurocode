clear all
close all
cd C:\WC_Germany\Cortical_analysis

load C:\WC_Germany\Cortical_analysis\cortical_dir
load C:\WC_Germany\Cortical_analysis\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Cortical_analysis\UDS_synch_state_dur\UDS_synch_state_dur_data
load C:\WC_Germany\Cortical_analysis\sigmoid_fit\sig_fit_alldata

Fs = 2016;

maxlag = 0.5;
xrange = linspace(-maxlag,maxlag,100);

for d = 1:length(sess_data)
    
    
    cd(sess_data(d).directory)
    disp(num2str(d))
    
    load spike_time_jmm
    t_50_wup{d} = rlid_wup{d};
    t_50_lup{d} = rlid_lup{d};
    
           
    %initialize
    tL50_tM50{d} = zeros(1,length(synch_ups{d}));
    tFS_tM50{d} = zeros(1,length(synch_ups{d}));
    tL50_tFS{d} = zeros(1,length(synch_ups{d}));
    tL50_tM10{d} = zeros(1,length(synch_ups{d}));
    
    for i = 1:length(t_50_wup{d})
        
       if ~isnan(t_50_wup{d}(i))
            [dummy,nearest_lfp_up] = min(abs(t_50_wup{d}(i) - t_50_lup{d}));
            next_spike = find(spkid > t_50_wup{d}(i),1,'first');
            tL50_tM50{d}(i) = t_50_lup{d}(nearest_lfp_up) - t_50_wup{d}(i);
            tL50_tM10{d}(i) = t_50_lup{d}(nearest_lfp_up) - t_10_wup{d}(i);
            if ~isempty(next_spike)
                tFS_tM50{d}(i) = spkid(next_spike) - t_50_wup{d}(i); 
                tL50_tFS{d}(i) = t_50_lup{d}(nearest_lfp_up) - spkid(next_spike);
            else
                tFS_tM50{d}(i) = nan;
                tL50_tFS{d}(i) = nan;
            end
       else
            tL50_tM50{d}(i) = nan;
            tFS_tM50{d}(i) = nan;
            tL50_tFS{d}(i) = nan;
            tL50_M10{d}(i) = nan;
       end
        
    end
    
    tL50_tM50{d} = tL50_tM50{d}/Fs;
    tFS_tM50{d} = tFS_tM50{d}/Fs;
    tL50_tFS{d} = tL50_tFS{d}/Fs;
    tL50_tM10{d} = tL50_tM10{d}/Fs;
    
%     figure
%     plot(tL50_tM50{d},tFS_tM50{d},'.')
%     figure
%     plot(tL50_tFS{d},tL50_tM50{d},'.')
    
    n_tL50_tM50(d,:) = histc(tL50_tM50{d},xrange);
    n_tFS_tM50(d,:) = histc(tFS_tM50{d},xrange);
    n_tL50_tFS(d,:) = histc(tL50_tFS{d},xrange);
    n_tL50_tM10(d,:) = histc(tL50_tM10{d},xrange);
    
    n_tL50_tM50(d,:) = n_tL50_tM50(d,:)/sum(n_tL50_tM50(d,:));
     n_tFS_tM50(d,:) = n_tFS_tM50(d,:)/sum(n_tFS_tM50(d,:));
       n_tL50_tFS(d,:) = n_tL50_tFS(d,:)/sum(n_tL50_tFS(d,:));
    n_tL50_tM10(d,:) = n_tL50_tM10(d,:)/sum(n_tL50_tM10(d,:));

    med_tL50_tM50(d) = nanmedian(tL50_tM50{d});
    med_tFS_tM50(d) = nanmedian(tFS_tM50{d});
    med_tL50_tFS(d) = nanmedian(tL50_tFS{d});
    med_tL50_tM10(d) = nanmedian(tL50_tM10{d});
%     figure
%     bar(xrange,n_tL50_tM50(d,:))
%     xlabel('Time (s)')
%     title('L50 - M50')
%     cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
%     t_names = ['C:\WC_Germany\Cortical_analysis\spike_timing\tL50_tM50_' cell_name];
%     print('-dpng',t_names);
%     close
% % 
%     figure
%     bar(xrange,n_tFS_tM50(d,:))
%     xlabel('Time (s)')
%     title('FS - M50')
%     cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
%     t_names = ['C:\WC_Germany\Cortical_analysis\spike_timing\tFS_tM50_' cell_name];
%     print('-dpng',t_names);
%     close
% % 
%     figure
%     bar(xrange,n_tL50_tFS(d,:))
%     xlabel('Time (s)')
%     title('L50 - FS')
%     cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
%     t_names = ['C:\WC_Germany\Cortical_analysis\spike_timing\tL50_tFS_' cell_name];
%     print('-dpng',t_names);
%     close
% 
%         figure
%     bar(xrange,n_tL50_tM10(d,:))
%     xlabel('Time (s)')
%     title('L50 - M10')
%     cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
%     t_names = ['C:\WC_Germany\Cortical_analysis\spike_timing\tL50_M10_' cell_name];
%     print('-dpng',t_names);
%     close

end

cd C:\WC_Germany\Cortical_analysis\spike_timing
save relative_timing_data tL50* tFS* n_t*