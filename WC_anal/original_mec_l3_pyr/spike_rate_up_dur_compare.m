clear all

load C:\WC_Germany\JMM_analysis_pyr\dir_tree_update
load C:\WC_Germany\JMM_analysis_pyr\UDS_dur_run_hist_v2\up_per_data_10_28
load C:\WC_Germany\JMM_analysis_pyr\sig_fit_data
load C:\WC_Germany\JMM_analysis_pyr\sig_down_fit_data

for d = 1:length(dir_array)
    
    cd(dir_array{d})
    pwd
    
    if ~exist('spike_time.mat')
        load spike_time_br
    else
        load spike_time
    end
    
    spike_id = round(spkid/8);
    
    mean_isi = zeros(length(up_trans{d}),1);
    cv_isi = mean_isi;
    mean_rate = mean_isi;
    median_isi = mean_isi;
    spike_count = mean_isi;
   for i = 1:length(up_trans{d})

       cur_spikes = spike_id(find(spike_id > up_trans{d}(i) & spike_id < down_trans{d}(i)));
       if length(cur_spikes) > 1
           cur_isis = diff(cur_spikes);
           cur_isis(1) = [];
           mean_isi(i) = mean(cur_isis);
           cv_isi(i) = mean_isi(i)/std(cur_isis);
           mean_rate(i) = length(cur_spikes)/up_state_dur{d}(i);
           median_isi(i) = median(cur_isis);
           spike_count(i) = length(cur_spikes);
       else
           mean_isi(i) = nan;
           cv_isi(i) = nan;
           mean_rate(i) = nan;
           median_isi(i) = nan;
           spike_count(i) = nan;
       end
       
   end
   
%    subplot(3,1,1)
   good_pts = ~isnan(mean_isi);
%    [a,b] = corrcoef(up_state_dur{d}(good_pts),mean_isi(good_pts));
%    plot(up_state_dur{d},mean_isi,'o')
%    title(sprintf('up dur v mean isi C = %0.2g  P = %0.2g',a(2,1),b(2,1)))

% subplot(2,1,1)
% [a,b] = corrcoef(up_state_dur{d}(good_pts),spike_count(good_pts));
% plot(up_state_dur{d},spike_count,'o')
%    title(sprintf('up dur v spike count C = %0.2g  P = %0.2g',a(2,1),b(2,1)))


%    subplot(2,1,2)
   [a,b] = corrcoef(up_state_dur{d}(good_pts),mean_rate(good_pts));
   plot(up_state_dur{d},mean_rate,'o')
   title(sprintf('up dur v mean rate C = %0.2g  P = %0.2g',a(2,1),b(2,1)))
    t_names = ['C:\WC_Germany\JMM_analysis_pyr\mean_rate_v_up_dur\' f_names{d}];
    print('-dpng',t_names)
    close all

%    subplot(3,1,3)
%    [a,b] = corrcoef(up_state_dur{d}(good_pts),median_isi(good_pts))
%    plot(up_state_dur{d},median_isi,'o')
%    title(sprintf('up dur v median isi C = %0.2g  P = %0.2g',a(2,1),b(2,1)))
    
end