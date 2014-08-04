clear all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data


%% set parameters
Fs = 2016;
niqf = Fs/2;
lcf = 10/niqf;
[b,a] = butter(2,lcf,'high');
look_back_bins = 4;
look_forward_bins = 15;
total_bins = look_back_bins+look_forward_bins+1;
width_range = linspace(0.002,0.4,200);
isi_curve = zeros(17,length(width_range));
isi_bins = linspace(0,0.5,500);

for d = 1:length(dir_array)
d=5
    cd(dir_array{d})
    pwd

    load used_data wcv

%     wcv_f = wavefilter(wcv',7);
wcv_f = filtfilt(b,a,wcv);
    wcv_f = zscore(wcv_f);
    wcv_z = zscore(wcv);
    
%     wcv = filtfilt(b2,a2,wcv);
    
    %locate spike times in filtered recording
    [spike_amp,spike_loc] = findpeaks(wcv_f,'minpeakheight',3);
    
%     plot(wcv_z)
%     hold on
%     plot(spike_loc,ones(size(spike_loc)),'r*')
    
    isis = [0 diff(spike_loc)]/Fs;
    spike_hist_mat = zeros(length(up_trans{d}),length(isi_bins));
    
    for i = 1:length(up_trans{d})
       
        cur_spikes = find(spike_loc/8 > up_trans{d}(i) & spike_loc/8 < down_trans{d}(i));
        burst_fract(i) = length(find(isis(cur_spikes) < .015))/length(cur_spikes);
        spike_hist_mat(i,:) = hist(isis(cur_spikes),isi_bins);
        spike_hist_mat(i,:) = spike_hist_mat(i,:)/sum(spike_hist_mat(i,:));
    end
    
    [dummy,up_order] = sort(up_state_dur{d});
    spike_hist_mat(:,end) = [];
    for i = 1:length(up_trans{d})
        spike_hist_mat(i,:) = jmm_smooth_1d(spike_hist_mat(i,:),2);
    end
    spike_hist_mat(isnan(spike_hist_mat)) = 0;
    pcolor(spike_hist_mat(up_order,:));shading flat
    caxis([0 0.1]);colorbar
    
    small_isis = isis(isis < 1);
    hist(small_isis,1000)
    xlim([0 0.5])
    
    
end