clear all

load C:\WC_Germany\JMM_analysis_ste\dir_tree_ste
load C:\WC_Germany\JMM_analysis_ste\UDS_dur_run_hist\data


%% set parameters
Fs = 2016;
niqf = Fs/2;
lcf = 80/niqf;
[b,a] = butter(2,lcf,'high');
look_back_bins = 4;
look_forward_bins = 15;
total_bins = look_back_bins+look_forward_bins+1;
width_range = linspace(0.002,0.4,200);
isi_curve = zeros(17,length(width_range));
% lcf_2 = 1/niqf;
% [b2,a2] = butter(2,lcf_2,'high');

for d = 1:length(dir_array)

    cd(dir_array{d})
pwd

    load used_data wcv

%     wcv_f = filtfilt(b,a,wcv);
wcv_f = wavefilter(wcv',7);
    wcv_f = zscore(wcv_f);
    
%     wcv = filtfilt(b2,a2,wcv);
    
    %locate spike times in filtered recording
    [spike_amp,spike_loc] = findpeaks(wcv_f,'minpeakheight',5);
    bad_spikes = find(spike_loc > length(wcv)-look_forward_bins);
    spike_amp(bad_spikes) = [];
    spike_loc(bad_spikes) = [];
    
    spike_amp = zscore(spike_amp);
    
    %initialize spike matrix
    num_spikes = length(spike_loc);

    spkMat = zeros(num_spikes,total_bins);
    for i = 1:num_spikes
        spkMat(i,:) = wcv(spike_loc(i)-look_back_bins:spike_loc(i)+look_forward_bins);
    end

    int_spkMat = spline([1:total_bins],spkMat,[1:.1:total_bins]);

    %% calculate spike widths
    spk_widths = zeros(1,num_spikes);
    for i = 1:num_spikes

        spk_peak = max(int_spkMat(i,:));
        spk_min = min(int_spkMat(i,:));
        thresh = (spk_peak-spk_min)/2+spk_min;
        spk_start = find(int_spkMat(i,:)>thresh,1,'first');
        spk_end = spk_start+find(int_spkMat(i,spk_start:end)<thresh,1,'first');
        if ~isempty(spk_end)
            spk_widths(i) = (spk_end-spk_start)*.1/Fs;
        else
            spk_widths(i) = nan;
        end

    end
    
mean_spk_width(d) = nanmean(spk_widths);
median_spk_width(d) = nanmedian(spk_widths);
std_spk_width(d) = nanstd(spk_widths);


%     z_spk_widths = spk_widths - nanmean(spk_widths);
%     z_spk_widths = z_spk_widths/nanstd(z_spk_widths);
% 
%     %% calculate isis
%     isis = [0 diff(spike_loc)]/Fs;
%     
%     burst_rat(d) = length(find(isis < 0.05))/length(isis);
%     log_isis = log(isis);
%     log_isis(isinf(log_isis)) = nan;
%     burst_isis = isis;
%     burst_isis(burst_isis > 0.1) = nan;

% calculate width isi relationship
% for r = 1:length(width_range)-1
%     cur_spikes = find(isis >=width_range(r) & isis < width_range(r+1));
%     if ~isempty(cur_spikes)
%         isi_curve(d,r+1) = nanmean(spk_widths(cur_spikes));
%     else
%         isi_curve(d,r+1) = nan;
%     end
% end

%     %% cycle through up states and look at trends within up states
%     num_ups = length(up_trans{d});
%     up_transitions = up_trans{d}*8;
%     down_transitions = down_trans{d}*8;
%     for i = 1:num_ups
%        
%         %find spikes in that up state
%         cur_spikes = find(spike_loc > up_transitions(i) & spike_loc < down_transitions(i));
%         if ~isempty(cur_spikes)
%             
%             avg_spike{d}(i,:) = mean(int_spkMat(cur_spikes,:));
%             
%         beg_ind = spike_loc(cur_spikes(1));
%         cur_locs = spike_loc(cur_spikes)-beg_ind;
%         
%         cur_locs = cur_locs';
%         cur_widths = spk_widths(cur_spikes)';
%         
%         B = [ones(length(cur_locs),1) cur_locs]\cur_widths;
%         
%         width_slope{d}(i) = B(2);
%         
%         cur_amps = spike_amp(cur_spikes)';
%         
%         B = [ones(length(cur_locs),1) cur_locs]\cur_amps;
%         
%         amp_slope{d}(i) = B(2);
%         
%         mean_width{d}(i) = nanmean(cur_widths);
%         mean_amp{d}(i) = mean(cur_amps);
%         
%         B = [ones(length(cur_locs),1) cur_amps]\cur_widths;
%         
%         width_amp_slope(i) = B(2);
%         
%                     mean_isi{d}(i) = mean(isis(cur_spikes(2:end)));
%                     median_isi{d}(i) = median(isis(cur_spikes(2:end)));
%                     mean_logisi(i) = nanmean(log(isis(cur_spikes(2:end))));
%                     mean_rate{d}(i) = length(cur_spikes)/up_state_dur{d}(i);
%                     
%         else
%             
%             width_slope{d}(i) = nan;
%             width_amp_slope(i) = nan;
%             mean_amp{d}(i) = nan;
%             mean_width{d}(i) = nan;
%             amp_slope{d}(i) = nan;
%             mean_rate{d}(i) = nan;
%             mean_isi{d}(i) = nan;
%             
%         end
%         
%         if length(cur_spikes) > 5
%              min_isi_5{d}(i) = min(isis(cur_spikes(2:5)));
%              beg_spk_width{d}(i) = mean(cur_widths(1:5));
%              mean_isi_5{d}(i) = mean(isis(cur_spikes(2:end)));
% 
%         else
%              min_isi_5{d}(i) = nan;
%              beg_spk_width{d}(i) = nan;
%              mean_isi_5{d}(i) = nan;
% 
%         end
%     end
%     
%    burst2 = find(isis < 0.05);
%     
%    burst2(1) = []; %get rid of initial isi artifact
%    burst3 = find(diff(burst2) == 1);
%    burst2rat(d) = length(burst2)/length(spike_loc);
%    burst3rat(d) = length(burst3)/length(spike_loc);
%    burst4rat(d) = length(find(diff(burst3)==1))/length(spike_loc);
%     
    
%     burst_isis = burst_isis - nanmean(burst_isis);
%     burst_isis = burst_isis/nanstd(burst_isis);
%     
% z_w = zscore(wcv);
% z_spk_widths = spk_widths - nanmean(spk_widths);
% z_spk_widths = z_spk_widths/nanstd(z_spk_widths);
% 
%     plot(z_w);
%     hold on
% % %     plot(spike_loc,burst_isis,'r.')
%     plot(spike_loc,z_spk_widths,'r.')
%     plot(spike_loc,spike_amp,'k.')
%     plot(spike_loc,z_isis,'g.')


%     good_pts = find(~isnan(width_slope{d}));
%     good_pts_beg = find(~isnan(beg_spk_width{d}));
%     
% 	plot(up_state_dur{d},mean_width{d},'.')
%     [ct,pt] = corrcoef(up_state_dur{d}(good_pts),mean_width{d}(good_pts));
%     title(sprintf('duration v width  c = %0.2g  p = %0.2g',ct(2,1),pt(2,1)))
%     t_names = ['C:\WC_Germany\JMM_analysis_pyr\spike_analysis\dur_width_cor_' f_names{d}];
%     print('-dpng',t_names)
%     
%    	plot(up_state_dur{d},mean_rate{d},'.')
%     [ct,pt] = corrcoef(up_state_dur{d}(good_pts),mean_rate{d}(good_pts));
%     title(sprintf('duration v rate  c = %0.2g  p = %0.2g',ct(2,1),pt(2,1)))
%     t_names = ['C:\WC_Germany\JMM_analysis_pyr\spike_analysis\dur_rate_cor_' f_names{d}];
%     print('-dpng',t_names)
%  
%     plot(up_state_dur{d},mean_isi{d},'.')
%     [ct,pt] = corrcoef(up_state_dur{d}(good_pts),mean_isi{d}(good_pts));
%     title(sprintf('duration v isi  c = %0.2g  p = %0.2g',ct(2,1),pt(2,1)))
%     t_names = ['C:\WC_Germany\JMM_analysis_pyr\spike_analysis\dur_isi_cor_' f_names{d}];
%     print('-dpng',t_names)
% 
%     plot(up_state_dur{d},beg_spk_width{d},'.')
%     [ct,pt] = corrcoef(up_state_dur{d}(good_pts_beg),beg_spk_width{d}(good_pts_beg));
%     title(sprintf('duration v beg width  c = %0.2g  p = %0.2g',ct(2,1),pt(2,1)))
%     t_names = ['C:\WC_Germany\JMM_analysis_pyr\spike_analysis\dur_beg_width_cor_' f_names{d}];
%     print('-dpng',t_names)
% 
%     
%     plot(up_state_dur{d},mean_amp{d},'.')
%     [ct,pt] = corrcoef(up_state_dur{d}(good_pts),mean_amp{d}(good_pts));
%     title(sprintf('duration v amp  c = %0.2g  p = %0.2g',ct(2,1),pt(2,1)))
%     t_names = ['C:\WC_Germany\JMM_analysis_pyr\spike_analysis\dur_amp_cor_' f_names{d}];
%     print('-dpng',t_names)
% 
%     plot(up_state_dur{d},width_slope{d},'.')
%     [ct,pt] = corrcoef(up_state_dur{d}(good_pts),width_slope{d}(good_pts));
%     title(sprintf('duration v widthslope  c = %0.2g  p = %0.2g',ct(2,1),pt(2,1)))
%     t_names = ['C:\WC_Germany\JMM_analysis_pyr\spike_analysis\dur_widthslope_cor_' f_names{d}];
%     print('-dpng',t_names)
% 
%     
%     plot(up_state_dur{d},fan_isi{d},'.')
%     [ct,pt] = corrcoef(up_state_dur{d}(good_pts_beg),fan_isi{d}(good_pts_beg));
%     title(sprintf('duration v fanisi  c = %0.2g  p = %0.2g',ct(2,1),pt(2,1)))
%     t_names = ['C:\WC_Germany\JMM_analysis_pyr\spike_analysis\dur_fanisi_cor_' f_names{d}];
%     print('-dpng',t_names)

% isi_range = linspace(-2,1,200);
% lisis = log10(isis);
% width_fun(d,:) = zeros(1,length(isi_range)-1);
% for i = 1:length(isi_range)-1
%     cur_s = find(lisis>=isi_range(i) & lisis < isi_range(i+1));
%     if ~isempty(cur_s)
%         width_fun(d,i) = nanmean(spk_widths(cur_s));
%     end
% end
% plot(lisis,spk_widths,'r.')
% hold on
% plot(isi_range(2:end),width_fun(d,:),'k.','MarkerSize',20)
%     t_names = ['C:\WC_Germany\JMM_analysis_pyr\spike_analysis\isi_width_fun' f_names{d}];
%     title('spike widths vs log isis')
%     xlim([-2.5 1])
%     ylim([1 7]*1e-3)
%     print('-dpng',t_names)
% close all
%     t_names = ['C:\WC_Germany\JMM_analysis_pyr\spike_analysis\isi_width_fun' f_names{d}];
%     title('spike widths vs log isis')
%     xlim([-2.5 1])
%     ylim([0 7]*1e-3)
%     print('-dpng',t_names)

% plot(log10(isis),spk_widths,'.')
%     t_names = ['C:\WC_Germany\JMM_analysis_pyr\spike_analysis\isi_width_' f_names{d}];
%     title('spike widths vs log isis')
%     xlim([-2.5 1])
%     ylim([0 7]*1e-3)
%     print('-dpng',t_names)
    
end


% cd C:\WC_Germany\Persistent_activity\
% save spk_data burst_rat mean_spk_width median_spk_width