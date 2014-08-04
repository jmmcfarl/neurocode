clear all
% close all

load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Persistent_activity\UDS_dist\synch_UDS_inds
load C:\WC_Germany\Persistent_activity\lf8_period_f_data2

Fs = 2016;
dsf = 8;
Fsd = Fs/dsf;
forwardlag = 10*Fsd;
backlag = 0;
binsize = 40;
lags = 0:binsize:forwardlag;

lf8lags = 0:forwardlag;

state_dur_hist = 0.3:1:10;
state_hist_binwidth = 0.25;
avg_rate_dur = nan(17,length(state_dur_hist));
for d = 1:length(dir_array)
    cd(dir_array{d})
    disp(['session ' num2str(d)])

    load spike_time_jmm
    load used_data wcv lf8
    wcv_d = downsample(wcv,dsf);
    datalen = length(wcv_d);

    spkids = round(spkid/dsf);

    %create vector of lfp state at all downsampled times
    lf8_state = zeros(datalen,1);
    for i = 1:length(synch_ups8{d})
       
        cur_up_trans = up_trans8{d}(synch_ups8{d}(i));
        cur_down_trans = down_trans8{d}(synch_ups8{d}(i));
        lf8_state(cur_up_trans:cur_down_trans) = 1;
        
    end
    
    lf8_state = logical(lf8_state);
    
    spkmat = nan(length(synch_ups{d}),length(lags));
    lf8statemat = nan(length(synch_ups{d}),length(lf8lags));
    tot_spikes = 0;
    tot_up_time = 0;
    net_spks_lup = 0;
    net_spks_ldown = 0;
    net_time_lup = 0;
    net_time_ldown = 0;
    
    rel_lf8_rate = nan(length(synch_ups{d}),1);
    cur_lup_rate = rel_lf8_rate;
    cur_ldown_rate = rel_lf8_rate;
    state_rate = rel_lf8_rate;
    state_dur = rel_lf8_rate;
    for i = 1:length(synch_ups{d})
        
        cur_up_trans = up_trans{d}(synch_ups{d}(i));
        cur_down_trans = down_trans{d}(synch_ups{d}(i));
        cur_dur = cur_down_trans-cur_up_trans;
        if cur_up_trans > backlag 
           cur_spikes = spkids(find(spkids > cur_up_trans-backlag & spkids <= cur_down_trans));
           lup_spikes = lf8_state(cur_spikes);
           ldown_spikes = ~lf8_state(cur_spikes);
           net_spks_lup = net_spks_lup + sum(lup_spikes);
           net_spks_ldown = net_spks_ldown + sum(ldown_spikes);
           net_time_lup = net_time_lup + sum(lf8_state(cur_up_trans:cur_down_trans))/Fsd;
           net_time_ldown = net_time_ldown + sum(~lf8_state(cur_up_trans:cur_down_trans))/Fsd;
           
           cur_lup_rate(i) = sum(lup_spikes)/sum(lf8_state(cur_up_trans:cur_down_trans))*Fsd;
           cur_ldown_rate(i) = sum(ldown_spikes)/sum(~lf8_state(cur_up_trans:cur_down_trans))*Fsd;
           if ~isinf(cur_lup_rate(i)) & ~isnan(cur_lup_rate(i)) & ~isinf(cur_ldown_rate(i)) & ...
                   ~isnan(cur_ldown_rate(i)) & cur_ldown_rate(i) > 0
                rel_lf8_rate(i) = cur_lup_rate(i)/cur_ldown_rate(i);
           else
               rel_lf8_rate(i) = nan;
           end
           
           cur_spikes = cur_spikes - cur_up_trans;
           tot_spikes = tot_spikes + length(cur_spikes);
           tot_up_time = tot_up_time + cur_dur/Fsd;
           state_rate(i) = length(cur_spikes)/cur_dur*Fsd;
           state_dur(i) = cur_dur/Fsd;
           spike_hist = hist(cur_spikes,lags);
           spike_hist(end) = nan;
           
           if cur_down_trans-cur_up_trans < forwardlag
                lf8statemat(i,1:cur_down_trans-cur_up_trans+1) = lf8_state(cur_up_trans:cur_down_trans);
           else
               lf8statemat(i,:) = lf8_state(cur_up_trans:cur_up_trans+forwardlag);
           end
           if cur_dur < forwardlag
               endpt = find(lags > cur_dur,1,'first');
               if isempty(endpt)
                   endpt = length(lags);
               end
               spike_hist(endpt:end) = nan;
           end
           spkmat(i,:) = spike_hist;
        else
            rel_lf8_rate(i) = nan;
            cur_lup_rate(i) = nan;
            cur_ldown_rate(i) = nan;
        end
        
    end
    
%     plot(lf8lags/Fsd,nanmean(lf8statemat))
%     hold on
%     plot(lf8lags/Fsd,jmm_smooth_1d_cor(nanmean(lf8statemat),30),'r')
%     pause
%     clf
    
    mean_lup_rate(d) = net_spks_lup/net_time_lup;
    mean_ldown_rate(d) = net_spks_ldown/net_time_ldown;
    
    mean_lf8_state(d,:) = jmm_smooth_1d_cor(nanmean(lf8statemat),50);
    clear lf8statemat
    
    cell_avg_up_rate(d) = tot_spikes/tot_up_time;
    avg_rate_tslu(d,:) = nanmean(spkmat)*Fsd/binsize;
    norm_rate_tslu(d,:) = avg_rate_tslu(d,:)/cell_avg_up_rate(d);
%     norm_rate_tslu(d,:) = jmm_smooth_1d_cor(norm_rate_tslu(d,:),2);
%     plot(lags/Fsd,jmm_smooth_1d_cor(avg_rate_tslu(d,:),10))
%     pause
%     clf
    
% scatter(cur_lup_rate,cur_ldown_rate)
% xlabel('LFP UP Rate')
% ylabel('LFP Down Rate')
% line([0 20],[0 20],'Color','k')
% pause
% clf

%create estimate of firing rate as a function of state dur
for s = 1:length(state_dur_hist)
    %find relevent MP up states
    cur_states = find(state_dur > state_dur_hist(s)-state_hist_binwidth/2 & ...
        state_dur <= state_dur_hist(s)+state_hist_binwidth/2);
    if ~isempty(cur_states)
        avg_rate_dur(d,s) = nanmean(state_rate(cur_states));
    end
end

avg_rate_dur(d,:) = avg_rate_dur(d,:)/cell_avg_up_rate(d);
% plot(state_dur_hist,avg_rate_dur(d,:))
% pause
% clf

end

% %view dependence of rate on state duration
% num_con_cells = sum(~isnan(avg_rate_dur));
% se_rate_dur = nanstd(avg_rate_dur)./sqrt(num_con_cells);
% errorbar(state_dur_hist,nanmean(avg_rate_dur),se_rate_dur)
% xlabel('Up state duration (s)')
% ylabel('Mean % of average firing rate')


% %number of contributing cells at each lag
num_cont = sum(~isnan(norm_rate_tslu));
se_norm_rate = nanstd(norm_rate_tslu)./sqrt(num_cont);
% se_norm_rate = nanstd(norm_rate_tslu);
mean_norm_rate = nanmean(norm_rate_tslu);
% 
% %smooth mean and se
% se_norm_rate = jmm_smooth_1d_cor(se_norm_rate,5);
% mean_norm_rate = jmm_smooth_1d_cor(mean_norm_rate,5);

u_norm_rate = mean_norm_rate + se_norm_rate;
d_norm_rate = mean_norm_rate - se_norm_rate;


mean_lf8_state_prob = nanmean(mean_lf8_state);
num_cont = sum(~isnan(mean_lf8_state));
se_lf8_state_prob = nanstd(mean_lf8_state)./sqrt(num_cont);
u_lf8_state = mean_lf8_state_prob+se_lf8_state_prob;
d_lf8_state = mean_lf8_state_prob-se_lf8_state_prob;

figure
plot(lags/Fsd,mean_norm_rate,'k')
hold on
X = [lags(1:end-1)/Fsd fliplr(lags(1:end-1)/Fsd)];
Y = [d_norm_rate(1:end-1) fliplr(u_norm_rate(1:end-1))];
fill(X,Y,'k')

X = [lf8lags/Fsd fliplr(lf8lags/Fsd)];
Y = [u_lf8_state fliplr(d_lf8_state)];
plot(lf8lags/Fsd,mean_lf8_state_prob,'r')
fill(X,Y,'r')

xlim([0 10])
ylim([0 1.5])
xlabel('Time since last MP up transition (s)')
ylabel('Fraction of mean firing rate')
% 
scatter(mean_lup_rate,mean_ldown_rate,'.')
xlabel('Average rate during LFP up state')
ylabel('Average rate during LFP down state')
line([0 12],[0 12],'Color','k')

