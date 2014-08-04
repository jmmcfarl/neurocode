%% depolarization steps
close all
clear all

load C:\WC_Germany\current_injection\stellate\ste_current_dir

for d = 1:length(dir_array)
    cd(dir_array{d})
    load used_data
    load lf8_period_f
    load spike_time_jmm
    
    first_up = first_dep{d};
    first_down = first_hyp{d};
    wait_dur = 2260;
    
    %process lfps
    niqf = 2016/2;
    [b,a] = butter(2,[0.1/niqf 40/niqf]);
    lf8_f = filtfilt(b,a,lf8);
    lf8_f = zscore(lf8_f);
    lf3_f = filtfilt(b,a,lf3);
    lf3_f = zscore(lf3_f);
    wcv_z = zscore(wcv_minus_spike);

    lf8_phase = lf8_period_p{1};

    t = (1:length(wcv))/2016;

    %get a derivative signal without any spikes
    wcv_deriv = [0;diff(wcv)];
    wcv_deriv = zscore(wcv_deriv);
    %find all instances when absolute derivative exceeds 3 std
    [dummy,deriv_peaks] = findpeaks(abs(wcv_deriv),'minpeakheight',3);
    inter_peak_int = diff(deriv_peaks);
    %if two peaks are less than 2 samples apart get rid of both of them
    bad_peaks = find(inter_peak_int <= 5);
    bad_peaks = [bad_peaks (bad_peaks+1)];
    deriv_peaks(bad_peaks) = [];

    %find any spikes detected within 10 samples of a transition and get rid of
    %it
    bad_spikes = [];
    for i = 1:length(spkid)
       cur_dist = min(abs(deriv_peaks-spkid(i)));
       if cur_dist <= 10
           bad_spikes = [bad_spikes i];
       end
    end
    spkid(bad_spikes) = [];

    %only take derivative peaks where derivative is negative
    deriv_peaks(wcv_deriv(deriv_peaks) < 0) = [];


num_ups = 60;
up_dur = 2016;

up_start = [first_up];
down_start = [first_up+up_dur];
for i = 2:num_ups
    up_start(i) = up_start(end)+wait_dur+up_dur;
    down_start(i) = down_start(end)+wait_dur+up_dur;
end


%now find closest derivative peak to each up start and down start
for i = 1:num_ups
   [dummy,nearest_peak] = min(abs(deriv_peaks-up_start(i)));
   up_start(i) = deriv_peaks(nearest_peak);
   [dummy,nearest_peak] = min(abs(deriv_peaks-down_start(i)));
   down_start(i) = deriv_peaks(nearest_peak); 
end


maxlag = 3*2016;
lags = -maxlag:maxlag;

depol_mat = zeros(length(up_start),length(lags));
depol_smat = nan(length(up_start),length(lags));
up_lf8_phase = zeros(size(up_start));
depol_8mat = zeros(length(up_start),length(lags));
depol_3mat = zeros(length(up_start),length(lags));
for i = 1:length(up_start)
   depol_mat(i,:) = wcv_z(up_start(i)-maxlag:up_start(i)+maxlag);
   depol_8mat(i,:) = lf8_f(up_start(i)-maxlag:up_start(i)+maxlag);
   depol_3mat(i,:) = lf3_f(up_start(i)-maxlag:up_start(i)+maxlag);
   cur_spikes = find(spkid > up_start(i)-maxlag & spkid < up_start(i)+maxlag);
   cur_hist_spikes = hist(spkid(cur_spikes)-up_start(i),lags);
   spk_locs{i} = find(cur_hist_spikes > 0);
   depol_smat(i,:) = jmm_smooth_1d_cor(cur_hist_spikes,100);
   up_lf8_phase(i) = lf8_phase(up_start(i));
end

depol_smat = depol_smat*2016;

av_wcv = nanmean(depol_mat);
se_wcv = nanstd(depol_mat)/sqrt(60);
av_lf8 = nanmean(depol_8mat);
se_lf8 = nanstd(depol_8mat)/sqrt(60);
av_lf3 = nanmean(depol_3mat);
se_lf3 = nanstd(depol_3mat)/sqrt(60);

av_rate = nanmean(depol_smat);

[dummy,phase_order] = sort(up_lf8_phase);

imagesc(lags/2016,1:60,depol_mat(phase_order,:));shading flat
hold on
for i = 1:length(up_start)
   cur_loc = find(phase_order==i);
    plot(lags(spk_locs{i})/2016,ones(size(spk_locs{i}))*cur_loc,'k.')
    
end
xlim([-2 3])
t_names = ['C:\WC_Germany\current_injection\stellate\trig_mat\sorted_up_' f_names{d}];
print('-dpng',t_names)
close all

figure
imagesc(lags/2016,1:60,depol_mat);shading flat
hold on
for i = 1:length(up_start)
   
    plot(lags(spk_locs{i})/2016,ones(size(spk_locs{i}))*i,'k.')
    
end
xlim([-2 3])
t_names = ['C:\WC_Germany\current_injection\stellate\trig_mat\raw_up_' f_names{d}];
print('-dpng',t_names)
close all

%now plot depol triggered averages
figure
errorbar(lags/2016,av_wcv,se_wcv)
hold on
plot(lags/2016,av_lf8,'r')
plot(lags/2016,av_lf3,'k')
xlim([-1 2])
xlabel('Time (s)')
ylabel('Amplitude (z)')
t_names = ['C:\WC_Germany\current_injection\stellate\trig_avg\amp_' f_names{d}];
print('-dpng',t_names)
close all

figure
plot(lags/2016,av_rate,'k')
xlabel('Time (s)')
ylabel('Firing rate (Hz)')
xlim([-1 2])
t_names = ['C:\WC_Germany\current_injection\stellate\trig_avg\rate_' f_names{d}];
print('-dpng',t_names)
close all

end