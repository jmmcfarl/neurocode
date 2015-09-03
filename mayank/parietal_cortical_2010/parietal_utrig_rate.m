clear all
close all

addpath('G:\Code\smoothing\software')
addpath('G:\Code\FullBNT-1.0.4\KPMstats\')
addpath('G:\Code\FullBNT-1.0.4\netlab3.3')
addpath('G:\WC_Germany\new_stellate_analysis\')
addpath('G:\WC_Germany\hsmm_state_detection')
addpath('G:\WC_Germany\parietal_cortical_2010\')

cd G:\WC_Germany\parietal_cortical_2010\
load parietal_cortical_2010
load G:\WC_Germany\parietal_cortical_2010\desynch_times_mp_lf8

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');

sess_data = sess_data(parietal);
desynch_start_times = desynch_start_times(parietal);
desynch_stop_times = desynch_stop_times(parietal);

%get rid of interneurons
interneurons = find_struct_field_vals(sess_data,'cell_type','interneuron');
sess_data(interneurons) = [];
desynch_start_times(interneurons) = [];
desynch_stop_times(interneurons) = [];

frontal = find_struct_field_vals(sess_data,'region','frontal');
prefrontal = find_struct_field_vals(sess_data,'region','prefrontal');
parietal = find_struct_field_vals(sess_data,'region','parietal');
n = length(sess_data);

raw_Fs = 2016;
dsf = 8;
lcf = 0.05;
hcf = 40;
Fsd = raw_Fs/dsf;
niqf = raw_Fs/2;
backlag = round(0.5*Fsd);
forwardlag = round(2*Fsd);

lags = -backlag:forwardlag;
[b,a] = butter(2,[lcf/niqf hcf/niqf]);

rate_smooth = round(0.05*Fsd);

for d = 1:n
    to_dir = sess_data(d).directory;
    to_dir(1) = 'G';
    cd(to_dir)
    pwd
    
    load used_data wcv_minus_spike lf5 lf8
    load spike_time_jmm
    load hsmm_state_seq_lf_pert15
%     load hsmm_state_seq8_lf_pert15
        
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    wcv_f = zscore(downsample(wcv_f,dsf));
    lf8_f = filtfilt(b,a,lf8);
    lf8_f = zscore(downsample(lf8_f,dsf));
    lf5_f = filtfilt(b,a,lf5);
    lf5_f = zscore(downsample(lf5_f,dsf));

    %use adjusted lf signals
    hsmm_state_seq = hsmm_bbstate_seq;
    Fs = Fs_bb;
    
    spkid = round(spkid/dsf);
        
    up_trans = find(hsmm_state_seq(1:end-1) == 1 & hsmm_state_seq(2:end)==2);
    down_trans = find(hsmm_state_seq(1:end-1)==2 & hsmm_state_seq(2:end)==1);
    up_trans = up_trans*Fsd/Fs;
    down_trans = down_trans*Fsd/Fs;
    up_trans(up_trans < backlag | up_trans > length(wcv_f) - forwardlag) = [];
    down_trans(down_trans < backlag | down_trans > length(wcv_f) - forwardlag) = [];
    up_trans(up_trans > down_trans(end)) = [];
    down_trans(down_trans < up_trans(1)) = [];
    
%     %use down transitions
%     up_trans = down_trans;
%     up_trans8 = down_trans8;
    
    up_trig_mp = zeros(length(up_trans),length(lags));
    up_trig_lf8 = zeros(length(up_trans),length(lags));
     up_trig_lf5 = zeros(length(up_trans),length(lags));
   up_trig_rate = zeros(length(up_trans),length(lags));
    up_time = 0;
    up_spikes = 0;
    %% cycle through all MP up transitions
    for i = 1:length(up_trans)
        cur_seg = wcv_f(up_trans(i)-backlag:up_trans(i)+forwardlag);
        cur_lf8 = lf8_f(up_trans(i)-backlag:up_trans(i)+forwardlag);
        cur_lf5 = lf5_f(up_trans(i)-backlag:up_trans(i)+forwardlag);
        cur_spikes = spkid(spkid > up_trans(i)-backlag & spkid < up_trans(i) + forwardlag);
        cur_spkrate = hist(cur_spikes - up_trans(i),lags)*Fsd;
        if i > 1
            prev_down = up_trans(i) - down_trans(i-1);
            if prev_down < backlag
                cur_spkrate(1:prev_down) = nan;
                cur_seg(1:prev_down) = nan;
                 cur_lf8(1:prev_down) = nan;
                cur_lf5(1:prev_down) = nan;
           end
        end
        next_down = down_trans(i)-up_trans(i);
        if next_down < forwardlag
           cur_spkrate(backlag+next_down:end) = nan; 
           cur_seg(backlag+next_down:end) = nan;
           cur_lf8(backlag+next_down:end) = nan;
           cur_lf5(backlag+next_down:end) = nan;
        end
        up_trig_rate(i,:) = cur_spkrate;
        up_trig_mp(i,:) = cur_seg;
        up_trig_lf8(i,:) = cur_lf8;
        up_trig_lf5(i,:) = cur_lf5;
        up_time = up_time + (down_trans(i)-up_trans(i))/Fsd;
        up_spikes = up_spikes + length(find(spkid > up_trans(i) & spkid < down_trans(i)));
    end
    
    avg_utrig_mp(d,:) = nanmean(up_trig_mp);
    avg_utrig_lf8(d,:) = nanmean(up_trig_lf8);
    avg_utrig_lf5(d,:) = nanmean(up_trig_lf5);
    avg_utrig_rate(d,:) = jmm_smooth_1d_cor(nanmean(up_trig_rate),rate_smooth);
    avg_uprate(d) = up_spikes/up_time;
end

cd G:\WC_Germany\parietal_cortical_2010
% save parietal_up_rates avg_uprate
figure
errorbar(lags/Fsd,mean(avg_utrig_mp(1:11,:)),std(avg_utrig_mp(1:11,:))/sqrt(11))
hold on
errorbar(lags/Fsd,mean(avg_utrig_mp(12:end,:)),std(avg_utrig_mp(12:end,:))/sqrt(10),'r')
xlim([0 1])

figure
errorbar(lags/Fsd,mean(avg_utrig_lf8),std(avg_utrig_lf8)/sqrt(21))
hold on
errorbar(lags/Fsd,mean(avg_utrig_lf5),std(avg_utrig_lf5)/sqrt(21),'k')
xlim([0 1])

% hold on
% errorbar(lags/Fsd,mean(avg_utrig_mp(prefrontal,:)),std(avg_utrig_mp(prefrontal,:))/sqrt(length(prefrontal)),'r')
% errorbar(lags/Fsd,mean(avg_utrig_mp(frontal,:)),std(avg_utrig_mp(frontal,:))/sqrt(length(frontal)),'k')

figure
errorbar(lags/Fsd,mean(avg_utrig_rate(1:11,:)),std(avg_utrig_rate(1:11,:))/sqrt(11))
hold on
errorbar(lags/Fsd,mean(avg_utrig_rate(12:end,:)),std(avg_utrig_rate(12:end,:))/sqrt(10),'r')

% hold on
% errorbar(lags/Fsd,mean(avg_utrig_rate(prefrontal,:)),std(avg_utrig_rate(prefrontal,:))/sqrt(length(prefrontal)),'r')
% errorbar(lags/Fsd,mean(avg_utrig_rate(frontal,:)),std(avg_utrig_rate(frontal,:))/sqrt(length(frontal)),'k')

% figure
% errorbar(lags/Fsd,mean(avg_utrig_mp(superficial,:)),std(avg_utrig_mp(superficial,:))/sqrt(length(superficial)))
% hold on
% errorbar(lags/Fsd,mean(avg_utrig_mp(deep,:)),std(avg_utrig_mp(deep,:))/sqrt(length(deep)),'r')
% 
% figure
% errorbar(lags/Fsd,mean(avg_utrig_rate(superficial,:)),std(avg_utrig_rate(superficial,:))/sqrt(length(superficial)))
% hold on
% errorbar(lags/Fsd,mean(avg_utrig_rate(deep,:)),std(avg_utrig_rate(deep,:))/sqrt(length(deep)),'r')


