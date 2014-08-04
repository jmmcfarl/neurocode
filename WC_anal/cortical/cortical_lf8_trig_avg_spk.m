clear all
close all
cd C:\WC_Germany\Cortical_analysis

load C:\WC_Germany\Cortical_analysis\cortical_dir
load C:\WC_Germany\Cortical_analysis\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Cortical_analysis\UDS_synch_state_dur\UDS_synch_state_dur_data
load C:\WC_Germany\Cortical_analysis\sigmoid_fit\sig_fit_alldata

Fs = 2016;
niqf = Fs/2;
lcf = 0.05/niqf;
hcf = 40/niqf;
[b,a] = butter(2,[lcf hcf]);
dsf = 8;
Fsd = Fs/dsf;

maxlag = 3*Fsd;
lags = -maxlag:maxlag;


nbins = 500;
nmin = -3;
nmax = 3;
amp_range = linspace(nmin,nmax,nbins);


for d = 1:length(sess_data)
    
    
    cd(sess_data(d).directory)
    disp(num2str(d))
    
    load used_data lf8 lf5 wcv_minus_spike
    load spike_time_jmm
    
    down_w = downsample(wcv_minus_spike,dsf);
    
    spkid = round(spkid/dsf);
    rlid_lup{d} = round(rlid_lup{d}/dsf);
    t_10_lup{d} = round(t_10_lup{d}/dsf);
    t_90_lup{d} = round(t_90_lup{d}/dsf);
       
    %initialize
    lf8_utrig_spk_mat_t10 = zeros(length(synch_ups8{d}),length(lags));
    
    lf8_utrig_spk_mat_t50 = zeros(length(synch_ups8{d}),length(lags));
    
    lf8_utrig_spk_mat_t90 = zeros(length(synch_ups8{d}),length(lags));

    
    %calculate mp utrigs
    for i = 1:length(synch_ups8{d})
        
        if t_10_lup{d}(i) > maxlag && ...
                length(down_w) - t_90_lup{d}(i) > maxlag && ~isnan(rlid_lup{d}(i))
           
            lf8_utrig_spk_mat_t50(i,:) = histc(spkid - rlid_lup{d}(i) ,lags);
            lf8_utrig_spk_mat_t10(i,:) = histc(spkid - t_10_lup{d}(i) ,lags);
            lf8_utrig_spk_mat_t90(i,:) = histc(spkid - t_90_lup{d}(i) ,lags);
 

        else
            
            lf8_utrig_spk_mat_t50(i,:) = nan;
            lf8_utrig_spk_mat_t10(i,:) = nan;
            lf8_utrig_spk_mat_t90(i,:) = nan;
            
        end
        
    end
     
   lf8_utrig_spk_t50(d,:) = jmm_smooth_1d_cor(nanmean(lf8_utrig_spk_mat_t50),4)*Fsd;
    lf8_utrig_spk_t10(d,:) = jmm_smooth_1d_cor(nanmean(lf8_utrig_spk_mat_t10),4)*Fsd;
    lf8_utrig_spk_t90(d,:) = jmm_smooth_1d_cor(nanmean(lf8_utrig_spk_mat_t90),4)*Fsd;

    
%% plot mp trig averages    
%    plot(lags/Fsd,lf8_utrig_spk_t50(d,:),'linewidth',2)
%    hold on
%    plot(lags/Fsd,lf8_utrig_spk_t10(d,:),'r','linewidth',2)
%    plot(lags/Fsd,lf8_utrig_spk_t90(d,:),'k','linewidth',2)
%    
%    legend('50','10','90')
%    title('LF8 Up-triggered firing rate')
%               cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\lf8_utrig_rate_' cell_name];
%    xlim([-0.5 0.5]); grid on
%    print('-dpng',t_names);
%    close
   
end

clear *mat
save C:\WC_Germany\Cortical_analysis\trig_avgs\lf8_utrig_avg_spk lags lf8* 