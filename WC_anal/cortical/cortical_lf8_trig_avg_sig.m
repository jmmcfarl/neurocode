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
    
    spkid = round(spkid/dsf);
    rlid_lup{d} = round(rlid_lup{d}/dsf);
    t_10_lup{d} = round(t_10_lup{d}/dsf);
    t_90_lup{d} = round(t_90_lup{d}/dsf);
    
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);
    lf5_f = filtfilt(b,a,lf5);
    
    down_w = downsample(wcv_f,dsf);
    down_8 = downsample(lf8_f,dsf);
%     down_5 = downsample(lf5_f,dsf);
  
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
%     down_5 = zscore(down_5);
   
    %initialize
    lf8_utrig_mp_mat_t10 = zeros(length(synch_ups{d}),length(lags));
    lf8_utrig_lf8_mat_t10 = zeros(length(synch_ups{d}),length(lags));
%     lf8_utrig_lf5_mat_t10 = zeros(length(synch_ups{d}),length(lags));
    lf8_utrig_spk_mat_t10 = zeros(length(synch_ups{d}),length(lags));
    
    lf8_utrig_mp_mat_t50 = zeros(length(synch_ups{d}),length(lags));
    lf8_utrig_lf8_mat_t50 = zeros(length(synch_ups{d}),length(lags));
%     lf8_utrig_lf5_mat_t50 = zeros(length(synch_ups{d}),length(lags));
    lf8_utrig_spk_mat_t50 = zeros(length(synch_ups{d}),length(lags));
    
    lf8_utrig_mp_mat_t90 = zeros(length(synch_ups{d}),length(lags));
    lf8_utrig_lf8_mat_t90 = zeros(length(synch_ups{d}),length(lags));
%     lf8_utrig_lf5_mat_t90 = zeros(length(synch_ups{d}),length(lags));
    lf8_utrig_spk_mat_t90 = zeros(length(synch_ups{d}),length(lags));

    
    %calculate lf8 utrigs
    for i = 1:length(synch_ups8{d})
        
        if t_10_lup{d}(i) > maxlag && ...
                length(down_w) - t_90_lup{d}(i) > maxlag && ~isnan(rlid_lup{d}(i))
           
            lf8_utrig_spk_mat_t50(i,:) = histc(spkid - rlid_lup{d}(i) ,lags);
            lf8_utrig_mp_mat_t50(i,:) = down_w(rlid_lup{d}(i)-maxlag:...
                rlid_lup{d}(i)+maxlag);
            lf8_utrig_lf8_mat_t50(i,:) = down_8(rlid_lup{d}(i)-maxlag:...
                rlid_lup{d}(i)+maxlag);
%             lf8_utrig_lf5_mat_t50(i,:) = down_5(rlid_lup{d}(i)-maxlag:...
%                 rlid_wup{d}(i)+maxlag);
            lf8_utrig_spk_mat_t10(i,:) = histc(spkid - t_10_lup{d}(i) ,lags);
            lf8_utrig_mp_mat_t10(i,:) = down_w(t_10_lup{d}(i)-maxlag:...
                t_10_lup{d}(i)+maxlag);
            lf8_utrig_lf8_mat_t10(i,:) = down_8(t_10_lup{d}(i)-maxlag:...
                t_10_lup{d}(i)+maxlag);
%             lf8_utrig_lf5_mat_t10(i,:) = down_5(t_10_lup{d}(i)-maxlag:...
%                 t_10_lup{d}(i)+maxlag);
            lf8_utrig_spk_mat_t90(i,:) = histc(spkid - t_90_lup{d}(i) ,lags);
            lf8_utrig_mp_mat_t90(i,:) = down_w(t_90_lup{d}(i)-maxlag:...
                t_90_lup{d}(i)+maxlag);
            lf8_utrig_lf8_mat_t90(i,:) = down_8(t_90_lup{d}(i)-maxlag:...
                t_90_lup{d}(i)+maxlag);
%             lf8_utrig_lf5_mat_t90(i,:) = down_5(t_90_wup{d}(i)-maxlag:...
%                 t_90_wup{d}(i)+maxlag);

        else
            
            lf8_utrig_spk_mat_t50(i,:) = nan;
            lf8_utrig_mp_mat_t50(i,:) = nan;
            lf8_utrig_lf8_mat_t50(i,:) = nan;
%             lf8_utrig_lf5_mat_t50(i,:) = nan;
            lf8_utrig_spk_mat_t10(i,:) = nan;
            lf8_utrig_mp_mat_t10(i,:) = nan;
            lf8_utrig_lf8_mat_t10(i,:) = nan;
%             lf8_utrig_lf5_mat_t10(i,:) = nan; 
            lf8_utrig_spk_mat_t90(i,:) = nan;
            lf8_utrig_mp_mat_t90(i,:) = nan;
            lf8_utrig_lf8_mat_t90(i,:) = nan;
%             lf8_utrig_lf5_mat_t90(i,:) = nan;
            
        end
        
    end
     
  bad_rows = logical(max(isnan(lf8_utrig_mp_mat_t10),[],2));
  lf8_utrig_mp_mat_t10(bad_rows,:) = [];
  lf8_utrig_lf8_mat_t10(bad_rows,:) = [];
%   lf8_utrig_lf5_mat_t10(bad_rows,:) = [];
  lf8_utrig_spk_mat_t10(bad_rows,:) = [];
  lf8_utrig_mp_mat_t90(bad_rows,:) = [];
  lf8_utrig_lf8_mat_t90(bad_rows,:) = [];
%   lf8_utrig_lf5_mat_t90(bad_rows,:) = [];
  lf8_utrig_spk_mat_t90(bad_rows,:) = [];
    lf8_utrig_mp_mat_t50(bad_rows,:) = [];
  lf8_utrig_lf8_mat_t50(bad_rows,:) = [];
%   lf8_utrig_lf5_mat_t50(bad_rows,:) = [];
  lf8_utrig_spk_mat_t50(bad_rows,:) = [];

  synch_ups8{d}(bad_rows) = [];
  t_10_lup{d}(bad_rows) = [];
  t_90_lup{d}(bad_rows) = [];
  rlid_lup{d}(bad_rows) = [];
    
%% calculate transition dependent distribution

% lf8_utrig_mp_dist_t10 = zeros(length(amp_range),length(lags));
% lf8_utrig_lf8_dist_t10 = zeros(length(amp_range),length(lags));
% % lf8_utrig_lf5_dist_t10 = zeros(length(amp_range),length(lags));
% lf8_utrig_mp_dist_t50 = zeros(length(amp_range),length(lags));
% lf8_utrig_lf8_dist_t50 = zeros(length(amp_range),length(lags));
% % lf8_utrig_lf5_dist_t50 = zeros(length(amp_range),length(lags));
% lf8_utrig_mp_dist_t90 = zeros(length(amp_range),length(lags));
% lf8_utrig_lf8_dist_t90 = zeros(length(amp_range),length(lags));
% % lf8_utrig_lf5_dist_t90 = zeros(length(amp_range),length(lags));
% 
% for i = 1:length(lags)
%    lf8_utrig_mp_dist_t10(:,i) = gpkde(lf8_utrig_mp_mat_t10(:,i),0.1,[nmin;nmax;nbins]); 
%    lf8_utrig_lf8_dist_t10(:,i) = gpkde(lf8_utrig_lf8_mat_t10(:,i),0.1,[nmin;nmax;nbins]); 
% %    lf8_utrig_lf5_dist_t10(:,i) = gpkde(lf8_utrig_lf5_mat_t10(:,i),0.1,[nmin;nmax;nbins]); 
%    lf8_utrig_mp_dist_t50(:,i) = gpkde(lf8_utrig_mp_mat_t50(:,i),0.1,[nmin;nmax;nbins]); 
%    lf8_utrig_lf8_dist_t50(:,i) = gpkde(lf8_utrig_lf8_mat_t50(:,i),0.1,[nmin;nmax;nbins]); 
% %    lf8_utrig_lf5_dist_t50(:,i) = gpkde(lf8_utrig_lf5_mat_t50(:,i),0.1,[nmin;nmax;nbins]); 
%    lf8_utrig_mp_dist_t90(:,i) = gpkde(lf8_utrig_mp_mat_t90(:,i),0.1,[nmin;nmax;nbins]); 
%    lf8_utrig_lf8_dist_t90(:,i) = gpkde(lf8_utrig_lf8_mat_t90(:,i),0.1,[nmin;nmax;nbins]); 
% %    lf8_utrig_lf5_dist_t90(:,i) = gpkde(lf8_utrig_lf5_mat_t90(:,i),0.1,[nmin;nmax;nbins]); 
% end
%   
% imagesc(lags/Fsd,amp_range,lf8_utrig_mp_dist_t50);shading flat
% caxis([0 1.5]);colorbar
% xlim([-0.5 0.5])
% line([0 0],[nmin nmax],'Color','w')
% cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_dist\t50_lf8_utrig_mp_' cell_name];
% print('-dpng',t_names);
% close
% imagesc(lags/Fsd,amp_range,lf8_utrig_lf8_dist_t50);shading flat
% caxis([0 1]);colorbar
% xlim([-0.5 0.5])
% line([0 0],[nmin nmax],'Color','w')
%         cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_dist\t50_lf8_utrig_lf8_' cell_name];
% print('-dpng',t_names);
% close
% % imagesc(lags/Fsd,amp_range,lf8_utrig_lf5_dist_t50);shading flat
% % caxis([0 1]);colorbar
% % xlim([-0.5 0.5])
% % line([0 0],[nmin nmax],'Color','w')
% % t_names = ['C:\WC_Germany\Cortical_analysis\trig_dist\t50_lf8_utrig_lf5_' sess_data(d).name];
% % print('-dpng',t_names);
% % close
% 
%    imagesc(lags/Fsd,amp_range,lf8_utrig_mp_dist_t10);shading flat
% caxis([0 1.5]);colorbar
% xlim([-0.5 0.5])
% line([0 0],[nmin nmax],'Color','w')
%         cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_dist\t10_lf8_utrig_mp_' cell_name];
% print('-dpng',t_names);
% close
% imagesc(lags/Fsd,amp_range,lf8_utrig_lf8_dist_t10);shading flat
% caxis([0 1]);colorbar
% xlim([-0.5 0.5])
% line([0 0],[nmin nmax],'Color','w')
%         cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_dist\t10_lf8_utrig_lf8_' cell_name];
% print('-dpng',t_names);
% close
% % imagesc(lags/Fsd,amp_range,lf8_utrig_lf5_dist_t10);shading flat
% % caxis([0 1]);colorbar
% % xlim([-0.5 0.5])
% % line([0 0],[nmin nmax],'Color','w')
% % t_names = ['C:\WC_Germany\Cortical_analysis\trig_dist\t10_lf8_utrig_lf5_' sess_data(d).name];
% % print('-dpng',t_names);
% % close
%  
% imagesc(lags/Fsd,amp_range,lf8_utrig_mp_dist_t90);shading flat
% caxis([0 1.5]);colorbar
% xlim([-0.5 0.5])
% line([0 0],[nmin nmax],'Color','w')
%         cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_dist\t90_lf8_utrig_mp_' cell_name];
% print('-dpng',t_names);
% close
% imagesc(lags/Fsd,amp_range,lf8_utrig_lf8_dist_t90);shading flat
% caxis([0 1]);colorbar
% xlim([-0.5 0.5])
% line([0 0],[nmin nmax],'Color','w')
%         cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_dist\t90_lf8_utrig_lf8_' cell_name];
% print('-dpng',t_names);
% close
% % imagesc(lags/Fsd,amp_range,lf8_utrig_lf5_dist_t90);shading flat
% % caxis([0 1]);colorbar
% % xlim([-0.5 0.5])
% % line([0 0],[nmin nmax],'Color','w')
% % t_names = ['C:\WC_Germany\Cortical_analysis\trig_dist\t90_lf8_utrig_lf5_' sess_data(d).name];
% % print('-dpng',t_names);
% % close


    lf8_utrig_mp_t50(d,:) = nanmean(lf8_utrig_mp_mat_t50);
    lf8_utrig_lf8_t50(d,:) = nanmean(lf8_utrig_lf8_mat_t50);
%     lf8_utrig_lf5_t50(d,:) = nanmean(lf8_utrig_lf5_mat_t50);
    lf8_utrig_spk_t50(d,:) = jmm_smooth_1d_cor(nanmean(lf8_utrig_spk_mat_t50),4)*Fsd;
      lf8_utrig_mp_t10(d,:) = nanmean(lf8_utrig_mp_mat_t10);
    lf8_utrig_lf8_t10(d,:) = nanmean(lf8_utrig_lf8_mat_t10);
%     lf8_utrig_lf5_t10(d,:) = nanmean(lf8_utrig_lf5_mat_t10);
    lf8_utrig_spk_t10(d,:) = jmm_smooth_1d_cor(nanmean(lf8_utrig_spk_mat_t10),4)*Fsd;
    lf8_utrig_mp_t90(d,:) = nanmean(lf8_utrig_mp_mat_t90);
    lf8_utrig_lf8_t90(d,:) = nanmean(lf8_utrig_lf8_mat_t90);
%     lf8_utrig_lf5_t90(d,:) = nanmean(lf8_utrig_lf5_mat_t90);
    lf8_utrig_spk_t90(d,:) = jmm_smooth_1d_cor(nanmean(lf8_utrig_spk_mat_t90),4)*Fsd;

    
%% plot mp trig averages    
%    plot(lags/Fsd,lf8_utrig_mp_t50(d,:),'linewidth',2)
%    hold on
%    plot(lags/Fsd,lf8_utrig_lf8_t50(d,:),'r','linewidth',2)
% %    plot(lags/Fsd,lf8_utrig_lf5_t50(d,:),'k','linewidth',2)
%    plot(lags/Fsd,lf8_utrig_spk_t50(d,:),'c','linewidth',2)
%    legend('MP','LF8','rate')
%    title('MP Up Triggered Avg')
%            cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\t50_lf8_utrig_' cell_name];
%    print('-dpng',t_names);
%    xlim([-0.5 0.5]); grid on
%               cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\t50_lf8_utrig_zoom_' cell_name];
%    print('-dpng',t_names);
%    close
%    
%    plot(lags/Fsd,lf8_utrig_mp_t10(d,:),'linewidth',2)
%    hold on
%    plot(lags/Fsd,lf8_utrig_lf8_t10(d,:),'r','linewidth',2)
% %    plot(lags/Fsd,lf8_utrig_lf5_t10(d,:),'k','linewidth',2)
%    plot(lags/Fsd,lf8_utrig_spk_t10(d,:),'c','linewidth',2)
%    legend('MP','LF8','rate')
%    title('MP Up Triggered Avg')
%            cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\t10_lf8_utrig_' cell_name];
%    print('-dpng',t_names);
%    xlim([-0.5 0.5]); grid on
%               cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\t10_lf8_utrig_zoom_' cell_name];
%    print('-dpng',t_names);
%    close
% 
%    plot(lags/Fsd,lf8_utrig_mp_t90(d,:),'linewidth',2)
%    hold on
%    plot(lags/Fsd,lf8_utrig_lf8_t90(d,:),'r','linewidth',2)
% %    plot(lags/Fsd,lf8_utrig_lf5_t90(d,:),'k','linewidth',2)
%    plot(lags/Fsd,lf8_utrig_spk_t90(d,:),'c','linewidth',2)
%    legend('MP','LF8','rate')
%    title('MP Up Triggered Avg')
%            cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\t90_lf8_utrig_' cell_name];
%    print('-dpng',t_names);
%    xlim([-0.5 0.5]); grid on
%               cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\t90_lf8_utrig_zoom_' cell_name];
%    print('-dpng',t_names);
%    close

end

clear *mat
save C:\WC_Germany\Cortical_analysis\trig_avgs\lf8_trig_avg_data_sigfit lags lf8* 