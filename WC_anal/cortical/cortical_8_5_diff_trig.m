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

smooth_win = 10;

for d = 1:length(sess_data)
    
    rlid_lup{d} = round(rlid_lup{d}/dsf);

    cd(sess_data(d).directory)
    disp(num2str(d))
    
    load used_data lf8 lf7 lf5 wcv_minus_spike
    load spike_time_jmm
    
    spkid = round(spkid/dsf);
    
%     wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);
    lf5_f = filtfilt(b,a,lf5);
%     lf7_f = filtfilt(b,a,lf7);
%     down_w = downsample(wcv_f,dsf);
    down_8 = downsample(lf8_f,dsf);
    down_5 = downsample(lf5_f,dsf);
%     down_7 = downsample(lf7_f,dsf);
  
%     down_w = zscore(down_w);
    down_8 = zscore(down_8);
    down_5 = zscore(down_5);
%     down_8_5 = zscore(down_8-down_5);
%     down_7 = zscore(down_7);
    down_8_5 = zscore(down_8-down_5);
    
    %initialize
%     mp_utrig_mp_mat = zeros(length(synch_ups{d}),length(lags));
%     mp_utrig_lf8_mat = zeros(length(synch_ups{d}),length(lags));
%     mp_utrig_lf5_mat = zeros(length(synch_ups{d}),length(lags));
    lf8_utrig_spk_mat = zeros(length(synch_ups8{d}),length(lags));
%     mp_utrig_lf8_5_mat = zeros(length(synch_ups{d}),length(lags));
    lf8_utrig_lf8_5_mat = zeros(length(synch_ups8{d}),length(lags));
    
%     mp_dtrig_mp_mat = zeros(length(synch_downs{d}),length(lags));
%     mp_dtrig_lf8_mat = zeros(length(synch_downs{d}),length(lags));
%     mp_dtrig_lf5_mat = zeros(length(synch_downs{d}),length(lags));
    lf8_dtrig_lf8_5_mat = zeros(length(synch_ups8{d}),length(lags));
     lf8_dtrig_spk_mat = zeros(length(synch_downs8{d}),length(lags));
%     mp_dtrig_lf7_5_mat = zeros(length(synch_ups{d}),length(lags));
    
    %calculate mp utrigs
    for i = 1:length(synch_ups8{d})
        
        if rlid_lup{d}(i) > maxlag && ...
                length(down_8_5) - rlid_lup{d}(i) > maxlag && ~isnan(rlid_lup{d}(i))
           
            lf8_utrig_spk_mat(i,:) = jmm_smooth_1d_cor(histc(spkid - rlid_lup{d}(i) ,lags),smooth_win);
%             mp_utrig_mp_mat(i,:) = down_w(up_trans{d}(synch_ups{d}(i))-maxlag:...
%                 up_trans{d}(synch_ups{d}(i))+maxlag);
%             mp_utrig_lf8_mat(i,:) = down_8(up_trans{d}(synch_ups{d}(i))-maxlag:...
%                 up_trans{d}(synch_ups{d}(i))+maxlag);
%             mp_utrig_lf5_mat(i,:) = down_5(up_trans{d}(synch_ups{d}(i))-maxlag:...
%                 up_trans{d}(synch_ups{d}(i))+maxlag);
            lf8_utrig_lf8_5_mat(i,:) = down_8_5(rlid_lup{d}(i)-maxlag:...
                rlid_lup{d}(i)+maxlag);

        else
            
            lf8_utrig_spk_mat(i,:) = nan;
%             mp_utrig_mp_mat(i,:) = nan;
%             mp_utrig_lf8_mat(i,:) = nan;
%               mp_utrig_lf5_mat(i,:) = nan;
             lf8_utrig_lf8_5_mat(i,:) = nan;
         
        end
        
    end
    
       %calculate mp dtrigs
    for i = 1:length(synch_downs8{d})
        
        if down_trans8{d}(synch_downs8{d}(i)) > maxlag && ...
                length(down_8_5) - down_trans8{d}(synch_downs8{d}(i)) > maxlag
            lf8_dtrig_spk_mat(i,:) = jmm_smooth_1d_cor(histc(spkid - down_trans8{d}(synch_downs8{d}(i)) ,lags), smooth_win);

%             mp_dtrig_mp_mat(i,:) = down_w(down_trans{d}(synch_downs{d}(i))-maxlag:...
%                 down_trans{d}(synch_downs{d}(i))+maxlag);
%             mp_dtrig_lf8_mat(i,:) = down_8(down_trans{d}(synch_downs{d}(i))-maxlag:...
%                 down_trans{d}(synch_downs{d}(i))+maxlag);
%             mp_dtrig_lf5_mat(i,:) = down_5(down_trans{d}(synch_downs{d}(i))-maxlag:...
%                 down_trans{d}(synch_downs{d}(i))+maxlag);
            lf8_dtrig_lf8_5_mat(i,:) = down_8_5(down_trans8{d}(synch_downs8{d}(i))-maxlag:...
                down_trans8{d}(synch_downs8{d}(i))+maxlag);

        else
            
%             mp_dtrig_mp_mat(i,:) = nan;
%             mp_dtrig_lf8_mat(i,:) = nan;
%             mp_dtrig_lf5_mat(i,:) = nan;
%             mp_dtrig_spk_mat(i,:) = nan;
                        lf8_dtrig_lf8_5_mat(i,:) = nan;

        end
        
    end
 
%   bad_rows = logical(max(isnan(mp_dtrig_mp_mat),[],2));
%   mp_dtrig_mp_mat(bad_rows,:) = [];
%   mp_dtrig_lf8_mat(bad_rows,:) = [];
%   mp_dtrig_lf5_mat(bad_rows,:) = [];
%   mp_dtrig_spk_mat(bad_rows,:) = [];
  
%   synch_ups{d}(bad_rows) = [];
  
lf8_utrig_lf8_5(d,:) = nanmean(lf8_utrig_lf8_5_mat);
lf8_dtrig_lf8_5(d,:) = nanmean(lf8_dtrig_lf8_5_mat);

%% plot mp trig averages    
   plot(lags/Fsd,lf8_utrig_lf8_5(d,:),'linewidth',2)
   title('MP Up Triggered Avg')
         cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\lf8_utrig_zoom_8_5_' cell_name];
   xlim([-1 1]); grid on
   print('-dpng',t_names);
   close
   
   plot(lags/Fsd,lf8_dtrig_lf8_5(d,:),'linewidth',2)
   title('MP Down Triggered Avg')
      xlim([-1 1]); grid on
      cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\lf8_dtrig_zoom_8_5_' cell_name];
   print('-dpng',t_names);
   close

[dummy, up_order8] = sort(synch_up_dur8{d});
%% plot mp up trig matrices
    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    imagesc(lags/Fsd,(1:length(synch_ups8{d}))/length(synch_ups8{d}) ...
        ,lf8_utrig_lf8_5_mat(up_order8,:));shading flat; colorbar; hold on
    plot(up_state_dur8{d}(synch_ups8{d}(up_order8)),(1:length(synch_ups8{d}))/length(synch_ups8{d}),'w','linewidth',2)
    caxis([-1.5 1.5]);colorbar
    xlim([-2 2])
    line([0 0],[0 1],'Color','k')
%     num_trans = size(mp_utrig_spk_mat,1);
          cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\lf8_utrig_7_5_' cell_name];
   print('-dpng',t_names);
   close

    Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    imagesc(lags/Fsd,(1:length(synch_ups8{d}))/length(synch_ups8{d}) ...
        ,lf8_dtrig_lf8_5_mat);shading flat; colorbar; hold on
%     plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))/length(synch_ups{d}),'w','linewidth',2)
    caxis([-1.5 1.5]);colorbar
    xlim([-2 2])
    line([0 0],[0 1],'Color','k')
%     num_trans = size(mp_utrig_spk_mat,1);
          cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\lf8_dtrig_8_5_' cell_name];
   print('-dpng',t_names);
   close

   
       Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    imagesc(lags/Fsd,(1:length(synch_ups8{d}))/length(synch_ups8{d}) ...
        ,lf8_utrig_spk_mat(up_order8,:));shading flat; colorbar; hold on
    plot(up_state_dur8{d}(synch_ups8{d}(up_order8)),(1:length(synch_ups8{d}))/length(synch_ups8{d}),'w','linewidth',2)
    colorbar
    xlim([-1 1])
    line([0 0],[0 1],'Color','w')
%     num_trans = size(mp_utrig_spk_mat,1);
          cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\lf8_utrig_spk_' cell_name];
   print('-dpng',t_names);
   close

        Fig = figure(1)
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    imagesc(lags/Fsd,(1:length(synch_ups8{d}))/length(synch_ups8{d}) ...
        ,lf8_dtrig_spk_mat);shading flat; colorbar; hold on
%     plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))/length(synch_ups{d}),'w','linewidth',2)
    colorbar
    xlim([-1 1])
    line([0 0],[0 1],'Color','w')
%     num_trans = size(mp_utrig_spk_mat,1);
          cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\lf8_dtrig_spk_' cell_name];
   print('-dpng',t_names);
   close
  
end

save C:\WC_Germany\Cortical_analysis\trig_avgs\lf8_trig_lf8_5_data lf8_*trig* lags Fsd