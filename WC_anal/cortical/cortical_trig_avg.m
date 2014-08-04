clear all
close all
cd C:\WC_Germany\Cortical_analysis

load C:\WC_Germany\Cortical_analysis\cortical_dir
load C:\WC_Germany\Cortical_analysis\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Cortical_analysis\UDS_synch_state_dur\UDS_synch_state_dur_data

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
    
    wcv_f = filtfilt(b,a,wcv_minus_spike);
    lf8_f = filtfilt(b,a,lf8);
    lf5_f = filtfilt(b,a,lf5);
    
    down_w = downsample(wcv_f,dsf);
    down_8 = downsample(lf8_f,dsf);
    down_5 = downsample(lf5_f,dsf);
  
    down_w = zscore(down_w);
    down_8 = zscore(down_8);
    down_5 = zscore(down_5);
   
    %initialize
    mp_utrig_mp_mat = zeros(length(synch_ups{d}),length(lags));
    mp_utrig_lf8_mat = zeros(length(synch_ups{d}),length(lags));
    mp_utrig_lf5_mat = zeros(length(synch_ups{d}),length(lags));
    mp_utrig_spk_mat = zeros(length(synch_ups{d}),length(lags));
    
    mp_dtrig_mp_mat = zeros(length(synch_downs{d}),length(lags));
    mp_dtrig_lf8_mat = zeros(length(synch_downs{d}),length(lags));
    mp_dtrig_lf5_mat = zeros(length(synch_downs{d}),length(lags));
     mp_dtrig_spk_mat = zeros(length(synch_downs{d}),length(lags));

    
    %calculate mp utrigs
    for i = 1:length(synch_ups{d})
        
        if up_trans{d}(synch_ups{d}(i)) > maxlag && ...
                length(down_w) - up_trans{d}(synch_ups{d}(i)) > maxlag
           
            mp_utrig_spk_mat(i,:) = histc(spkid - up_trans{d}(synch_ups{d}(i)) ,lags);
            mp_utrig_mp_mat(i,:) = down_w(up_trans{d}(synch_ups{d}(i))-maxlag:...
                up_trans{d}(synch_ups{d}(i))+maxlag);
            mp_utrig_lf8_mat(i,:) = down_8(up_trans{d}(synch_ups{d}(i))-maxlag:...
                up_trans{d}(synch_ups{d}(i))+maxlag);
            mp_utrig_lf5_mat(i,:) = down_5(up_trans{d}(synch_ups{d}(i))-maxlag:...
                up_trans{d}(synch_ups{d}(i))+maxlag);

        else
            
            mp_utrig_spk_mat(i,:) = nan;
            mp_utrig_mp_mat(i,:) = nan;
            mp_utrig_lf8_mat(i,:) = nan;
              mp_utrig_lf5_mat(i,:) = nan;
          
        end
        
    end
    
       %calculate mp dtrigs
    for i = 1:length(synch_downs{d})
        
        if down_trans{d}(synch_downs{d}(i)) > maxlag && ...
                length(down_w) - down_trans{d}(synch_downs{d}(i)) > maxlag
            mp_dtrig_spk_mat(i,:) = histc(spkid - down_trans{d}(synch_downs{d}(i)) ,lags);

            mp_dtrig_mp_mat(i,:) = down_w(down_trans{d}(synch_downs{d}(i))-maxlag:...
                down_trans{d}(synch_downs{d}(i))+maxlag);
            mp_dtrig_lf8_mat(i,:) = down_8(down_trans{d}(synch_downs{d}(i))-maxlag:...
                down_trans{d}(synch_downs{d}(i))+maxlag);
            mp_dtrig_lf5_mat(i,:) = down_5(down_trans{d}(synch_downs{d}(i))-maxlag:...
                down_trans{d}(synch_downs{d}(i))+maxlag);
        else
            
            mp_dtrig_mp_mat(i,:) = nan;
            mp_dtrig_lf8_mat(i,:) = nan;
            mp_dtrig_lf5_mat(i,:) = nan;
            mp_dtrig_spk_mat(i,:) = nan;
        end
        
    end
 
%   bad_rows = logical(max(isnan(mp_dtrig_mp_mat),[],2));
%   mp_dtrig_mp_mat(bad_rows,:) = [];
%   mp_dtrig_lf8_mat(bad_rows,:) = [];
%   mp_dtrig_lf5_mat(bad_rows,:) = [];
%   mp_dtrig_spk_mat(bad_rows,:) = [];
%   
%   synch_ups{d}(bad_rows) = [];
  
%% calculate transition dependent distribution

% % mp_dtrig_mp_dist = zeros(length(amp_range),length(lags));
% % mp_dtrig_lf8_dist = zeros(length(amp_range),length(lags));
% % mp_dtrig_lf5_dist = zeros(length(amp_range),length(lags));
% % 
% % for i = 1:length(lags)
% %    mp_dtrig_mp_dist(:,i) = gpkde(mp_dtrig_mp_mat(:,i),0.1,[nmin;nmax;nbins]); 
% %    mp_dtrig_lf8_dist(:,i) = gpkde(mp_dtrig_lf8_mat(:,i),0.1,[nmin;nmax;nbins]); 
% %    mp_dtrig_lf5_dist(:,i) = gpkde(mp_dtrig_lf5_mat(:,i),0.1,[nmin;nmax;nbins]); 
% % end
% %   
% % pcolor(lags/Fsd,amp_range,mp_dtrig_mp_dist);shading flat
% % caxis([0 1.5]);colorbar
% % xlim([-0.5 0.5])
% % line([0 0],[nmin nmax],'Color','w')
% % cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% % t_names = ['C:\WC_Germany\Cortical_analysis\trig_dist\mp_dtrig_mp_' cell_name];
% % print('-dpng',t_names);
% % close
% % pcolor(lags/Fsd,amp_range,mp_dtrig_lf8_dist);shading flat
% % caxis([0 1]);colorbar
% % xlim([-0.5 0.5])
% % line([0 0],[nmin nmax],'Color','w')
% % cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% % t_names = ['C:\WC_Germany\Cortical_analysis\trig_dist\mp_dtrig_lf8_' cell_name];
% % print('-dpng',t_names);
% % close
% % pcolor(lags/Fsd,amp_range,mp_dtrig_lf5_dist);shading flat
% % caxis([0 1]);colorbar
% % xlim([-0.5 0.5])
% % line([0 0],[nmin nmax],'Color','w')
% % cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% % t_names = ['C:\WC_Germany\Cortical_analysis\trig_dist\mp_dtrig_lf5_' cell_name];
% % print('-dpng',t_names);
% % close

    
    mp_utrig_mp(d,:) = nanmean(mp_utrig_mp_mat);
    mp_utrig_lf8(d,:) = nanmean(mp_utrig_lf8_mat);
    mp_utrig_lf5(d,:) = nanmean(mp_utrig_lf5_mat);

    mp_utrig_spk(d,:) = jmm_smooth_1d_cor(nanmean(mp_utrig_spk_mat),4)*Fsd;

    mp_dtrig_mp(d,:) = nanmean(mp_dtrig_mp_mat);
    mp_dtrig_lf8(d,:) = nanmean(mp_dtrig_lf8_mat);
      mp_dtrig_lf5(d,:) = nanmean(mp_dtrig_lf5_mat);
      mp_dtrig_spk(d,:) = jmm_smooth_1d_cor(nanmean(mp_dtrig_spk_mat),4)*Fsd;

    [dummy,up_order] = sort(up_state_dur{d}(synch_ups{d}));
%     [dummy,down_order] = sort(down_state_dur{d}(synch_downs{d}(1:end-1)));
%     [dummy,up_order8] = sort(up_state_dur8{d}(synch_ups8{d}));
%     [dummy,down_order8] = sort(down_state_dur8{d}(synch_downs8{d}));
    
% %% plot mp trig averages    
%    plot(lags/Fsd,mp_utrig_mp(d,:),'linewidth',2)
%    hold on
%    plot(lags/Fsd,mp_utrig_lf8(d,:),'r','linewidth',2)
%    plot(lags/Fsd,mp_utrig_lf5(d,:),'k','linewidth',2)
%    plot(lags/Fsd,mp_utrig_spk(d,:),'c','linewidth',2)
%    legend('MP','LF8','LF5','rate')
%    title('MP Up Triggered Avg')
%          cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\mp_utrig_' cell_name];
%    print('-dpng',t_names);
%    xlim([-1 1]); grid on
%             cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\mp_utrig_zoom_' cell_name];
%    print('-dpng',t_names);
%    close
%    
%    plot(lags/Fsd,mp_dtrig_mp(d,:),'linewidth',2)
%    hold on
%    plot(lags/Fsd,mp_dtrig_lf8(d,:),'r','linewidth',2)
%    plot(lags/Fsd,mp_dtrig_lf5(d,:),'k','linewidth',2)
%    plot(lags/Fsd,mp_dtrig_spk(d,:),'c','linewidth',2)
%    legend('MP','LF8','LF5','rate')
%    title('MP Down Triggered Avg')
%          cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\mp_dtrig_' cell_name];
%    print('-dpng',t_names);
%       xlim([-1 1]); grid on
%       cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\mp_dtrig_zoom_' cell_name];
%    print('-dpng',t_names);
%    close


%% plot mp up trig matrices
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     imagesc(lags/Fsd,(1:length(synch_ups{d}))/length(synch_ups{d}) ...
%         ,mp_utrig_mp_mat(up_order,:));shading flat; colorbar; hold on
%     plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))/length(synch_ups{d}),'w','linewidth',2)
%     caxis([-3 3]);colorbar
%     xlim([-2 3])
%     line([0 0],[0 1],'Color','k')
%     num_trans = size(mp_utrig_spk_mat,1);
%     for i = 1:num_trans
%        cur_spikes = find((spkid > up_trans{d}(synch_ups{d}(up_order(i)))-maxlag) & (spkid < up_trans{d}(synch_ups{d}(up_order(i)))+maxlag));
%        cur_spikes = spkid(cur_spikes) - up_trans{d}(synch_ups{d}(up_order(i)));
%        plot(cur_spikes/Fsd,ones(size(cur_spikes))*i/num_trans,'k.')
%     end
%     
%         cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\mp_utrig_mp_' cell_name];
%     print('-dpng',t_names);
%     close
%    
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     imagesc(lags/Fsd,(1:length(synch_ups{d}))/length(synch_ups{d}) ...
%         ,mp_utrig_lf8_mat(up_order,:));shading flat; colorbar; hold on
%     plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))/length(synch_ups{d}),'w','linewidth',2)
%         xlim([-2 3])
%     caxis([-3 3]);colorbar
%     line([0 0],[0 1],'Color','k')
%     num_trans = size(mp_utrig_spk_mat,1);
%     for i = 1:num_trans
%        cur_spikes = find((spkid > up_trans{d}(synch_ups{d}(up_order(i)))-maxlag) & (spkid < up_trans{d}(synch_ups{d}(up_order(i)))+maxlag));
%        cur_spikes = spkid(cur_spikes) - up_trans{d}(synch_ups{d}(up_order(i)));
%        plot(cur_spikes/Fsd,ones(size(cur_spikes))*i/num_trans,'k.')
%     end
%        cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
%  t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\mp_utrig_lf8_' cell_name];
%     print('-dpng',t_names);
%     close
% 
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     imagesc(lags/Fsd,(1:length(synch_ups{d}))/length(synch_ups{d}) ...
%         ,mp_utrig_lf5_mat(up_order,:));shading flat; colorbar; hold on
%     plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))/length(synch_ups{d}),'w','linewidth',2)
%         xlim([-2 3])
%     caxis([-3 3]);colorbar
%     line([0 0],[0 1],'Color','k')
%     num_trans = size(mp_utrig_spk_mat,1);
%     for i = 1:num_trans
%        cur_spikes = find((spkid > up_trans{d}(synch_ups{d}(up_order(i)))-maxlag) & (spkid < up_trans{d}(synch_ups{d}(up_order(i)))+maxlag));
%        cur_spikes = spkid(cur_spikes) - up_trans{d}(synch_ups{d}(up_order(i)));
%        plot(cur_spikes/Fsd,ones(size(cur_spikes))*i/num_trans,'k.')
%     end
%        cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
%  t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\mp_utrig_lf5_' cell_name];
%     print('-dpng',t_names);
%     close

     
%% plot mp down trig matrices
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     imagesc(lags/Fsd,(1:length(synch_downs{d}))/length(synch_downs{d}) ...
%         ,mp_dtrig_mp_mat);shading flat; colorbar; hold on
% %     plot(down_state_dur{d}(synch_downs{d}(down_order)),(1:length(synch_downs{d}))/length(synch_downs{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%         xlim([-2 3])
%     caxis([-3 3]);colorbar
%     num_trans = size(mp_dtrig_spk_mat,1);
%     for i = 1:num_trans
%        cur_spikes = find((spkid > down_trans{d}(synch_downs{d}(i))-maxlag) & (spkid < down_trans{d}(synch_downs{d}(i))+maxlag));
%        cur_spikes = spkid(cur_spikes) - down_trans{d}(synch_downs{d}(i));
%        plot(cur_spikes/Fsd,ones(size(cur_spikes))*i/num_trans,'k.')
%     end
% 
%     cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\mp_dtrig_mp_' cell_name];
%     print('-dpng',t_names);
%     close
% %    
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     imagesc(lags/Fsd,(1:length(synch_downs{d}))/length(synch_downs{d}) ...
%         ,mp_dtrig_lf8_mat);shading flat; colorbar; hold on
% %     plot(down_state_dur{d}(synch_downs{d}(down_order)),(1:length(synch_downs{d}))/length(synch_downs{d}),'w','linewidth',2)
%         xlim([-2 3])
%     caxis([-3 3]);colorbar
%     line([0 0],[0 1],'Color','k')
%         num_trans = size(mp_dtrig_spk_mat,1);
%     for i = 1:num_trans
%        cur_spikes = find((spkid > down_trans{d}(synch_downs{d}(i))-maxlag) & (spkid < down_trans{d}(synch_downs{d}(i))+maxlag));
%        cur_spikes = spkid(cur_spikes) - down_trans{d}(synch_downs{d}(i));
%        plot(cur_spikes/Fsd,ones(size(cur_spikes))*i/num_trans,'k.')
%     end
% 
%     cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
% t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\mp_dtrig_lf8_' cell_name];
%     print('-dpng',t_names);
%     close
%     
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     imagesc(lags/Fsd,(1:length(synch_downs{d}))/length(synch_downs{d}) ...
%         ,mp_dtrig_lf5_mat);shading flat; colorbar; hold on
% %     plot(down_state_dur{d}(synch_downs{d}(down_order)),(1:length(synch_downs{d}))/length(synch_downs{d}),'w','linewidth',2)
%         xlim([-2 3])
%     caxis([-3 3]);colorbar
%     line([0 0],[0 1],'Color','k')
%         num_trans = size(mp_dtrig_spk_mat,1);
%     for i = 1:num_trans
%        cur_spikes = find((spkid > down_trans{d}(synch_downs{d}(i))-maxlag) & (spkid < down_trans{d}(synch_downs{d}(i))+maxlag));
%        cur_spikes = spkid(cur_spikes) - down_trans{d}(synch_downs{d}(i));
%        plot(cur_spikes/Fsd,ones(size(cur_spikes))*i/num_trans,'k.')
%     end
% 
%    cell_name = ['L' sess_data(d).layer '_' sess_data(d).cell_type '_' sess_data(d).region '_' sess_data(d).name];
%  t_names = ['C:\WC_Germany\Cortical_analysis\trig_avgs\mp_dtrig_lf5_' cell_name];
%     print('-dpng',t_names);
%     close

end

clear *mat
save C:\WC_Germany\Cortical_analysis\trig_avgs\trig_avg_data_wideband lags mp* lf8* lf5*