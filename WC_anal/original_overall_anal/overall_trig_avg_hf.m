clear all
close all

load('C:\WC_Germany\overall_calcs\overall_dir.mat')
load C:\WC_Germany\overall_calcs\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\overall_calcs\UDS_synch_state_dur\UDS_synch_state_dur_data
load C:\WC_Germany\overall_calcs\overall_info_file_data

Fs = 2016;
niqf = Fs/2;
lcf = 5/niqf;
hcf = 80/niqf;
[b,a] = butter(1,[lcf hcf]);
lcf2 = 0.05/niqf;
hcf2 = 10/niqf;
[b2,a2] = butter(1,[lcf2 hcf2]);
dsf = 8;
Fsd = Fs/dsf;

backlag = 5*Fsd;
forwardlag = 5*Fsd;
lags = (-backlag:forwardlag)/Fsd;

for d = 1:length(over_dir)
    
    cd(over_dir{d})
    disp(num2str(d))
    
    load used_data lf8 lf2 lf3 wcv_minus_spike
    wcv_f = filter(b,a,wcv_minus_spike);
    lf8_f = filter(b,a,lf8);
%     lf2_f = filter(b,a,lf2);
    lf3_f = filter(b,a,lf3);
   down_w = downsample(wcv_f,dsf);
   down_8 = downsample(lf8_f,dsf);
   down_3 = downsample(lf3_f,dsf);
   
   down_w = down_w.^2;
   down_8 = down_8.^2;
%     down_w(down_w < 0) = 0;
%     down_8(down_8 < 0) = 0;
%     down_3(down_3 < 0) = 0;
    
    
    wcv_fl = filter(b2,a2,wcv_minus_spike);
    lf8_fl = filter(b2,a2,lf8);
%     lf2_fl = filter(b2,a2,lf2);
    lf3_fl = filter(b2,a2,lf3);
    
%     wcv_f = [0;diff(wcv_f)];
%     lf8_f = [0;diff(lf8_f)];
%     lf3_f = [0;diff(lf3_f)];
%     
%     wcv_f = wcv_f.^2;
%     lf8_f = lf8_f.^2;
%     lf3_f = lf3_f.^2;
    
%     down_w = downsample(wcv_f,dsf);
%     down_8 = downsample(lf8_f,dsf);
%     down_3 = downsample(lf3_f,dsf);
%     
%     down_w = log(jmm_smooth_1d(down_w,4));
%     down_3 = log(jmm_smooth_1d(down_3,4));
%     down_8 = log(jmm_smooth_1d(down_8,4));
    
%     down_w = zscore(down_w);
%     down_8 = zscore(down_8);
%     down_3 = zscore(down_3);
% down_w(down_w > 1) = 1;
% down_8(down_8 > 1) = 1;
% down_3(down_3 > 1) = 1;
    
    
    down_wl = downsample(wcv_fl,dsf);
    down_8l = downsample(lf8_fl,dsf);
%     down_2l = downsample(lf2_fl,dsf);
    down_3l = downsample(lf3_fl,dsf);
     
%     down_wl = down_wl-min(down_wl)+0.01;
%     down_8l = down_8l-min(down_8l)+0.01;
%     down_3l = down_3l-min(down_3l)+0.01;
    
%     down_wl = log(down_wl);
%     down_8l = log(down_8l);
%     down_3l = log(down_3l);
    
        down_wl = zscore(down_wl);
    down_8l = zscore(down_8l);
    down_3l = zscore(down_3l);

%     down_2l = zscore(down_2l);

%     temp_cor = corrcoef(down_3,down_8);
%     hc_cor = temp_cor(2,1);
%     down3_sub = down_3-hc_cor*down_8;
    
%% find LFP up transitions which occur during an MP down state
used_up8 = [];
corres_down = [];
for i = 1:length(synch_ups8{d})
    prev_down = find(down_trans{d} < up_trans8{d}(synch_ups8{d}(i)),1,'last');
    prev_up = find(up_trans{d} < up_trans8{d}(synch_ups8{d}(i)),1,'last');
    if ~isempty(prev_down) & ~isempty(prev_up)
        if prev_down >= prev_up
            used_up8 = [used_up8 synch_ups8{d}(i)];
            corres_down = [corres_down prev_down];
        end
    end
end


% 
% synch_ups8{d} = used_up8;
% used_down8 = [];
% for i = 1:length(synch_downs8{d})
%     prev_up = find(up_trans{d} < down_trans8{d}(synch_downs8{d}(i)),1,'last');
%     prev_down = find(down_trans{d} < down_trans8{d}(synch_downs8{d}(i)),1,'last');
%     if ~isempty(prev_down) & ~isempty(prev_down)
%         if prev_up > prev_down
%             used_down8 = [used_down8 synch_downs8{d}(i)];
%         end
%     end
% end

synch_ups8{d} = used_up8;


    %initialize
    mp_utrig_mp_mat = zeros(length(synch_ups{d}),length(lags));
    mp_utrig_lf8_mat = zeros(length(synch_ups{d}),length(lags));
    mp_utrig_lf3_mat = zeros(length(synch_ups{d}),length(lags));
%     mp_utrig_lf2_mat = zeros(length(synch_ups{d}),length(lags));

    mp_utrig_mpl_mat = zeros(length(synch_ups{d}),length(lags));
    mp_utrig_lf8l_mat = zeros(length(synch_ups{d}),length(lags));
    mp_utrig_lf3l_mat = zeros(length(synch_ups{d}),length(lags));
%     mp_utrig_lf2l_mat = zeros(length(synch_ups{d}),length(lags));

lf8_dtrig_mp_mat = zeros(length(synch_downs8{d}),length(lags));
lf8_dtrig_lf8_mat = zeros(length(synch_downs8{d}),length(lags));
lf8_dtrig_mpl_mat = zeros(length(synch_downs8{d}),length(lags));
lf8_dtrig_lf8l_mat = zeros(length(synch_downs8{d}),length(lags));

    lf8_utrig_mp_mat = zeros(length(synch_ups8{d}),length(lags));
    lf8_utrig_lf8_mat = zeros(length(synch_ups8{d}),length(lags));
    lf8_utrig_lf3_mat = zeros(length(synch_ups8{d}),length(lags));
%     lf8_utrig_lf2_mat = zeros(length(synch_ups8{d}),length(lags));
    
    lf8_utrig_mpl_mat = zeros(length(synch_ups8{d}),length(lags));
    lf8_utrig_lf8l_mat = zeros(length(synch_ups8{d}),length(lags));
    lf8_utrig_lf3l_mat = zeros(length(synch_ups8{d}),length(lags));
%     lf8_utrig_lf2l_mat = zeros(length(synch_ups8{d}),length(lags));    
    
down_w_uo = down_w;
down_8_uo = down_8;

% set filtered signals to nans during MP up states
for i = 1:length(synch_ups{d})
    
    down_w(up_trans{d}(synch_ups{d}(i))-10:down_trans{d}(synch_ups{d}(i))) = nan;
        down_8(up_trans{d}(synch_ups{d}(i))-10:down_trans{d}(synch_ups{d}(i))) = nan;
%     down_3(up_trans{d}(synch_ups{d}(i))-10:down_trans{d}(synch_ups{d}(i))) = nan;

end
% set filtered signals to nans during MP down states
for i = 1:(length(synch_ups{d})-1)
    
    down_w_uo(down_trans{d}(synch_ups{d}(i))-10:up_trans{d}(synch_ups{d}(i)+1)) = nan;
        down_8_uo(down_trans{d}(synch_ups{d}(i))-10:up_trans{d}(synch_ups{d}(i)+1)) = nan;
%     down_3(up_trans{d}(synch_ups{d}(i))-10:down_trans{d}(synch_ups{d}(i))) = nan;

end

%zscore down state signals
down_w = (down_w-nanmean(down_w))/nanstd(down_w);
down_8 = (down_8-nanmean(down_8))/nanstd(down_8);
down_3 = (down_3-nanmean(down_3))/nanstd(down_3);

%zscore down state signals
down_w_uo = (down_w_uo-nanmean(down_w_uo))/nanstd(down_w_uo);
down_8_uo = (down_8_uo-nanmean(down_8_uo))/nanstd(down_8_uo);
% down_3 = (down_3-nanmean(down_3))/nanstd(down_3);



    %calculate lf8 dtrigs
    for i = 1:length(synch_downs8{d})
        
        if down_trans8{d}(synch_downs8{d}(i)) > backlag && ...
                length(down_w) - down_trans8{d}(synch_downs8{d}(i)) > forwardlag
            
            lf8_dtrig_mp_mat(i,:) = down_w_uo(down_trans8{d}(synch_downs8{d}(i))-backlag:...
                down_trans8{d}(synch_downs8{d}(i))+forwardlag);
            lf8_dtrig_lf8_mat(i,:) = down_8_uo(down_trans8{d}(synch_downs8{d}(i))-backlag:...
                down_trans8{d}(synch_downs8{d}(i))+forwardlag);
%             lf8_utrig_lf3_mat(i,:) = down_3(down_trans{d}(synch_downs{d}(i))-backlag:...
%                 down_trans{d}(synch_downs{d}(i))+forwardlag);
%             mp_utrig_lf2_mat(i,:) = down_2(up_trans{d}(synch_ups{d}(i))-backlag:...
%                 up_trans{d}(synch_ups{d}(i))+forwardlag);

            lf8_dtrig_mpl_mat(i,:) = down_wl(down_trans8{d}(synch_downs8{d}(i))-backlag:...
                down_trans8{d}(synch_downs8{d}(i))+forwardlag);
            lf8_dtrig_lf8l_mat(i,:) = down_8l(down_trans8{d}(synch_downs8{d}(i))-backlag:...
                down_trans8{d}(synch_downs8{d}(i))+forwardlag);
%             lf8_utrig_lf3l_mat(i,:) = down_3l(down_trans{d}(synch_downs{d}(i))-backlag:...
%                 down_trans{d}(synch_downs{d}(i))+forwardlag);
%             mp_utrig_lf2l_mat(i,:) = down_2l(up_trans{d}(synch_ups{d}(i))-backlag:...
%                 up_trans{d}(synch_ups{d}(i))+forwardlag);

        else
            
            lf8_dtrig_mp_mat(i,:) = nan;
            lf8_dtrig_lf8_mat(i,:) = nan;
%             mp_utrig_lf3_mat(i,:) = nan;
%             mp_utrig_lf2_mat(i,:) = nan;
% 
%             lf8_dtrig_mpl_mat(i,:) = nan;
%             lf8_dtrig_lf8l_mat(i,:) = nan;
%             mp_utrig_lf3l_mat(i,:) = nan;
%             mp_utrig_lf2l_mat(i,:) = nan;

        end
        
    end
    
 
      %calculate lfp utrigs
    for i = 1:length(synch_ups8{d})
        
        if up_trans8{d}(synch_ups8{d}(i)) > backlag && ...
                length(down_w) - up_trans8{d}(synch_ups8{d}(i)) > forwardlag
            %
            prev_mp_down = find(down_trans{d} < up_trans8{d}(synch_ups8{d}(i)),1,'last');
            back_dist = up_trans8{d}(synch_ups8{d}(i))-down_trans{d}(prev_mp_down);
            lf8_utrig_mp_mat(i,:) = down_w(up_trans8{d}(synch_ups8{d}(i))-backlag:...
                up_trans8{d}(synch_ups8{d}(i))+forwardlag);
            lf8_utrig_lf8_mat(i,:) = down_8(up_trans8{d}(synch_ups8{d}(i))-backlag:...
                up_trans8{d}(synch_ups8{d}(i))+forwardlag);
            lf8_utrig_lf3_mat(i,:) = down_3(up_trans8{d}(synch_ups8{d}(i))-backlag:...
                up_trans8{d}(synch_ups8{d}(i))+forwardlag);
%             lf8_utrig_lf2_mat(i,:) = down_2(up_trans8{d}(synch_ups8{d}(i))-backlag:...
%                 up_trans8{d}(synch_ups8{d}(i))+forwardlag);

            lf8_utrig_mpl_mat(i,:) = down_wl(up_trans8{d}(synch_ups8{d}(i))-backlag:...
                up_trans8{d}(synch_ups8{d}(i))+forwardlag);
            lf8_utrig_lf8l_mat(i,:) = down_8l(up_trans8{d}(synch_ups8{d}(i))-backlag:...
                up_trans8{d}(synch_ups8{d}(i))+forwardlag);
            lf8_utrig_lf3l_mat(i,:) = down_3l(up_trans8{d}(synch_ups8{d}(i))-backlag:...
                up_trans8{d}(synch_ups8{d}(i))+forwardlag);
%             lf8_utrig_lf2l_mat(i,:) = down_2l(up_trans8{d}(synch_ups8{d}(i))-backlag:...
%                 up_trans8{d}(synch_ups8{d}(i))+forwardlag);
% 
            if back_dist < backlag
               lf8_utrig_mp_mat(i,1:backlag-back_dist) = nan; 
               lf8_utrig_mpl_mat(i,1:backlag-back_dist) = nan;
            end


        else
            
            lf8_utrig_mp_mat(i,:) = nan;
            lf8_utrig_lf8_mat(i,:) = nan;
            lf8_utrig_lf3_mat(i,:) = nan;
%             lf8_utrig_lf2_mat(i,:) = nan;
               
            lf8_utrig_mpl_mat(i,:) = nan;
            lf8_utrig_lf8l_mat(i,:) = nan;
            lf8_utrig_lf3l_mat(i,:) = nan;
            lf8_utrig_lf2l_mat(i,:) = nan;
         
        end
        
    end
    
  
%     mp_utrig_mp(d,:) = nanmean(mp_utrig_mp_mat);
%     mp_utrig_lf8(d,:) = nanmean(mp_utrig_lf8_mat);
%     mp_utrig_lf3(d,:) = nanmean(mp_utrig_lf3_mat);
%     mp_utrig_lf2(d,:) = nanmean(mp_utrig_lf2_mat);
 
%     mp_utrig_mpl(d,:) = nanmean(mp_utrig_mpl_mat);
%     mp_utrig_lf8l(d,:) = nanmean(mp_utrig_lf8l_mat);
%     mp_utrig_lf3l(d,:) = nanmean(mp_utrig_lf3l_mat);
%     mp_utrig_lf2l(d,:) = nanmean(mp_utrig_lf2l_mat);

% mp_utrig_mp(d,lags < 0) = nan;
% mp_utrig_lf8(d,lags < 0) = nan;
% mp_utrig_lf3(d,lags < 0) = nan;




    lf8_utrig_mp(d,:) = nanmean(lf8_utrig_mp_mat);
    lf8_utrig_lf8(d,:) = nanmean(lf8_utrig_lf8_mat);
    lf8_utrig_lf3(d,:) = nanmean(lf8_utrig_lf3_mat);
%     lf8_utrig_lf2(d,:) = nanmean(lf8_utrig_lf2_mat);

    lf8_utrig_mpl(d,:) = nanmean(lf8_utrig_mpl_mat);
    lf8_utrig_lf8l(d,:) = nanmean(lf8_utrig_lf8l_mat);
    lf8_utrig_lf3l(d,:) = nanmean(lf8_utrig_lf3l_mat);
%     lf8_utrig_lf2l(d,:) = nanmean(lf8_utrig_lf2l_mat);


bad_trans = find(corres_down > length(down_state_dur{d}));
synch_ups8{d}(bad_trans) = [];
corres_down(bad_trans) = [];

%     [dummy,down_order] = sort(down_state_dur{d}(synch_downs{d}));
    [dummy,up_order] = sort(up_state_dur{d}(synch_ups{d}));
%     [dummy,up_order8] = sort(down_state_dur{d}(corres_down));
    [dummy,up_order8] = sort(up_state_dur8{d}(synch_ups8{d}));
    
%% plot mp trig averages    
%    plot(lags,mp_utrig_mpl(d,:))
%    hold on
%    plot(lags,mp_utrig_lf8l(d,:),'r')
%    plot(lags,mp_utrig_lf3l(d,:),'g')
% %    plot(lags,mp_utrig_lf2l(d,:),'k','linewidth',2)
%    legend('MP','LF8','LF3')
%       plot(lags,mp_utrig_mp(d,:),'--')
%       plot(lags,mp_utrig_lf8(d,:),'r--')
%       plot(lags,mp_utrig_lf3(d,:),'g--')
% %       plot(lags,mp_utrig_lf2(d,:),'k--')
% 
%    title('MP Up Triggered Avg')
%    xlim([-2 2]); grid on
%       t_names = ['C:\WC_Germany\overall_calcs\trig_avgs_hf\mp_down_trig_' num2str(cell_type(d)) '_' over_names{d}];
%    print('-dpng',t_names);
% close
%     

% %% plot lf8 trig averages    
%    plot(lags,lf8_utrig_mpl(d,:))
%    hold on
%    plot(lags,lf8_utrig_lf8l(d,:),'r')
% %    plot(lags,lf8_utrig_lf3l(d,:),'g')
% %    plot(lags,lf8_utrig_lf2l(d,:),'k','linewidth',2)
%    legend('MP','LF8')
%    plot(lags,lf8_utrig_mp(d,:),'--')
%    plot(lags,lf8_utrig_lf8(d,:),'r--')
% %    plot(lags,lf8_utrig_lf3(d,:),'g--')
% %    plot(lags,lf8_utrig_lf2(d,:),'k--')
%    title('LF8 Up Triggered Avg')
%       xlim([-2 2]); grid on
%       t_names = ['C:\WC_Germany\overall_calcs\trig_avgs_hf\new_hf_lf8_utrig_' num2str(cell_type(d)) '_' over_names{d}];
%    print('-dpng',t_names);
% close
    
% %% plot mp down trig matrices
%     Fig = figure(1);
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     subplot(2,1,1)
%     pcolor(lags,(1:length(synch_downs{d}))/length(synch_downs{d}) ...
%         ,mp_utrig_mpl_mat(down_order,:));shading flat; colorbar; hold on
%     plot(down_state_dur{d}(synch_downs{d}(down_order)),(1:length(synch_downs{d}))/length(synch_downs{d}),'w','linewidth',2)
%     caxis([-2 2.5]);colorbar
%     xlim([-1 2])
%     line([0 0],[0 1],'Color','k')
%         subplot(2,1,2)
%     pcolor(lags,(1:length(synch_downs{d}))/length(synch_downs{d}) ...
%         ,mp_utrig_mp_mat(down_order,:));shading flat; colorbar; hold on
%     plot(down_state_dur{d}(synch_downs{d}(down_order)),(1:length(synch_downs{d}))/length(synch_downs{d}),'w','linewidth',2)
% caxis([-0.5 1.5]);colorbar
% xlim([-1 2])
%     line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats_hf\new_mp_dtrig_mp_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close
   
%     Fig = figure(1);
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     subplot(2,1,1)
%     pcolor(lags,(1:length(synch_downs{d}))/length(synch_downs{d}) ...
%         ,mp_utrig_lf8l_mat(down_order,:));shading flat; colorbar; hold on
%     plot(down_state_dur{d}(synch_downs{d}(down_order)),(1:length(synch_downs{d}))/length(synch_downs{d}),'w','linewidth',2)
%     xlim([-1 2])
%     caxis([-2 2.5]);colorbar
% line([0 0],[0 1],'Color','k')
%     subplot(2,1,2)
%     pcolor(lags,(1:length(synch_downs{d}))/length(synch_downs{d}) ...
%         ,mp_utrig_lf8_mat(down_order,:));shading flat; colorbar; hold on
%     plot(down_state_dur{d}(synch_downs{d}(down_order)),(1:length(synch_downs{d}))/length(synch_downs{d}),'w','linewidth',2)
%     xlim([-1 2])
% caxis([-1 1]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats_hf\new_mp_dtrig_lf8_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close
% 
%         Fig = figure(1);
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     subplot(2,1,1)
%     pcolor(lags,(1:length(synch_downs{d}))/length(synch_downs{d}) ...
%         ,mp_utrig_lf3l_mat(down_order,:));shading flat; colorbar; hold on
%     plot(down_state_dur{d}(synch_downs{d}(down_order)),(1:length(synch_downs{d}))/length(synch_downs{d}),'w','linewidth',2)
%     xlim([-1 2])
%     caxis([-2 2.5]);colorbar
% line([0 0],[0 1],'Color','k')
%     subplot(2,1,2)
%     pcolor(lags,(1:length(synch_downs{d}))/length(synch_downs{d}) ...
%         ,mp_utrig_lf3_mat(down_order,:));shading flat; colorbar; hold on
%     plot(down_state_dur{d}(synch_downs{d}(down_order)),(1:length(synch_downs{d}))/length(synch_downs{d}),'w','linewidth',2)
%     xlim([-1 2])
% caxis([-1 1]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats_hf\new_mp_dtrig_lf3_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close

%         Fig = figure(1);
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     subplot(2,1,1)
%     pcolor(lags,(1:length(synch_ups{d}))/length(synch_ups{d}) ...
%         ,mp_utrig_lf2l_mat(up_order,:));shading flat; colorbar; hold on
%     plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))/length(synch_ups{d}),'w','linewidth',2)
%     xlim([-1 1])
% caxis([-2 2.5]);colorbar
% line([0 0],[0 1],'Color','k')
%     subplot(2,1,2)
%     pcolor(lags,(1:length(synch_ups{d}))/length(synch_ups{d}) ...
%         ,mp_utrig_lf2_mat(up_order,:));shading flat; colorbar; hold on
%     plot(up_state_dur{d}(synch_ups{d}(up_order)),(1:length(synch_ups{d}))/length(synch_ups{d}),'w','linewidth',2)
%     xlim([-1 1])
%     caxis([-2 2]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats_hf\mp_utrig_lf2_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close

%% plot lfp up trig matrices
    Fig = figure(1);
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
    pcolor(lags,(1:length(synch_ups8{d}))/length(synch_ups8{d}) ...
        ,lf8_utrig_mpl_mat(up_order8,:));shading flat; colorbar; hold on
    plot(up_state_dur8{d}(synch_ups8{d}(up_order8)),(1:length(synch_ups8{d}))/length(synch_ups8{d}),'w','linewidth',2)
    xlim([-1 0.5])
    caxis([-2 2.5]);colorbar
line([0 0],[0 1],'Color','k')
subplot(2,1,2)    
pcolor(lags,(1:length(synch_ups8{d}))/length(synch_ups8{d}) ...
        ,lf8_utrig_mp_mat(up_order8,:));shading flat; colorbar; hold on
    plot(up_state_dur8{d}(synch_ups8{d}(up_order8)),(1:length(synch_ups8{d}))/length(synch_ups8{d}),'w','linewidth',2)
    xlim([-1 0.5])
caxis([-0.5 0.5]);colorbar
line([0 0],[0 1],'Color','k')
    t_names = ['C:\WC_Germany\overall_calcs\trig_mats_hf\full_up_lf8_utrig_mp_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close
   
%     Fig = figure(1);
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     subplot(2,1,1)
%     pcolor(lags,(1:length(synch_ups8{d}))/length(synch_ups8{d}) ...
%         ,lf8_utrig_lf8l_mat(up_order8,:));shading flat; colorbar; hold on
%     plot(up_state_dur8{d}(synch_ups8{d}(up_order8)),(1:length(synch_ups8{d}))/length(synch_ups8{d}),'w','linewidth',2)
%     xlim([-1 2])
%     caxis([-2 2.5]);colorbar
% line([0 0],[0 1],'Color','k')
%     subplot(2,1,2)
%     pcolor(lags,(1:length(synch_ups8{d}))/length(synch_ups8{d}) ...
%         ,lf8_utrig_lf8_mat(up_order8,:));shading flat; colorbar; hold on
%     plot(up_state_dur8{d}(synch_ups8{d}(up_order8)),(1:length(synch_ups8{d}))/length(synch_ups8{d}),'w','linewidth',2)
%     xlim([-1 2])
% caxis([-0.5 1]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats_hf\full_up_lf8_utrig_lf8_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close
    
    
        Fig = figure(1);
    clf
    set(Fig,'PaperUnits','centimeters');
    set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
    set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
    subplot(2,1,1)
    pcolor(lags,(1:length(synch_downs8{d}))/length(synch_downs8{d}) ...
        ,lf8_dtrig_mpl_mat);shading flat; colorbar; hold on
%     plot(up_state_dur8{d}(synch_ups8{d}(up_order8)),(1:length(synch_ups8{d}))/length(synch_ups8{d}),'w','linewidth',2)
    xlim([-1 1])
    caxis([-2 2.5]);colorbar
line([0 0],[0 1],'Color','k')
subplot(2,1,2)    
pcolor(lags,(1:length(synch_downs8{d}))/length(synch_downs8{d}) ...
        ,lf8_dtrig_mp_mat);shading flat; colorbar; hold on
%     plot(up_state_dur8{d}(synch_ups8{d}(up_order8)),(1:length(synch_ups8{d}))/length(synch_ups8{d}),'w','linewidth',2)
    xlim([-1 1])
caxis([-0.5 0.5]);colorbar
line([0 0],[0 1],'Color','k')
    t_names = ['C:\WC_Germany\overall_calcs\trig_mats_hf\full_up_lf8_dtrig_mp_' num2str(cell_type(d)) '_' over_names{d}];
    print('-dpng',t_names);
    close
   
%     Fig = figure(1);
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     subplot(2,1,1)
%     pcolor(lags,(1:length(synch_ups8{d}))/length(synch_ups8{d}) ...
%         ,lf8_utrig_lf8l_mat(up_order8,:));shading flat; colorbar; hold on
%     plot(up_state_dur8{d}(synch_ups8{d}(up_order8)),(1:length(synch_ups8{d}))/length(synch_ups8{d}),'w','linewidth',2)
%     xlim([-1 2])
%     caxis([-2 2.5]);colorbar
% line([0 0],[0 1],'Color','k')
%     subplot(2,1,2)
%     pcolor(lags,(1:length(synch_ups8{d}))/length(synch_ups8{d}) ...
%         ,lf8_utrig_lf8_mat(up_order8,:));shading flat; colorbar; hold on
%     plot(up_state_dur8{d}(synch_ups8{d}(up_order8)),(1:length(synch_ups8{d}))/length(synch_ups8{d}),'w','linewidth',2)
%     xlim([-1 2])
% caxis([-0.5 1]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats_hf\full_up_lf8_dtrig_lf8_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close

% 
%         Fig = figure(1);
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     subplot(2,1,1)
%     pcolor(lags,(1:length(synch_ups8{d}))/length(synch_ups8{d}) ...
%         ,lf8_utrig_lf3l_mat);shading flat; colorbar; hold on
%     plot(up_state_dur8{d}(synch_ups8{d}(up_order8)),(1:length(synch_ups8{d}))/length(synch_ups8{d}),'w','linewidth',2)
%     xlim([-1 2])
%     caxis([-2 2.5]);colorbar
% line([0 0],[0 1],'Color','k')
%     subplot(2,1,2)
%     pcolor(lags,(1:length(synch_ups8{d}))/length(synch_ups8{d}) ...
%         ,lf8_utrig_lf3_mat(up_order8,:));shading flat; colorbar; hold on
%     plot(up_state_dur8{d}(synch_ups8{d}(up_order8)),(1:length(synch_ups8{d}))/length(synch_ups8{d}),'w','linewidth',2)
%     xlim([-1 2])
% caxis([-0.5 1]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats_hf\full_up_lf8_utrig_lf3_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close
% 
%         Fig = figure(1);
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     subplot(2,1,1)
%     pcolor(lags,(1:length(synch_ups8{d}))/length(synch_ups8{d}) ...
%         ,lf8_utrig_lf2l_mat(up_order8,:));shading flat; colorbar; hold on
% %     plot(up_state_dur8{d}(synch_ups8{d}(up_order8)),(1:length(synch_ups8{d}))/length(synch_ups8{d}),'w','linewidth',2)
%     xlim([-1 1])
% caxis([-2 2.5]);colorbar
% line([0 0],[0 1],'Color','k')
%     subplot(2,1,2)
%     pcolor(lags,(1:length(synch_ups8{d}))/length(synch_ups8{d}) ...
%         ,lf8_utrig_lf2_mat(up_order8,:));shading flat; colorbar; hold on
% %     plot(up_state_dur8{d}(synch_ups8{d}(up_order8)),(1:length(synch_ups8{d}))/length(synch_ups8{d}),'w','linewidth',2)
%     xlim([-1 1])
%     caxis([-2 2]);colorbar
% line([0 0],[0 1],'Color','k')
%     t_names = ['C:\WC_Germany\overall_calcs\trig_mats_hf\lf8_utrig_lf2_' num2str(cell_type(d)) '_' over_names{d}];
%     print('-dpng',t_names);
%     close
clear lf8_utrig*
end

clear *mat
% save C:\WC_Germany\overall_calcs\trig_avgs_comb\trig_avg_data lags mp* lf8*