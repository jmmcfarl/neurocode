clear all
close all
load C:\WC_Germany\Persistent_activity\dir_tree_update
load C:\WC_Germany\Persistent_activity\UDS_dur_raw\UDS_raw_data
load C:\WC_Germany\Persistent_activity\UDS_dist\synch_UDS_inds

Fs = 2016;
niqf = Fs/2;
dsf = 8;
Fsd = Fs/dsf;

lcf = 0.05/niqf;
hcf = 4/niqf;
[blow,alow] = butter(2,[lcf hcf]);
maxlag = 1*Fsd;
lags = -maxlag:maxlag;

backlag = 3*Fsd;
forwardlag = 8*Fsd
lags = -backlag:forwardlag;

% lcf = 2/niqf;
% hcf = 10/niqf;
% [bmid,amid] = butter(2,[lcf hcf]);
%
% lcf = 100/niqf;
% hcf = 250/niqf;
% [bhigh,ahigh] = butter(2,[lcf hcf]);
flip_sign = [1 1 1 1 1 1 -1 1 -1 1 -1 -1 -1 1 1 -1 -1];

win = 20;
% scale_fact = [.05 0 -0.01 0.15 0.1 0.15 0.52 0.28 0.53 0.08 0.75 0.81 0.62 0.08 0.14 0.7 0.35];

for d = 1:length(dir_array)
    

    cd(dir_array{d})
    pwd

    load used_data lf2 lf3 lf5 lf6 lf7 lf8 wcv_minus_spike
    wcv_f = filtfilt(blow,alow,wcv_minus_spike);
    lf3_uds = filtfilt(blow,alow,lf3);
    lf5_uds = filtfilt(blow,alow,lf5);
    lf8_uds = filtfilt(blow,alow,lf8);
    wcv_d = downsample(wcv_f,dsf);

    lf3_uds = downsample(lf3_uds,dsf);
    lf8_uds = downsample(lf8_uds,dsf);
    lf5_uds = downsample(lf5_uds,dsf);

    lf8_uds = zscore(lf8_uds);
    lf3_uds = zscore(lf3_uds);
    lf5_uds = zscore(lf5_uds);
    
    a = corrcoef(lf3_uds,lf5_uds);
    scale_fact = a(2,1);
    lf3_t = lf3_uds-scale_fact*lf5_uds;
    t_axis = (1:length(lf3_uds))/Fsd;
    wcv_d = zscore(wcv_d);
    
    clear *_f clear lf8 lf2 lf3 lf5 lf7
    
    
    wcv_up_ctrig_mat = zeros(length(up_trans{d}),length(lags));
    lf8_up_ctrig_mat = wcv_up_ctrig_mat;
%     lf2_up_ctrig_mat = wcv_up_ctrig_mat;
    lf3_up_ctrig_mat = wcv_up_ctrig_mat;
        lf3t_up_ctrig_mat = wcv_up_ctrig_mat;

%     lf5_up_ctrig_mat = wcv_up_ctrig_mat;
%     lf7_up_ctrig_mat = wcv_up_ctrig_mat;
    

    for i = 1:length(up_trans{d})

        
        %find closest lfp up transition
        %        [dummy,close_lf8up_ind] = min(abs(up_trans8{d}-up_trans{d}(i)));
        %        close_lf8up = up_trans8{d}(close_lf8up_ind);
%         lfp_period_dur{d}(i) = lf8_period_f{d}(down_trans{d}(i))-lf8_period_f{d}(up_trans{d}(i));
        
%         first_in_lfp_down = find(down_trans8{d} > up_trans{d}(i),1,'first');
%         if ~isempty(first_in_lfp_down)
%             %find all lfp up transitions that occur before mp comes back down
%             secondary_lfp_states{i} = find(down_trans8{d} < down_trans{d}(i) & up_trans8{d} > down_trans8{d}(first_in_lfp_down));
%             if ~isempty(secondary_lfp_states{i})
%                 for q = 1:length(secondary_lfp_states{i})
%                     lf8_d(up_trans8{d}(secondary_lfp_states{i}(q)):down_trans8{d}(secondary_lfp_states{i}(q))) = nan;
%                 end
%             end
%         end
%         %does mp down occur during an lfp down
%         prev_lfp_up = find(up_trans8{d}<up_trans{d}(i),1,'last');
%         prev_lfp_down = find(down_trans8{d}<down_trans{d}(i),1,'last');
%         next_lfp_up = find(up_trans8{d}>down_trans{d}(i),1,'first');
%         if ~isempty(prev_lfp_up)
%             lfp_pup(i) = up_trans8{d}(prev_lfp_up);
%             lfp_pupdur(i) = up_state_dur8{d}(prev_lfp_up);
%         else
%             lfp_pup(i) = nan;
%             lfp_pupdur(i) = nan;
%         end
%         if ~isempty(prev_lfp_down)
%             lfp_pdown(i) = down_trans8{d}(prev_lfp_down);
%         else
%             lfp_pdown(i) = nan;
%         end
%         if ~isempty(next_lfp_up)
%             lfp_nup(i) = up_trans8{d}(next_lfp_up);
%         else
%             lfp_nup(i) = nan;
%         end
%         
%         if ~isempty(prev_lfp_down) & prev_lfp_down < length(up_trans8{d})
%             mp_down_state(i) = up_trans8{d}(prev_lfp_down+1) > down_trans{d}(i);
%         else
%             mp_down_state(i) = nan;
%         end
        if up_trans{d}(i) > backlag & up_trans{d}(i) < length(wcv_d)-forwardlag
%             wcv_up_ctrig_mat(i,:) = wcv_d(up_trans{d}(i)-backlag:up_trans{d}(i)+forwardlag);
            lf8_up_ctrig_mat(i,:) = lf8_uds(up_trans{d}(i)-backlag:up_trans{d}(i)+forwardlag);
%             lf2_up_ctrig_mat(i,:) = lf2_d(up_trans{d}(i)-backlag:up_trans{d}(i)+forwardlag);
            lf3_up_ctrig_mat(i,:) = lf3_uds(up_trans{d}(i)-backlag:up_trans{d}(i)+forwardlag);
            lf3t_up_ctrig_mat(i,:) = lf3_t(up_trans{d}(i)-backlag:up_trans{d}(i)+forwardlag);

%             lf5_up_ctrig_mat(i,:) = lf5_d(up_trans{d}(i)-backlag:up_trans{d}(i)+forwardlag);
%             lf7_up_ctrig_mat(i,:) = lf7_d(up_trans{d}(i)-backlag:up_trans{d}(i)+forwardlag);


        else
%             wcv_up_ctrig_mat(i,:) = nan;
            lf8_up_ctrig_mat(i,:) = nan;
%             lf2_up_ctrig_mat(i,:) = nan;
            lf3_up_ctrig_mat(i,:) = nan;
                        lf3t_up_ctrig_mat(i,:) = nan;

%             lf5_up_ctrig_mat(i,:) = nan;
%             lf7_up_ctrig_mat(i,:) = nan;

        end

    end

%         wcv_down_ctrig_mat = zeros(length(down_trans{d}),length(dlags));
%         lf8_down_ctrig_mat = wcv_down_ctrig_mat;
%          lf2_down_ctrig_mat = wcv_down_ctrig_mat;
%          lf3_down_ctrig_mat = wcv_down_ctrig_mat;
%           lf5_down_ctrig_mat = wcv_down_ctrig_mat;
%         lf7_down_ctrig_mat = wcv_down_ctrig_mat;
% 
%             for i = 1:length(down_trans{d})
%     
%     
%     %             prev_lfp_down(i) = down_trans8{d}(find(down_trans8{d} < down_trans{d}(i),1,'first'));
%     
%                if down_trans{d}(i) > forwardlag & down_trans{d}(i) < length(wcv_d)-backlag
%                    wcv_down_ctrig_mat(i,:) = wcv_d(down_trans{d}(i)-backlag:down_trans{d}(i)+backlag);
%                    lf8_down_ctrig_mat(i,:) = lf8_d(down_trans{d}(i)-backlag:down_trans{d}(i)+backlag);
%                    lf2_down_ctrig_mat(i,:) = lf2_d(down_trans{d}(i)-backlag:down_trans{d}(i)+backlag);
%                    lf3_down_ctrig_mat(i,:) = lf3_d(down_trans{d}(i)-backlag:down_trans{d}(i)+backlag);
%                    lf5_down_ctrig_mat(i,:) = lf5_d(down_trans{d}(i)-backlag:down_trans{d}(i)+backlag);
%                    lf7_down_ctrig_mat(i,:) = lf7_d(down_trans{d}(i)-backlag:down_trans{d}(i)+backlag);
% 
%                else
%                    wcv_down_ctrig_mat(i,:) = nan;
%                    lf8_down_ctrig_mat(i,:) = nan;
%                    lf2_down_ctrig_mat(i,:) = nan;
%                    lf3_down_ctrig_mat(i,:) = nan;
%                    lf5_down_ctrig_mat(i,:) = nan;
%                    lf7_down_ctrig_mat(i,:) = nan;
% 
%                end
%     
%             end
% 
%         std_down_diff(d) = std(lfp_down_diff);
%         std_up_diff(d) = std(lfp_up_diff);
%        cv_down_diff(d) = std(lfp_down_diff)/mean(lfp_down_diff);
%        cv_up_diff(d) = std(lfp_up_diff)/mean(lfp_up_diff);

% wcv_down_ctrig_avg(d,:) = nanmean(wcv_down_ctrig_mat);
% lf2_down_ctrig_avg(d,:) = nanmean(lf2_down_ctrig_mat);
% lf3_down_ctrig_avg(d,:) = nanmean(lf3_down_ctrig_mat);
% lf5_down_ctrig_avg(d,:) = nanmean(lf5_down_ctrig_mat);
% lf7_down_ctrig_avg(d,:) = nanmean(lf7_down_ctrig_mat);
% lf8_down_ctrig_avg(d,:) = nanmean(lf8_down_ctrig_mat);
% wcv_up_ctrig_avg(d,:) = nanmean(wcv_up_ctrig_mat);
% lf2_up_ctrig_avg(d,:) = nanmean(lf2_up_ctrig_mat);
% lf3_up_ctrig_avg(d,:) = nanmean(lf3_up_ctrig_mat);
% lf5_up_ctrig_avg(d,:) = nanmean(lf5_up_ctrig_mat);
% lf7_up_ctrig_avg(d,:) = nanmean(lf7_up_ctrig_mat);
% lf8_up_ctrig_avg(d,:) = nanmean(lf8_up_ctrig_mat);

%     lfp_dur_calc = (lfp_pdown-lfp_pup)/Fsd;
%     tsld = (down_trans{d}-lfp_pdown)/Fsd;
%     ttnu = (lfp_nup-down_trans{d})/Fsd;
%     tslu = (up_trans{d}-lfp_pup)/Fsd;
%     
%     up_dur_ratio{d} = up_state_dur{d}./lfp_dur_calc;
%   
% 
    [dummy,up_order] = sort(up_state_dur{d});
%     Fig = figure(1)
%     clf
%     set(Fig,'PaperUnits','centimeters');
%     set(gcf, 'PaperSize', [30 20]);% paper size is in [width height] format
%     set(Fig,'PaperPosition',[0,0,(get(Fig,'PaperSize'))])
%     %
%     subplot(2,1,1)
%     pcolor(lags/Fsd,(1:length(up_trans{d}))/length(up_trans{d}) ...
%         ,wcv_up_ctrig_mat(up_order,:));shading flat;
%     colorbar
%     hold on
%     plot(up_state_dur{d}(up_order),(1:length(up_trans{d}))/length(up_trans{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%     subplot(2,1,2)

%     pcolor(lags/Fsd,(1:length(up_trans{d}))/length(up_trans{d})...
%         ,lf8_up_ctrig_mat(up_order,:));shading flat;
%     caxis([-1 4]);colorbar
%     hold on
%     plot(up_state_dur{d}(up_order),(1:length(up_trans{d}))/length(up_trans{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%     grid on
%     cd 'C:\WC_Germany\JMM_Analysis_Pyr\lf8_up_ctrig_all'
%     tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_up_ctrig_all\lf8_' f_names{d}];
%     print('-dpng',tname);
%     tname = ['lf8_' f_names{d} '.fig'];
% 
% saveas(gcf,tname)
%     close
% 
%         pcolor(lags/Fsd,(1:length(up_trans{d}))/length(up_trans{d})...
%         ,lf2_up_ctrig_mat(up_order,:));shading flat;
%     caxis([-1 4]);colorbar
%     hold on
%     grid on
%     plot(up_state_dur{d}(up_order),(1:length(up_trans{d}))/length(up_trans{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
% 
%     tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_up_ctrig_all\lf2_' f_names{d}];
%     print('-dpng',tname);
% tname = ['lf2_' f_names{d} '.fig'];
% 
% saveas(gcf,tname)
% 
%     close
% 
%         pcolor(lags/Fsd,(1:length(up_trans{d}))/length(up_trans{d})...
%         ,lf3_up_ctrig_mat(up_order,:));shading flat;
%     caxis([-1 4]);colorbar
%     hold on
%     grid on
%     plot(up_state_dur{d}(up_order),(1:length(up_trans{d}))/length(up_trans{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%     tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_up_ctrig_all\lf3_' f_names{d}];
%     print('-dpng',tname);
% tname = ['lf3_' f_names{d} '.fig'];
% 
% saveas(gcf,tname)
% close
% 
%         pcolor(lags/Fsd,(1:length(up_trans{d}))/length(up_trans{d})...
%         ,lf5_up_ctrig_mat(up_order,:));shading flat;
%     caxis([-1 4]);colorbar
%     hold on
%     grid on
%     plot(up_state_dur{d}(up_order),(1:length(up_trans{d}))/length(up_trans{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
% tname = ['lf5_' f_names{d} '.fig'];
%     
% saveas(gcf,tname)
%     tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_up_ctrig_all\lf5_' f_names{d}];
% 
%     print('-dpng',tname);
%     close
% 
%         pcolor(lags/Fsd,(1:length(up_trans{d}))/length(up_trans{d})...
%         ,lf7_up_ctrig_mat(up_order,:));shading flat;
%     caxis([-1 4]);colorbar
%     hold on
%     grid on
%     plot(up_state_dur{d}(up_order),(1:length(up_trans{d}))/length(up_trans{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
% tname = ['lf7_' f_names{d} '.fig'];
% saveas(gcf,tname)
%     tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_up_ctrig_all\lf7_' f_names{d}];
%     print('-dpng',tname);
%     close

    %calculate fraction of MP up states with secondary states
%     sec_lfp_check{d} = zeros(length(up_trans{d}),1);
%     num_sec_states{d} = zeros(length(up_trans{d}),1);
% 
%     for i = 1:length(up_trans{d})
% 
%         sec_lfp_check{d}(i) = ~isempty(secondary_lfp_states{i});
%         num_sec_states{d}(i) = length(secondary_lfp_states{i});
% 
%     end
% 
%     pers_fract(d) = nansum(sec_lfp_check{d})/length(sec_lfp_check{d});
%     mp_down_frac(d) = nansum(mp_down_state)/length(mp_down_state);
% 
% hrange = linspace(0,3,50);

% no_pers_states = find(num_sec_states{d}==0);
% pers_states = find(num_sec_states{d} > 0);
% plot(up_state_dur{d}(no_pers_states),lfp_pupdur(no_pers_states),'o')
% hold on
% plot(up_state_dur{d}(pers_states),lfp_pupdur(pers_states),'ro')
% ylim([0 3])
% tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_up_ctrig\up_dur_cor_PS_' f_names{d}];
% print('-dpng',tname);
% close
% 
% plot(tslu(no_pers_states),tsld(no_pers_states),'o')
% hold on
% plot(tslu(pers_states),tsld(pers_states),'ro')
% xlim([0 3])
% ylim([0 3])
% tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_up_ctrig\tsld_tslu_PS_' f_names{d}];
% print('-dpng',tname);
% close

% up_locking_cv(d) = nanstd(tslu)/nanmean(tslu);
% down_locking_cv(d) = nanstd(tsld)/nanmean(tsld);

% hist(lfp_period_dur{d},100)
% tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_up_ctrig\lfp_fract_period_' f_names{d}];
% print('-dpng',tname);
% close


    d
    %
    %   figure
    %       subplot(2,1,1)
    %     pcolor(dlags/Fsd,1:length(down_trans{d}),wcv_down_ctrig_mat(up_order,:));shading flat
    %     subplot(2,1,2)
%     close all
        pcolor(lags/Fsd,(1:length(up_trans{d}))/length(up_trans{d}),lf8_up_ctrig_mat(up_order,:));shading flat
        hold on
        caxis([-1 3]);colorbar
    plot(up_state_dur{d}(up_order),(1:length(up_trans{d}))/length(up_trans{d}),'w','linewidth',2)
    line([0 0],[0 1],'Color','k')
        tname = ['C:\WC_Germany\JMM_Analysis_Pyr\newest_mp_up_trig_lf8_lf3_lf2\low_freq_lf8_' f_names{d}];
        print('-dpng',tname);
        close
% %     %
%         pcolor(lags/Fsd,(1:length(up_trans{d}))/length(up_trans{d}),lf2_up_ctrig_mat(up_order,:));shading flat
%         hold on
%                 caxis([-1 3]);colorbar
%     plot(up_state_dur{d}(up_order),(1:length(up_trans{d}))/length(up_trans{d}),'w','linewidth',2)
%     line([0 0],[0 1],'Color','k')
%         tname = ['C:\WC_Germany\JMM_Analysis_Pyr\newest_mp_up_trig_lf8_lf3_lf2\lf2_' f_names{d}];
%         print('-dpng',tname);
%         close
        
        pcolor(lags/Fsd,(1:length(up_trans{d}))/length(up_trans{d}),lf3_up_ctrig_mat(up_order,:));shading flat
        hold on
                caxis([-1 3]);colorbar
    plot(up_state_dur{d}(up_order),(1:length(up_trans{d}))/length(up_trans{d}),'w','linewidth',2)
    line([0 0],[0 1],'Color','k')
        tname = ['C:\WC_Germany\JMM_Analysis_Pyr\newest_mp_up_trig_lf8_lf3_lf2\low_freq_lf3_' f_names{d}];
        print('-dpng',tname);
        close

        pcolor(lags/Fsd,(1:length(up_trans{d}))/length(up_trans{d}),lf3t_up_ctrig_mat(up_order,:));shading flat
        hold on
                caxis([-1 3]);colorbar
    plot(up_state_dur{d}(up_order),(1:length(up_trans{d}))/length(up_trans{d}),'w','linewidth',2)
    line([0 0],[0 1],'Color','k')
    caxis([-1 2]);colorbar
        tname = ['C:\WC_Germany\JMM_Analysis_Pyr\newest_mp_up_trig_lf8_lf3_lf2\low_freq_lf3_corrected_' f_names{d}];
        print('-dpng',tname);
        close

%         
%         pcolor(dlags/Fsd,1:length(down_trans{d}),lf5_down_ctrig_mat(up_order,:));shading flat
%         tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_down_ctrig_all\lf5_' f_names{d}];
%         print('-dpng',tname);
%         close
% 
%                 pcolor(dlags/Fsd,1:length(down_trans{d}),lf7_down_ctrig_mat(up_order,:));shading flat
%         tname = ['C:\WC_Germany\JMM_Analysis_Pyr\lf8_down_ctrig_all\lf7_' f_names{d}];
%         print('-dpng',tname);
%         close

%     clear near_lfp_ind
%     clear prev_lfp_down

%     clear wcv_d lf2_d lf3_d lf5_d lf7_d lf8_d
    

end