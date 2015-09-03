clear all
close all

load G:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')

load G:\WC_Germany\persistent_9_27_2010\pa_state_dur_stats

cd G:\WC_Germany\persistent_9_27_2010\
load ./depthofanesth_data lf8_sep mp_sep

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

drive_letter = 'G';
dsf = 8;
Fsd = 2016/dsf;
lf8_up_fract = 0.5;

hist_range = linspace(0,7,1000);
binsize = hist_range(2)-hist_range(1);
x = linspace(0,7,100);
xbinsize = x(2)-x(1);

drange = linspace(-5,7,1000);

% up_range = [0.1 10];
% down_range = [0.1 10];
% % up_range = [0.1 3];
% % down_range = [0.1 3];
% numBins = 50;
% lin_grid = linspace(up_range(1),up_range(2),numBins+1);
up_range = [0.3 7];
down_range = [0.3 7];
numBins = 25;
lin_grid = linspace(up_range(1),up_range(2),numBins+1);
log_grid = logspace(log10(up_range(1)),log10(up_range(2)),numBins+1);

for d = 1:length(sess_data)
    
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    
    load ./used_data wcv
    wcv_d = downsample(wcv,dsf);
    datalen = length(wcv_d);
    load ./pa_hsmm_state_seq
    mp_state_seq = hsmm_bbstate_seq;
    load ./pa_hsmm_state_seq8
    lfp_state_seq = hsmm_bbstate_seq8;
    
    load ./pa_lf8_period_data
    
    [new_seg_inds] = resample_uds_seg_inds(hmm.UDS_segs,hmm.Fs,Fsd,length(wcv_d));
    seg_durs = new_seg_inds(:,2)-new_seg_inds(:,1)+1;
    [mp_up_inds{d},mp_down_inds{d}] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
    [lf8_up_inds{d},lf8_down_inds{d}] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);
    
    mp_updur{d} = nan(size(mp_up_inds{d})); %MP up state durations
    mp_updur_lfpc_delay{d} = nan(size(mp_up_inds{d})); %MP up state duration (in LFP cycles)
    mp_updur_lfpc{d} = nan(size(mp_up_inds{d})); %MP up state duration (in LFP cycles)
    mp_downdur{d} = nan(size(mp_up_inds{d})); %MP up state durations
    mp_downdur_lfpc_delay{d} = nan(size(mp_up_inds{d})); %MP up state duration (in LFP cycles)
    mp_downdur_lfpc{d} = nan(size(mp_up_inds{d})); %MP up state duration (in LFP cycles)
    corresp_lfp_updur{d} = nan(size(mp_up_inds{d}));
    corresp_lfp_downdur{d} = nan(size(mp_up_inds{d}));
    corresp_lfp_maxskippeddowndur{d} = nan(size(mp_up_inds{d}));
    corresp_lfp_mindowndur{d} = nan(size(mp_up_inds{d}));
    corresp_lfp_dmaxskippeddowndur{d} = nan(size(mp_up_inds{d}));
    corresp_lfp_dmindowndur{d} = nan(size(mp_up_inds{d}));
    corresp_lfp_maxskippedupdur{d} = nan(size(mp_up_inds{d}));
    corresp_lfp_minupdur{d} = nan(size(mp_up_inds{d}));
    corresp_lfp_dmaxskippedupdur{d} = nan(size(mp_up_inds{d}));
    corresp_lfp_dminupdur{d} = nan(size(mp_up_inds{d}));
    lfp_skipped_downs{d} = [];
    lfp_dskipped_downs{d} = [];
    lfp_skipped_ups{d} = [];
    lfp_dskipped_ups{d} = [];
    cnt = 1;
    
    mp_state_vec = nan(size(wcv_d));
    
    [lfp_state_durations] = compute_state_durations_seg(lfp_state_seq,Fsd);
    
    for ns = 1:hmm.Nsegs
        
        cur_up_trans = 1+find(mp_state_seq{ns}(1:end-1) == 1 & mp_state_seq{ns}(2:end) == 2);
        cur_down_trans = 1+find(mp_state_seq{ns}(1:end-1) == 2 & mp_state_seq{ns}(2:end) == 1);
        cur_up_trans(cur_up_trans > cur_down_trans(end)) = [];
        cur_down_trans(cur_down_trans < cur_up_trans(1)) = [];
        
        cur_up_trans8 = 1+find(lfp_state_seq{ns}(1:end-1) == 1 & lfp_state_seq{ns}(2:end) == 2);
        cur_down_trans8 = 1+find(lfp_state_seq{ns}(1:end-1) == 2 & lfp_state_seq{ns}(2:end) == 1);
        cur_up_trans8(cur_up_trans8 > cur_down_trans8(end)) = [];
        cur_down_trans8(cur_down_trans8 < cur_up_trans8(1)) = [];
        
        n_lfp_ups = length(cur_up_trans8);
        n_mp_ups = length(cur_up_trans);
        
        for i = 1:length(cur_up_trans)
            mp_state_vec(cur_up_trans(i):cur_down_trans(i)) = 1;
        end
        for i = 1:length(cur_up_trans)-1
            mp_state_vec(cur_down_trans(i):cur_up_trans(i+1)) = 0;
        end
        
        for i = 1:length(cur_up_trans)-1
            
            %find nearest lfp up transition
            [dummy,near_lfp_up] = min(abs(cur_up_trans(i) - cur_up_trans8));
            
            %mp up state dur
            mp_updur{d}(cnt) = (cur_down_trans(i)-cur_up_trans(i))/Fsd;
            mp_downdur{d}(cnt) = (cur_up_trans(i+1)-cur_down_trans(i))/Fsd;
            
            cur_lfp_up_lag = (cur_up_trans(i)-cur_up_trans8(near_lfp_up));
            
            %find corresponding down transition
            next_lfp_down = cur_down_trans8(near_lfp_up);
            
            %if the MP down transition happens before the down-transition
            %following the nearest LFP up-transition than call this the
            %corresponding down transition
            if cur_down_trans(i) < next_lfp_down
                corresp_lfp_down = near_lfp_up;
            else
                %otherwise use the LFP down-transition which precedes the
                %MP down transition most closely
                corresp_lfp_down = find(cur_down_trans8 < cur_down_trans(i),1,'last');
            end
            if ~isempty(corresp_lfp_down)
                cur_lfp_down_lag = cur_down_trans(i)-cur_down_trans8(corresp_lfp_down);
            end
            
            %store mp up state dur in lfp cycles
            cur_lfpc = lf8_period_f{ns}(cur_down_trans(i) - cur_lfp_up_lag) ...
                - lf8_period_f{ns}(cur_up_trans(i)-cur_lfp_up_lag);
            mp_updur_lfpc_delay{d}(cnt) = cur_lfpc;
            cur_lfpc = lf8_period_f{ns}(cur_down_trans(i)) - lf8_period_f{ns}(cur_up_trans(i));
            mp_updur_lfpc{d}(cnt) = cur_lfpc;
            
            cur_lfpc = lf8_period_f{ns}(cur_up_trans(i+1)-cur_lfp_up_lag) ...
                - lf8_period_f{ns}(cur_down_trans(i)-cur_lfp_down_lag);
            mp_downdur_lfpc_delay{d}(cnt) = cur_lfpc;
            cur_lfpc = lf8_period_f{ns}(cur_up_trans(i+1)) - lf8_period_f{ns}(cur_down_trans(i));
            mp_downdur_lfpc{d}(cnt) = cur_lfpc;
            
            corresp_lfp_updur{d}(cnt) = (cur_down_trans8(near_lfp_up)-cur_up_trans8(near_lfp_up))/Fsd;
            temp = find(lf8_up_inds{d} == cur_up_trans8(near_lfp_up) + new_seg_inds(ns)-1);
            if ~isempty(temp)
                corresp_lfp_upid{d}(cnt) = temp;
            end
            if n_lfp_ups > corresp_lfp_down+1 %if the nearest up transition is not the final LFP up-transition
                corresp_lfp_downdur{d}(cnt) = (cur_up_trans8(corresp_lfp_down+1)-cur_down_trans8(corresp_lfp_down))/Fsd;
            end
            
            %% find skipped down states
            skipped_lfp_downs = find(cur_down_trans8(1:end-1) > cur_up_trans8(near_lfp_up) & cur_down_trans8(1:end-1) < cur_down_trans(i) ...
                & cur_down_trans8(1:end-1) > cur_up_trans(i) & cur_up_trans8(2:end) < cur_down_trans(i));
            dskipped_lfp_downs = find(cur_down_trans8(1:end-1) > cur_up_trans8(near_lfp_up) & cur_down_trans8(1:end-1) < cur_down_trans(i)-cur_lfp_up_lag ...
                & cur_down_trans8(1:end-1) > cur_up_trans(i)-cur_lfp_up_lag & cur_up_trans8(2:end) < cur_down_trans(i)-cur_lfp_up_lag);
            lfp_skipped_downs{d} = [lfp_skipped_downs{d}; find(ismember(lf8_down_inds{d},cur_down_trans8(skipped_lfp_downs)+new_seg_inds(ns)-1))];
            lfp_dskipped_downs{d} = [lfp_dskipped_downs{d}; find(ismember(lf8_down_inds{d},cur_down_trans8(dskipped_lfp_downs)+new_seg_inds(ns)-1))];
            if ~isempty(skipped_lfp_downs)
                corresp_lfp_maxskippeddowndur{d}(cnt) = nanmax(lfp_state_durations{1}(skipped_lfp_downs));
                corresp_lfp_mindowndur{d}(cnt) = nanmin(lfp_state_durations{1}(skipped_lfp_downs));
            else
                corresp_lfp_mindowndur{d}(cnt) = corresp_lfp_downdur{d}(cnt);
            end
            if ~isempty(dskipped_lfp_downs)
                corresp_lfp_dmaxskippeddowndur{d}(cnt) = nanmax(lfp_state_durations{1}(dskipped_lfp_downs));
                corresp_lfp_dmindowndur{d}(cnt) = nanmin(lfp_state_durations{1}(dskipped_lfp_downs));
            else
                corresp_lfp_dmindowndur{d}(cnt) = corresp_lfp_downdur{d}(cnt);
            end
            
            %% find skipped up states
            if n_lfp_ups > i
                skipped_lfp_ups = find(cur_up_trans8(1:end-1) > cur_down_trans8(corresp_lfp_down) & cur_up_trans8(1:end-1) < cur_up_trans(i+1) ...
                    & cur_up_trans8(1:end-1) > cur_down_trans(i) & cur_down_trans8(1:end-1) < cur_up_trans(i+1));
                dskipped_lfp_ups = find(cur_up_trans8(1:end-1) > cur_down_trans8(corresp_lfp_down) & cur_up_trans8(1:end-1) < cur_up_trans(i+1) - cur_lfp_down_lag ...
                    & cur_up_trans8(1:end-1) > cur_down_trans(i)-cur_lfp_down_lag & cur_down_trans8(1:end-1) < cur_up_trans(i+1) - cur_lfp_down_lag);
                lfp_skipped_ups{d} = [lfp_skipped_ups{d}; find(ismember(lf8_up_inds{d},cur_up_trans8(skipped_lfp_ups)+new_seg_inds(ns)-1))];
                lfp_dskipped_ups{d} = [lfp_dskipped_ups{d}; find(ismember(lf8_up_inds{d},cur_up_trans8(dskipped_lfp_ups)+new_seg_inds(ns)-1))];
                if ~isempty(skipped_lfp_ups)
                    corresp_lfp_maxskippedupdur{d}(cnt) = nanmax(lfp_state_durations{2}(skipped_lfp_ups));
                    corresp_lfp_mindowndur{d}(cnt) = nanmin(lfp_state_durations{2}(skipped_lfp_ups));
                else
                    corresp_lfp_mindowndur{d}(cnt) = corresp_lfp_updur{d}(cnt);
                end
                if ~isempty(dskipped_lfp_ups)
                    corresp_lfp_dmaxskippedupdur{d}(cnt) = nanmax(lfp_state_durations{2}(dskipped_lfp_ups));
                    corresp_lfp_dminupdur{d}(cnt) = nanmin(lfp_state_durations{2}(dskipped_lfp_ups));
                else
                    corresp_lfp_dminupdur{d}(cnt) = corresp_lfp_updur{d}(cnt);
                end
            end
            %%
            
            cnt = cnt + 1;
        end
    end
    
    rel_down_dur{d} = mp_downdur{d}./corresp_lfp_downdur{d};
    diff_down_dur{d} = mp_downdur{d} - corresp_lfp_downdur{d};
    rel_up_dur{d} = mp_updur{d}./corresp_lfp_updur{d};
    diff_up_dur{d} = mp_updur{d} - corresp_lfp_updur{d};

    med_down_dur(d) = nanmedian(mp_downdur{d});
    uset = find(~isnan(mp_downdur{d}) & ~isnan(corresp_lfp_downdur{d}));
    down_diff_test(d) = signrank(mp_downdur{d}(uset),corresp_lfp_downdur{d}(uset));
    down_med_diff(d) = nanmedian(diff_down_dur{d});
    down_med_rel(d) = nanmedian(rel_down_dur{d});
    
    med_up_dur(d) = nanmedian(mp_updur{d});
    uset = find(~isnan(mp_updur{d}) & ~isnan(corresp_lfp_updur{d}));
    up_diff_test(d) = signrank(mp_updur{d}(uset),corresp_lfp_updur{d}(uset));
    up_med_diff(d) = nanmedian(diff_up_dur{d});
    up_med_rel(d) = nanmedian(rel_up_dur{d});
    
    %these are the MP up states which were found to be shorter than one UDS
    %cycle
    npersd_across_ids{d} = find(mp_updur_lfpc_delay{d} < 1);
    npers_across_ids{d} = find(mp_updur_lfpc{d} < 1);
    
    %compute the fraction of MP UP states which are last more than the
    %corresponding LFP up state plus at least one LFP down state lasting >
    %0.5s (MP up states which get excluded for having too short cortical
    %DOWN states get excluded from the fraction calculation)
    rpersd_across_ids{d} = find(corresp_lfp_dmaxskippeddowndur{d} > 0.5);
    rpersd_across(d) = length(rpersd_across_ids{d})/(length(rpersd_across_ids{d}) + length(npersd_across_ids{d}));
    rpers_across_ids{d} = find(corresp_lfp_maxskippeddowndur{d} > 0.5);
    rpers_across(d) = length(rpers_across_ids{d})/(length(rpers_across_ids{d}) + length(npers_across_ids{d}));
    
    persd_within(d) = sum(mp_updur_lfpc_delay{d} > lf8_up_fract)/sum(~isnan(mp_updur_lfpc_delay{d}));
    persd_across_ids{d} = find(mp_updur_lfpc_delay{d} > 1);
    persd_across(d) = length(persd_across_ids{d})/sum(~isnan(mp_updur_lfpc_delay{d}));
    persd_within_np(d) = sum(mp_updur_lfpc_delay{d} > lf8_up_fract & mp_updur_lfpc_delay{d} < 1)...
        /sum(mp_updur_lfpc_delay{d} < 1);
    pers_within(d) = sum(mp_updur_lfpc{d} > lf8_up_fract)/sum(~isnan(mp_updur_lfpc{d}));
    pers_across(d) = sum(mp_updur_lfpc{d} > 1)/sum(~isnan(mp_updur_lfpc{d}));
    pers_within_np(d) = sum(mp_updur_lfpc{d} > lf8_up_fract & mp_updur_lfpc{d} < 1)...
        /sum(mp_updur_lfpc{d} < 1);

    
%% for up state skipping
    ndown_persd_across_ids{d} = find(mp_downdur_lfpc_delay{d} < 1);
    ndown_pers_across_ids{d} = find(mp_downdur_lfpc{d} < 1);

    rdown_persd_across_ids{d} = find(corresp_lfp_dmaxskippedupdur{d} > 0.5);
    rdown_persd_across(d) = length(rdown_persd_across_ids{d})/(length(rdown_persd_across_ids{d}) + length(ndown_persd_across_ids{d}));
    rdown_pers_across_ids{d} = find(corresp_lfp_maxskippedupdur{d} > 0.5);
    rdown_pers_across(d) = length(rdown_pers_across_ids{d})/(length(rdown_pers_across_ids{d}) + length(ndown_pers_across_ids{d}));
    
    down_persd_within(d) = sum(mp_downdur_lfpc_delay{d} > 0.5)/sum(~isnan(mp_downdur_lfpc_delay{d}));
    down_persd_across_ids{d} = find(mp_downdur_lfpc_delay{d} > 1);
    down_persd_across(d) = length(down_persd_across_ids{d})/sum(~isnan(mp_downdur_lfpc_delay{d}));
    down_pers_within(d) = sum(mp_downdur_lfpc{d} > 0.5)/sum(~isnan(mp_downdur_lfpc{d}));
    down_pers_across(d) = sum(mp_downdur_lfpc{d} > 1)/sum(~isnan(mp_downdur_lfpc{d}));

    
    clear lf8_period*
    
    cell_uhistd(d,:) = hist(mp_updur_lfpc_delay{d},hist_range);
    cell_uhistd(d,:) = cell_uhistd(d,:)/sum(cell_uhistd(d,:))/binsize;
    cell_uhistd(d,:) = jmm_smooth_1d_cor(cell_uhistd(d,:),25);
    cell_uhist(d,:) = hist(mp_updur_lfpc{d},hist_range);
    cell_uhist(d,:) = cell_uhist(d,:)/sum(cell_uhist(d,:))/binsize;
    cell_uhist(d,:) = jmm_smooth_1d_cor(cell_uhist(d,:),25);
    
    cell_dhistd(d,:) = hist(mp_downdur_lfpc_delay{d},hist_range);
    cell_dhistd(d,:) = cell_dhistd(d,:)/sum(cell_dhistd(d,:))/binsize;
    cell_dhistd(d,:) = jmm_smooth_1d_cor(cell_dhistd(d,:),25);
    cell_dhist(d,:) = hist(mp_downdur_lfpc{d},hist_range);
    cell_dhist(d,:) = cell_dhist(d,:)/sum(cell_dhist(d,:))/binsize;
    cell_dhist(d,:) = jmm_smooth_1d_cor(cell_dhist(d,:),25);
    
    cor_uhist(d,:) = hist(rel_up_dur{d},hist_range);
    cor_uhist(d,:) = cor_uhist(d,:)/sum(cor_uhist(d,:))/binsize;
    cor_uhist(d,:) = jmm_smooth_1d_cor(cor_uhist(d,:),25);
    cor_dhist(d,:) = hist(rel_down_dur{d},hist_range);
    cor_dhist(d,:) = cor_dhist(d,:)/sum(cor_dhist(d,:))/binsize;
    cor_dhist(d,:) = jmm_smooth_1d_cor(cor_dhist(d,:),25);
    
    diff_uhist(d,:) = hist(diff_up_dur{d},drange);
    diff_uhist(d,:) = diff_uhist(d,:)/sum(diff_uhist(d,:))/binsize;
    diff_uhist(d,:) = jmm_smooth_1d_cor(diff_uhist(d,:),25);
    diff_dhist(d,:) = hist(diff_down_dur{d},drange);
    diff_dhist(d,:) = diff_dhist(d,:)/sum(diff_dhist(d,:))/binsize;
    diff_dhist(d,:) = jmm_smooth_1d_cor(diff_dhist(d,:),25);
    
    %     y = hist(mp_updur_lfpc{d},x);
    %     sy = sum(y);
    
    %     eps = 1e-4;
    %     figure
    %     subplot(2,1,1)
    %     bar(x,(y/sy/xbinsize),1,'k')
    %     hold on
    %     temp = jmm_smooth_1d_cor(y,2);
    %     temp = temp/sy;
    %     plot(hist_range,cell_hist(d,:),'r','linewidth',2)
    %     xlim([0 7])
    %     subplot(2,1,2)
    %     plot(hist_range,cell_hist(d,:)+eps,'r','linewidth',2)
    %     xlim([0 7])
    %     set(gca,'yscale','log')
    %     ylim([1e-4 5])
    %     t_names = ['F:\WC_Germany\persistent_9_27_2010\quantization\' s_name];
    %     print('-dpng',t_names);
    %     close
    
    %     plot(lf8_sep{d},mp_updur_lfpc_delay{d},'.')
    %     pause
    %     clf
    
end


cd G:\WC_Germany\persistent_9_27_2010\
save pa_corresponding_lfp_state_data_downs pers* mp_up* mp_down* npers* lfp_* cell*hist* rpers* rdown* down* med*

%%
figure
mec_cells = 1:length(l3mec_p);
lec_cells = 23:36;
h = errorbar(hist_range,nanmean(cell_uhistd(mec_cells,:)),nanstd(cell_uhistd(mec_cells,:))/sqrt(length(mec_cells)));
errorbar_tick(h,.001,'units');
hold on
h = errorbar(hist_range,nanmean(cell_dhistd(mec_cells,:)),nanstd(cell_dhistd(mec_cells,:))/sqrt(length(mec_cells)),'r');
errorbar_tick(h,.001,'units');
h = errorbar(hist_range,nanmean(cell_uhistd(lec_cells,:)),nanstd(cell_uhistd(lec_cells,:))/sqrt(length(lec_cells)),'k');
errorbar_tick(h,.001,'units');
h = errorbar(hist_range,nanmean(cell_dhistd(lec_cells,:)),nanstd(cell_dhistd(lec_cells,:))/sqrt(length(lec_cells)),'g');
errorbar_tick(h,.001,'units');
xlim([0 4])
% ylim([0 1.6])
% line([0.5 0.5],[0 1.6],'color','k')
% line([1 1],[0 1.6],'color','k')
set(gca,'yscale','log')
line([0.5 0.5],[1e-3 2],'color','k')
line([1 1],[1e-3 2],'color','k')
ylim([1e-3 2])
legend('MEC UP','MEC DOWN','LEC UP','LEC DOWN')
xlabel('Duration (Ncx UDS Cycles)','Fontsize',14)
ylabel('Probability density','Fontsize',14)

figure
mec_cells = 1:length(l3mec_p);
lec_cells = 23:36;
h = errorbar(hist_range,nanmean(cor_uhist(mec_cells,:)),nanstd(cor_uhist(mec_cells,:))/sqrt(length(mec_cells)));
errorbar_tick(h,.001,'units');
hold on
h = errorbar(hist_range,nanmean(cor_dhist(mec_cells,:)),nanstd(cor_dhist(mec_cells,:))/sqrt(length(mec_cells)),'r');
errorbar_tick(h,.001,'units');
h = errorbar(hist_range,nanmean(cor_uhist(lec_cells,:)),nanstd(cor_uhist(lec_cells,:))/sqrt(length(lec_cells)),'k');
errorbar_tick(h,.001,'units');
h = errorbar(hist_range,nanmean(cor_dhist(lec_cells,:)),nanstd(cor_dhist(lec_cells,:))/sqrt(length(lec_cells)),'g');
errorbar_tick(h,.001,'units');
xlim([0 4])
ylim([0 1.2])
line([1 1],[0 1.6],'color','k')
legend('MEC UP','MEC DOWN','LEC UP','LEC DOWN')
xlabel('Duration (Ncx UDS Cycles)','Fontsize',14)
ylabel('Probability density','Fontsize',14)

figure
h = errorbar(drange,nanmean(diff_uhist(mec_cells,:)),nanstd(diff_uhist(mec_cells,:))/sqrt(length(mec_cells)));
errorbar_tick(h,.001,'units');
hold on
h = errorbar(drange,nanmean(diff_dhist(mec_cells,:)),nanstd(diff_dhist(mec_cells,:))/sqrt(length(mec_cells)),'r');
errorbar_tick(h,.001,'units');
h = errorbar(drange,nanmean(diff_uhist(lec_cells,:)),nanstd(diff_uhist(lec_cells,:))/sqrt(length(lec_cells)),'k');
errorbar_tick(h,.001,'units');
h = errorbar(drange,nanmean(diff_dhist(lec_cells,:)),nanstd(diff_dhist(lec_cells,:))/sqrt(length(lec_cells)),'g');
errorbar_tick(h,.001,'units');
xlim([-5 5])
ylim([0 1.6])
line([0 0],[0 1.6],'color','k')
legend('MEC UP','MEC DOWN','LEC UP','LEC DOWN')
xlabel('Duration (Ncx UDS Cycles)','Fontsize',14)
ylabel('Probability density','Fontsize',14)