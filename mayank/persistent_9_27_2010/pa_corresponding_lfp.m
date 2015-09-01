clear all

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
    lfp_up_lag{d} = nan(size(mp_up_inds{d})); %MP up-transition lag (relative to LFP up-transition)
    lfp_down_lag{d} = nan(size(mp_up_inds{d})); %MP down-transition lag (rel LFP down-trans)
    mp_updur_lfpc_delay{d} = nan(size(mp_up_inds{d})); %MP up state duration (in LFP cycles)
    mp_updur_lfpc{d} = nan(size(mp_up_inds{d})); %MP up state duration (in LFP cycles)
    corresp_lfp_updur{d} = nan(size(mp_up_inds{d})); %duration of the corresponding LFP up state
    corresp_lfp_downdur{d} = nan(size(mp_up_inds{d})); %duration of the corresponding LFP down state
    corresp_lfp_cycledur{d} = nan(size(mp_up_inds{d})); %duration of the corresponding LFP cycle
    corresp_lfp_dutycycle{d} = nan(size(mp_up_inds{d})); %corresponding LFP duty cycle
    corresp_lfp_upid{d} = nan(size(mp_up_inds{d}));
    corresp_lfp_maxskippeddowndur{d} = nan(size(mp_up_inds{d})); %minimum duration of any LFP down states skipped by the MP Up state
    corresp_lfp_mindowndur{d} = nan(size(mp_up_inds{d})); %minimum duration of any LFP down states skipped by the MP Up state
    corresp_lfp_dmaxskippeddowndur{d} = nan(size(mp_up_inds{d})); %minimum duration of any LFP down states skipped by the MP Up state
    corresp_lfp_dmindowndur{d} = nan(size(mp_up_inds{d})); %minimum duration of any LFP down states skipped by the MP Up state
    lfp_skipped_downs{d} = [];
    lfp_dskipped_downs{d} = [];
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
        
        for i = 1:length(cur_up_trans)
            
            %find nearest lfp up transition
            [dummy,near_lfp_up] = min(abs(cur_up_trans(i) - cur_up_trans8));
            
            %mp up state dur
            mp_updur{d}(cnt) = (cur_down_trans(i)-cur_up_trans(i))/Fsd;
            
            %store the LFP up-transition lag in sec
            cur_lfp_up_lag = (cur_up_trans(i)-cur_up_trans8(near_lfp_up));
            lfp_up_lag{d}(cnt) = cur_lfp_up_lag/Fsd;
            
            %store mp up state dur in lfp cycles
            cur_lfpc = lf8_period_f{ns}(cur_down_trans(i) - cur_lfp_up_lag) ...
                - lf8_period_f{ns}(cur_up_trans(i)-cur_lfp_up_lag);
            mp_updur_lfpc_delay{d}(cnt) = cur_lfpc;
            cur_lfpc = lf8_period_f{ns}(cur_down_trans(i)) - lf8_period_f{ns}(cur_up_trans(i));
            mp_updur_lfpc{d}(cnt) = cur_lfpc;
            
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
                lfp_down_lag{d}(cnt) = cur_lfp_down_lag/Fsd;
            end
            
            corresp_lfp_updur{d}(cnt) = (cur_down_trans8(near_lfp_up)-cur_up_trans8(near_lfp_up))/Fsd;
            temp = find(lf8_up_inds{d} == cur_up_trans8(near_lfp_up) + new_seg_inds(ns)-1);
            if ~isempty(temp)
                corresp_lfp_upid{d}(cnt) = temp;
            end
            if n_lfp_ups > corresp_lfp_down+1 %if the nearest up transition is not the final LFP up-transition
                corresp_lfp_downdur{d}(cnt) = (cur_up_trans8(corresp_lfp_down+1)-cur_down_trans8(corresp_lfp_down))/Fsd;
%                 corresp_lfp_cycledur{d}(cnt) = (cur_up_trans8(near_lfp_up+1)-cur_up_trans8(near_lfp_up))/Fsd;
%                 corresp_lfp_dutycycle{d}(cnt) = corresp_lfp_updur{d}(cnt)/corresp_lfp_cycledur{d}(cnt);
            end
            skipped_lfp_downs = find(cur_down_trans8(1:end-1) > cur_up_trans8(near_lfp_up) & cur_down_trans8(1:end-1) < cur_down_trans(i) ...
                & cur_up_trans8(2:end) < cur_down_trans(i));
            dskipped_lfp_downs = find(cur_down_trans8(1:end-1) > cur_up_trans8(near_lfp_up) & cur_down_trans8(1:end-1) < cur_down_trans(i)-cur_lfp_up_lag ...
                & cur_up_trans8(2:end) < cur_down_trans(i)-cur_lfp_up_lag);
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
            
            t1(d) = sum(mp_updur{d} > corresp_lfp_updur{d})/length(mp_updur{d});
            t2(d) = sum(mp_updur{d} > corresp_lfp_updur{d} + corresp_lfp_downdur{d})/length(mp_updur{d});
            
            cnt = cnt + 1;
        end
    end
    
%     
%     is_lf8_down_skipped = zeros(size(lfp_state_durations{1}));
%     for i = 1:length(is_lf8_down_skipped)-1
%         curset = lf8_down_inds{d}(i):lf8_up_inds{d}(i+1);
%         if min(mp_state_vec(curset)) == 1
%             is_lf8_down_skipped(i) = 1;
%         elseif max(mp_state_vec(curset)) == 0
%             is_lf8_down_skipped(i) = 0;
%         end
%     end


    lfp_dskipped_downs{d}(lfp_state_durations{1}(lfp_dskipped_downs{d}) < 0.5) = [];
    lfp_ndskipped_downs{d} = setdiff(1:length(lf8_down_inds{d}),lfp_dskipped_downs{d});
    lfp_ndskipped_downs{d}(lfp_state_durations{1}(lfp_ndskipped_downs{d}) < 0.5) = [];

    lfp_skipped_downs{d}(lfp_state_durations{1}(lfp_skipped_downs{d}) < 0.5) = [];
    lfp_nskipped_downs{d} = setdiff(1:length(lf8_down_inds{d}),lfp_skipped_downs{d});
    lfp_nskipped_downs{d}(lfp_state_durations{1}(lfp_nskipped_downs{d}) < 0.5) = [];

    skipped_down_durs(d,:) = histc(lfp_state_durations{1}(lfp_dskipped_downs{d}),log_grid);
    skipped_down_durs(d,:) = skipped_down_durs(d,:)/sum(skipped_down_durs(d,:));
    nonskipped_down_durs(d,:) = histc(lfp_state_durations{1}(lfp_ndskipped_downs{d}),log_grid);
    nonskipped_down_durs(d,:) = nonskipped_down_durs(d,:)/sum(nonskipped_down_durs(d,:));
    ov_lfp_down_durs(d,:) = histc(lfp_state_durations{1},log_grid);
    ov_lfp_down_durs(d,:) = ov_lfp_down_durs(d,:)/sum(ov_lfp_down_durs(d,:));
    
    
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
    
    
    rpersd_lfp_updurs(d,:) = histc(corresp_lfp_updur{d}(rpersd_across_ids{d}),log_grid);
    rpersd_lfp_updurs(d,:) = rpersd_lfp_updurs(d,:)/sum(rpersd_lfp_updurs(d,:));
    persd_lfp_updurs(d,:) = histc(corresp_lfp_updur{d}(persd_across_ids{d}),log_grid);
    persd_lfp_updurs(d,:) = persd_lfp_updurs(d,:)/sum(persd_lfp_updurs(d,:));
    npersd_lfp_updurs(d,:) = histc(corresp_lfp_updur{d}(npersd_across_ids{d}),log_grid);
    npersd_lfp_updurs(d,:) = npersd_lfp_updurs(d,:)/sum(npersd_lfp_updurs(d,:));
    
    avg_persd_dur(d) = nanmean(mp_updur{d}(mp_updur_lfpc_delay{d} > 1));
    avg_pers_dur(d) = nanmean(mp_updur{d}(mp_updur_lfpc{d} > 1));
    avg_rpersd_dur(d) = nanmean(mp_updur{d}(rpersd_across_ids{d}));
    avg_rpers_dur(d) = nanmean(mp_updur{d}(rpersd_across_ids{d}));
    
    mean_up_lag(d) = nanmean(lfp_up_lag{d});
    mean_down_lag(d) = nanmean(lfp_down_lag{d});
    median_up_lag(d) = nanmedian(lfp_up_lag{d});
    median_down_lag(d) = nanmedian(lfp_down_lag{d});
    
    clear lf8_period*
    
    cell_histd(d,:) = hist(mp_updur_lfpc_delay{d},hist_range);
    cell_histd(d,:) = cell_histd(d,:)/sum(cell_histd(d,:))/binsize;
    cell_histd(d,:) = jmm_smooth_1d_cor(cell_histd(d,:),25);
    cell_hist(d,:) = hist(mp_updur_lfpc{d},hist_range);
    cell_hist(d,:) = cell_hist(d,:)/sum(cell_hist(d,:))/binsize;
    cell_hist(d,:) = jmm_smooth_1d_cor(cell_hist(d,:),25);
    
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
save pa_corresponding_lfp_state_data_rtest pers* mean*lag median*lag lfp*lag mp_up* ...
    mp_down* corresp* rpers* skipped* nonskipped* ov_lfp* npers* avg_* lfp_*

%%
figure
mec_cells = 1:length(l3mec_p);
h = errorbar(log_grid,nanmean(skipped_down_durs(mec_cells,:)),nanstd(skipped_down_durs(mec_cells,:))/sqrt(length(mec_cells)));
hold on
h = errorbar(log_grid,nanmean(nonskipped_down_durs(mec_cells,:)),nanstd(nonskipped_down_durs(mec_cells,:))/sqrt(length(mec_cells)),'r');
h = errorbar(log_grid,nanmean(ov_lfp_down_durs(mec_cells,:)),nanstd(ov_lfp_down_durs(mec_cells,:))/sqrt(length(mec_cells)),'k');
set(gca,'yscale','log')
ylim([1e-3 0.2])
xlim([0 7])
line([0.5 0.5],[1e-3 0.2])

figure
h = errorbar(log_grid,nanmean(rpersd_lfp_updurs(mec_cells,:)),nanstd(rpersd_lfp_updurs(mec_cells,:))/sqrt(length(mec_cells)));
hold on
h = errorbar(log_grid,nanmean(persd_lfp_updurs(mec_cells,:)),nanstd(persd_lfp_updurs(mec_cells,:))/sqrt(length(mec_cells)),'r');
h = errorbar(log_grid,nanmean(npersd_lfp_updurs(mec_cells,:)),nanstd(npersd_lfp_updurs(mec_cells,:))/sqrt(length(mec_cells)),'k');
set(gca,'yscale','log')
ylim([1e-3 0.2])
xlim([0 7])
line([0.5 0.5],[1e-3 0.2],'color','k')
