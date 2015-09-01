clear all

load C:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('C:\WC_Germany\parietal_cortical_2010\')
addpath('C:\WC_Germany\persistent_9_27_2010\')
addpath('C:\WC_Germany\hsmm_state_detection\')

cd C:\WC_Germany\persistent_9_27_2010\
used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

drive_letter = 'C';
dsf = 8;
Fsd = 2016/dsf;

thresh_lf8_updur = 0.5;
thresh_lf8_downdur = 0.5;

for d = 1:length(sess_data)
%     d=36
    fprintf('session %d\n',d)
    cdir = sess_data(d).directory;
    cdir(1) = drive_letter;
    cd(cdir)
    s_name = strcat(sess_data(d).region,'_l',sess_data(d).layer,'_',sess_data(d).name);
    
    load ./used_data wcv
    wcv_d = downsample(wcv,dsf);
    datalen = length(wcv_d);
    load ./pa_hsmm_state_seq_new2.mat
    mp_state_seq = hsmm_bbstate_seq;
    load ./pa_hsmm_state_seq8_new2.mat
    lfp_state_seq = hsmm_bbstate_seq8;
    
    t_axis = (1:length(wcv_d))/Fsd;
    rec_dur(d) = t_axis(end);
    
    load ./pa_lf8_period_data_new2.mat
    
%     [new_seg_inds] = resample_uds_seg_inds(hmm.UDS_segs,hmm.Fs,Fsd,length(wcv_d));
    [new_seg_inds] = resample_uds_seg_inds_v2(hmm.UDS_segs,hmm.Fs,Fsd,mp_state_seq);
    seg_durs = new_seg_inds(:,2)-new_seg_inds(:,1)+1;
    
%     if you want to flip the roles of the MP and LFP
%         temp = lfp_state_seq;
%         lfp_state_seq = mp_state_seq;
%         mp_state_seq = temp;
    
    [mp_uptrans{d},mp_downtrans{d},mp_upsegnums{d},mp_downsegnums{d}] = compute_state_transitions_seg(new_seg_inds,mp_state_seq);
    [lf8_uptrans{d},lf8_downtrans{d},lf8_upsegnums{d},lf8_downsegnums{d}] = compute_state_transitions_seg(new_seg_inds,lfp_state_seq);
    
    %     %if you want to switch the up and down states
    %     temp = mp_uptrans{d};
    %     temp2 = lf8_uptrans{d};
    %     mp_uptrans{d} = mp_downtrans{d};
    %     lf8_uptrans{d} = lf8_downtrans{d};
    %     mp_downtrans{d} = temp;
    %     lf8_downtrans{d} = temp2;
    %     mp_downtrans{d}(1) = [];
    %     lf8_downtrans{d}(1) = [];
    %     mp_uptrans{d}(end) = [];
    %     lf8_uptrans{d}(end) = [];
    
    lf8_period_vec = nan(size(wcv_d));
    for i = 1:size(new_seg_inds,1)
        lf8_period_vec(new_seg_inds(i,1):new_seg_inds(i,2)) = lf8_period_f{i};
    end
    
    mp_state_vec = zeros(size(wcv_d));
    lf8_state_vec = zeros(size(wcv_d));
    for i = 1:length(mp_uptrans{d})
        mp_state_vec(mp_uptrans{d}(i):mp_downtrans{d}(i)) = 1;
    end
    for i = 1:length(lf8_uptrans{d})
        lf8_state_vec(lf8_uptrans{d}(i):lf8_downtrans{d}(i)) = 1;
    end
    uds_vec = zeros(size(wcv_d));
    for i = 1:size(new_seg_inds,1)
        uds_vec(new_seg_inds(i,1):new_seg_inds(i,2)) = 1;
    end
    mp_state_vec(uds_vec==0) = nan;
    lf8_state_vec(uds_vec==0) = nan;
    
    %     % compute state durations
    [lfp_state_durations] = compute_state_durations_seg(lfp_state_seq,Fsd);
    [mp_state_durations] = compute_state_durations_seg(mp_state_seq,Fsd);
    %
    mp_updurs{d} = mp_state_durations{2};
    mp_downdurs{d} = mp_state_durations{1};
    lf8_updurs{d} = lfp_state_durations{2};
    lf8_downdurs{d} = lfp_state_durations{1};

    mp_dutycycs{d} = mp_updurs{d}./(mp_updurs{d}+mp_downdurs{d});
    lf8_dutycycs{d} = lf8_updurs{d}./(lf8_updurs{d}+lf8_downdurs{d});
    
    %     if length(lfp_state_durations{1}) ~= length(lf8_uptrans{d})
    %         error('state mismatch')
    %     end
    %     if length(lfp_state_durations{2}) ~= length(lf8_uptrans{d})
    %         error('state mismatch')
    %     end
    
    %assign 1-1 corresponding state transitions
    % %     [corresp_lf8_upinds{d},corresp_lf8_downinds{d}] = find_corresponding_ncx_state_transitions(...
    % %         mp_uptrans{d},mp_upsegnums{d},mp_downtrans{d},mp_downsegnums{d},lf8_uptrans{d},...
    % %         lf8_upsegnums{d},lf8_downtrans{d},lf8_downsegnums{d});
    % [corresp_lf8_upinds{d},corresp_lf8_downinds{d}] = greedy_find_corresponding_ncx_state_transitions(...
    %     mp_uptrans{d},mp_downtrans{d},lf8_uptrans{d},lf8_downtrans{d});
    [corresp_lf8_upinds{d},corresp_lf8_downinds{d}] = greedy_find_corresponding_ncx_state_transitions_simp(...
        mp_uptrans{d},mp_downtrans{d},lf8_uptrans{d},lf8_downtrans{d});
    
    % %     [corresp_mp_upinds{d},corresp_mp_downinds{d}] = find_corresponding_ncx_state_transitions(...
    % %         lf8_uptrans{d},lf8_upsegnums{d},lf8_downtrans{d},lf8_downsegnums{d},mp_uptrans{d},...
    % %         mp_upsegnums{d},mp_downtrans{d},mp_downsegnums{d});
    % [corresp_mp_upinds{d},corresp_mp_downinds{d}] = greedy_find_corresponding_ncx_state_transitions(...
    %     lf8_uptrans{d},lf8_downtrans{d},mp_uptrans{d},mp_downtrans{d});
    %%
    %find non-skipped states
    non_skipped_mp_up_states{d} = find(~isnan(corresp_lf8_upinds{d}));
    non_skipped_mp_down_states{d} = find(~isnan(corresp_lf8_downinds{d}));
        
    lf8_uplags{d} = nan(size(mp_uptrans{d}));
    lf8_downlags{d} = nan(size(mp_uptrans{d}));
    
    %compute transition lags for non-skipped states
    lf8_uplags{d}(non_skipped_mp_up_states{d}) = mp_uptrans{d}(non_skipped_mp_up_states{d}) - lf8_uptrans{d}(corresp_lf8_upinds{d}(non_skipped_mp_up_states{d}));
    lf8_downlags{d}(non_skipped_mp_down_states{d}) = mp_downtrans{d}(non_skipped_mp_down_states{d}) - lf8_downtrans{d}(corresp_lf8_downinds{d}(non_skipped_mp_down_states{d}));
    
    lf8_reldownlags{d} = nan(size(mp_uptrans{d}));
    lf8_reluplags{d} = nan(size(mp_uptrans{d}));
    for i = 1:length(mp_uptrans{d})
        if ~isnan(corresp_lf8_downinds{d}(i))
            lf8_reldownlags{d}(i) = lf8_downlags{d}(i)/lf8_downdurs{d}(corresp_lf8_downinds{d}(i))/Fsd;
        end
        if ~isnan(corresp_lf8_upinds{d}(i))
            lf8_reluplags{d}(i) = lf8_uplags{d}(i)/lf8_updurs{d}(corresp_lf8_upinds{d}(i))/Fsd;
        end
    end
    
    %compute durations in units of Ncx UDS cycles
    [mp_updurs_lfpc{d},mp_downdurs_lfpc{d}] = find_duration_ncx_uds_cycles(mp_uptrans{d},mp_downtrans{d},...
        lf8_uplags{d},lf8_downlags{d},lf8_period_vec);
    
    % %     [mp_upskipped{d},mp_downskipped{d}] = find_skipped_ncx_states(mp_uptrans{d},...
    % %         mp_downtrans{d},lf8_uptrans{d},lf8_downtrans{d},corresp_lf8_upinds{d},corresp_lf8_downinds{d},...
    % %         lf8_uplags{d},lf8_downlags{d},lfp_state_durations{2},lfp_state_durations{1},thresh_lf8_updur,thresh_lf8_downdur);
    [mp_upskipped{d},mp_downskipped{d}] = greedy_find_skipped_ncx_states(...
        corresp_lf8_upinds{d},corresp_lf8_downinds{d},lfp_state_durations{2},lfp_state_durations{1},thresh_lf8_downdur,thresh_lf8_updur);
            
    rt2_ups{d} = find(mp_upskipped{d}.rnum_skipped > 0);
    t2_ups{d} = find(mp_upskipped{d}.num_skipped > 0);
    nrt2_ups{d} = find(mp_upskipped{d}.num_skipped == 0);
    fract_rt2_ups(d) = length(rt2_ups{d})/(length(rt2_ups{d}) + length(nrt2_ups{d}));
    fract_t2_ups(d) = length(t2_ups{d})/(length(t2_ups{d}) + length(nrt2_ups{d}));
      
    mp_corresp_lf8_mindur{d} = [mp_upskipped{d}.min_dur];
    mp_corresp_lf8_mindur{d}(nrt2_ups{d}) = lf8_updurs{d}(corresp_lf8_upinds{d}(nrt2_ups{d}));
    
    rt2_downs{d} = find(mp_downskipped{d}.rnum_skipped > 0);
    nrt2_downs{d} = find(mp_downskipped{d}.num_skipped==0); %number of down states that didn't skip any LFP up states
    fract_rt2_downs(d) = length(rt2_downs{d})/(length(rt2_downs{d}) + length(nrt2_downs{d}));
    
    lf8_uplags_t2{d} = lf8_uplags{d}(rt2_ups{d});
    lf8_uplags_nt2{d} = lf8_uplags{d}(nrt2_ups{d});
    lf8_downlags_t2{d} = lf8_downlags{d}(rt2_downs{d});
    lf8_downlags_nt2{d} = lf8_downlags{d}(nrt2_downs{d});
    
    robust_non_skipped_mp_ups = [rt2_ups{d}; nrt2_ups{d}];
    robust_non_skipped_mp_downs = [rt2_downs{d}; nrt2_downs{d}];
    
    mp_rel_updurs{d} = nan(size(mp_uptrans{d}));
    mp_rel_downdurs{d} = nan(size(mp_uptrans{d}));
    mp_rel_updurs{d}(robust_non_skipped_mp_ups) = mp_state_durations{2}(robust_non_skipped_mp_ups) - ...
        lfp_state_durations{2}(corresp_lf8_upinds{d}(robust_non_skipped_mp_ups));
    mp_rel_downdurs{d}(robust_non_skipped_mp_downs) = mp_state_durations{1}(robust_non_skipped_mp_downs) - ...
        lfp_state_durations{1}(corresp_lf8_downinds{d}(robust_non_skipped_mp_downs));
    mp_rel_updurs_t2{d} = mp_rel_updurs{d}(rt2_ups{d});
    mp_rel_updurs_nt2{d} = mp_rel_updurs{d}(nrt2_ups{d});
    mp_rel_downdurs_t2{d} = mp_rel_downdurs{d}(rt2_downs{d});
    mp_rel_downdurs_nt2{d} = mp_rel_downdurs{d}(nrt2_downs{d});
    
    fract_rt1_ups(d) = nansum(mp_rel_updurs{d} > 0)/sum(~isnan(mp_rel_updurs{d}));
    fract_rt1_ups_nt2(d) = nansum(mp_rel_updurs_nt2{d} > 0)/sum(~isnan(mp_rel_updurs_nt2{d}));
    fract_rt1_downs(d) = nansum(mp_rel_downdurs{d} > 0)/sum(~isnan(mp_rel_downdurs{d}));
    
    mp_rel_updur_skew(d) = skewness(mp_rel_updurs{d});
    mp_rel_downdur_skew(d) = skewness(mp_rel_downdurs{d});
    mp_rel_updur_nt2_skew(d) = skewness(mp_rel_updurs_nt2{d});
    mp_rel_downdur_nt2_skew(d) = skewness(mp_rel_downdurs_nt2{d});
    
    median_uplag(d) = nanmedian(lf8_uplags{d});
    median_uplag_t2(d) = nanmedian(lf8_uplags{d}(rt2_ups{d}));
    median_uplag_nt2(d) = nanmedian(lf8_uplags{d}(nrt2_ups{d}));
    median_downlag(d) = nanmedian(lf8_downlags{d});
    median_downlag_t2(d) = nanmedian(lf8_downlags{d}(rt2_ups{d}));
    median_downlag_nt2(d) = nanmedian(lf8_downlags{d}(nrt2_ups{d}));
    mean_uplag(d) = nanmean(lf8_uplags{d});
    mean_downlag(d) = nanmean(lf8_downlags{d});
    median_reluplag(d) = nanmedian(lf8_reluplags{d});
    median_reldownlag(d) = nanmedian(lf8_reldownlags{d});
    mean_reluplag(d) = nanmean(lf8_reluplags{d});
    mean_reldownlag(d) = nanmean(lf8_reldownlags{d});

    median_updur(d) = nanmedian(mp_updurs{d});
    median_downdur(d) = nanmedian(mp_downdurs{d});
    median_updur8(d) = nanmedian(lf8_updurs{d});
    median_downdur8(d) = nanmedian(lf8_downdurs{d});
    median_dc(d) = nanmedian(mp_dutycycs{d});
    mean_dc(d) = nanmean(mp_dutycycs{d});
    mean_dc8(d) = nanmean(lf8_dutycycs{d});
    median_dc8(d) = nanmedian(lf8_dutycycs{d});
    
    mean_rel_updur_nt2(d) = nanmean(mp_rel_updurs_nt2{d});
    median_rel_updur_nt2(d) = nanmedian(mp_rel_updurs_nt2{d});
    mean_rel_downdur_nt2(d) = nanmean(mp_rel_downdurs_nt2{d});
    median_rel_downdur_nt2(d) = nanmedian(mp_rel_downdurs_nt2{d});
    
%     hist(lf8_downlags{d}/Fsd,50)
%     median_downlag(d)/Fsd
%     pause
%     clf
    
    
    %%
% %     load G:\WC_Germany\persistent_9_27_2010\pa_corresponding_lfp_state_data_downs
%     win = 10;
%     plot(t_axis,mp_state_vec)
%     hold on
%     plot(t_axis,lf8_state_vec+1.1,'r')
%     ylim([-0.2 2.2])
%     plot(t_axis(mp_uptrans{d}(rt2_ups{d})),0.2*ones(size(rt2_ups{d})),'k*')
%     plot(t_axis(mp_downtrans{d}(rt2_downs{d})),0.2*ones(size(rt2_downs{d})),'g*')
%     
%     plot(t_axis(mp_uptrans{d}(robust_non_skipped_mp_ups)),0.4*ones(size(robust_non_skipped_mp_ups)),'ko')
%     plot(t_axis(mp_downtrans{d}(robust_non_skipped_mp_downs)),0.4*ones(size(robust_non_skipped_mp_downs)),'go')
%     plot(t_axis(mp_uptrans{d}(rt2_ups{d})),0.4*ones(size(rt2_ups{d})),'k*')
%     plot(t_axis(mp_downtrans{d}(rt2_downs{d})),0.4*ones(size(rt2_downs{d})),'g*')
%     plot(t_axis(mp_uptrans{d}(nrt2_ups{d})),0.1*ones(size(nrt2_ups{d})),'k.')
%     plot(t_axis(mp_downtrans{d}(nrt2_downs{d})),0.1*ones(size(nrt2_downs{d})),'g.')
%     for i = 1:length(corresp_lf8_upinds{d})
%         plot(t_axis(mp_uptrans{d}(i)),0.5,'ko','linewidth',2)
%         plot(t_axis(mp_downtrans{d}(i)),0.5,'go','linewidth',2)
%         if ~isnan(corresp_lf8_upinds{d}(i))
%             plot(t_axis(lf8_uptrans{d}(corresp_lf8_upinds{d}(i))),1.5,'ko')
%             line([t_axis(mp_uptrans{d}(i)) t_axis(lf8_uptrans{d}(corresp_lf8_upinds{d}(i)))],[0.5 1.5],'color','k','linewidth',2)
%         end
%         if ~isnan(corresp_lf8_downinds{d}(i))
%             plot(t_axis(lf8_downtrans{d}(corresp_lf8_downinds{d}(i))),1.5,'go','linewidth',2)
%             line([t_axis(mp_downtrans{d}(i)) t_axis(lf8_downtrans{d}(corresp_lf8_downinds{d}(i)))],[0.5 1.5],'color','g','linewidth',2)
%         end
%         xlim([t_axis(mp_uptrans{d}(i))-win t_axis(mp_uptrans{d}(i))+win])
%         mp_updurs_lfpc{d}(i)
%         pause
%         plot(t_axis(mp_uptrans{d}(i)),0.5,'co')
%         plot(t_axis(mp_downtrans{d}(i)),0.5,'co')
%         if ~isnan(corresp_lf8_upinds{d}(i))
%         plot(t_axis(lf8_uptrans{d}(corresp_lf8_upinds{d}(i))),1.5,'co')
%         end
%         if ~isnan(corresp_lf8_downinds{d}(i))
%         plot(t_axis(lf8_downtrans{d}(corresp_lf8_downinds{d}(i))),1.5,'co')
%         end
%     end
    
end

mec_cells = 1:length(l3mec_p);
lec_cells = 23:36;
cd C:\WC_Germany\persistent_9_27_2010\
save pa_corresponding_lfp_revised_12_2_2011 lf8_* mp_* fract_* corresp_* mean_* median_* rt2* t2* nrt2*

%%
[dlag_mec_c,dlag_mec_p] = corrcoef(median_downlag(mec_cells),fract_rt2_ups(mec_cells));
[dlag_lec_c,dlag_lec_p] = corrcoef(median_downlag(lec_cells),fract_rt2_ups(lec_cells));
uset = find(fract_rt2_ups(mec_cells) > 0);
[dlag_mec_cl,dlag_mec_pl] = corrcoef(median_downlag(mec_cells(uset)),log(fract_rt2_ups(mec_cells(uset))));
uset = find(fract_rt2_ups(lec_cells) > 0);
[dlag_lec_cl,dlag_lec_pl] = corrcoef(median_downlag(lec_cells(uset)),log(fract_rt2_ups(lec_cells(uset))));
[ulag_mec_c,ulag_mec_p] = corrcoef(median_uplag(mec_cells),fract_rt2_downs(mec_cells));
[ulag_lec_c,ulag_lec_p] = corrcoef(median_uplag(lec_cells),fract_rt2_downs(lec_cells));
uset = find(fract_rt2_downs(mec_cells) > 0);
[ulag_mec_cl,ulag_mec_pl] = corrcoef(median_uplag(mec_cells(uset)),log(fract_rt2_downs(mec_cells(uset))));
uset = find(fract_rt2_downs(lec_cells) > 0);
[ulag_lec_cl,ulag_lec_pl] = corrcoef(median_uplag(lec_cells(uset)),log(fract_rt2_downs(lec_cells(uset))));
   
eps = 0;
figure
plot(median_downlag(mec_cells)/Fsd,fract_rt2_ups(mec_cells)+eps,'o','markersize',8)
hold on
plot(median_downlag(lec_cells)/Fsd,fract_rt2_ups(lec_cells)+eps,'ro','markersize',8)
used = find(fract_rt2_ups + eps > 0);
temp_all = polyfit(median_downlag(used)/Fsd,log(fract_rt2_ups(used)+eps),1);
temp_mec = polyfit(median_downlag(mec_cells)/Fsd,log(fract_rt2_ups(mec_cells)+eps),1);
temp_lec = polyfit(median_downlag(lec_cells)/Fsd,log(fract_rt2_ups(lec_cells)+eps),1);
x_ax = linspace(-0.2,1.25,100);
% plot(x_ax,polyval(temp_mec,x_ax),'b')
% plot(x_ax,polyval(temp_lec,x_ax),'r')
set(gca,'fontsize',14,'fontname','arial')
plot(x_ax,exp(polyval(temp_all,x_ax)),'k')
xlabel('Median down-transition lag (s)','fontsize',14,'fontname','arial')
ylabel('Log probability persistent up state','fontsize',14,'fontname','arial')
set(gca,'yscale','log')
xlim([-0.2 1.25])
ylim([1e-3 1])

% figure
% plot(median_uplag(mec_cells)/Fsd,log(fract_rt2_ups(mec_cells)+eps),'o')
% hold on
% plot(median_uplag(lec_cells)/Fsd,log(fract_rt2_ups(lec_cells)+eps),'ro')
% temp_mec = polyfit(median_uplag(mec_cells)/Fsd,log(fract_rt2_ups(mec_cells)+eps),1);
% temp_lec = polyfit(median_uplag(lec_cells)/Fsd,log(fract_rt2_ups(lec_cells)+eps),1);
% x_ax = linspace(0,0.5,100);
% plot(x_ax,polyval(temp_mec,x_ax),'b')
% plot(x_ax,polyval(temp_lec,x_ax),'r')
% xlabel('Median up-transition lag (s)','fontsize',14)
% ylabel('Log probability persistent up state','fontsize',14)
% ylim([-7.5 0])

figure
set(gca,'fontsize',14,'fontname','arial')
plot(median_uplag(mec_cells)/Fsd,fract_rt2_downs(mec_cells)+eps,'o','markersize',8)
hold on
plot(median_uplag(lec_cells)/Fsd,fract_rt2_downs(lec_cells)+eps,'ro','markersize',8)
used = find(fract_rt2_downs+eps > 0);
temp_all = polyfit(median_uplag(used)/Fsd,log(fract_rt2_downs(used)+eps),1);
temp_mec = polyfit(median_uplag(mec_cells)/Fsd,log(fract_rt2_downs(mec_cells)+eps),1);
temp_lec = polyfit(median_uplag(lec_cells)/Fsd,log(fract_rt2_downs(lec_cells)+eps),1);
x_ax = linspace(0,0.45,100);
% plot(x_ax,polyval(temp_mec,x_ax),'b')
% plot(x_ax,polyval(temp_lec,x_ax),'r')
plot(x_ax,exp(polyval(temp_all,x_ax)),'k')
xlabel('Median up-transition lag (s)','fontsize',14,'fontname','arial')
ylabel('Log probability persistent down state','fontsize',14,'fontname','arial')
set(gca,'yscale','log')
xlim([0 0.45])
ylim([1e-3 1])

% figure
% plot(median_downlag(mec_cells)/Fsd,log(fract_rt2_downs(mec_cells)+eps),'o')
% hold on
% plot(median_downlag(lec_cells)/Fsd,log(fract_rt2_downs(lec_cells)+eps),'ro')
% temp_mec = polyfit(median_downlag(mec_cells)/Fsd,log(fract_rt2_downs(mec_cells)+eps),1);
% temp_lec = polyfit(median_downlag(lec_cells)/Fsd,log(fract_rt2_downs(lec_cells)+eps),1);
% x_ax = linspace(-0.2,1.2,100);
% plot(x_ax,polyval(temp_mec,x_ax),'b')
% plot(x_ax,polyval(temp_lec,x_ax),'r')
% xlabel('Median down-transition lag (s)','fontsize',14)
% ylabel('Log probability persistent down state','fontsize',14)
% ylim([-7.5 0])

%%
eps = 0;
um = mec_cells(fract_rt2_ups(mec_cells) > 0);
ul = lec_cells(fract_rt2_ups(lec_cells) > 0);
ua = fract_rt2_ups > 0;

sf = log(fract_rt2_ups./(1-fract_rt2_ups));
temp_m = polyfit(median_downlag(um)/Fsd,sf(um),1);
temp_l = polyfit(median_downlag(ul)/Fsd,sf(ul),1);
temp_o = polyfit(median_downlag(ua)/Fsd,sf(ua),1);
x_ax = linspace(-0.2,1.3,100);
ev_m = 1./(1+exp(-polyval(temp_m,x_ax)));
ev_l = 1./(1+exp(-polyval(temp_l,x_ax)));
ev_o = 1./(1+exp(-polyval(temp_o,x_ax)));

figure
set(gca,'fontname','arial')
plot(median_downlag(um)/Fsd,fract_rt2_ups(um),'ro','markersize',8)
hold on
plot(median_downlag(ul)/Fsd,fract_rt2_ups(ul),'o','markersize',8)
% plot(x_ax,ev_m,'r')
% plot(x_ax,ev_l,'b')
plot(x_ax,ev_o,'k')
ylim([4e-3 0.6])
xlabel('Median down lag (s)','fontname','arial')
ylabel('Prob T2 up','fontname','arial')
xlim([-0.2 1.3])

set(gca,'yscale','log')
shg

%%
temp_m = polyfit(median_downlag(mec_cells)/Fsd,fract_rt2_ups(mec_cells),1);
temp_l = polyfit(median_downlag(lec_cells)/Fsd,fract_rt2_ups(lec_cells),1);
x_ax1 = linspace(-0.2,0.25,100);
x_ax2 = linspace(0.25,1.25,100);
ev_m = polyval(temp_m,x_ax2);
ev_l = polyval(temp_l,x_ax1);
figure
set(gca,'fontname','arial')
plot(median_downlag(mec_cells)/Fsd,fract_rt2_ups(mec_cells),'ro','markersize',8)
hold on
plot(median_downlag(lec_cells)/Fsd,fract_rt2_ups(lec_cells),'o','markersize',8)
plot(x_ax2,ev_m,'r')
plot(x_ax1,ev_l,'b')
xlim([-0.2 1.3])
ylim([0 0.5])
xlabel('Median down lag (s)','fontname','arial')
ylabel('Prob T2 up','fontname','arial')

%%
figure
set(gca,'fontsize',14)
set(gca,'fontname','arial')
plot(100*fract_rt2_ups(mec_cells),100*fract_rt2_downs(mec_cells),'ro','markersize',8)
hold on
plot(100*fract_rt2_ups(lec_cells),100*fract_rt2_downs(lec_cells),'o','markersize',8)
xlim([0 0.6]*100)
ylim([0 0.6]*100)
xlabel('Percent persistent Ups','fontsize',16,'fontname','arial')
ylabel('Percent persistent Downs','fontsize',16,'fontname','arial')
line([0 0.6]*100,[0 0.6]*100,'color','k')

%% UP transition lag
close all
clear lag_* sm_*
% lag_range = linspace(-0.5,1.,1000);
lag_range = linspace(-0.5,1.5,1000);
dlag = lag_range(2)-lag_range(1);
for d = 1:length(sess_data)
    lag_hist(d,:) = histc(lf8_uplags{d}/Fsd,lag_range);
    lag_hist(d,:) = lag_hist(d,:)/sum(lag_hist(d,:))/dlag;
    lag_hist_nt2(d,:) = histc(lf8_uplags_nt2{d}/Fsd,lag_range);
    lag_hist_nt2(d,:) = lag_hist_nt2(d,:)/sum(lag_hist_nt2(d,:))/dlag;
    sm_lag_hist(d,:) = jmm_smooth_1d_cor(lag_hist(d,:),20);
    sm_lag_hist_nt2(d,:) = jmm_smooth_1d_cor(lag_hist_nt2(d,:),20);
    %     phase_hist(d,:) = fgsmooth(phase_hist(d,:),4);
end

mean_lag_mec = mean(sm_lag_hist(mec_cells,:));
u_lag_mec = mean_lag_mec + 1*std(sm_lag_hist(mec_cells,:))/sqrt(length(mec_cells));
l_lag_mec = mean_lag_mec - 1*std(sm_lag_hist(mec_cells,:))/sqrt(length(mec_cells));
mean_lag_lec = mean(sm_lag_hist(lec_cells,:));
u_lag_lec = mean_lag_lec + 1*std(sm_lag_hist(lec_cells,:))/sqrt(length(lec_cells));
l_lag_lec = mean_lag_lec - 1*std(sm_lag_hist(lec_cells,:))/sqrt(length(lec_cells));
mean_lag_np_mec = nanmean(sm_lag_hist_nt2(mec_cells,:));
u_lag_np_mec = mean_lag_np_mec + 1*nanstd(sm_lag_hist_nt2(mec_cells,:))/sqrt(length(mec_cells));
l_lag_np_mec = mean_lag_np_mec - 1*nanstd(sm_lag_hist_nt2(mec_cells,:))/sqrt(length(mec_cells));

figure
set(gca,'fontsize',14,'fontname','arial')
plot(lag_range,mean_lag_mec,'linewidth',2)
hold on
plot(lag_range,mean_lag_lec,'g','linewidth',2)
% plot(lag_range,mean_lag_np_mec,'c','linewidth',2)
% legend('MEC','LEC')

X = [lag_range fliplr(lag_range)];
Y = [u_lag_lec fliplr(l_lag_lec)];
fill(X,Y,'g')
X = [lag_range fliplr(lag_range)];
Y = [u_lag_mec fliplr(l_lag_mec)];
fill(X,Y,'b')
% X = [lag_range fliplr(lag_range)];
% Y = [u_lag_np_mec fliplr(l_lag_np_mec)];
% fill(X,Y,'c')
% xlim([-0.5 1.])
xlim([-0.4 1.2])
ylim([0 3.])
line([0 0],[0 3],'Color','k')
xlabel('Lag (s)','fontsize',16,'fontname','arial')
ylabel('Probability','fontsize',16,'fontname','arial')

%% DOWN transition lag
close all
clear lag_* sm_*
lag_range = linspace(-1,2.3,1000);
% lag_range = linspace(-4,4,4000);
dlag = lag_range(2)-lag_range(1);
for d = 1:length(sess_data)
    lag_hist(d,:) = histc(lf8_downlags{d}/Fsd,lag_range);
    lag_hist(d,:) = lag_hist(d,:)/sum(lag_hist(d,:))/dlag;
    lag_hist_nt2(d,:) = histc(lf8_downlags_nt2{d}/Fsd,lag_range);
    lag_hist_nt2(d,:) = lag_hist_nt2(d,:)/sum(lag_hist_nt2(d,:))/dlag;
    sm_lag_hist(d,:) = jmm_smooth_1d_cor(lag_hist(d,:),20);
    sm_lag_hist_nt2(d,:) = jmm_smooth_1d_cor(lag_hist_nt2(d,:),20);
    %     phase_hist(d,:) = fgsmooth(phase_hist(d,:),4);
end

mean_lag_mec = mean(sm_lag_hist(mec_cells,:));
u_lag_mec = mean_lag_mec + 1*std(sm_lag_hist(mec_cells,:))/sqrt(length(mec_cells));
l_lag_mec = mean_lag_mec - 1*std(sm_lag_hist(mec_cells,:))/sqrt(length(mec_cells));
mean_lag_lec = mean(sm_lag_hist(lec_cells,:));
u_lag_lec = mean_lag_lec + 1*std(sm_lag_hist(lec_cells,:))/sqrt(length(lec_cells));
l_lag_lec = mean_lag_lec - 1*std(sm_lag_hist(lec_cells,:))/sqrt(length(lec_cells));
mean_lag_np_mec = mean(sm_lag_hist_nt2(mec_cells,:));
u_lag_np_mec = mean_lag_np_mec + 1*std(sm_lag_hist_nt2(mec_cells,:))/sqrt(length(mec_cells));
l_lag_np_mec = mean_lag_np_mec - 1*std(sm_lag_hist_nt2(mec_cells,:))/sqrt(length(mec_cells));

figure
set(gca,'fontsize',14,'fontname','arial')
plot(lag_range,mean_lag_mec,'linewidth',2)
hold on
plot(lag_range,mean_lag_lec,'g','linewidth',2)
% plot(lag_range,mean_lag_np_mec,'y','linewidth',2)
% legend('MEC','LEC')

X = [lag_range fliplr(lag_range)];
Y = [u_lag_lec fliplr(l_lag_lec)];
fill(X,Y,'g')
X = [lag_range fliplr(lag_range)];
Y = [u_lag_mec fliplr(l_lag_mec)];
fill(X,Y,'b')
% X = [lag_range fliplr(lag_range)];
% Y = [u_lag_np_mec fliplr(l_lag_np_mec)];
% fill(X,Y,'c')
xlim([-0.7 2.])
% xlim([-4 2])
ylim([0 2.])
line([0 0],[0 2],'Color','k')
xlabel('Lag (s)','fontsize',16,'fontname','arial')
ylabel('Probability','fontsize',16,'fontname','arial')

%% RELATIVE UP STATE DURATION
close all
clear dur_* sm_*

dur_range = linspace(-3,6,1000);
ddur = dur_range(2)-dur_range(1);
for d = 1:length(sess_data)
    dur_hist(d,:) = histc(mp_rel_updurs{d},dur_range);
    dur_hist_nt2(d,:) = histc(mp_rel_updurs{d}(nrt2_ups{d}),dur_range);
    dur_hist_nt2(d,:) = dur_hist_nt2(d,:)/sum(dur_hist(d,:))/ddur;
    dur_hist(d,:) = dur_hist(d,:)/sum(dur_hist(d,:))/ddur;
    %     dur_hist_t2(d,:) = histc(mp_rel_updurs_t2{d},dur_range);
    %     dur_hist_t2(d,:) = dur_hist_t2(d,:)/sum(dur_hist_t2(d,:))/ddur;
    %     dur_hist_s2(d,:) = histc(mp_rel_updurs_s2{d},dur_range);
    %     dur_hist_s2(d,:) = dur_hist_s2(d,:)/sum(dur_hist_s2(d,:))/ddur;
    %     dur_hist_s3(d,:) = histc(mp_rel_updurs_s3{d},dur_range);
    %     dur_hist_s3(d,:) = dur_hist_s3(d,:)/sum(dur_hist_s3(d,:))/ddur;
    %     dur_hist_s4(d,:) = histc(mp_rel_updurs_s4{d},dur_range);
    %     dur_hist_s4(d,:) = dur_hist_s4(d,:)/sum(dur_hist_s4(d,:))/ddur;
    sm_dur_hist(d,:) = jmm_smooth_1d_cor(dur_hist(d,:),20);
    sm_dur_hist_nt2(d,:) = jmm_smooth_1d_cor(dur_hist_nt2(d,:),20);
    %     sm_dur_hist_t2(d,:) = jmm_smooth_1d_cor(dur_hist_t2(d,:),20);
    %     sm_dur_hist_s2(d,:) = jmm_smooth_1d_cor(dur_hist_s2(d,:),25);
    %     sm_dur_hist_s3(d,:) = jmm_smooth_1d_cor(dur_hist_s3(d,:),25);
    %     sm_dur_hist_s4(d,:) = jmm_smooth_1d_cor(dur_hist_s4(d,:),25);
end
eps = 1e-5;
mean_dur_mec = mean(sm_dur_hist(mec_cells,:));
u_dur_mec = mean_dur_mec + 1*std(sm_dur_hist(mec_cells,:))/sqrt(length(mec_cells));
l_dur_mec = mean_dur_mec - 1*std(sm_dur_hist(mec_cells,:))/sqrt(length(mec_cells));
mean_dur_mec(mean_dur_mec <= eps) = eps;
u_dur_mec(u_dur_mec <= eps) = eps;
l_dur_mec(l_dur_mec <= eps) = eps;
mean_dur_lec = mean(sm_dur_hist(lec_cells,:));
u_dur_lec = mean_dur_lec + 1*std(sm_dur_hist(lec_cells,:))/sqrt(length(lec_cells));
l_dur_lec = mean_dur_lec - 1*std(sm_dur_hist(lec_cells,:))/sqrt(length(lec_cells));
mean_dur_lec(mean_dur_lec <= eps) = eps;
u_dur_lec(u_dur_lec <= eps) = eps;
l_dur_lec(l_dur_lec <= eps) = eps;
mean_dur_np_mec = mean(sm_dur_hist_nt2(mec_cells,:));
u_dur_np_mec = mean_dur_np_mec + 1*std(sm_dur_hist_nt2(mec_cells,:))/sqrt(length(mec_cells));
l_dur_np_mec = mean_dur_np_mec - 1*std(sm_dur_hist_nt2(mec_cells,:))/sqrt(length(mec_cells));
mean_dur_np_mec(mean_dur_np_mec <= eps) = eps;
u_dur_np_mec(u_dur_np_mec <= eps) = eps;
l_dur_np_mec(l_dur_np_mec <= eps) = eps;
mean_dur_np_lec = mean(sm_dur_hist_nt2(lec_cells,:));
u_dur_np_lec = mean_dur_np_lec + 1*std(sm_dur_hist_nt2(lec_cells,:))/sqrt(length(lec_cells));
l_dur_np_lec = mean_dur_np_lec - 1*std(sm_dur_hist_nt2(lec_cells,:))/sqrt(length(lec_cells));
mean_dur_np_lec(mean_dur_np_lec <= eps) = eps;
u_dur_np_lec(u_dur_np_lec <= eps) = eps;
l_dur_np_lec(l_dur_np_lec <= eps) = eps;
% mean_dur_p_mec = mean(sm_dur_hist_t2(mec_cells,:));
% u_dur_p_mec = mean_dur_p_mec + 1*std(sm_dur_hist_t2(mec_cells,:))/sqrt(length(mec_cells));
% l_dur_p_mec = mean_dur_p_mec - 1*std(sm_dur_hist_t2(mec_cells,:))/sqrt(length(mec_cells));
% mean_dur_s2_mec = nanmean(sm_dur_hist_s2(mec_cells,:));
% u_dur_s2_mec = mean_dur_s2_mec + 1*nanstd(sm_dur_hist_s2(mec_cells,:))/sqrt(length(mec_cells));
% l_dur_s2_mec = mean_dur_s2_mec - 1*nanstd(sm_dur_hist_s2(mec_cells,:))/sqrt(length(mec_cells));
% mean_dur_s3_mec = nanmean(sm_dur_hist_s3(mec_cells,:));
% u_dur_s3_mec = mean_dur_s3_mec + 1*nanstd(sm_dur_hist_s3(mec_cells,:))/sqrt(length(mec_cells));
% l_dur_s3_mec = mean_dur_s3_mec - 1*nanstd(sm_dur_hist_s3(mec_cells,:))/sqrt(length(mec_cells));
% mean_dur_s4_mec = nanmean(sm_dur_hist_s4(mec_cells,:));
% u_dur_s4_mec = mean_dur_s4_mec + 1*nanstd(sm_dur_hist_s4(mec_cells,:))/sqrt(length(mec_cells));
% l_dur_s4_mec = mean_dur_s4_mec - 1*nanstd(sm_dur_hist_s4(mec_cells,:))/sqrt(length(mec_cells));

figure
set(gca,'fontsize',14,'fontname','arial')
plot(dur_range,mean_dur_mec,'linewidth',2)
hold on
plot(dur_range,mean_dur_lec,'g','linewidth',2)
plot(dur_range,mean_dur_np_mec,'c','linewidth',2)
plot(dur_range,mean_dur_np_lec,'r','linewidth',2)
% plot(dur_range,mean_dur_p_mec,'k','linewidth',2)
% plot(dur_range,mean_dur_s2_mec,'k','linewidth',2)
% plot(dur_range,mean_dur_s3_mec,'r','linewidth',2)
% plot(dur_range,mean_dur_s3_mec,'y','linewidth',2)
% legend('MEC','LEC','MEC non-persistent','MEC persistent')
% h = errorbar(dur_range,nanmean(sm_dur_hist(mec_cells,:)),nanstd(sm_dur_hist(mec_cells,:))/sqrt(length(mec_cells)));
% errorbar_tick(h,.001,'units')
% hold on
% h = errorbar(dur_range,nanmean(sm_dur_hist(lec_cells,:)),nanstd(sm_dur_hist(lec_cells,:))/sqrt(length(lec_cells)),'g');
% errorbar_tick(h,.001,'units')
% h = errorbar(dur_range,nanmean(sm_dur_hist_nt2(mec_cells,:)),nanstd(sm_dur_hist_nt2(mec_cells,:))/sqrt(length(mec_cells)),'c');
% errorbar_tick(h,.001,'units')
% h = errorbar(dur_range,nanmean(sm_dur_hist_nt2(lec_cells,:)),nanstd(sm_dur_hist_nt2(lec_cells,:))/sqrt(length(mec_cells)),'r');
% errorbar_tick(h,.001,'units')
% xlim([-3 6])
% xlabel('Dur (s)')
% ylabel('Probability')
% ylim([0 1.3])
% line([0 0],[1e-3 1.3],'color','k')
% set(gca,'yscale','log')
% line([1e-3 0],[0 1.3],'color','k')
% ylim([2e-3 1.3])


X = [dur_range fliplr(dur_range)];
Y = [u_dur_lec fliplr(l_dur_lec)];
fill(X,Y,'g')
X = [dur_range fliplr(dur_range)];
Y = [u_dur_mec fliplr(l_dur_mec)];
fill(X,Y,'b')
X = [dur_range fliplr(dur_range)];
Y = [u_dur_np_mec fliplr(l_dur_np_mec)];
fill(X,Y,'c')
X = [dur_range fliplr(dur_range)];
Y = [u_dur_np_lec fliplr(l_dur_np_lec)];
fill(X,Y,'r')
% X = [dur_range fliplr(dur_range)];
% Y = [u_dur_s2_mec fliplr(l_dur_s2_mec)];
% fill(X,Y,'k')
% X = [dur_range fliplr(dur_range)];
% Y = [u_dur_s3_mec fliplr(l_dur_s3_mec)];
% fill(X,Y,'r')
% X = [dur_range fliplr(dur_range)];
% Y = [u_dur_s4_mec fliplr(l_dur_s4_mec)];
% fill(X,Y,'y')
xlim([-3 5])
xlabel('Relative Duration (s)','fontsize',16,'fontname','arial')
ylabel('Probability Density','fontsize',16,'fontname','arial')
% ylim([0 1.3])
% line([0 0],[0 1.3],'color','k')
set(gca,'yscale','log')
line([0 0],[1e-3 1.3],'color','k')
ylim([3e-3 1.3])

%% RELATIVE DOWN STATE DURATION
% close all
dur_range = linspace(-4,5,1000);
ddur = dur_range(2)-dur_range(1);
for d = 1:length(sess_data)
    dur_hist(d,:) = histc(mp_rel_downdurs{d},dur_range);
    dur_hist_nt2(d,:) = histc(mp_rel_downdurs{d}(nrt2_downs{d}),dur_range);
    dur_hist_nt2(d,:) = dur_hist_nt2(d,:)/sum(dur_hist(d,:))/ddur;
    dur_hist(d,:) = dur_hist(d,:)/sum(dur_hist(d,:))/ddur;
    
%     dur_hist_t2(d,:) = histc(mp_rel_downdurs_t2{d},dur_range);
%     dur_hist_t2(d,:) = dur_hist_t2(d,:)/sum(dur_hist_t2(d,:))/ddur;
    sm_dur_hist(d,:) = jmm_smooth_1d_cor(dur_hist(d,:),20);
    sm_dur_hist_nt2(d,:) = jmm_smooth_1d_cor(dur_hist_nt2(d,:),20);
%     sm_dur_hist_t2(d,:) = jmm_smooth_1d_cor(dur_hist_t2(d,:),20);
end

eps = 1e-5;
mean_dur_mec = mean(sm_dur_hist(mec_cells,:));
u_dur_mec = mean_dur_mec + 1*std(sm_dur_hist(mec_cells,:))/sqrt(length(mec_cells));
l_dur_mec = mean_dur_mec - 1*std(sm_dur_hist(mec_cells,:))/sqrt(length(mec_cells));
mean_dur_mec(mean_dur_mec <= eps) = eps;
u_dur_mec(u_dur_mec <= eps) = eps;
l_dur_mec(l_dur_mec <= eps) = eps;
mean_dur_lec = mean(sm_dur_hist(lec_cells,:));
u_dur_lec = mean_dur_lec + 1*std(sm_dur_hist(lec_cells,:))/sqrt(length(lec_cells));
l_dur_lec = mean_dur_lec - 1*std(sm_dur_hist(lec_cells,:))/sqrt(length(lec_cells));
mean_dur_lec(mean_dur_lec <= eps) = eps;
u_dur_lec(u_dur_lec <= eps) = eps;
l_dur_lec(l_dur_lec <= eps) = eps;
mean_dur_np_mec = nanmean(sm_dur_hist_nt2(mec_cells,:));
u_dur_np_mec = mean_dur_np_mec + 1*nanstd(sm_dur_hist_nt2(mec_cells,:))/sqrt(length(mec_cells));
l_dur_np_mec = mean_dur_np_mec - 1*nanstd(sm_dur_hist_nt2(mec_cells,:))/sqrt(length(mec_cells));
mean_dur_np_mec(mean_dur_np_mec <= eps) = eps;
u_dur_np_mec(u_dur_np_mec <= eps) = eps;
l_dur_np_mec(l_dur_np_mec <= eps) = eps;
mean_dur_np_mec(mean_dur_np_mec <= 0) = eps;
% mean_dur_p_mec = nanmean(sm_dur_hist_t2(mec_cells,:));
% u_dur_p_mec = mean_dur_p_mec + 1*nanstd(sm_dur_hist_t2(mec_cells,:))/sqrt(length(mec_cells));
% l_dur_p_mec = mean_dur_p_mec - 1*nanstd(sm_dur_hist_t2(mec_cells,:))/sqrt(length(mec_cells));
mean_dur_np_lec = nanmean(sm_dur_hist_nt2(lec_cells,:));
u_dur_np_lec = mean_dur_np_lec + 1*nanstd(sm_dur_hist_nt2(lec_cells,:))/sqrt(length(lec_cells));
l_dur_np_lec = mean_dur_np_lec - 1*nanstd(sm_dur_hist_nt2(lec_cells,:))/sqrt(length(lec_cells));
mean_dur_np_lec(mean_dur_np_lec <= eps) = eps;
u_dur_np_lec(u_dur_np_lec <= eps) = eps;
l_dur_np_lec(l_dur_np_lec <= eps) = eps;
% mean_dur_p_lec = nanmean(sm_dur_hist_t2(lec_cells,:));
% u_dur_p_lec = mean_dur_p_lec + 1*nanstd(sm_dur_hist_t2(lec_cells,:))/sqrt(length(lec_cells));
% l_dur_p_lec = mean_dur_p_lec - 1*nanstd(sm_dur_hist_t2(lec_cells,:))/sqrt(length(lec_cells));

% figure
% h = errorbar(dur_range,nanmean(sm_dur_hist(mec_cells,:)),nanstd(sm_dur_hist(mec_cells,:))/sqrt(length(mec_cells)));
% errorbar_tick(h,.001,'units')
% hold on
% h = errorbar(dur_range,nanmean(sm_dur_hist(lec_cells,:)),nanstd(sm_dur_hist(lec_cells,:))/sqrt(length(lec_cells)),'g');
% errorbar_tick(h,.001,'units')
% h = errorbar(dur_range,nanmean(sm_dur_hist_nt2(mec_cells,:)),nanstd(sm_dur_hist_nt2(mec_cells,:))/sqrt(length(mec_cells)),'c');
% errorbar_tick(h,.001,'units')
% h = errorbar(dur_range,nanmean(sm_dur_hist_nt2(lec_cells,:)),nanstd(sm_dur_hist_nt2(lec_cells,:))/sqrt(length(mec_cells)),'r');
% errorbar_tick(h,.001,'units')
% xlim([-4 5])
% xlabel('Dur (s)')
% ylabel('Probability')
% % ylim([0 1.3])
% line([0 0],[1e-3 1.3],'color','k')
% set(gca,'yscale','log')
% ylim([1e-3 1.3])

figure
set(gca,'fontsize',14,'fontname','arial')
hold on
plot(dur_range,mean_dur_mec,'linewidth',2)
hold on
plot(dur_range,mean_dur_lec,'g','linewidth',2)
plot(dur_range,mean_dur_np_mec,'c','linewidth',2)
% plot(dur_range,mean_dur_p_mec,'k','linewidth',2)
plot(dur_range,mean_dur_np_lec,'r','linewidth',2)
% plot(dur_range,mean_dur_p_lec,'y','linewidth',2)
% legend('MEC','LEC','MEC non-persistent','MEC persistent')
%
X = [dur_range fliplr(dur_range)];
Y = [u_dur_lec fliplr(l_dur_lec)];
fill(X,Y,'g')
X = [dur_range fliplr(dur_range)];
Y = [u_dur_mec fliplr(l_dur_mec)];
fill(X,Y,'b')
X = [dur_range fliplr(dur_range)];
Y = [u_dur_np_mec fliplr(l_dur_np_mec)];
fill(X,Y,'c')
% X = [dur_range fliplr(dur_range)];
% Y = [u_dur_p_mec fliplr(l_dur_p_mec)];
% fill(X,Y,'k')
X = [dur_range fliplr(dur_range)];
Y = [u_dur_np_lec fliplr(l_dur_np_lec)];
fill(X,Y,'r')
% X = [dur_range fliplr(dur_range)];
% Y = [u_dur_p_lec fliplr(l_dur_p_lec)];
% fill(X,Y,'y')


% xlim([-2 5])
% xlabel('Dur (s)')
% ylabel('Probability')
% ylim([0 1.4])
% line([0 0],[0 1.4],'color','k')
line([0 0],[1e-3 1.7],'color','k')
set(gca,'yscale','log')
ylim([3e-3 1.7])
xlabel('Relative Duration (s)','fontsize',16,'fontname','arial')
ylabel('Probability Density','fontsize',16,'fontname','arial')
xlim([-3. 5])

%% UP quantization
% hist_range = linspace(0,6,1000);
hist_range = linspace(0,8,1500);
binsize = hist_range(2)-hist_range(1);

for i = 1:36
    cell_hist(i,:) = hist(rmp_updurs_lfpc{i},hist_range);
    cell_hist(i,:) = cell_hist(i,:)/sum(cell_hist(i,:))/binsize;
    cell_hist(i,:) = jmm_smooth_1d_cor(cell_hist(i,:),20);
end

% clear cell_hist
% chi_square = [];
% for i = mec_cells
%     temp = ceil(rmp_updurs_lfpc{i});
%     temp(temp > 10) = [];
%     cell_freq(i,:) = hist(temp,1:10);
%     cell_rfreq(i,:) = cell_freq(i,:)/sum(cell_freq(i,:));
%     p(i) = 1/nanmean(temp);
%     n_rups(i) = length(temp);
%     pred_freq(i,:) = n_rups(i)*p(i).*(1-p(i)).^(0:9);
%     good_pts = find(pred_freq(i,:) >= 5);
%     chi_square = [chi_square (pred_freq(i,good_pts) - cell_freq(i,good_pts)).^2./pred_freq(i,good_pts)];
% end
% sum_chi_square = sum(chi_square);
% dof = length(chi_square) - length(mec_cells) - 1;
% chi2f_p = 1 - chi2cdf(sum_chi_square,dof);
% errorbar(1:6,nanmean(cell_hist(mec_cells,:)),nanstd(cell_hist(mec_cells,:))/sqrt(22));
% hold on
% errorbar(1:6,nanmean(cell_hist(lec_cells,:)),nanstd(cell_hist(lec_cells,:))/sqrt(14),'r')
% set(gca,'yscale','log')

mean_hist_mec = nanmean(cell_hist(mec_cells,:));
u_mec = mean_hist_mec+1*nanstd(cell_hist(mec_cells,:))./sqrt(sum(~isnan(cell_hist(mec_cells,:))));
l_mec = mean_hist_mec-1*nanstd(cell_hist(mec_cells,:))./sqrt(sum(~isnan(cell_hist(mec_cells,:))));

mean_hist_lec = nanmean(cell_hist(lec_cells,:));
u_lec = mean_hist_lec+1*nanstd(cell_hist(lec_cells,:))./sqrt(sum(~isnan(cell_hist(lec_cells,:))));
l_lec = mean_hist_lec-1*nanstd(cell_hist(lec_cells,:))./sqrt(sum(~isnan(cell_hist(lec_cells,:))));

figure
set(gca,'fontsize',14,'fontname','arial')
plot(hist_range,mean_hist_mec,'linewidth',2)
hold on
plot(hist_range,mean_hist_lec,'g','linewidth',2)
legend('MEC','LEC')

epsilon = 1e-10;
cell_hist(cell_hist < epsilon) = epsilon;
u_mec(u_mec < epsilon) = epsilon;
l_mec(l_mec < epsilon) = epsilon;
u_lec(u_lec < epsilon) = epsilon;
l_lec(l_lec < epsilon) = epsilon;

X = [hist_range fliplr(hist_range)];
Y = [u_mec fliplr(l_mec)];
fill(X,Y,'b')
Y = [u_lec fliplr(l_lec)];
fill(X,Y,'g')
%
set(gca,'yscale','log')
% xlim([0 5.3])
xlim([0 7.3])
ylim([3e-4 4])
xlabel('Duration (Ncx UDS Cycles)','fontsize',16,'fontname','arial')
ylabel('Probability Density','fontsize',16,'fontname','arial')
[peak_amps,b] = findpeaks(mean_hist_mec);
peak_locs = hist_range(b);
% plot(peak_locs,peak_amps,'ro')
% for i = 1:length(peak_locs)
%     line([peak_locs(i) peak_locs(i)],[5e-4 2],'Color','k')
% end
% temp_p = polyfit(peak_locs,log10(peak_amps),1);
% pred_p = 10.^polyval(temp_p,hist_range);
% plot(hist_range,pred_p,'k','linewidth',2)
% xlabel('MP Up state duration (LFP cycles)')
% ylabel('Probability')

%% DOWN quantization
hist_range = linspace(0,6,1000);
binsize = hist_range(2)-hist_range(1);

for i = 1:36
    cell_hist(i,:) = hist(rmp_downdurs_lfpc{i},hist_range);
    cell_hist(i,:) = cell_hist(i,:)/sum(cell_hist(i,:))/binsize;
    cell_hist(i,:) = jmm_smooth_1d_cor(cell_hist(i,:),20);
end

% % clear cell_hist
% chi_square = [];
% for i = lec_cells
%     temp = ceil(rmp_downdurs_lfpc{i});
%     temp(temp > 10) = [];
%     cell_freq(i,:) = hist(temp,1:10);
%     cell_rfreq(i,:) = cell_freq(i,:)/sum(cell_freq(i,:));
%     p(i) = 1/nanmean(temp);
%     n_rups(i) = length(temp);
%     pred_freq(i,:) = n_rups(i)*p(i).*(1-p(i)).^(0:9);
%     good_pts = find(pred_freq(i,:) >= 5);
%     chi_square = [chi_square (pred_freq(i,good_pts) - cell_freq(i,good_pts)).^2./pred_freq(i,good_pts)];
% end
% sum_chi_square = sum(chi_square);
% dof = length(lec_cells) + 1;
% chi2f_p = 1 - chi2cdf(sum_chi_square,dof);
% errorbar(1:6,nanmean(cell_hist(mec_cells,:)),nanstd(cell_hist(mec_cells,:))/sqrt(22));
% hold on
% errorbar(1:6,nanmean(cell_hist(lec_cells,:)),nanstd(cell_hist(lec_cells,:))/sqrt(14),'r')
% set(gca,'yscale','log')

mean_hist_mec = nanmean(cell_hist(mec_cells,:));
u_mec = mean_hist_mec+1*nanstd(cell_hist(mec_cells,:))./sqrt(sum(~isnan(cell_hist(mec_cells,:))));
l_mec = mean_hist_mec-1*nanstd(cell_hist(mec_cells,:))./sqrt(sum(~isnan(cell_hist(mec_cells,:))));

mean_hist_lec = nanmean(cell_hist(lec_cells,:));
u_lec = mean_hist_lec+1*nanstd(cell_hist(lec_cells,:))./sqrt(sum(~isnan(cell_hist(lec_cells,:))));
l_lec = mean_hist_lec-1*nanstd(cell_hist(lec_cells,:))./sqrt(sum(~isnan(cell_hist(lec_cells,:))));

figure
set(gca,'fontsize',14,'fontname','arial')
plot(hist_range,mean_hist_mec,'linewidth',2)
hold on
plot(hist_range,mean_hist_lec,'g','linewidth',2)
legend('MEC','LEC')

epsilon = 1e-10;
cell_hist(cell_hist < epsilon) = epsilon;
u_mec(u_mec < epsilon) = epsilon;
l_mec(l_mec < epsilon) = epsilon;
u_lec(u_lec < epsilon) = epsilon;
l_lec(l_lec < epsilon) = epsilon;

X = [hist_range fliplr(hist_range)];
Y = [u_mec fliplr(l_mec)];
fill(X,Y,'b')
Y = [u_lec fliplr(l_lec)];
fill(X,Y,'g')

set(gca,'yscale','log')
xlim([0 5.3])
ylim([3e-4 3])
xlabel('Duration (Ncx UDS Cycles)','fontsize',16,'fontname','arial')
ylabel('Probability Density','fontsize',16,'fontname','arial')

[peak_amps,b] = findpeaks(mean_hist_mec);
peak_locs = hist_range(b);
% plot(peak_locs,peak_amps,'ro')
% for i = 1:length(peak_locs)
%     line([peak_locs(i) peak_locs(i)],[5e-4 2],'Color','k')
% end
% temp_p = polyfit(peak_locs,log10(peak_amps),1);
% pred_p = 10.^polyval(temp_p,hist_range);
% plot(hist_range,pred_p,'k','linewidth',2)
% xlabel('MP Up state duration (LFP cycles)')
% ylabel('Probability')
