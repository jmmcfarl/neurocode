clear all
close all

load G:\WC_Germany\overall_EC\overall_EC_dir.mat
addpath('G:\WC_Germany\parietal_cortical_2010\')
addpath('G:\WC_Germany\hsmm_state_detection\')
addpath('G:\WC_Germany\persistent_revised\')

used_data = [l3mec_p l3lec_p];
sess_data = sess_data(used_data);

cd G:\WC_Germany\persistent_9_27_2010\
load ./pa_corresponding_lfp_revised_simp_new2
load ./depthofanesth_data_new2
% load ./pa_heka_UDS_data_allcells_sm500
% load ./spike_rate_data_sm500

% dep_variable = 'Frequency';
dep_variable = 'Duty Cycle';

% xl = [0.2 0.6];
xl = [0.3 0.6];
nbins = 10;
% binedges = linspace(0.4,2,nbins); %up duration
% binedges = linspace(1,6,nbins); %cycle duration
% binedges = linspace(0.3,0.6,nbins); %duty cycle
binedges = linspace(0.2,0.6,nbins+1); %for freq
% binedges = linspace(0.4,5,nbins); %down dur
% binedges = linspace(-5,5,nbins); %lf8 pow

bincenters = (binedges(1:end-1)+binedges(2:end))/2;

lf8_up_fract = 0.5;

% temp = [];
for d = 1:length(sess_data)
  
%     if length(lf8_dutycyc{d}) ~= length(mp_updurs_lfpc{d}) ||... 
%         length(lf8_udsfreq{d}) ~= length(mp_updurs_lfpc{d}) ||...
%         length(mp_upstate_heka_meanamp{d}) ~= length(mp_updurs_lfpc{d})
%         error('LENGTH MISMATCH')
%     end
    
%     temp = [temp; lf8_udsfreq{d}];
    
%     [counts,binlocs] = histc(lf8_medupdur{d},binedges);
%     [counts,binlocs] = histc(lf8_medcycdur{d},binedges);
%     [counts,binlocs] = histc(lf8_dutycyc{d},binedges);    
    [counts,binlocs] = histc(lf8_udsfreq{d},binedges);
%     [counts,binlocs] = histc(lf8_meddowndur{d},binedges);
%     [counts,binlocs] = histc(lf8_udspow{d},binedges);
    
%     cond_dist_dur(d,:,:) = nan(length(lfp_updur_bincenters),length(mp_updur_bins));
%     cond_dist_cdur(d,:,:) = nan(length(lfp_updur_bincenters),length(mp_upcdur_bins));
    cond_prob_p2(d,:) = nan(1,length(bincenters));
    cond_prob_p1(d,:) = nan(1,length(bincenters));
    cond_mean_updur(d,:) = nan(1,length(bincenters));
    cond_mean_upcdur(d,:) = nan(1,length(bincenters));
    cond_mean_uplag(d,:) = nan(1,length(bincenters));
    cond_mean_downlag(d,:) = nan(1,length(bincenters));
%     cond_mean_uprate(d,:) = nan(1,length(bincenters)); 
%     cond_mean_upamp(d,:) = nan(1,length(bincenters));
%     cond_mean_downamp(d,:) = nan(1,length(bincenters));
%     cond_mean_upstd(d,:) = nan(1,length(bincenters));
%     cond_mean_downstd(d,:) = nan(1,length(bincenters));
    cond_mean_relupdur(d,:) = nan(1,length(bincenters));
    
    rt2_up_ind = nan(size(lf8_udsfreq{d}));
    rt2_up_ind(rt2_ups{d}) = 1;
    rt2_up_ind(nrt2_ups{d}) = 0;
    
    use_pts = find(~isnan(lf8_udsfreq{d}) & ~isnan(lf8_downlags{d}));
    [a,b] = corr(lf8_udsfreq{d}(use_pts),lf8_downlags{d}(use_pts),'type','spearman');
    udsfreq_downlag_corr(d) = a;
    udsfreq_downlag_p(d) = b;
    use_pts = find(~isnan(lf8_udsfreq{d}) & ~isnan(lf8_uplags{d}));
    [a,b] = corr(lf8_udsfreq{d}(use_pts),lf8_uplags{d}(use_pts),'type','spearman');
    udsfreq_uplag_corr(d) = a;
    udsfreq_uplag_p(d) = b;
    use_pts = find(~isnan(lf8_udsfreq{d}) & ~isnan(rt2_up_ind));
    [a,b] = corr(lf8_udsfreq{d}(use_pts),rt2_up_ind(use_pts),'type','spearman');
    udsfreq_rt2_corr(d) = a;
    udsfreq_rt2_p(d) = b;
    
    use_pts = find(~isnan(lf8_dutycyc{d}) & ~isnan(lf8_downlags{d}));
    [a,b] = corr(lf8_dutycyc{d}(use_pts),lf8_downlags{d}(use_pts),'type','spearman');
    udscyc_downlag_corr(d) = a;
    udscyc_downlag_p(d) = b;
    use_pts = find(~isnan(lf8_dutycyc{d}) & ~isnan(lf8_uplags{d}));
    [a,b] = corr(lf8_dutycyc{d}(use_pts),lf8_uplags{d}(use_pts),'type','spearman');
    udscyc_uplag_corr(d) = a;
    udscyc_uplag_p(d) = b;
    use_pts = find(~isnan(lf8_dutycyc{d}) & ~isnan(rt2_up_ind));
    [a,b] = corr(lf8_dutycyc{d}(use_pts),rt2_up_ind(use_pts),'type','spearman');
    udscyc_rt2_corr(d) = a;
    udscyc_rt2_p(d) = b;

    use_pts = find(~isnan(lf8_meddowndur{d}) & ~isnan(lf8_downlags{d}));
    [a,b] = corr(lf8_meddowndur{d}(use_pts),lf8_downlags{d}(use_pts),'type','spearman');
    downdur_downlag_corr(d) = a;
    downdur_downlag_p(d) = b;
    use_pts = find(~isnan(lf8_meddowndur{d}) & ~isnan(lf8_uplags{d}));
    [a,b] = corr(lf8_meddowndur{d}(use_pts),lf8_uplags{d}(use_pts),'type','spearman');
    downdur_uplag_corr(d) = a;
    downdur_uplag_p(d) = b;
    use_pts = find(~isnan(lf8_meddowndur{d}) & ~isnan(rt2_up_ind));
    [a,b] = corr(lf8_meddowndur{d}(use_pts),rt2_up_ind(use_pts),'type','spearman');
    downdur_rt2_corr(d) = a;
    downdur_rt2_p(d) = b;

    for i = 1:length(bincenters)
        cur_points = find(binlocs==i);
        num_counts(d,i) = length(cur_points);
        if length(cur_points) > 2
%             n_rpers_ups = sum(ismember(t2_ups{d},cur_points));
            n_rpers_ups = sum(ismember(rt2_ups{d},cur_points));
            n_npers_ups = sum(ismember(nrt2_ups{d},cur_points));
%             cond_dist_dur(d,i,:) = ksdensity(mp_updur{d}(cur_points),mp_updur_bins);
%             cond_dist_cdur(d,i,:) = ksdensity(mp_updur_corresp_lfp{d}(cur_points),mp_upcdur_bins);
%             cond_mean_uprate(d,i) = nanmean(mp_upstate_rate{d}(cur_points));
            cond_prob_p2(d,i) = n_rpers_ups/(n_rpers_ups+n_npers_ups);
%             cond_prob_p2(d,i) = n_rpers_ups/length(cur_points);

%             cond_prob_p2(d,i) = nansum(mp_updur_lfpc_delay{d}(cur_points) > 1)/length(cur_points);
            cond_prob_p1(d,i) = nansum(mp_rel_updurs{d}(cur_points) > 0)/length(cur_points);
            cond_mean_updur(d,i) = nanmean(mp_updurs{d}(cur_points));
            cond_mean_upcdur(d,i) = nanmean(mp_updurs_lfpc{d}(cur_points));
            cond_mean_uplag(d,i) = nanmean(lf8_uplags{d}(cur_points));
            cond_mean_downlag(d,i) = nanmean(lf8_downlags{d}(cur_points));
%             cond_mean_upamp(d,i) = nanmean(mp_upstate_heka_meanamp_sig{d}(cur_points));
%             cond_mean_downamp(d,i) = nanmean(mp_downstate_heka_meanamp_sig{d}(cur_points));
%             cond_mean_upstd(d,i) = nanmean(mp_upstate_nlx_stdamp_sig{d}(cur_points));
%             cond_mean_downstd(d,i) = nanmean(mp_downstate_nlx_stdamp_sig{d}(cur_points));
            cur_nt2s = find(ismember(nrt2_ups{d},cur_points));
            if ~isempty(cur_nt2s)
            cond_mean_relupdur(d,i) = nanmean(mp_rel_updurs_nt2{d}(cur_nt2s));
            end
        end
    end
    fract_rt1_ups(d) = nansum(mp_rel_updurs{d} > 0)/sum(~isnan(mp_rel_updurs{d}));
    cond_nprob_p2(d,:) = cond_prob_p2(d,:)/fract_rt2_ups(d);
    cond_nprob_p1(d,:) = cond_prob_p1(d,:)/fract_rt1_ups(d);
    cond_nmean_updur(d,:) = cond_mean_updur(d,:)/nanmean(mp_updurs{d});
    cond_nmean_upcdur(d,:) = cond_mean_upcdur(d,:)/nanmean(mp_updurs_lfpc{d});
    cond_nmean_uplag(d,:) = cond_mean_uplag(d,:) - nanmean(lf8_uplags{d});
    cond_nmean_downlag(d,:) = cond_mean_downlag(d,:) - nanmean(lf8_downlags{d});
%     cond_nmean_upamp(d,:) = cond_mean_upamp(d,:) - upstate_mean_spksub(d);
%     cond_nmean_downamp(d,:) = cond_mean_downamp(d,:) - downstate_mean_spksub(d);
%     cond_nmean_upstd(d,:) = cond_mean_upstd(d,:)/sqrt(upstate_var_hfnlx(d));
%     cond_nmean_downstd(d,:) = cond_mean_downstd(d,:)/sqrt(downstate_var_hfnlx(d));
%     cond_nmean_uprate(d,:) = cond_mean_uprate(d,:)/mp_up_rate(d);
end

mec_cells = 1:length(l3mec_p);
lec_cells = (mec_cells(end)+1):(mec_cells(end)+length(l3lec_p));

%% Unnormalized plots    
% plot_mec_lec_comparison(cond_prob_p2,bincenters,mec_cells,lec_cells)
% set(gca,'fontsize',14,'fontname','arial')
% xlabel(dep_variable,'fontsize',16,'fontname','arial')
% ylabel('Prob Type 2 Persistence','fontsize',16,'fontname','arial')
% legend('MEC','LEC')
% xlim(xl)
% ylim([0 0.35])
% 
% plot_mec_lec_comparison(cond_prob_p1,bincenters,mec_cells,lec_cells)
% set(gca,'fontsize',14,'fontname','arial')
% xlabel(dep_variable,'fontsize',16,'fontname','arial')
% ylabel('Prob Type 1 Persistence','fontsize',16,'fontname','arial')
% legend('MEC','LEC')
% xlim(xl)
% ylim([0.2 1])
% 
% plot_mec_lec_comparison(cond_mean_updur,bincenters,mec_cells,lec_cells)
% set(gca,'fontsize',14,'fontname','arial')
% xlabel(dep_variable,'fontsize',16,'fontname','arial')
% ylabel('Mean MP up state duration (s)','fontsize',16,'fontname','arial')
% legend('MEC','LEC')
% xlim(xl)
% ylim([0.5 3])
% 
% % plot_mec_lec_comparison(cond_mean_relupdur,bincenters,mec_cells,lec_cells)
% % set(gca,'fontsize',14,'fontname','arial')
% % xlabel(dep_variable,'fontsize',16,'fontname','arial')
% % ylabel('Mean relative up state duration (s)','fontsize',16,'fontname','arial')
% % legend('MEC','LEC')
% % xlim([0.25 0.6])
% 
% % plot_mec_lec_comparison(cond_mean_upcdur,bincenters,mec_cells,lec_cells)
% % xlabel(dep_variable)
% % ylabel('Mean MP up state duration (LFP cycles)')
% % legend('MEC','LEC')
% % xlim([0.2 0.6])
% 
% % plot_mec_lec_comparison(cond_mean_uplag,bincenters,mec_cells,lec_cells)
% % xlabel(dep_variable)
% % ylabel('Mean up lag (s)')
% % legend('MEC','LEC')
% % 
% plot_mec_lec_comparison(cond_mean_downlag/252,bincenters,mec_cells,lec_cells)
% xlabel(dep_variable)
% ylabel('Mean MP down lag (s)')
% legend('MEC','LEC')
% 
% % plot_mec_lec_comparison(cond_mean_upamp,bincenters,mec_cells,lec_cells)
% % xlabel(dep_variable)
% % ylabel('Mean MP up-state amplitude')
% % legend('MEC','LEC')
% % 
% % plot_mec_lec_comparison(cond_mean_downamp,bincenters,mec_cells,lec_cells)
% % xlabel(dep_variable)
% % ylabel('Mean MP down-state amplitude')
% % legend('MEC','LEC')
% % 
% % plot_mec_lec_comparison(cond_mean_upstd,bincenters,mec_cells,lec_cells)
% % xlabel(dep_variable)
% % ylabel('Mean MP up-state std')
% % legend('MEC','LEC')
% % 
% % plot_mec_lec_comparison(cond_mean_downstd,bincenters,mec_cells,lec_cells)
% % xlabel(dep_variable)
% % ylabel('Mean MP down-state std')
% % legend('MEC','LEC')
% % 
% % plot_mec_lec_comparison(cond_mean_uprate,bincenters,mec_cells,lec_cells)
% % xlabel(dep_variable)
% % ylabel('Mean MP up rate (Hz)')
% % legend('MEC','LEC')
% 
% %% for normalized plots
% % plot_mec_lec_comparison(cond_nprob_p2,bincenters,mec_cells,lec_cells)
% % xlabel(dep_variable)
% % ylabel('relative Prob Type 2 Persistence')
% % legend('MEC','LEC')
% % 
% % plot_mec_lec_comparison(cond_nprob_p1,bincenters,mec_cells,lec_cells)
% % xlabel(dep_variable)
% % ylabel('relative Prob Type 1 Persistence')
% % legend('MEC','LEC')
% % 
% % plot_mec_lec_comparison(cond_nmean_updur,bincenters,mec_cells,lec_cells)    
% % xlabel(dep_variable)
% % ylabel('relative MP up state duration (s)')
% % legend('MEC','LEC')
% % 
% % 
% % plot_mec_lec_comparison(cond_nmean_upcdur,bincenters,mec_cells,lec_cells)
% % xlabel(dep_variable)
% % ylabel('relative MP up state duration (LFP cycles)')
% % legend('MEC','LEC')
% % 
% % plot_mec_lec_comparison(cond_nmean_uplag,bincenters,mec_cells,lec_cells)
% % xlabel(dep_variable)
% % ylabel('relative up lag')
% % legend('MEC','LEC')
% % 
% % plot_mec_lec_comparison(cond_nmean_downlag,bincenters,mec_cells,lec_cells)
% % xlabel(dep_variable)
% % ylabel('relative down lag')
% % legend('MEC','LEC')
% % 
% % plot_mec_lec_comparison(cond_nmean_upamp,bincenters,mec_cells,lec_cells)
% % xlabel(dep_variable)
% % ylabel('relative up amplitude')
% % legend('MEC','LEC')
% % 
% % plot_mec_lec_comparison(cond_nmean_downamp,bincenters,mec_cells,lec_cells)
% % xlabel(dep_variable)
% % ylabel('relative down amplitude')
% % legend('MEC','LEC')
% % 
% % plot_mec_lec_comparison(cond_nmean_upstd,bincenters,mec_cells,lec_cells)
% % xlabel(dep_variable)
% % ylabel('relative up std')
% % legend('MEC','LEC')
% % 
% % plot_mec_lec_comparison(cond_nmean_downstd,bincenters,mec_cells,lec_cells)
% % xlabel(dep_variable)
% % ylabel('relative down std')
% % legend('MEC','LEC')
% % 
% % plot_mec_lec_comparison(cond_nmean_uprate,bincenters,mec_cells,lec_cells)
% % xlabel(dep_variable)
% % ylabel('relative up rate')
% % legend('MEC','LEC')