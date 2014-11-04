clear all
close all
addpath('C:\WC_Germany\parietal_cortical_2010\')
addpath('C:\WC_Germany\hsmm_state_detection\')
addpath('C:\WC_Germany\persistent_revised\')

cd C:\WC_Germany\sven_thomas_combined\
load ./combined_dir_nd.mat
load ./combined_core_analysis_fin_nd.mat
load ./combined_depthofanesth_data_fin_nd_v2.mat
% load ./combined_depthofanesth_data_fin_nd.mat

uset = sort([l3mec l3lec]);
all_cells = 1:length(combined_dir);
l3mec = find(ismember(all_cells(uset),l3mec));
l3lec = find(ismember(all_cells(uset),l3lec));
combined_dir = combined_dir(uset);
hpc_mua = hpc_mua(uset);
hpc_lfp = hpc_lfp(uset);
ctx_lfp = ctx_lfp(uset);

dep_variable = 'Frequency';
% dep_variable = 'Duty Cycle';
% dep_variable = 'Down Duration';
% dep_variable = 'Power';

nbins = 10;
if strcmp(dep_variable,'Duty Cycle')
    binedges = linspace(0.3,0.6,nbins+1); %duty cycle
    xl = [0.3 0.6];
elseif strcmp(dep_variable,'Frequency')
    binedges = linspace(0.2,0.6,nbins+1); %for freq
    xl = [0.2 0.6];
elseif strcmp(dep_variable,'Down Duration')
    binedges = linspace(0.7,2.5,nbins+1);
    xl = [0.7 2.5];
elseif strcmp(dep_variable,'Power')
%     binedges = linspace(-2.5,2.5,nbins+1);
    binedges = linspace(-3,2,nbins+1);
    xl = [-3 2];
end
bincenters = (binedges(1:end-1)+binedges(2:end))/2;

lf8_up_fract = 0.5;

for d = 1:length(combined_dir)
  
    mp_updurs = mp_state_durations{d}{2};
    mp_downdurs = mp_state_durations{d}{1};
    lf8_updurs = lfp_state_durations{d}{2};
    lf8_downdurs = lfp_state_durations{d}{1};
    
    if strcmp(dep_variable,'Duty Cycle')
        pred_var = lf8_dutycyc{d};
    elseif strcmp(dep_variable,'Frequency')
        pred_var = lf8_udsfreq{d};
    elseif strcmp(dep_variable,'Down Duration')
        pred_var = lf8_meddowndur{d};
    elseif strcmp(dep_variable,'Power')
        pred_var = zscore(lf8_udspow{d});
    end
    
    [counts,binlocs] = histc(pred_var,binedges);
    
    cond_prob_p2(d,:) = nan(1,length(bincenters));
    cond_prob_p1(d,:) = nan(1,length(bincenters));
    cond_mean_updur(d,:) = nan(1,length(bincenters));
    cond_mean_upcdur(d,:) = nan(1,length(bincenters));
    cond_mean_uplag(d,:) = nan(1,length(bincenters));
    cond_mean_downlag(d,:) = nan(1,length(bincenters));
    cond_mean_reluplag(d,:) = nan(1,length(bincenters));
    cond_mean_reldownlag(d,:) = nan(1,length(bincenters));
    cond_mean_relupdur(d,:) = nan(1,length(bincenters));
    
    rt2_up_ind = nan(size(lf8_udsfreq{d}));
    rt2_up_ind(rt2_ups{d}) = 1;
    rt2_up_ind(nrt2_ups{d}) = 0;
    
    cor_type = 'pearson';
    mean_uds_pow(d) = nanmean(lf8_udspow{d});
    mean_uds_freq(d) = nanmean(lf8_udsfreq{d});
    
        use_pts = find(~isnan(lf8_udspow{d}) & ~isnan(mp_downlags{d}));
    [a,b] = corr(lf8_udspow{d}(use_pts),mp_downlags{d}(use_pts),'type',cor_type);
    udspow_downlag_corr(d) = a;
    udspow_downlag_p(d) = b;
    [a,b] = corr(lf8_udspow{d}(use_pts),mp_reldownlags{d}(use_pts),'type',cor_type);
    udspow_reldownlag_corr(d) = a;
    udspow_reldownlag_p(d) = b;
    use_pts = find(~isnan(lf8_udspow{d}) & ~isnan(mp_uplags{d}));
    [a,b] = corr(lf8_udspow{d}(use_pts),mp_uplags{d}(use_pts),'type',cor_type);
    udspow_uplag_corr(d) = a;
    udspow_uplag_p(d) = b;
    [a,b] = corr(lf8_udspow{d}(use_pts),mp_reluplags{d}(use_pts),'type',cor_type);
    udspow_reluplag_corr(d) = a;
    udspow_reluplag_p(d) = b;
    use_pts = find(~isnan(lf8_udspow{d}) & ~isnan(rt2_up_ind));
    [a,b] = corr(lf8_udspow{d}(use_pts),rt2_up_ind(use_pts),'type',cor_type);
% [B,dev,stats] = glmfit(lf8_udspow{d}(use_pts),rt2_up_ind(use_pts),'binomial');
    udspow_rt2_corr(d) = a;
    udspow_rt2_p(d) = b;

    use_pts = find(~isnan(lf8_udsfreq{d}) & ~isnan(mp_downlags{d}));
    a = corr(lf8_udsfreq{d}(use_pts),mp_downlags{d}(use_pts),'type','spearman');
    udsfreq_downlag_corr(d) = a;
%     udsfreq_downlag_p(d) = b;
    a = corr(lf8_udsfreq{d}(use_pts),mp_reldownlags{d}(use_pts),'type','spearman');
    udsfreq_reldownlag_corr(d) = a;
%     udsfreq_reldownlag_p(d) = b;
    use_pts = find(~isnan(lf8_udsfreq{d}) & ~isnan(mp_uplags{d}));
    a = corr(lf8_udsfreq{d}(use_pts),mp_uplags{d}(use_pts),'type','spearman');
    udsfreq_uplag_corr(d) = a;
%     udsfreq_uplag_p(d) = b;
    a = corr(lf8_udsfreq{d}(use_pts),mp_reluplags{d}(use_pts),'type','spearman');
    udsfreq_reluplag_corr(d) = a;
%     udsfreq_reluplag_p(d) = b;
    use_pts = find(~isnan(lf8_udsfreq{d}) & ~isnan(rt2_up_ind));
    a = corr(lf8_udsfreq{d}(use_pts),rt2_up_ind(use_pts),'type','spearman');
% [B,dev,stats] = glmfit(lf8_udsfreq{d}(use_pts),rt2_up_ind(use_pts),'binomial');
    udsfreq_rt2_corr(d) = a;
%     udsfreq_rt2_p(d) = b;
%     
    use_pts = find(~isnan(lf8_dutycyc{d}) & ~isnan(mp_downlags{d}));
    a = corr(lf8_dutycyc{d}(use_pts),mp_downlags{d}(use_pts),'type','spearman');
    udsdc_downlag_corr(d) = a;
%     udsdc_downlag_p(d) = b;
    use_pts = find(~isnan(lf8_dutycyc{d}) & ~isnan(mp_reldownlags{d}));
    a = corr(lf8_dutycyc{d}(use_pts),mp_reldownlags{d}(use_pts),'type','spearman');
    udsdc_reldownlag_corr(d) = a;
%     udsdc_reldownlag_p(d) = b;
    use_pts = find(~isnan(lf8_dutycyc{d}) & ~isnan(mp_uplags{d}));
    a = corr(lf8_dutycyc{d}(use_pts),mp_uplags{d}(use_pts),'type','spearman');
    udsdc_uplag_corr(d) = a;
%     udsdc_uplag_p(d) = b;
    a = corr(lf8_dutycyc{d}(use_pts),mp_reluplags{d}(use_pts),'type','spearman');
    udsdc_reluplag_corr(d) = a;
%     udsdc_reluplag_p(d) = b;
    use_pts = find(~isnan(lf8_dutycyc{d}) & ~isnan(rt2_up_ind));
    a = corr(lf8_dutycyc{d}(use_pts),rt2_up_ind(use_pts),'type','spearman');
% [B,dev,stats] = glmfit(lf8_dutycyc{d}(use_pts),rt2_up_ind(use_pts),'binomial');
    udsdc_rt2_corr(d) = a;
%     udsdc_rt2_p(d) = b;

    use_pts = find(~isnan(lf8_meddowndur{d}) & ~isnan(mp_downlags{d}));
    a = corr(lf8_meddowndur{d}(use_pts),mp_downlags{d}(use_pts),'type','spearman');
    downdur_downlag_corr(d) = a;
%     downdur_downlag_p(d) = b;
    use_pts = find(~isnan(lf8_meddowndur{d}) & ~isnan(mp_reldownlags{d}));
    a = corr(lf8_meddowndur{d}(use_pts),mp_reldownlags{d}(use_pts),'type','spearman');
    downdur_reldownlag_corr(d) = a;
%     downdur_reldownlag_p(d) = b;
    use_pts = find(~isnan(lf8_meddowndur{d}) & ~isnan(mp_uplags{d}));
    a = corr(lf8_meddowndur{d}(use_pts),mp_uplags{d}(use_pts),'type','spearman');
    downdur_uplag_corr(d) = a;
%     downdur_uplag_p(d) = b;
    a = corr(lf8_meddowndur{d}(use_pts),mp_reluplags{d}(use_pts),'type','spearman');
    downdur_reluplag_corr(d) = a;
%     downdur_reluplag_p(d) = b;
    use_pts = find(~isnan(lf8_meddowndur{d}) & ~isnan(rt2_up_ind));
    a = corr(lf8_meddowndur{d}(use_pts),rt2_up_ind(use_pts),'type','spearman');
% [B,dev,stats] = glmfit(lf8_meddowndur{d}(use_pts),rt2_up_ind(use_pts),'binomial');
    downdur_rt2_corr(d) = a;
%     downdur_rt2_p(d) = b;

    for i = 1:length(bincenters)
        cur_points = find(binlocs==i);
        num_counts(d,i) = length(cur_points);
        if length(cur_points) > 2
            n_rpers_ups = sum(ismember(rt2_ups{d},cur_points));
            n_npers_ups = sum(ismember(nrt2_ups{d},cur_points));
            cond_prob_p2(d,i) = n_rpers_ups/(n_rpers_ups+n_npers_ups);

            cond_prob_p1(d,i) = nansum(mp_rel_updurs{d}(cur_points) > 0)/length(cur_points);
            cond_mean_updur(d,i) = nanmean(mp_updurs(cur_points));
%             cond_mean_upcdur(d,i) = nanmean(mp_updurs_lfpc{d}(cur_points));
            cond_mean_uplag(d,i) = nanmean(mp_uplags{d}(cur_points));
            cond_mean_downlag(d,i) = nanmean(mp_downlags{d}(cur_points));
           cond_mean_reluplag(d,i) = nanmean(mp_reluplags{d}(cur_points));
            cond_mean_reldownlag(d,i) = nanmean(mp_reldownlags{d}(cur_points));
%             cur_nt2s = find(ismember(nrt2_ups{d},cur_points));
%             if ~isempty(cur_nt2s)
%             cond_mean_relupdur(d,i) = nanmean(mp_rel_updurs_nt2{d}(cur_nt2s));
%             end
        end
    end
    fract_rt1_ups(d) = nansum(mp_rel_updurs{d} > 0)/sum(~isnan(mp_rel_updurs{d}));
    cond_nprob_p2(d,:) = cond_prob_p2(d,:)/fract_rt2_ups(d);
    cond_nprob_p1(d,:) = cond_prob_p1(d,:)/fract_rt1_ups(d);
    cond_nmean_updur(d,:) = cond_mean_updur(d,:)/nanmean(mp_updurs);
%     cond_nmean_upcdur(d,:) = cond_mean_upcdur(d,:)/nanmean(mp_updurs_lfpc{d});
    cond_nmean_uplag(d,:) = cond_mean_uplag(d,:) - nanmean(mp_uplags{d});
    cond_nmean_downlag(d,:) = cond_mean_downlag(d,:) - nanmean(mp_downlags{d});
end

%% Unnormalized plots    
f = plot_mec_lec_comparison(100*cond_prob_p2,bincenters,l3mec,l3lec)
set(gca,'fontsize',14,'fontname','arial')
xlabel(dep_variable,'fontsize',16,'fontname','arial')
ylabel('Percentage Type 2 Persistence','fontsize',16,'fontname','arial')
legend('MEC','LEC')
xlim(xl)
ylim([0 25])

% plot_mec_lec_comparison(cond_prob_p1,bincenters,l3mec,l3lec)
% set(gca,'fontsize',14,'fontname','arial')
% xlabel(dep_variable,'fontsize',16,'fontname','arial')
% ylabel('Prob Type 1 Persistence','fontsize',16,'fontname','arial')
% legend('MEC','LEC')
% xlim(xl)
% ylim([0.2 1])

% plot_mec_lec_comparison(cond_mean_downlag/Fsd,bincenters,l3mec,l3lec)
% set(gca,'fontsize',14,'fontname','arial')
% xlabel(dep_variable,'fontsize',16,'fontname','arial')
% ylabel('Prob Type 1 Persistence','fontsize',16,'fontname','arial')
% legend('MEC','LEC')
% xlim(xl)
% ylim([0.2 1])
% 
% plot_mec_lec_comparison(cond_mean_reldownlag,bincenters,l3mec,l3lec)
% set(gca,'fontsize',14,'fontname','arial')
% xlabel(dep_variable,'fontsize',16,'fontname','arial')
% ylabel('Prob Type 1 Persistence','fontsize',16,'fontname','arial')
% legend('MEC','LEC')
% xlim(xl)
% ylim([0.2 1])
% 
% plot_mec_lec_comparison(cond_mean_updur,bincenters,l3mec,l3lec)
% set(gca,'fontsize',14,'fontname','arial')
% xlabel(dep_variable,'fontsize',16,'fontname','arial')
% ylabel('Mean MP up state duration (s)','fontsize',16,'fontname','arial')
% legend('MEC','LEC')
% xlim(xl)
% ylim([0.5 3])
% 
% % plot_mec_lec_comparison(cond_mean_relupdur,bincenters,l3mec,l3lec)
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