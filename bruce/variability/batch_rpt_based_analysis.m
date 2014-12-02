close all
clear all
clc

% fig_dir = '/Users/james/Analysis/bruce/variability/figures/';
fig_dir = '/home/james/Analysis/bruce/variability/figures/';
base_sname = 'rpt_variability_analysis4';
base_tname = 'model_variability_analysis';

all_SU_data = [];
n_probes = 24;

%% LOAD JBE
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296','M297'};
ori_list = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 160 nan; 40 nan; 45 nan; 0 90];
rmfield_list = {};

for ee = 1:length(Expt_list)
% for ee = 1:7
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
    cur_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
    
    clear ori_data
    for ii =  1:2
        if ~isnan(ori_list(ee,ii))
            %load trig avg data
            tname = strcat(cur_dir,base_sname,sprintf('_ori%d',ori_list(ee,ii)));
            load(tname);
            tname = strcat(cur_dir,base_tname,sprintf('_ori%d',ori_list(ee,ii)));
            load(tname);

                        for jj = 1:length(Rpt_Data)
                Rpt_Data(jj).Mod_EP = EP_data(jj);
            end

            [Rpt_Data.expt_num] = deal(Expt_num);
            [Rpt_Data.bar_ori] = deal(ori_list(ee,ii));
            [Rpt_Data.animal] = deal('jbe');
            
            ori_data(ii,:) = Rpt_Data((n_probes+1):end);
        end
    end
    
    [cur_Noris,cur_Nunits] = size(ori_data);
    used = arrayfun(@(x) isstruct(x.ModData),ori_data);
    no_data = find(arrayfun(@(x) isempty(x.ModData.unit_data),ori_data(used)));
    uset = find(used);
    used(uset(no_data)) = false;
    xvLLimps = nan(cur_Noris,cur_Nunits);
    xvLLimps(used) = arrayfun(@(x)x.ModData.rectGQM.xvLLimp,ori_data(used));
    avg_rates = nan(cur_Noris,cur_Nunits);
    avg_rates(used) = arrayfun(@(x)x.ModData.unit_data.avg_rate,ori_data(used));
    xvLLrates = xvLLimps.*avg_rates;
    
    [mvals,mlocs] = max(xvLLrates,[],1);
    for ii = 1:cur_Nunits
        if ~isnan(xvLLrates(mlocs(ii),ii))
            all_SU_data = cat(1,all_SU_data,ori_data(mlocs(ii),ii));
        end
    end
end

%% SELECT USABLE CELLS
dt = 0.01;

%selection criteria
min_rate = 5; % min avg rate in Hz (5)
min_xvLLimp = 0.0; %(0.05);
min_rpt_trials = 30;

tot_Nunits = length(all_SU_data);
avg_rates = arrayfun(@(x) x.ModData.unit_data.avg_rate,all_SU_data);
tot_spikes = arrayfun(@(x) x.ModData.unit_data.tot_spikes,all_SU_data);
rec_dur = arrayfun(@(x) x.ModData.unit_data.N_used_samps,all_SU_data)*dt/60;
mod_xvLLimps = arrayfun(@(x) x.ModData.rectGQM.xvLLimp,all_SU_data);
expt_nums = [all_SU_data(:).expt_num];
expt_oris = [all_SU_data(:).bar_ori];

clust_iso_dist = arrayfun(@(x) x.ModData.unit_data.SU_isodist,all_SU_data);
clust_Lratio = arrayfun(@(x) x.ModData.unit_data.SU_Lratio,all_SU_data);
clust_refract = arrayfun(@(x) x.ModData.unit_data.SU_refract,all_SU_data);
clust_dprime = arrayfun(@(x) x.ModData.unit_data.SU_dprime,all_SU_data);
rate_stability_cv = arrayfun(@(x) x.ModData.unit_data.rate_stability_cv,all_SU_data);
dprime_stability_cv = arrayfun(@(x) x.ModData.unit_data.dprime_stability_cv,all_SU_data);

jbe_SUs = find(strcmp('jbe',{all_SU_data(:).animal}));
lem_SUs = find(strcmp('lem',{all_SU_data(:).animal}));

n_rpt_trials = arrayfun(@(x) x.n_utrials,all_SU_data);

%%
cur_SUs = find(avg_rates >= min_rate & mod_xvLLimps > min_xvLLimp & n_rpt_trials >= min_rpt_trials);

actual_EP_SD = arrayfun(@(x) x.Mod_EP.actual_EP_SD,all_SU_data(cur_SUs));
alpha_funs = cell2mat(arrayfun(@(x) x.Mod_EP.alpha_funs, all_SU_data(cur_SUs),'uniformoutput',0));
spline_sig_var = cell2mat(arrayfun(@(x) x.spline_resp_ZPT, all_SU_data(cur_SUs),'uniformoutput',0));
rpt_psth_var = cell2mat(arrayfun(@(x) x.rand_psth_var, all_SU_data(cur_SUs),'uniformoutput',0));
spline_alpha = rpt_psth_var./spline_sig_var;
ms_resp_var = cell2mat(arrayfun(@(x) x.ms_resp_var, all_SU_data(cur_SUs),'uniformoutput',0));

SNR = spline_sig_var./ms_resp_var;

mod_alphas = nan(length(cur_SUs),1);
for ss = 1:length(cur_SUs)
    mod_alphas(ss) = interp1(poss_SDs,alpha_funs(ss,:),actual_EP_SD(ss));
end

%%
% min_psth_var = 0;
% used = find(rpt_psth_var >= min_psth_var);
used = 1:length(cur_SUs);

msize = 10;
f1 = figure();
plot(mod_alphas(used),spline_alpha(used),'.','markersize',msize);
line([0 1],[0 1],'color','k');
xlabel('Model-predicted alpha');
ylabel('Estimated alpha');
set(gca,'xtick',[0:0.2:1],'ytick',[0:0.2:1]);
% 
% %%PRINT FIGURE
% fig_width = 3.5; rel_height = 0.9;
% figufy(f1);
% fname = [fig_dir 'Mod_spline_alphas.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%%
cur_SUs = find(avg_rates >= min_rate & mod_xvLLimps > min_xvLLimp & n_rpt_trials >= min_rpt_trials);

maxtlag = 10;

ulags = find(abs(tlags) <= 5);

all_psth_noise_corrs = [];
all_spline_noise_corrs = [];
all_sig_corrs = [];
all_psth_corrs = [];
all_pair_IDS = [];
all_expt_IDS = [];
for ss = 1:length(cur_SUs)
    cur_sig_covs = squeeze(all_SU_data(cur_SUs(ss)).spline_xcov(:,maxtlag+1));
    cur_psth_covs = squeeze(all_SU_data(cur_SUs(ss)).psth_xcov(:,maxtlag+1));
    cur_emp_covs = squeeze(all_SU_data(cur_SUs(ss)).emp_xcov(:,maxtlag+1));
%     cur_sig_covs = squeeze(mean(all_SU_data(cur_SUs(ss)).spline_xcov(:,ulags),2));
% %     cur_psth_covs = squeeze(mean(all_SU_data(cur_SUs(ss)).rand_xcov(:,ulags),2));
%     cur_psth_covs = squeeze(mean(all_SU_data(cur_SUs(ss)).psth_xcov(:,ulags),2));
%     cur_emp_covs = squeeze(mean(all_SU_data(cur_SUs(ss)).emp_xcov(:,ulags),2));
    
    psth_noise_covs = cur_emp_covs - cur_psth_covs;
    spline_noise_covs = cur_emp_covs - cur_sig_covs;
    
    cur_noise_norms = all_SU_data(cur_SUs(ss)).noise_varnorm;
    cur_sig_norms = all_SU_data(cur_SUs(ss)).psth_varnorm;
    cur_snoise_norms = all_SU_data(cur_SUs(ss)).spline_noise_varnorm;
    cur_ssig_norms = all_SU_data(cur_SUs(ss)).spline_varnorm;
    cur_tot_norms = all_SU_data(cur_SUs(ss)).tot_varnorm;
    
%     cur_sig_norms = all_SU_data(cur_SUs(ss)).tot_varnorm;
    
    cur_ssig_norms(cur_ssig_norms <= 0) = nan;
    cur_snoise_norms(cur_snoise_norms <= 0) = nan;
    cur_sig_norms(cur_sig_norms <= 0) = nan;
    cur_noise_norms(cur_noise_norms <= 0) = nan;

    cur_pair_IDS = [repmat(all_SU_data(cur_SUs(ss)).unit_num,length(cur_sig_norms),1) (1:length(cur_sig_norms))'];
    
    all_psth_noise_corrs = cat(1,all_psth_noise_corrs,bsxfun(@rdivide,psth_noise_covs,cur_noise_norms));
    all_spline_noise_corrs = cat(1,all_spline_noise_corrs,bsxfun(@rdivide,spline_noise_covs,cur_noise_norms));
    all_sig_corrs = cat(1,all_sig_corrs,bsxfun(@rdivide,cur_sig_covs,cur_sig_norms));
    all_psth_corrs = cat(1,all_psth_corrs,bsxfun(@rdivide,cur_psth_covs,cur_sig_norms));
%     all_psth_noise_corrs = cat(1,all_psth_noise_corrs,bsxfun(@rdivide,psth_noise_covs,cur_tot_norms));
%     all_spline_noise_corrs = cat(1,all_spline_noise_corrs,bsxfun(@rdivide,spline_noise_covs,cur_tot_norms));
%     all_sig_corrs = cat(1,all_sig_corrs,bsxfun(@rdivide,cur_sig_covs,cur_tot_norms));
%     all_psth_corrs = cat(1,all_psth_corrs,bsxfun(@rdivide,cur_psth_covs,cur_tot_norms));
    
    all_pair_IDS = cat(1,all_pair_IDS,cur_pair_IDS);
    all_expt_IDS = cat(1,all_expt_IDS,repmat(all_SU_data(cur_SUs(ss)).expt_num,length(cur_sig_norms),1));
end
bad = find(all_sig_corrs == 1);
all_sig_corrs(bad) = nan;
all_psth_corrs(bad) = nan;
all_spline_noise_corrs(bad) = nan;
all_psth_noise_corrs(bad) = nan;


% close all
uset = find(all_pair_IDS(:,1) > all_pair_IDS(:,2));
% uset = find(all_pair_IDS(:,1) ~= all_pair_IDS(:,2) & min(all_pair_IDS,[],2) > 24 & ~isnan(all_sig_corrs) & ~isnan(all_spline_noise_corrs) & ~isnan(all_psth_corrs));

% poss_expts = [289];
poss_expts = [266 270 275 277 281 287 289 294 296 297];
uset(~ismember(all_expt_IDS(uset),poss_expts)) = [];

all_sig_corrs = all_sig_corrs(uset);
all_psth_corrs = all_psth_corrs(uset);
all_spline_noise_corrs = all_spline_noise_corrs(uset);
all_psth_noise_corrs = all_psth_noise_corrs(uset);

xx = linspace(-0.3,1,100);
msize = 6;
f1 = figure(); hold on
plot(all_psth_corrs,all_psth_noise_corrs,'k.','markersize',msize); 
r = robustfit(all_psth_corrs,all_psth_noise_corrs);
plot(xx,r(1)+r(2)*xx,'b','linewidth',2)
% xlim([-0.3 1]);
% ylim([-0.05 0.2]);
% set(gca,'xtick',[-0.2:0.2:1]);
% 

xx = linspace(-0.3,1,100);
msize = 6;
f1 = figure(); hold on
plot(all_sig_corrs,all_spline_noise_corrs,'r.','markersize',msize); 
r2 = robustfit(all_sig_corrs,all_spline_noise_corrs);
% r2 = regress(all_spline_noise_corrs,[ones(length(all_sig_corrs),1) all_sig_corrs]);
plot(xx,r2(1)+r2(2)*xx,'g','linewidth',2)
% xlim([-0.3 1]);
% ylim([-0.05 0.2]);
% set(gca,'xtick',[-0.2:0.2:1]);


% % %PRINT FIGURE
% fig_width = 3.5; rel_height = 0.9;
% figufy(f1);
% fname = [fig_dir 'Sigc_Noisec_SU_scatter.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%%
pbins = 0:10:100;
bin_edges = prctile(all_sig_corrs,pbins);
spline_cond_noisecorrs = nan(length(bin_edges)-1,3);
spline_cond_sigcorrs = nan(length(bin_edges)-1,3);
for ii = 1:length(bin_edges)-1
    curset = find(all_sig_corrs >= bin_edges(ii) & all_sig_corrs < bin_edges(ii+1));
    spline_cond_noisecorrs(ii,:) = prctile(all_spline_noise_corrs(curset),[25 50 75]);
    spline_cond_sigcorrs(ii,:) = prctile(all_sig_corrs(curset),[25 50 75]);
end
x = spline_cond_sigcorrs(:,2); lerrx = spline_cond_sigcorrs(:,1); rerrx = spline_cond_sigcorrs(:,3);
y = spline_cond_noisecorrs(:,2); lerry = y - spline_cond_noisecorrs(:,1); rerry = spline_cond_noisecorrs(:,3) - y;


bin_edges = prctile(all_psth_corrs,pbins);
psth_cond_noisecorrs = nan(length(bin_edges)-1,3);
psth_cond_sigcorrs = nan(length(bin_edges)-1,3);
for ii = 1:length(bin_edges)-1
    curset = find(all_psth_corrs >= bin_edges(ii) & all_psth_corrs < bin_edges(ii+1));
    psth_cond_noisecorrs(ii,:) = prctile(all_psth_noise_corrs(curset),[25 50 75]);
    psth_cond_sigcorrs(ii,:) = prctile(all_psth_corrs(curset),[25 50 75]);
end
x2 = psth_cond_sigcorrs(:,2); 
y2 = psth_cond_noisecorrs(:,2); lerry2 = y2 - psth_cond_noisecorrs(:,1); rerry2 = psth_cond_noisecorrs(:,3) - y2;

f1 = figure();
errorbar(x,y,lerry,rerry)
hold on
errorbar(x2,y2,lerry2,rerry2,'r')
xlim([-0.2 1]);
ylim([-0.02 0.2])
xlabel('Signal Correlation');
ylabel('Noise Correlation');
set(gca,'xtick',[-0.2:0.2:1]);

% % %PRINT FIGURE
% fig_width = 3.5; rel_height = 0.9;
% figufy(f1);
% fname = [fig_dir 'Sigc_Noisec_MUs.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%%
cur_SUs = find(avg_rates >= min_rate & mod_xvLLimps > min_xvLLimp & n_rpt_trials >= min_rpt_trials);

maxtlag = 10;

all_pnoise_corrs = [];
all_snoise_corrs = [];
all_sig_corrs = [];
all_psth_corrs = [];
all_trial_corrs = [];
all_pair_IDS = [];
all_pair_avgrates = [];
all_pair_ntrials = [];
all_expt_IDS = [];
for ss = 1:length(cur_SUs)
    cur_sig_covs = all_SU_data(cur_SUs(ss)).spline_cnt_cov;
    cur_psth_covs = all_SU_data(cur_SUs(ss)).psth_cnt_cov;
    cur_emp_covs = all_SU_data(cur_SUs(ss)).obs_cnt_cov;
    
    cur_snoise_covs = cur_emp_covs - cur_sig_covs;
    cur_pnoise_covs = cur_emp_covs - cur_psth_covs;
    
    cur_trial_corrs = all_SU_data(cur_SUs(ss)).trial_corrs;
    all_trial_corrs = cat(1,all_trial_corrs,cur_trial_corrs(:));
    
    psth_norm = all_SU_data(cur_SUs(ss)).psth_cnt_varnorm;
    pnoise_norm = all_SU_data(cur_SUs(ss)).psthnoise_cnt_varnorm;
    
    psth_norm1 = repmat(psth_norm(:,1),1,length(poss_cnt_wins));
    pnoise_norm1 = repmat(pnoise_norm(:,1),1,length(poss_cnt_wins));
    psth_norm(psth_norm < psth_norm1) = psth_norm1(psth_norm < psth_norm1);
    pnoise_norm(pnoise_norm < pnoise_norm1) = pnoise_norm1(pnoise_norm < pnoise_norm1);
    
    if any(~isreal(psth_norm(:)))
        fprintf('Warning, nonreal psth norms\n');
    end
    if any(~isreal(pnoise_norm(:)))
        fprintf('Warning, nonreal pnoise norms\n');
    end
    
    cur_psth_corrs = cur_psth_covs./psth_norm;
    cur_sig_corrs = cur_sig_covs./psth_norm;
    cur_pnoise_corrs = cur_pnoise_covs./pnoise_norm;
    cur_snoise_corrs = cur_snoise_covs./pnoise_norm;
%     cur_psth_corrs = bsxfun(@rdivide,cur_psth_covs,psth_norm(:,1));
%     cur_sig_corrs = bsxfun(@rdivide,cur_sig_covs,psth_norm(:,1));
%     cur_pnoise_corrs = bsxfun(@rdivide,cur_pnoise_covs,pnoise_norm(:,1));
%     cur_snoise_corrs = bsxfun(@rdivide,cur_snoise_covs,pnoise_norm(:,1));
    
    cur_pair_IDS = [repmat(all_SU_data(cur_SUs(ss)).unit_num,size(cur_psth_corrs,1),1) (1:size(cur_psth_corrs,1))'];
    
    all_pnoise_corrs = cat(1,all_pnoise_corrs,cur_pnoise_corrs);
    all_snoise_corrs = cat(1,all_snoise_corrs,cur_snoise_corrs);
    all_sig_corrs = cat(1,all_sig_corrs,cur_sig_corrs);
    all_psth_corrs = cat(1,all_psth_corrs,cur_psth_corrs);
    
    all_pair_IDS = cat(1,all_pair_IDS,cur_pair_IDS);
    all_expt_IDS = cat(1,all_expt_IDS,repmat(all_SU_data(cur_SUs(ss)).expt_num,size(cur_psth_corrs,1),1));
    all_pair_ntrials = cat(1,all_pair_ntrials,all_SU_data(cur_SUs(ss)).otherunit_nutrials);
    all_pair_avgrates = cat(1,all_pair_avgrates,all_SU_data(cur_SUs(ss)).otherunit_avgrates);
end

pair_min_rate = 1;

% uset = find(all_pair_IDS(:,1) > all_pair_IDS(:,2));
uset = find(all_pair_IDS(:,1) > all_pair_IDS(:,2) & all_pair_ntrials >= min_rpt_trials & ...
    all_pair_avgrates/dt >= pair_min_rate);

% uset = find(all_pair_IDS(:,1) > all_pair_IDS(:,2) & min(all_pair_IDS,[],2) > 24);
% uset = find(all_pair_IDS(:,1) > all_pair_IDS(:,2) & min(all_pair_IDS,[],2) > 24 & ...
%     all_pair_ntrials >= min_rpt_trials & all_pair_avgrates/dt >= pair_min_rate);

uu = uset(~isnan(all_psth_corrs(uset,1)) & ~isnan(all_sig_corrs(uset,1)) & ~isnan(all_pnoise_corrs(uset,1)) & ~isnan(all_snoise_corrs(uset,1)));



%%
% yl = [-0.1 0.5];
xl1 = [-0.5 1];
% xl2 = [-0.5 1.5];
xl2 = [-0.5 2];
% yl = [-0.05 0.3];
yl = [-0.2 0.3];
% xl1 = [-0.4 1];
% xl2 = [-0.5 1.5];
sm_span = 201;
% sm_span = 51;
% sm_span = 101;
rmpath('~/James_scripts/bruce/bruce_code/')

use_win = 1;

% close all
xx = all_psth_corrs(uu,use_win); yy = all_pnoise_corrs(uu,use_win);
[psth_spear,psth_p] = corr(xx,yy,'type','spearman');
[psth_pear] = corr(xx,yy,'type','pearson');
f1 = figure();
subplot(1,2,1);
plot(xx,yy,'r.')
hold on
[xs,ord] = sort(xx); yy = yy(ord);
ysm = smooth(xs,yy,sm_span,'rlowess');
plot(xs,ysm,'k','linewidth',2);
xlim(xl1);
ylim(yl);
% axis tight

xx2 = all_sig_corrs(uu,use_win); yy2 = all_snoise_corrs(uu,use_win);
[sig_spear,sig_p] = corr(xx2,yy2,'type','spearman');
[sig_pear] = corr(xx2,yy2,'type','pearson');
subplot(1,2,2);
plot(xx2,yy2,'.')
hold on
[xs2,ord] = sort(xx2); yy2 = yy2(ord);
ysm2 = smooth(xs2,yy2,sm_span,'rlowess');
plot(xs2,ysm2,'k','linewidth',2);
xlim(xl2);
ylim(yl);
% axis tight

% f2 = figure();
% nc_range = linspace(-0.2,0.3,50);
% psth_h = histc(all_pnoise_corrs(uu,use_win),nc_range);
% psth_h = psth_h/length(uu);
% sig_h = histc(all_snoise_corrs(uu,use_win),nc_range);
% sig_h = sig_h/length(uu);
% stairs(nc_range,psth_h); hold on
% stairs(nc_range,sig_h,'r');

% f3 = figure(); hold on
% plot(xs,ysm,'r--');
% plot(xs2,ysm2,'b--');

% % %PRINT FIGURE
% fig_width = 7; rel_height = 0.4;
% figufy(f1);
% fname = [fig_dir 'Sigc_Noisec_MUs2.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);

%%
sm_span = 101;

un_expt_ids = unique(all_expt_IDS);
expt_psth_ss = nan(length(un_expt_ids),length(poss_cnt_wins));
expt_sig_ss = nan(length(un_expt_ids),length(poss_cnt_wins));
for cm = 1:length(poss_cnt_wins)
%     psth_ss(cm) = corr(all_psth_corrs(uu,cm),all_pnoise_corrs(uu,cm),'type','spearman');
%     sig_ss(cm) = corr(all_sig_corrs(uu,cm),all_snoise_corrs(uu,cm),'type','spearman');
    psth_ss(cm) = corr(all_psth_corrs(uu,cm),all_pnoise_corrs(uu,cm),'type','pearson');
    sig_ss(cm) = corr(all_sig_corrs(uu,cm),all_snoise_corrs(uu,cm),'type','pearson');

%     xx = all_psth_corrs(uu,cm); yy = all_pnoise_corrs(uu,cm);
%     [xs,ord] = sort(xx); yy = yy(ord);
%     ysm = smooth(xs,yy,sm_span,'rlowess');
%     ncc(cm) = range(ysm);
%     
%     xx = all_sig_corrs(uu,cm); yy = all_snoise_corrs(uu,cm);
%     [xs,ord] = sort(xx); yy = yy(ord);
%     ysm = smooth(xs,yy,sm_span,'rlowess');
%     ncc2(cm) = range(ysm);
%     for jj = 1:length(un_expt_ids)
%         curset = uu(all_expt_IDS(uu) == un_expt_ids(jj));
%         expt_psth_ss(jj,cm) = corr(all_psth_corrs(curset,cm),all_pnoise_corrs(curset,cm),'type','spearman');
%         expt_sig_ss(jj,cm) = corr(all_sig_corrs(curset,cm),all_snoise_corrs(curset,cm),'type','spearman');
%     end
end


msize = 6;
f1 = figure();
plot((poss_cnt_wins+1)*dt*1e3,psth_ss,'ko-','markersize',msize)
hold on
plot((poss_cnt_wins+1)*dt*1e3,sig_ss,'ro-','markersize',msize);
xlabel('Time window (ms)');
ylabel('Correlation between sig and noise corrs');
% xlim([10 100]);

% % %PRINT FIGURE
% fig_width = 3.5; rel_height = 0.8;
% figufy(f1);
% fname = [fig_dir 'Signoisecorr_cntwin.pdf'];
% exportfig(f1,fname,'width',fig_width,'height',rel_height*fig_width,'fontmode','scaled','fontsize',1);
% close(f1);
% 

%%
cur_SUs = find(avg_rates >= min_rate & mod_xvLLimps > min_xvLLimp & n_rpt_trials >= min_rpt_trials);

maxtlag = 10;

all_sig_corrs = [];
all_psth_corrs = [];
all_emp_corrs = [];
all_pair_IDS = [];
all_pair_avgrates = [];
all_pair_ntrials = [];
all_expt_IDS = [];
for ss = 1:length(cur_SUs)
    cur_sig_covs = all_SU_data(cur_SUs(ss)).spline_xcov;
%     cur_psth_covs = all_SU_data(cur_SUs(ss)).psth_xcov;
    cur_psth_covs = all_SU_data(cur_SUs(ss)).rand_xcov;
    cur_emp_covs = all_SU_data(cur_SUs(ss)).emp_xcov;
        
    psth_norm = all_SU_data(cur_SUs(ss)).psth_varnorm;
        
    cur_psth_corrs = bsxfun(@rdivide,cur_psth_covs,psth_norm);
    cur_sig_corrs = bsxfun(@rdivide,cur_sig_covs,psth_norm);
    cur_emp_corrs = bsxfun(@rdivide,cur_emp_covs,psth_norm);
    
    cur_pair_IDS = [repmat(all_SU_data(cur_SUs(ss)).unit_num,size(cur_psth_corrs,1),1) (1:size(cur_psth_corrs,1))'];
    
    all_sig_corrs = cat(1,all_sig_corrs,cur_sig_corrs);
    all_psth_corrs = cat(1,all_psth_corrs,cur_psth_corrs);
    all_emp_corrs = cat(1,all_emp_corrs,cur_emp_corrs);
    
    all_pair_IDS = cat(1,all_pair_IDS,cur_pair_IDS);
    all_expt_IDS = cat(1,all_expt_IDS,repmat(all_SU_data(cur_SUs(ss)).expt_num,size(cur_psth_corrs,1),1));
    all_pair_ntrials = cat(1,all_pair_ntrials,all_SU_data(cur_SUs(ss)).otherunit_nutrials);
    all_pair_avgrates = cat(1,all_pair_avgrates,all_SU_data(cur_SUs(ss)).otherunit_avgrates);
end

pair_min_rate = 1;

uset = find(all_pair_IDS(:,1) > all_pair_IDS(:,2) & all_pair_ntrials >= min_rpt_trials & ...
    all_pair_avgrates/dt >= pair_min_rate);

% uset = find(all_pair_IDS(:,1) > all_pair_IDS(:,2) & min(all_pair_IDS,[],2) > 24 & ...
%     all_pair_ntrials >= min_rpt_trials & all_pair_avgrates/dt >= pair_min_rate);

uu = uset(~isnan(all_psth_corrs(uset,1)) & ~isnan(all_sig_corrs(uset,1)));

