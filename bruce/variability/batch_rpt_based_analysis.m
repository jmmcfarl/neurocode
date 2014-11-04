close all
clear all
clc

fig_dir = '/Users/james/Analysis/bruce/variability/figures/';
base_sname = 'rpt_variability_analysis';
base_tname = 'model_variability_analysis';

all_SU_data = [];
n_probes = 24;

%% LOAD JBE
Expt_list = {'M266','M270','M275','M277','M281','M287','M289','M294','M296','M297'};
ori_list = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 160 nan; 40 nan; 45 nan; 0 90];
rmfield_list = {};

for ee = 1:length(Expt_list)
% for ee = 1:3
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
min_rpt_trials = 50;

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
msize = 10;
f1 = figure();
plot(mod_alphas,spline_alpha,'.','markersize',msize);
line([0 1],[0 1],'color','k');
%%
cur_SUs = find(avg_rates >= min_rate & mod_xvLLimps > min_xvLLimp & n_rpt_trials >= min_rpt_trials);

maxtlag = 10;

all_psth_noise_corrs = [];
all_spline_noise_corrs = [];
all_sig_corrs = [];
all_psth_corrs = [];
all_pair_IDS = [];
all_expt_IDS = [];
for ss = 1:length(cur_SUs)
    cur_sig_covs = squeeze(all_SU_data(cur_SUs(ss)).spline_xcov(:,maxtlag+1));
    cur_psth_covs = squeeze(all_SU_data(cur_SUs(ss)).rand_xcov(:,maxtlag+1));
    cur_emp_covs = squeeze(all_SU_data(cur_SUs(ss)).emp_xcov(:,maxtlag+1));
    
    psth_noise_covs = cur_emp_covs - cur_psth_covs;
    spline_noise_covs = cur_emp_covs - cur_sig_covs;
    
%     cur_norms = all_SU_data(cur_SUs(ss)).varnorm_mat(:,maxtlag+1);
      cur_norms = all_SU_data(cur_SUs(ss)).noise_varnorm_mat(:,maxtlag+1);
  
    cur_pair_IDS = [repmat(all_SU_data(cur_SUs(ss)).unit_num,length(cur_norms),1) (1:length(cur_norms))'];
    
    all_psth_noise_corrs = cat(1,all_psth_noise_corrs,psth_noise_covs./cur_norms);
    all_spline_noise_corrs = cat(1,all_spline_noise_corrs,spline_noise_covs./cur_norms);
    all_sig_corrs = cat(1,all_sig_corrs,cur_sig_covs./cur_norms);
    all_psth_corrs = cat(1,all_psth_corrs,cur_psth_covs./cur_norms);
    all_pair_IDS = cat(1,all_pair_IDS,cur_pair_IDS);
    all_expt_IDS = cat(1,all_expt_IDS,repmat(all_SU_data(cur_SUs(ss)).expt_num,length(cur_norms),1));
end
% bad = find(all_sig_corrs == 1);
% all_sig_corrs(bad) = nan;
% all_noise_corrs(bad) = nan;
% 

close all
% uset = find(all_pair_IDS(:,1) ~= all_pair_IDS(:,2));
uset = find(all_pair_IDS(:,1) ~= all_pair_IDS(:,2) & min(all_pair_IDS,[],2) > 24);

% poss_expts = [289];
poss_expts = [266 270 275 277 281 287 289 294 296 297];
uset(~ismember(all_expt_IDS(uset),poss_expts)) = [];

xx = linspace(-0.2,1,100);
uu = uset(~isnan(all_psth_noise_corrs(uset)));

f1 = figure(); hold on
plot(all_sig_corrs(uset),all_psth_noise_corrs(uset),'k.'); 
r = robustfit(all_sig_corrs(uset),all_psth_noise_corrs(uset));
plot(xx,r(1)+r(2)*xx,'r')

plot(all_sig_corrs(uset),all_spline_noise_corrs(uset),'b.'); 
r = robustfit(all_sig_corrs(uset),all_spline_noise_corrs(uset));
plot(xx,r(1)+r(2)*xx,'g')

xlim([-0.1 0.5]);
ylim([-0.1 0.5]);


f2 = figure(); hold on
plot(all_psth_corrs(uset),all_psth_noise_corrs(uset),'k.'); 
r = robustfit(all_psth_corrs(uset),all_psth_noise_corrs(uset));
plot(xx,r(1)+r(2)*xx,'r')

plot(all_sig_corrs(uset),all_spline_noise_corrs(uset),'b.'); 
r = robustfit(all_sig_corrs(uset),all_spline_noise_corrs(uset));
plot(xx,r(1)+r(2)*xx,'g')
xlim([-0.1 0.5]);
ylim([-0.1 0.5]);

