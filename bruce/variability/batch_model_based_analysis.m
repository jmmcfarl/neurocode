close all
clear all
clc

% fig_dir = '/Users/james/Analysis/bruce/variability/figures/';
fig_dir = '/home/james/Analysis/bruce/variability/figures/';
base_sname = 'model_variability_analysis';
% base_sname = 'model_variability_analysis_unCor';

all_SU_data = [];

%% LOAD LEM
Expt_list = {'M266','M270','M275','M277','M281','M287','M294','M296','M297'};%NOTE: Excluding M289 because fixation point jumps in and out of RFs, could refine analysis to handle this
n_probes = 24;
ori_list = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 40 nan; 45 nan; 0 90];
rmfield_list = {};

for ee = 1:length(Expt_list)
% for ee = 1:8
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
    cur_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
    
    clear ori_data
    for ii =  1:2
        if ~isnan(ori_list(ee,ii))
            %load trig avg data
            tname = strcat(cur_dir,base_sname,sprintf('_ori%d',ori_list(ee,ii)));
            load(tname);
            
            [EP_data.expt_num] = deal(Expt_num);
            [EP_data.bar_ori] = deal(ori_list(ee,ii));
            [EP_data.animal] = deal('lem');
            
            ori_data(ii,:) = EP_data;
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

%% LOAD JBE
Expt_list = {'G085','G086','G087','G088','G089','G091','G093','G095'};
n_probes = 96;
ori_list = [0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 90; 0 nan];
% ori_list = [0 nan; 0 nan; 0 nan; 0 nan; 0 nan; 0 nan; 0 nan; 0 nan];
rmfield_list = {};

for ee = 1:length(Expt_list)
% for ee = 1
    Expt_name = Expt_list{ee};
    Expt_num = str2num(Expt_name(2:end));
    cur_dir = ['~/Analysis/bruce/' Expt_name '/variability/'];
    
    clear ori_data
    for ii =  1:2
        if ~isnan(ori_list(ee,ii))
            %load trig avg data
            tname = strcat(cur_dir,base_sname,sprintf('_ori%d',ori_list(ee,ii)));
            load(tname);
            
            [EP_data.expt_num] = deal(Expt_num);
            [EP_data.bar_ori] = deal(ori_list(ee,ii));
            [EP_data.animal] = deal('jbe');
            
            ori_data(ii,:) = EP_data;
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

RF_ecc = arrayfun(@(x) x.ModData.tune_props.RF_ecc,all_SU_data);
RF_FSF = arrayfun(@(x) x.ModData.tune_props.RF_FSF,all_SU_data);
RF_gSF = arrayfun(@(x) x.ModData.tune_props.RF_gSF,all_SU_data);
RF_sigma = arrayfun(@(x) x.ModData.tune_props.RF_sigma,all_SU_data);
RF_PRM = arrayfun(@(x) x.ModData.tune_props.PRM,all_SU_data);
%%
cur_SUs = find(avg_rates >= min_rate & mod_xvLLimps > min_xvLLimp);
actual_EP_SD = arrayfun(@(x) x.actual_EP_SD,all_SU_data(cur_SUs));
alpha_funs = cell2mat(arrayfun(@(x) x.alpha_funs, all_SU_data(cur_SUs),'uniformoutput',0));

median_ep_SD = median(actual_EP_SD);
mm = minmax(actual_EP_SD);
min_ep_SD = mm(1); max_ep_SD = mm(2);

figure;hold on
plot(poss_SDs,alpha_funs,'r');
shadedErrorBar(poss_SDs,nanmean(alpha_funs),nanstd(alpha_funs));
yl = ylim();
line([0 0]+median_ep_SD,yl,'color','k')
line([0 0]+min_ep_SD,yl,'color','k');
line([0 0]+max_ep_SD,yl,'color','k');
target_alphas_med = interp1(poss_SDs,alpha_funs',median_ep_SD);
target_alphas_min = interp1(poss_SDs,alpha_funs',min_ep_SD);
target_alphas_max = interp1(poss_SDs,alpha_funs',max_ep_SD);

minmax_alpha_folddiff = target_alphas_min./target_alphas_max;

%%

figure;
subplot(2,1,1)
plot(RF_ecc(cur_SUs),target_alphas_med,'o');
subplot(2,1,2);
plot(RF_sigma(cur_SUs)*2,target_alphas_med,'o');
xlim([0 0.6])


%%
close all

maxtlag = 10;
sig_corr_bin_edges = [-0.99 -0.75 -0.5 -0.25 -0.1 -0.05 0 0.05 0.1 0.25 0.5 0.75 0.99];
sig_corr_bin_cents = 0.5*(sig_corr_bin_edges(1:end-1) + sig_corr_bin_edges(2:end));

all_noisecorr_avgs = zeros(length(sig_corr_bin_cents),maxtlag*2+1,length(poss_SDs));
all_noisecorr_n = zeros(length(sig_corr_bin_cents),1);

all_noise_corrs = [];
all_sig_corrs = [];
all_psth_corrs = [];

for ss = 1:length(cur_SUs)
    cur_sig_corrs = squeeze(all_SU_data(cur_SUs(ss)).sig_corr_mat(:,maxtlag+1,1));
    [n,bins] = histc(cur_sig_corrs,sig_corr_bin_edges);
    for ii = 1:length(sig_corr_bin_cents)
        if n(ii) > 0
        all_noisecorr_avgs(ii,:,:) = all_noisecorr_avgs(ii,:,:) + (sum(all_SU_data(cur_SUs(ss)).noise_corr_mat(bins==ii,:,:),1));
        all_noisecorr_n(ii) = all_noisecorr_n(ii) + n(ii);
        end
    end
    all_sig_corrs = cat(1,all_sig_corrs,squeeze(all_SU_data(cur_SUs(ss)).sig_corr_mat(:,maxtlag+1,:)));
    all_noise_corrs = cat(1,all_noise_corrs,squeeze(all_SU_data(cur_SUs(ss)).noise_corr_mat(:,maxtlag+1,:)));
    all_psth_corrs = cat(1,all_psth_corrs,squeeze(all_SU_data(cur_SUs(ss)).psth_corr_mat(:,maxtlag+1,:)));
end
bad = find(all_sig_corrs == 1);
all_sig_corrs(bad) = nan;
all_noise_corrs(bad) = nan;
all_psth_corrs(bad) = nan;

all_noisecorr_avgs = bsxfun(@rdivide,all_noisecorr_avgs,all_noisecorr_n);

all_sig_corrs_med = interp1(poss_SDs,all_sig_corrs',median_ep_SD);
all_psth_corrs_med = interp1(poss_SDs,all_psth_corrs',median_ep_SD);
all_noise_corrs_med = interp1(poss_SDs,all_noise_corrs',median_ep_SD);

xx = linspace(-0.3,0.7,100);
f = figure(); hold on
plot(all_psth_corrs_med,all_noise_corrs_med,'.');
r = regress(all_noise_corrs_med',[ones(size(all_noise_corrs_med)); all_psth_corrs_med]');
plot(xx,xx*r(2)+r(1),'r');
xlim([-0.3 0.6]);
ylim([-0.3 0.6])

noise_psth_corr_slope = nan(9,1);
f2 = figure; hold on
cmap = jet(9);
for cc = 1:9
    plot(all_psth_corrs(:,cc),all_noise_corrs(:,cc),'.','markersize',5,'color',cmap(cc,:));
    r = regress(all_noise_corrs(:,cc),[ones(size(all_noise_corrs_med)); all_psth_corrs(:,cc)']');
    plot(xx,xx*r(2)+r(1),'color',cmap(cc,:));
    noise_psth_corr_slope(cc) = r(2);
end
xlim([-0.3 0.7]);
ylim([-0.3 0.7])


noise_sig_corr_slope = nan(9,1);
f4 = figure; hold on
for cc = 1:9
    plot(all_sig_corrs(:,cc),all_noise_corrs(:,cc),'.','markersize',5,'color',cmap(cc,:));
    r = regress(all_noise_corrs(:,cc),[ones(size(all_noise_corrs_med)); all_sig_corrs(:,cc)']');
    plot(xx,xx*r(2)+r(1),'color',cmap(cc,:));
    noise_sig_corr_slope(cc) = r(2);
end
xlim([-0.3 0.7]);
ylim([-0.3 0.7])

f3 = figure();hold on
plot(poss_SDs,noise_psth_corr_slope,'k--');
plot(poss_SDs,noise_sig_corr_slope,'r--');
for cc = 1:9
    plot(poss_SDs(cc),noise_psth_corr_slope(cc),'o','color',cmap(cc,:),'linewidth',2,'markersize',8);
    plot(poss_SDs(cc),noise_sig_corr_slope(cc),'o','color',cmap(cc,:),'linewidth',2,'markersize',8);
end

%%
tlags = -maxtlag:maxtlag;
target_SD_ind = 5;
figure;
pcolor(tlags*dt,sig_corr_bin_cents,squeeze(all_noisecorr_avgs(:,:,target_SD_ind)));shading flat
caxis([-0.2 0.2]);

%%
cc = 91;
sd = 5;
sig_xcovs = squeeze(all_SU_data(cc).true_rate_cov(:,:,sd));
