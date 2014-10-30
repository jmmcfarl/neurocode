close all
clear all
clc

fig_dir = '/Users/james/Analysis/bruce/variability/figures/';
base_sname = 'model_variability_analysis';

all_SU_data = [];

%% LOAD JBE
Expt_list = {'M266','M270','M275','M277','M281','M287','M294','M296','M297'};%NOTE: Excluding M289 because fixation point jumps in and out of RFs, could refine analysis to handle this
n_probes = 24;
ori_list = [80 nan; 60 nan; 135 nan; 70 nan; 140 nan; 90 nan; 40 nan; 45 nan; 0 90];
rmfield_list = {};

% for ee = 1:length(Expt_list)
for ee = 1:8
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

figure;
shadedErrorBar(poss_SDs,nanmean(alpha_funs),nanstd(alpha_funs));
hold on
plot(poss_SDs,alpha_funs,'r');

target_SD = 0.1;
target_alphas = interp1(poss_SDs,alpha_funs',target_SD);

target_SD1 = 0.08;
target_alphas1 = interp1(poss_SDs,alpha_funs',target_SD1);
target_SD2 = 0.13;
target_alphas2 = interp1(poss_SDs,alpha_funs',target_SD2);

figure;
subplot(2,1,1)
plot(RF_ecc(cur_SUs),target_alphas,'o');subplot(2,1,2);
plot(RF_sigma(cur_SUs),target_alphas,'o');

