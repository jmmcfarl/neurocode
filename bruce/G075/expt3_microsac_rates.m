clear all
% close all
cd ~/Data/bruce/G075/
load jbeG075Expts.mat

%%
load ./Expt3_fixbased_data.mat
Pix2Deg = 0.018837;
NT = length(full_stim_ids);
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

load ./expt3_eyedata_NS.mat

%%
sim_sac_blocks = [14 18 24 28 37 40 50] - 6;
% sim_sac_blocks = [18 28 40] - 6; %sim sac
% % sim_sac_blocks = [10 22 32 42] - 6; %sim sac a
% sim_sac_blocks = [14 24 37 50] - 6; %sim sac b

dt = median(diff(full_t));
backlag = round(0.1/dt);
forwardlag = round(0.6/dt);
interp_eyespeed = interp1(all_t,all_eyespeed,full_t);

sac_bins = round(interp1(full_t,1:length(full_t),all_t(all_sac_start_inds)));
sac_bins(isnan(sac_bins)) = [];
sac_binned = zeros(size(full_t));
sac_binned(sac_bins) = 1;
sac_binned = jmm_smooth_1d_cor(sac_binned,1);
for ll = 1:length(sim_sac_blocks)
    cur_set = find(full_expt_vec == sim_sac_blocks(ll));
    trial_durs = trial_stop_inds - trial_start_inds;
    used_start_inds = trial_start_inds(trial_start_inds > backlag & trial_start_inds < length(full_t)-forwardlag & trial_durs >= 75);
    used_start_inds(~ismember(used_start_inds,cur_set)) = [];
    [stim_trg_avg_eyespeed(ll,:),lags] = get_event_trig_avg(interp_eyespeed,used_start_inds,backlag,forwardlag);
    [stim_trg_avg_sac(ll,:),lags] = get_event_trig_avg(sac_binned,used_start_inds,backlag,forwardlag);
end

%%
figure
shadedErrorBar(lags*dt,mean(stim_trg_avg_sac)/dt,std(stim_trg_avg_sac)/dt/sqrt(length(sim_sac_blocks)),{'color','r'});
ylim([0 4]); xlim([0 0.47]);

figure
shadedErrorBar(lags*dt,mean(stim_trg_avg_eyespeed),std(stim_trg_avg_eyespeed)/sqrt(length(sim_sac_blocks)),{'color','r'});
% ylim([0 4]); 
xlim([0 0.47]);

%%
% load ./expt3_eyedata_NS_sensitive.mat
% 
% sac_bins = round(interp1(full_t,1:length(full_t),all_t(all_sac_start_inds)));
% sac_bins(isnan(sac_bins)) = [];
% sac_binned = zeros(size(full_t));
% sac_binned(sac_bins) = 1;
% sac_binned = jmm_smooth_1d_cor(sac_binned,1);
% for ll = 1:length(sim_sac_blocks)
%     cur_set = find(full_expt_vec == sim_sac_blocks(ll));
%     trial_durs = trial_stop_inds - trial_start_inds;
%     used_start_inds = trial_start_inds(trial_start_inds > backlag & trial_start_inds < length(full_t)-forwardlag & trial_durs >= 75);
%     used_start_inds(~ismember(used_start_inds,cur_set)) = [];
%     [stim_trg_avg_eyespeed(ll,:),lags] = get_event_trig_avg(interp_eyespeed,used_start_inds,backlag,forwardlag);
%     [stim_trg_avg_sac(ll,:),lags] = get_event_trig_avg(sac_binned,used_start_inds,backlag,forwardlag);
% end
% hold on
% shadedErrorBar(lags*dt,mean(stim_trg_avg_sac)/dt,std(stim_trg_avg_sac)/dt/sqrt(length(sim_sac_blocks)),{'color','b'});
