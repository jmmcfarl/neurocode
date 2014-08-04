%% Load Data
clear all;
cd ~/Data/bruce/7_15_12/G034/

obj_info_dir = '~/James_scripts/data_processing/Images/object_im_info/';
nat_set = 1:685;
white_set = 686:913;
sim_set = 914:1141;
obj_set = 1142:1369;
noise_set = [white_set sim_set];

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt3_fixbased_data.mat

Pix2Deg = 0.018837;
Nyp = 1024;
Nxp = 1280;
use_win = [-7 7;-7 7];

min_fix_dur = 0.15;
use_lfps = [1:2:96];

dts = median(diff(full_t));
spk_smwin = round(.006/dts);

forwardlag = round(0.6/dts);
backlag = round(0.6/dts);
lags = -backlag:forwardlag;

x_pix = xax(xpatch_inds);
y_pix = yax(ypatch_inds);
[X_pix,Y_pix] = meshgrid(x_pix,y_pix);

% load ./expt2_lfp_amp_models_full_v2
% % fft_pred_out = fft_pred_out(1:2:end,:,:);
% use_gamma_freq = 20;
% gamma_mod_out = squeeze(fft_pred_out(:,use_gamma_freq,:));
% gamma_mod_out = zscore(gamma_mod_out')';
% gamma_prctile_mat = prctile(gamma_mod_out',[25 75]);

%%
n_trials = length(trial_start_inds);
trial_stop_sacinds = trial_stop_inds;
for i = 1:n_trials
    cur_inds = trial_start_inds(i):trial_stop_inds(i);
    cur_sac_stops = find(full_insac(cur_inds) == 1,1,'first');
    if ~isempty(cur_sac_stops)
        cur_sac_stops = cur_inds(cur_sac_stops);
        trial_stop_sacinds(i) = cur_sac_stops;
    end
end
trial_stop_inds = trial_stop_sacinds;
trial_start_times = full_t(trial_start_inds);
trial_stop_times = full_t(trial_stop_inds);
trial_durs = trial_stop_times - trial_start_times;
expt_start_inds = 1 + find(diff(full_trial_vec) > 0);
expt_start_inds = expt_start_inds(ismember(expt_start_inds,trial_start_inds));

% bad_inds = find(ismember(trial_start_inds,expt_start_inds));
% trial_start_inds(bad_inds) = [];
% trial_stop_inds(bad_inds) = [];
n_trials = length(trial_start_inds);
in_trial = zeros(size(full_t));
for i =1 :n_trials
        cur_inds = trial_start_inds(i):trial_stop_inds(i);
in_trial(cur_inds) = 1;
end
%%
full_avg_rate = mean(full_binned_spks(in_trial==1,:));
full_binned_nspks = bsxfun(@rdivide,full_binned_spks,full_avg_rate);
full_smbinned_spks = nan(size(full_binned_spks));
for c = 1:96
    full_smbinned_spks(:,c) = jmm_smooth_1d_cor(full_binned_nspks(:,c),spk_smwin);
end
%%

use_trials = find(trial_start_inds > backlag & trial_stop_inds < length(full_t) - forwardlag);
% use_trials = find(trial_start_inds > backlag & trial_stop_inds < length(full_t) - forwardlag & ~ismember(trial_start_inds,expt_start_inds));
binned_spk_mat = zeros(96,length(lags));
nbinned_spk_mat = zeros(96,length(lags));
smbinned_spk_mat = zeros(96,length(lags));
cnt = zeros(1,length(lags));
for i = 1:length(use_trials)
    cur_inds = (trial_start_inds(use_trials(i))-backlag):(trial_start_inds(use_trials(i))+forwardlag);
    cur_inds(cur_inds > trial_stop_inds(use_trials(i))) = [];
    cur_inds(length(lags)+1:end) = [];
    cl = length(cur_inds);
    binned_spk_mat(:,1:cl) = binned_spk_mat(:,1:cl) + full_binned_spks(cur_inds,:)';
    nbinned_spk_mat(:,1:cl) = nbinned_spk_mat(:,1:cl) + full_binned_nspks(cur_inds,:)';
    smbinned_spk_mat(:,1:cl) = smbinned_spk_mat(:,1:cl) + full_smbinned_spks(cur_inds,:)';
cnt(1:cl) = cnt(1:cl) + 1;
end
binned_spk_mat = bsxfun(@rdivide,binned_spk_mat,cnt);
nbinned_spk_mat = bsxfun(@rdivide,nbinned_spk_mat,cnt);
smbinned_spk_mat = bsxfun(@rdivide,smbinned_spk_mat,cnt);
clear full_binned_spks
%% Load Data
% clear all
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12/G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt2_compiled_windata_d1p5

min_fix_dur = 0.15;
%%
n_fixs = length(full_fix_wends);
dte = median(diff(full_t));
spk_smwin = round(.006/dte);

forwardlage = round(0.6/dte);
backlage = round(0.4/dte);
lagse = -backlage:forwardlage;

sac_inds = find(~isnan(fix_sac_amps));

use_inds = zeros(size(full_binned_spks,1),1);
for i = 1:length(use_trials)
    cur_inds = (full_fix_starts(i)):(full_fix_ends(i));
    use_inds(cur_inds) = 1;
end
%%
full_avg_rate = mean(full_binned_spks(use_inds==1,:));
full_binned_nspks = bsxfun(@rdivide,full_binned_spks,full_avg_rate);
full_smbinned_spks = nan(size(full_binned_nspks));
for c = 1:96
    full_smbinned_spks(:,c) = jmm_smooth_1d_cor(full_binned_nspks(:,c),spk_smwin);
end
%%
use_trials = sac_inds(full_fix_starts(sac_inds) > backlage & full_fix_ends(sac_inds) < length(full_t) - forwardlage);
binned_spk_mate = zeros(96,length(lagse));
nbinned_spk_mate = zeros(96,length(lagse));
smbinned_spk_mate = zeros(96,length(lagse));
cnt = zeros(1,length(lagse));
for i = 1:length(use_trials)
    cur_inds = (full_fix_starts(use_trials(i))-backlage):(full_fix_starts(use_trials(i))+forwardlage);
%     cur_inds(cur_inds > full_fix_ends(use_trials(i))) = [];
    cur_inds(length(lagse)+1:end) = [];
    cl = length(cur_inds);
    binned_spk_mate(:,1:cl) = binned_spk_mate(:,1:cl) + full_binned_spks(cur_inds,:)';
    nbinned_spk_mate(:,1:cl) = nbinned_spk_mate(:,1:cl) + full_binned_nspks(cur_inds,:)';
    smbinned_spk_mate(:,1:cl) = smbinned_spk_mate(:,1:cl) + full_smbinned_spks(cur_inds,:)';
cnt(1:cl) = cnt(1:cl) + 1;
end
binned_spk_mate = bsxfun(@rdivide,binned_spk_mate,cnt);
nbinned_spk_mate = bsxfun(@rdivide,nbinned_spk_mate,cnt);
smbinned_spk_mate = bsxfun(@rdivide,smbinned_spk_mate,cnt);

%%
figure
shadedErrorBar(lagse*dte,mean(smbinned_spk_mate),std(smbinned_spk_mate)/sqrt(96),{'color','r'})
hold on
shadedErrorBar(lags*dts,mean(smbinned_spk_mat),std(smbinned_spk_mat)/sqrt(96))
