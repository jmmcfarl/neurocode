%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
% cd ~/Data/bruce/7_15_12/G029/
cd ~/Data/bruce/7_15_12/G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt3_fixbased_data.mat
resh_all_stims = resh_all_stims/std(resh_all_stims(:));

fix_win_dur = 0.15;

Pix2Deg = 0.018837;
NT = length(full_stim_ids);
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

Fs = 3e4;
dsf = 40;Fsd = Fs/dsf;
scales = logspace(log10(1.5),log10(30),40);
scales = [scales 35 40 50 60 70 80 90];
scales = scales*100/dsf;
% scales = logspace(log10(2.5),log10(40),25);
% scales = [5:10 12 14 16 20 25 30 35 40 45 50 55 60];
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);
use_lfps = [1:96];

%%
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;
cur_dists = squareform(pdist([X_pos' Y_pos']));
cur_dists(logical(eye(96))) = Inf;
[min_dist,min_loc] = min(cur_dists);
cur_dists(logical(eye(96))) = -Inf;
[max_dist,max_loc] = max(cur_dists);

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

min_trial_dur = 0.4;
used_trials = find(trial_durs >= min_trial_dur);
trial_start_inds = trial_start_inds(used_trials);
trial_stop_inds = trial_stop_inds(used_trials);
trial_imnum = trial_imnum(used_trials);
trial_stimnum = trial_stimnum(used_trials);
trial_start_times = trial_start_times(used_trials);
trial_stop_times = trial_stop_times(used_trials);
trial_stop_winds = trial_start_inds + round(min_trial_dur/dt);
trial_stop_wtimes = full_t(trial_stop_winds);

n_trials = length(trial_start_inds);
fix_expt_num = nan(n_trials,1);
full_intrial = zeros(size(full_t));
full_islate = zeros(size(full_t));
full_t_sincefix = nan(size(full_t));
for i = 1:n_trials
    cur_inds = trial_start_inds(i):trial_stop_winds(i);
    fix_expt_num(i) = unique(full_expt_vec(cur_inds));
    cur_inds = (trial_start_inds(i)+round(0.15/dt)):trial_stop_inds(i);
    full_islate(cur_inds) = 1;
    cur_inds = trial_start_inds(i):trial_stop_inds(i);
    full_intrial(cur_inds) = 1;
    full_t_sincefix(cur_inds) = full_t(cur_inds) - full_t(cur_inds(1));
end
%%
% Expt_nu = [13 14 15 16 25 28 29];
Expt_nu = [7:12]; %expt 3 34
n_expts = length(Expt_nu);
fixavg_lfp_amps = nan(length(trial_start_inds),length(wfreqs),length(use_lfps));

for cc =1:length(use_lfps)
    spike_bins{cc} = [];
    spike_phases{cc} = [];
    spike_phases_d{cc} = [];
    spike_phases_s{cc} = [];
    spike_bins{cc} = [];
end
% pool_ampgrams = [];
pool_phasegrams = [];
pool_intrial = [];
pool_tsince = [];
for ee = 1:n_expts
    fprintf('Expt %d of %d\n',ee,n_expts);
    clear ampgram phasegram 
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d\n',ll);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = decimate(V,dsf);
        V = V(:);
        cur_cwt = cwt(V,scales,'cmor1-1');
%         ampgram(:,:,ll) = abs(cur_cwt)';
        phasegram(:,:,ll) = angle(cur_cwt)';
    end
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
    
    fix_inds = find(full_expt_vec == Expt_nu(ee));
    interp_intrial = round(interp1(full_t(fix_inds),full_intrial(fix_inds),t_ax));
    interp_tsince = interp1(full_t(fix_inds),full_t_sincefix(fix_inds),t_ax);
    
    offset_cnt = length(pool_tsince);
    load(sprintf('Expt%dClusterTimes.mat',Expt_nu(ee)));
    for cc = 1:length(use_lfps)
       cur_spike_bins = round(interp1(t_ax,1:length(t_ax),Clusters{use_lfps(cc)}.times)); 
       spike_bins{cc} = [spike_bins{cc} (cur_spike_bins + offset_cnt)];
       spike_phases{cc} = [spike_phases{cc}; phasegram(cur_spike_bins,:,min_loc(cc))];
       spike_phases_d{cc} = [spike_phases_d{cc}; phasegram(cur_spike_bins,:,max_loc(cc))];
       spike_phases_s{cc} = [spike_phases_s{cc}; phasegram(cur_spike_bins,:,cc)];
    end
    
    pool_intrial = [pool_intrial; interp_intrial'];
    pool_tsince = [pool_tsince; interp_tsince'];
% %     pool_ampgrams = [pool_ampgrams; ampgram];
%     pool_phasegrams = [pool_phasegrams; phasegram];    
end
%%
phase_lock = nan(length(use_lfps),length(wfreqs));
phase_lock_early = nan(length(use_lfps),length(wfreqs));
phase_lock_late = nan(length(use_lfps),length(wfreqs));
phase_lockd = nan(length(use_lfps),length(wfreqs));
phase_lock_earlyd = nan(length(use_lfps),length(wfreqs));
phase_lock_lated = nan(length(use_lfps),length(wfreqs));
phase_locks = nan(length(use_lfps),length(wfreqs));
phase_lock_earlys = nan(length(use_lfps),length(wfreqs));
phase_lock_lates = nan(length(use_lfps),length(wfreqs));
% phase_lock_set2 = nan(length(use_lfps),length(wfreqs));
% phase_lock_set3 = nan(length(use_lfps),length(wfreqs));
% phase_lock_set4 = nan(length(use_lfps),length(wfreqs));

pref_phase = nan(length(use_lfps),length(wfreqs));
pref_phase_early = nan(length(use_lfps),length(wfreqs));
pref_phase_late = nan(length(use_lfps),length(wfreqs));

usable_set = find(~isnan(pool_tsince));
early_set = find(pool_tsince < 0.15);
% late_set = find(pool_tsince > 0.15);
late_set = find(pool_tsince > 0.2);
set2 = find(pool_tsince > 0.1 & pool_tsince < 0.2);
set3 = find(pool_tsince > 0.2 & pool_tsince < 0.3);
set4 = find(pool_tsince > 0.3 & pool_tsince < 0.4);

for ww = 1:length(wfreqs)
    for cc = 1:length(use_lfps)
%         cur_spk_phases = squeeze(pool_phasegrams(spike_bins{cc},ww,cc));
%         phase_lock(cc,ww) = circ_kappa(cur_spk_phases);
        uset = find(ismember(spike_bins{cc},usable_set));
        phase_lock(cc,ww) = circ_kappa(squeeze(spike_phases{cc}(uset,ww)));
        pref_phase(cc,ww) = circ_mean(squeeze(spike_phases{cc}(uset,ww)));
        uset = find(ismember(spike_bins{cc},early_set));
        phase_lock_early(cc,ww) = circ_kappa(squeeze(spike_phases{cc}(uset,ww)));
        pref_phase_early(cc,ww) = circ_mean(squeeze(spike_phases{cc}(uset,ww)));
        uset = find(ismember(spike_bins{cc},late_set));
        phase_lock_late(cc,ww) = circ_kappa(squeeze(spike_phases{cc}(uset,ww)));
        pref_phase_late(cc,ww) = circ_mean(squeeze(spike_phases{cc}(uset,ww)));

                uset = find(ismember(spike_bins{cc},usable_set));
        phase_lockd(cc,ww) = circ_kappa(squeeze(spike_phases_d{cc}(uset,ww)));
        uset = find(ismember(spike_bins{cc},early_set));
        phase_lock_earlyd(cc,ww) = circ_kappa(squeeze(spike_phases_d{cc}(uset,ww)));
        uset = find(ismember(spike_bins{cc},late_set));
        phase_lock_lated(cc,ww) = circ_kappa(squeeze(spike_phases_d{cc}(uset,ww)));

                        uset = find(ismember(spike_bins{cc},usable_set));
        phase_locks(cc,ww) = circ_kappa(squeeze(spike_phases_s{cc}(uset,ww)));
        uset = find(ismember(spike_bins{cc},early_set));
        phase_lock_earlys(cc,ww) = circ_kappa(squeeze(spike_phases_s{cc}(uset,ww)));
        uset = find(ismember(spike_bins{cc},late_set));
        phase_lock_lates(cc,ww) = circ_kappa(squeeze(spike_phases_s{cc}(uset,ww)));
%         uset = find(ismember(spike_bins{cc},set2));
%         phase_lock_set2(cc,ww) = circ_kappa(squeeze(spike_phases{cc}(uset,ww)));
%         uset = find(ismember(spike_bins{cc},set3));
%         phase_lock_set3(cc,ww) = circ_kappa(squeeze(spike_phases{cc}(uset,ww)));
%         uset = find(ismember(spike_bins{cc},set4));
%         phase_lock_set4(cc,ww) = circ_kappa(squeeze(spike_phases{cc}(uset,ww)));
    end
end

%%
save phase_lock_analysis_v3 wfreqs phase_lock* pref_phase*
%%
close all
for i = 1:96
    i
    plot(wfreqs,phase_lock_late(i,:))
    set(gca,'xscale','log')
    xlim([5 100])
    pause
    clf
end