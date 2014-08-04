cclear all
% close all
addpath(genpath('~/James_scripts'));

cd ~/Data/bruce/2_27_12/
load Blocks.mat
accept_window = [-5 5;-5 5];
sac_eyespeed = 10;
thresh_eyespeed = 2.5;
sac_buffer = 0;
min_fix_dur = 0.2;
blink_thresh = 5;
min_blink_dur = 0.05;

Pix2Deg = 0.018837;

% down-sampling fraction for image
dsfrac = 2;
Fsd = 1/Pix2Deg/dsfrac;

Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;

RF_patch = [2.75 6.25; -3.5 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

cd ~/Data/bruce/2_27_12/stimrecon
load fixation_data_v3

cd ~/Data/bruce/2_27_12/stimrecon
load sing_eye_pgabortrack_d2
load sing_eye_pgabortrack_fin_est_d2
cd ~/Data/bruce/2_27_12/ExptA

SDIM = length(xpatch_inds);
cd ~/Data/bruce/2_27_12/stimrecon
load fixation_image_patches_corrected_raw

cd ~/James_scripts/bruce/modelfits
load pref_oris
cd ~/James_scripts/bruce/modelfits/
% load fixation_gabor_models
load gabor_tempmodfits_en_lin


NT = length(used_fixs);
X_resh = reshape(X,NT,SDIM^2);
dt = .005;

cellids = [1 2 3 4 5 6 7 8 10];
muaids = [1 2 4 5 6 11 12 13 14];
su_probes = [Blocks{1}.suprobes(cellids)];
mu_probes = [Blocks{1}.muprobes(muaids)];
all_probes = [su_probes mu_probes];
spk_cnts = [spk_cnts(:,cellids) mua_cnts(:,muaids)];

n_used_cells = size(spk_cnts,2);

%%
Fs = 1000;
dsf = 3;
Fsd = Fs/dsf;
niqf = Fs/2;

scales = logspace(log10(5),log10(60),30);
scales = [scales 70 80 90];
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

new_dt = .0025;

%%
new_model_time_axis = [];
new_model_blockids = [];
new_model_fixids = [];
spikes_binned = [];
time_since_fix = [];
rel_spk_times = cell(n_used_cells,1);
spk_fix_inds = cell(n_used_cells,1);
for blockid = 1:4;
    fprintf('Block %d of %d\n',blockid,4);
    cur_set = find(blockids(used_fixs) == blockid);
    
    for i = 1:length(cur_set)
        start_T = all_fix_start_times(used_fixs(cur_set(i)));
        end_T = all_fix_stop_times(used_fixs(cur_set(i)));
        cur_tedges = start_T:new_dt:end_T;
        cur_tcents = 0.5*cur_tedges(1:end-1)+0.5*cur_tedges(2:end);
        new_model_time_axis = [new_model_time_axis cur_tcents];
        new_model_blockids = [new_model_blockids blockid*ones(1,length(cur_tcents))];
        new_model_fixids = [new_model_fixids cur_set(i)*ones(1,length(cur_tcents))];
        temp_binned = zeros(10,length(cur_tcents));
        temp_modouts = zeros(10,length(cur_tcents));
        for c = 1:n_used_cells
            if c <= 9
                cur_spk_times = find(Blocks{blockid}.spktimes{cellids(c)} > start_T & ...
                    Blocks{blockid}.spktimes{cellids(c)} < end_T);
                temp = hist(Blocks{blockid}.spktimes{cellids(c)}(cur_spk_times),cur_tcents);
            else
                cur_spk_times = find(Blocks{blockid}.mutimes{muaids(c-9)} > start_T & ...
                    Blocks{blockid}.mutimes{muaids(c-9)} < end_T);
                temp = hist(Blocks{blockid}.mutimes{muaids(c-9)}(cur_spk_times),cur_tcents);
            end
            temp_binned(c,:) = temp;
        end
        spikes_binned = [spikes_binned; temp_binned'];
        
        time_since_fix = [time_since_fix cur_tcents-start_T];
        
    end
end

%%
new_phasegram = [];
new_ampgram = [];
% all_t = [];
for blockid = 1:4;
    fprintf('Block %d of %d\n',blockid,4);
    
    %%
    cd ~/Data/bruce/2_27_12
    load(sprintf('lemM232A.5%d.lfp.mat',blockid));
    %get start times of each LFP trial
    if blockid == 4
        LFP.Trials = LFP.Trials(1:5);
    end
    n_lfp_trials = length(LFP.Trials);
    lfp_trial_start = nan(n_lfp_trials,1);
    lfp_trial_stop = nan(n_lfp_trials,1);
    for i = 1:n_lfp_trials
        lfp_trial_start(i) = LFP.Trials(i).Start/1e4; %originally in tenths of ms
        lfp_trial_stop(i) = LFP.Trials(i).End/1e4;
        lfp_dur(i) = size(LFP.Trials(i).LFP,1)/1000;
    end
    
    lfp_time = [];
    lfp_samps = [];
    for i = 1:n_lfp_trials
        cur_len(i) = size(LFP.Trials(i).LFP,1);
        if i < n_lfp_trials
            next_start = lfp_trial_start(i+1);
            start_len = length(lfp_trial_start(i):1/Fs:next_start);
        else
            next_start = Inf;
            start_len = Inf;
        end
        cur_end(i) = min(cur_len(i),start_len);
        cur_t = lfp_trial_start(i):1/Fs:(lfp_trial_start(i)+cur_end(i)/Fs);
        cur_t(cur_end(i)+1:end) = [];
        lfp_time = [lfp_time cur_t];
        lfp_samps = [lfp_samps; LFP.Trials(i).LFP(1:cur_end(i),:)];
    end
    
    lfp_timed = downsample(lfp_time,dsf);
    
    cur_phasegram = nan(length(lfp_timed),length(wfreqs),24);
    cur_ampgram = nan(length(lfp_timed),length(wfreqs),24);
    for ll = 1:24
        fprintf('Processing channel %d of %d\n',ll,24);
        lfp_sampsd = decimate(lfp_samps(:,ll),dsf);
        temp = cwt(lfp_sampsd,scales,'cmor1-1');
        cur_phasegram(:,:,ll) = angle(temp)';
        cur_ampgram(:,:,ll) = abs(temp)';
    end
    
    cur_set = find(new_model_blockids==blockid);
    %     interp_phasegram = interp1(lfp_timed,cur_phasegram,new_model_time_axis(cur_set));
    unwr_phasegram = unwrap(cur_phasegram);
    interp_phasegram = interp1(lfp_timed,unwr_phasegram,new_model_time_axis(cur_set));
    interp_phasegram = mod(interp_phasegram+pi,2*pi)-pi;

    interp_ampgram = interp1(lfp_timed,cur_ampgram,new_model_time_axis(cur_set));
    
    new_phasegram = cat(1,new_phasegram,interp_phasegram);
    new_ampgram = cat(1,new_ampgram,interp_ampgram);
    %     all_t = [all_t; lfp_timed(:)];
end


%%
% save('all_lfp_phase_data.mat','-v7.3','new_phasegram','new_ampgram','new_*')

%%

cur_used_fixs = unique(new_model_fixids);
n_fixs = length(cur_used_fixs);
trial_start_inds = 1 + find(diff(new_model_fixids) ~= 0);
trial_stop_inds = find(diff(new_model_fixids) ~= 0);
trial_start_inds = [1 trial_start_inds];
trial_stop_inds = [trial_stop_inds length(new_model_fixids)];
%%
cur_dt = 0.006;
flen_t = 0.5;
tent_centers = [0:cur_dt:0.15];
cur_sp = dt;
while max(tent_centers) < flen_t
    tent_centers = [tent_centers (tent_centers(end)+cur_sp)];
    cur_sp = cur_sp + dt/5;
end

tent_centers = round(tent_centers/new_dt);
tbmat = construct_tent_bases(tent_centers,1);
[ntents,tblen] = size(tbmat);

tbmat = [zeros(ntents,tblen-1) tbmat];

is_bad = max(max(isnan(new_phasegram),[],3),[],2);
uset = find(time_since_fix <= max(tent_centers) & is_bad' == 0);

trial_inds = zeros(size(new_model_time_axis));
trial_inds(trial_start_inds) = 1;
trial_Tmat = zeros(length(new_model_time_axis),ntents);
for i = 1:ntents
    trial_Tmat(:,i) = conv(trial_inds,tbmat(i,:),'same');
end

beg_dur = round(0.15/new_dt);
late_indicator = zeros(size(new_model_time_axis));
for i = 1:length(trial_start_inds)
    cur_inds = (beg_dur+trial_start_inds(i)):trial_stop_inds(i);
    late_indicator(cur_inds) = 1;
end

NT = length(new_model_time_axis);
%%
cd /home/james/Data/bruce/2_27_12/M232/
load ./bar_phase_models_ftime_v2
bar_se_mods = stim_ind_phase_pfilt;
bar_full_mods = [stim_ind_sinphase_cfilt(all_probes,:) stim_ind_sinphase_sfilt(all_probes,:)];

cd /home/james/Data/bruce/2_27_12
load ./free_viewing_phase_models
new_phase_set = [reshape(cos(new_phasegram),length(trial_inds),length(wfreqs)*24) reshape(sin(new_phasegram),length(trial_inds),length(wfreqs)*24)];
phase_elec_set = [repmat(1:24,1,length(wfreqs)) repmat(1:24,1,length(wfreqs))];
phase_elec_set = phase_elec_set(:);

%%
for cur_cell = 1:length(all_probes);
fprintf('Cell %d of %d\n',cur_cell,length(all_probes));

    Robs = spikes_binned(:,cur_cell);
    if all_probes(cur_cell) < 24
        nearest_probe = all_probes(cur_cell)+1;
    else
        nearest_probe = all_probes(cur_cell)-1;
    end
    % nearest_probe = cc;
cur_alpha_phase = squeeze(new_phasegram(:,:,nearest_probe));
Pmat = nan(NT,length(wfreqs)*nbins);
for ww = 1:length(wfreqs)
    cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
    cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
    cur_tb(:,end) = [];
    cur_set = ((ww-1)*nbins+1):ww*nbins;
    Pmat(:,cur_set) = cur_tb;
end

se_phasemod_out(cur_cell,:) = Pmat*stim_ind_phase_pfilt(cur_cell,:)';
se_phasemod_out(cur_cell,:) = se_phasemod_out(cur_cell,:) - nanmean(se_phasemod_out(cur_cell,:));
bse_phasemod_out(cur_cell,:) = Pmat*bar_se_mods(cur_cell,:)';
bse_phasemod_out(cur_cell,:) = bse_phasemod_out(cur_cell,:) - nanmean(bse_phasemod_out(cur_cell,:));

use_elecs = 1:24;
use_elecs(cur_cell) = [];
use_set = find(ismember(phase_elec_set,use_elecs));
full_phase_filt = [stim_ind_sinphase_cfilt(cur_cell,:) stim_ind_sinphase_sfilt(cur_cell,:)];
full_phasemod_out(cur_cell,:) = new_phase_set(:,use_set)*full_phase_filt';
full_phasemod_out(cur_cell,:) = full_phasemod_out(cur_cell,:) - nanmean(full_phasemod_out(cur_cell,:));
bfull_phasemod_out(cur_cell,:) = new_phase_set(:,use_set)*bar_full_mods(cur_cell,:)';
bfull_phasemod_out(cur_cell,:) = bfull_phasemod_out(cur_cell,:) - nanmean(bfull_phasemod_out(cur_cell,:));

% glm_kern = get_k_mat(glm_fit(cur_cell));
% stim_out = stim_Xmat*glm_kern;
% stim_out_interp(cur_cell,:) = stim_out(full_old_inds(tr_inds_new));
% 
% rate_sm = 1;
% sm_rate(cur_cell,:) = jmm_smooth_1d_cor(Robs,rate_sm);
end

%%
[~,ord] = sort(all_probes);
for i = 1:length(trial_start_inds)
    cur_set = trial_start_inds(i):trial_stop_inds(i);
    subplot(3,1,1)
    imagesc(new_model_time_axis(cur_set)-new_model_time_axis(cur_set(1)),1:length(all_probes),se_phasemod_out(ord,cur_set));
    caxis([-2 1.5])
        xlim([0 0.5])
    subplot(3,1,2)
    imagesc(new_model_time_axis(cur_set)-new_model_time_axis(cur_set(1)),1:length(all_probes),full_phasemod_out(ord,cur_set));
    caxis([-2 1.5])
        xlim([0 0.5])

    subplot(3,1,3)
    imagesc(new_model_time_axis(cur_set)-new_model_time_axis(cur_set(1)),1:length(all_probes),bfull_phasemod_out(ord,cur_set));
    caxis([-2 1.5])
         xlim([0 0.5])
   pause
    clf
end

%%
maxlag = round(0.5/new_dt);
trg_avg_full = zeros(n_used_cells,maxlag);
trg_avg_bfull = zeros(n_used_cells,maxlag);
trg_avg_se = zeros(n_used_cells,maxlag);
trg_avg_bse = zeros(n_used_cells,maxlag);
cnt = zeros(maxlag,1);
for i = 1:length(trial_start_inds)
    cur_set = trial_start_inds(i):trial_stop_inds(i);
    cur_set(maxlag+1:end) = [];
    cl = length(cur_set);
    
    trg_avg_se(:,1:cl) = trg_avg_se(:,1:cl) + se_phasemod_out(:,cur_set);
    trg_avg_bse(:,1:cl) = trg_avg_bse(:,1:cl) + bse_phasemod_out(:,cur_set);
    trg_avg_full(:,1:cl) = trg_avg_full(:,1:cl) + full_phasemod_out(:,cur_set);
    trg_avg_bfull(:,1:cl) = trg_avg_bfull(:,1:cl) + bfull_phasemod_out(:,cur_set);
    
    cnt(1:cl) = cnt(1:cl) + 1;
end
trg_avg_se = bsxfun(@rdivide,trg_avg_se,cnt');
trg_avg_bse = bsxfun(@rdivide,trg_avg_bse,cnt');
trg_avg_full = bsxfun(@rdivide,trg_avg_full,cnt');
trg_avg_bfull = bsxfun(@rdivide,trg_avg_bfull,cnt');