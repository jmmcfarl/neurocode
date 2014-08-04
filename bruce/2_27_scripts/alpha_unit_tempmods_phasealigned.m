clear all
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
% nlags = 4;

Pix2Deg = 0.018837;

% down-sampling fraction for image
dsfrac = 2;
Fsd = 1/Pix2Deg/dsfrac;

Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
% [XX,YY] = meshgrid(xax,yax);

% RF_patch = [3.5 6; -3.5 -1]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [2.5 6.5; -4 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [2.75 6.25; -3.5 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

% patch_inds = find(XX >= RF_patch(1,1) & XX <= RF_patch(1,2) & YY >= RF_patch(2,1) & YY <= RF_patch(2,2));
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

NT = length(used_fixs);
X_resh = reshape(X,NT,SDIM^2);
dt = .005;

cd ~/James_scripts/bruce/modelfits
load pref_oris
cd ~/James_scripts/bruce/modelfits/
% load fixation_gabor_models
load gabor_tempmodfits_en_lin


cellids = [1 2 3 4 5 6 7 8 10];
muaids = [1 2 4 5 6 11 12 13 14];
su_probes = [Blocks{1}.suprobes];
mu_probes = [Blocks{1}.muprobes];
all_probes = [su_probes mu_probes];
spk_cnts = [spk_cnts(:,cellids) mua_cnts(:,muaids)];

n_used_cells = size(spk_cnts,2);

%%
Fs = 1000;
dsf = 4;
Fsd = Fs/dsf;
niqf = Fs/2;
[b,a] = butter(2,120/niqf,'low');
[b_alpha,a_alpha] = butter(2,[5 12]/niqf);
% alpha_ch = 18;
%%
all_lfp_pow = [];
within_fix_avgs = [];
fix_nums = [];
all_used_inds = [];
all_alpha_phase = [];
all_alpha_uphase = [];
spikes_binned = [];

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
    
    lfp_samps = filtfilt(b,a,lfp_samps);
    lfp_sampsd = downsample(lfp_samps,dsf);
    lfp_timed = downsample(lfp_time,dsf);
    lfp_alpha = filtfilt(b_alpha,a_alpha,lfp_samps);
    lfp_alpha = zscore(downsample(lfp_alpha,dsf));
    
    cur_all_model = find(all_model_blockids==blockid);
    %     alpha_phase = angle(hilbert(lfp_alpha(:,alpha_ch)));
    %     alpha_phase = unwrap_phase_monotonic(alpha_phase)';
    clear alpha_phase alpha_uphase
    for c = 1:24
        alpha_phase(:,c) = angle(hilbert(lfp_alpha(:,c)));
        alpha_uphase(:,c) = unwrap_phase_monotonic(alpha_phase(:,c));
    end
    interp_alpha_phase = interp1(lfp_timed,alpha_phase,all_model_time_axis(cur_all_model));
    all_alpha_phase = [all_alpha_phase; interp_alpha_phase];
     interp_alpha_uphase = interp1(lfp_timed,alpha_uphase,all_model_time_axis(cur_all_model));
    all_alpha_uphase = [all_alpha_uphase; interp_alpha_uphase];
   
end

%%
all_model_time_axis = [];
all_model_blockids = [];
all_model_fixids = [];
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
        cur_tedges = start_T:dt:end_T;
        cur_tcents = 0.5*cur_tedges(1:end-1)+0.5*cur_tedges(2:end);
        all_model_time_axis = [all_model_time_axis cur_tcents];
        all_model_blockids = [all_model_blockids blockid*ones(1,length(cur_tcents))];
        all_model_fixids = [all_model_fixids cur_set(i)*ones(1,length(cur_tcents))];
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
cur_used_fixs = unique(all_model_fixids);
n_fixs = length(cur_used_fixs);
trial_start_inds = 1 + find(diff(all_model_fixids) ~= 0);
trial_stop_inds = find(diff(all_model_fixids) ~= 0);
trial_start_inds = [1 trial_start_inds];
trial_stop_inds = [trial_stop_inds length(all_model_fixids)];

% phase_since_trial = nan(size(all_alpha_phase));
% for i = 1:n_fixs
%     cur_set = find(all_model_fixids == cur_used_fixs(i));
% %    phase_since_trial(cur_set) = all_alpha_phase(cur_set) - all_alpha_phase(cur_set(1));
%    phase_since_trial(cur_set,:) = bsxfun(@minus,all_alpha_phase(cur_set,:),all_alpha_phase(cur_set(1),:));
% 
% end
% phase_since_trial = phase_since_trial(:,2);

% for ll = 1:length(use_lfps)
%     zero_crossing_inds{ll} = find(full_alpha_phase(1:end-1,ll) < 0 & full_alpha_phase(2:end,ll) > 0);
% end
phase_deriv = [zeros(1,24); diff(all_alpha_phase)];
for ll = 1:24
    pi_crossing_inds{ll} = find(phase_deriv(1:end-1,ll) > -1 & phase_deriv(2:end,ll) < -1);
end

% buffer_t = 0.2;
% buffer_inds = round(buffer_t/dt);
ref_elec = 16;
phase_start_inds = zeros(length(trial_start_inds),1);
phase_since_trial = nan(length(trial_start_inds),24);
for i = 1:n_fixs
    %     ll = 14;
    %     for ll = 1:24
    %         cur_zdiffs = abs(zero_crossing_inds{ll} - trial_start_inds(i));
    %         [~,b] = min(cur_zdiffs);
    %         phase_start_inds(i,ll) = zero_crossing_inds{ll}(b);
    next_pcross = find(pi_crossing_inds{ref_elec} >= trial_start_inds(i),1,'first');
    if pi_crossing_inds{ref_elec}(next_pcross) < trial_stop_inds(i)
        phase_start_inds(i) = pi_crossing_inds{ref_elec}(next_pcross);
    else
        prev_pcross = find(pi_crossing_inds{ref_elec} < trial_start_inds(i),1,'last');
        phase_start_inds(i) = pi_crossing_inds{ref_elec}(prev_pcross);
    end
    cur_set = phase_start_inds(i):trial_stop_inds(i);
    phase_since_trial(cur_set,:) = bsxfun(@minus,all_alpha_uphase(cur_set,:),all_alpha_uphase(cur_set(1),:));
    %     end
end

% full_unique_trial = nan(length(full_t),length(use_lfps));
% for i = 1:length(trial_start_inds)
%     for ll = 1:length(use_lfps)
%         cur_set = phase_start_inds(i,ll):trial_stop_inds(i);
%         full_unique_trial(cur_set,ll) = i;
%     end
% end
% 
%% Compute TBR time-since fix onset
max_phase = 2*pi*5;
min_phase = 0;
nbins = 60;
% tax = linspace(0,max_phase,nbins);
uset = find(phase_since_trial(:,ref_elec) <= max_phase & phase_since_trial(:,ref_elec) >= min_phase);
tax = linspace(min_phase,max_phase,nbins);

Tmat = [];
for ll = 1:24
Tmat = [Tmat tbrep(phase_since_trial(:,ll),tax)];
end

%%
l2_ind = 500;
l2_dep = 1000;
l1 = 2;
silent = 1;
for cc = 1:n_used_cells
%     uset = find(phase_since_trial(:,cc) <= max_phase);
% Tmat = tbrep(phase_since_trial(:,cc),tax);

    disp('Fitting total stim model')
    
    Robs = spikes_binned(uset,cc);
    spkbs = convert_to_spikebins(Robs);


    X_resh = reshape(X,size(X,1),SDIM^2);
    cur_mask1 = get_pgabor_mask(gabor_params_fin(cc,1:6),0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(gabor_params_fin(cc,1:6),pi/2,[SDIM SDIM]);
    mask1_out = X_resh*cur_mask1(:);
    mask2_out = X_resh*cur_mask2(:);
    
    energy_out = gabor_params_fin(cc,7)*sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
    lin_out = gabor_params_fin(cc,8)*mask1_out(all_model_fixids) + gabor_params_fin(cc,9)*mask1_out(all_model_fixids);
    
    total_out = lin_out + energy_out;
    total_out = zscore(total_out);
    
%     Xmat = [bsxfun(@times,Tmat,total_out) Tmat];
     Xmat = [Tmat];
       
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    lamrange2 = [];
    for ll = 1:24
        lamrange2 = [lamrange2; l2_dep (ll-1)*nbins+1 ll*nbins 0];
    end
%     for ll = 1:24
%         lamrange2 = [lamrange2; l2_ind 24*nbins+(ll-1)*nbins+1 24*nbins+ll*nbins 0];
%     end
    llist = [l1 1:klen];
    [fitp,grad] = GLMsolve_jmm(Xmat(uset,:), spkbs, K0, silent, [], lamrange2, [], [],llist, [], 0);
        k(cc,:) = fitp.k(1:end-1);
% k_ind = k(24*nbins+1:end);
% k_dep = k(1:24*nbins);

end

%%
cd ~/Data/bruce/2_27_12
save unit_alpha_stimmod_ch18_peakaligned stim_* ov_const tax

%%
close all
for c = 1:18
    c
    subplot(2,1,1)
    plot(tax/(2*pi),stim_dep_kern(c,:))
    subplot(2,1,2)
    plot(tax/(2*pi),stim_ind_kern(c,:))
    pause
    clf
end