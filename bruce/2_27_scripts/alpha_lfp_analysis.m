clear all
% close all

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

cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_data_v3

% cd /Users/James/Data/bruce/2_27_12/stimrecon
% load sing_eye_pgabortrack_d2
% load sing_eye_pgabortrack_fin_est_d2
cd ~/Data/bruce/2_27_12/ExptA

% SDIM = length(xpatch_inds);
% cd /Users/James/Data/bruce/2_27_12/stimrecon
% load fixation_image_patches_corrected
% cd /Users/James/James_scripts/bruce/modelfits/
% % load fixation_gabor_models
% load gabor_tempmodfits_en_lin


Fs = 1000;
niqf = Fs/2;
[b,a] = butter(2,[1 450]/niqf);
scales = logspace(log10(3.25),log10(50),25);
wfreqs = scal2frq(scales,'cmor1-1',1/Fs);
nwfreqs = length(wfreqs);

n_phase_bins = 15;
phase_bins = linspace(-pi,pi,n_phase_bins+1);
%%
cd ~/Data/bruce/2_27_12

for blockid = 1:4;
    
    fprintf('Processing block %d of %d...\n',blockid,4);
    
    block_times = Blocks{blockid}.blocktimes;
    stim_times = Blocks{blockid}.stimtime;
    stimID = Blocks{blockid}.stimids;
    
    n_stims = length(stim_times);
    
    cd ~/Data/bruce/2_27_12/saccades/
    load(sprintf('lemM232.5%d.em.sac.mat',blockid))
    %     load(sprintf('lemM232.5%d.em.hor.mat',blockid))
    
    % identify saccade start and stop times
    EyeStartT = Expt.Trials.Start/10000; % time of first eye sample
    EyeEndT = Expt.Trials.End/10000; % time of last eye sample
    Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
    eyets = EyeStartT:Eyedt:EyeEndT; %eye tracking time axis (sec)
    sac_buffer_inds = round(sac_buffer/Eyedt);
    
    reye_pos = [Expt.Trials.Eyevals.rh Expt.Trials.Eyevals.rv];
    leye_pos = [Expt.Trials.Eyevals.lh Expt.Trials.Eyevals.lv];
    
%     avg_eyepos = (reye_pos + leye_pos)/2;
%     clear sm_avg_eyepos eye_vel
%     sm_avg_eyepos(:,1) = smooth(avg_eyepos(:,1),3);
%     sm_avg_eyepos(:,2) = smooth(avg_eyepos(:,2),3);
%     eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/Eyedt;
%     eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/Eyedt;
%     eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
%     vel_thresh = 5;
%     eye_vel(eye_vel(:,1) > vel_thresh,1) = vel_thresh;
%     eye_vel(eye_vel(:,1) < -vel_thresh,1) = -vel_thresh;
%     eye_vel(eye_vel(:,2) > vel_thresh,2) = vel_thresh;
%     eye_vel(eye_vel(:,2) < -vel_thresh,2) = -vel_thresh;
    
    %%
    cd /Users/James/Data/bruce/2_27_12
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
    bad_t = 1+find(diff(lfp_time)==0);
    lfp_time(bad_t) = [];
    lfp_samps(bad_t,:) = [];
    % lfp_samps = zscore(lfp_samps);
    
    [b,a] = butter(2,[1 450]/niqf);
    lfp_samps = filtfilt(b,a,lfp_samps);
    [b,a] = butter(2,[4 12]/niqf);
    alpha_ch = 2;
    lfp_samps_alpha = filtfilt(b,a,lfp_samps);
    lfp_samps_alpha_phase = angle(hilbert(lfp_samps_alpha));
    phase_revs = 1+find(diff(lfp_samps_alpha_phase) < -pi);
    [b,a] = butter(2,[25 60]/niqf);
    lfp_samps_gamma = filtfilt(b,a,lfp_samps);
    lfp_samps_gamma_phase = angle(hilbert(lfp_samps_gamma));
    
    used_ch = 1:1:24;
    NT = size(lfp_samps,1);
    ampgram = nan(NT,length(wfreqs),length(used_ch));
    for t = 1:length(used_ch)
        t
        ampgram(:,:,t) = zscore(abs(cwt(lfp_samps(:,used_ch(t)),scales,'cmor1-1'))');
    end
    %%
    spikes_binned = zeros(24,length(lfp_time));
    for c = 1:10
        c
        temp = round(interp1(lfp_time,1:length(lfp_time),Blocks{blockid}.spktimes{c}));
        temp(isnan(temp)) = [];
        all_spike_lfp_inds{c} = temp;
        
        cur_hist = hist(Blocks{blockid}.spktimes{c},lfp_time);
        cur_hist(end) = 0;
        spikes_binned(c,:) = smooth(cur_hist,10);
    end
     for c = 1:14
         c
        temp = round(interp1(lfp_time,1:length(lfp_time),Blocks{blockid}.mutimes{c}));
        temp(isnan(temp)) = [];
        all_mua_lfp_inds{c} = temp;
        
        cur_hist = hist(Blocks{blockid}.mutimes{c},lfp_time);
        cur_hist(end) = 0;
        spikes_binned(c+10,:) = smooth(cur_hist,10);
     end
     
    spike_rates = spikes_binned*Fs;
    spike_rates = zscore(spike_rates')';
    net_rate = mean(spike_rates);
    
    suprobes = Blocks{1}.suprobes;
    muprobes = Blocks{1}.muprobes;
    allprobes = [suprobes muprobes];
    
    %%
    min_fix_durs = 0.3;
    cur_fixs = find(blockids == blockid);
    cur_fix_start_inds = round(interp1(lfp_time,1:length(lfp_time),all_fix_start_times(cur_fixs)));
    cur_fix_stop_inds = round(interp1(lfp_time,1:length(lfp_time),all_fix_stop_times(cur_fixs)));
    cur_nfixs = length(cur_fixs);
    fix_durs = all_fix_stop_times(cur_fixs)-all_fix_start_times(cur_fixs);
    used_fixs = find(fix_durs > min_fix_durs & ~isnan(cur_fix_start_inds) & ~isnan(cur_fix_stop_inds));
    
    in_beg = zeros(length(lfp_time),1);
    in_mid = zeros(length(lfp_time),1);
    in_late = zeros(length(lfp_time),1);
    for i = 1:length(used_fixs)
        cur_set = cur_fix_start_inds(used_fixs(i)) + (1:round(0.15*Fs));
        in_beg(cur_set) = 1;
        cur_set = cur_fix_start_inds(used_fixs(i)) + round(0.15*Fs) + (1:round(0.15*Fs));
        in_mid(cur_set) = 1;
        cur_set = cur_fix_start_inds(used_fixs(i)):cur_fix_stop_inds(used_fixs(i));
        cur_set(1:round(0.15*Fs)) = [];
        in_late(cur_set) = 1;
    end
    
    %% COMPUTE SPIKE PHASE-LOCKING STRENGTHS 
    su_alpha_lock(blockid,:,:) = nan(10,24);
    mu_alpha_lock(blockid,:,:) = nan(14,24);
    su_alpha_lock_beg(blockid,:,:) = nan(10,24);
    mu_alpha_lock_beg(blockid,:,:) = nan(14,24);
    su_alpha_lock_late(blockid,:,:) = nan(10,24);
    mu_alpha_lock_late(blockid,:,:) = nan(14,24);
    su_gamma_lock(blockid,:,:) = nan(10,24);
    mu_gamma_lock(blockid,:,:) = nan(14,24);
    su_gamma_lock_beg(blockid,:,:) = nan(10,24);
    mu_gamma_lock_beg(blockid,:,:) = nan(14,24);
    su_gamma_lock_late(blockid,:,:) = nan(10,24);
    mu_gamma_lock_late(blockid,:,:) = nan(14,24);
    for c = 1:10
        cur_set = all_spike_lfp_inds{c};
        cur_beg_set = cur_set(in_beg(cur_set) ==1);
        cur_late_set = cur_set(in_late(cur_set) == 1);
        for ch = 1:24
            su_alpha_lock(blockid,c,ch) = circ_kappa(lfp_samps_alpha_phase(cur_set,ch));
            su_alpha_lock_beg(blockid,c,ch) = circ_kappa(lfp_samps_alpha_phase(cur_beg_set,ch));
            su_alpha_lock_late(blockid,c,ch) = circ_kappa(lfp_samps_alpha_phase(cur_late_set,ch));
            su_gamma_lock(blockid,c,ch) = circ_kappa(lfp_samps_gamma_phase(cur_set,ch));
            su_gamma_lock_beg(blockid,c,ch) = circ_kappa(lfp_samps_gamma_phase(cur_beg_set,ch));
            su_gamma_lock_late(blockid,c,ch) = circ_kappa(lfp_samps_gamma_phase(cur_late_set,ch));
        end
    end
    for c = 1:14
        cur_set = all_mua_lfp_inds{c};
        cur_beg_set = cur_set(in_beg(cur_set) ==1);
        cur_late_set = cur_set(in_late(cur_set) == 1);
        for ch = 1:24
            mu_alpha_lock(blockid,c,ch) = circ_kappa(lfp_samps_alpha_phase(cur_set,ch));
            mu_alpha_lock_beg(blockid,c,ch) = circ_kappa(lfp_samps_alpha_phase(cur_beg_set,ch));
            mu_alpha_lock_late(blockid,c,ch) = circ_kappa(lfp_samps_alpha_phase(cur_late_set,ch));
            mu_gamma_lock(blockid,c,ch) = circ_kappa(lfp_samps_gamma_phase(cur_set,ch));
            mu_gamma_lock_beg(blockid,c,ch) = circ_kappa(lfp_samps_gamma_phase(cur_beg_set,ch));
            mu_gamma_lock_late(blockid,c,ch) = circ_kappa(lfp_samps_gamma_phase(cur_late_set,ch));
        end
    end
    
    
    %%
    phase_amp_corr(blockid,:,:) = nan(24,length(wfreqs));
    phase_amp_p(blockid,:,:) = nan(24,length(wfreqs));
    phase_amp_corr_beg(blockid,:,:) = nan(24,length(wfreqs));
    phase_amp_p_beg(blockid,:,:) = nan(24,length(wfreqs));
    phase_amp_corr_late(blockid,:,:) = nan(24,length(wfreqs));
    phase_amp_p_late(blockid,:,:) = nan(24,length(wfreqs));
    
    gphase_amp_corr(blockid,:,:) = nan(24,length(wfreqs));
    gphase_amp_p(blockid,:,:) = nan(24,length(wfreqs));
    gphase_amp_corr_beg(blockid,:,:) = nan(24,length(wfreqs));
    gphase_amp_p_beg(blockid,:,:) = nan(24,length(wfreqs));
    gphase_amp_corr_late(blockid,:,:) = nan(24,length(wfreqs));
    gphase_amp_p_late(blockid,:,:) = nan(24,length(wfreqs));
    
%     bphase_amp_corr(blockid,:,:) = nan(24,length(wfreqs));
%     bphase_amp_p(blockid,:,:) = nan(24,length(wfreqs));
    bphase_amp_corr_beg(blockid,:,:) = nan(24,length(wfreqs));
    bphase_amp_p_beg(blockid,:,:) = nan(24,length(wfreqs));
    bphase_amp_corr_late(blockid,:,:) = nan(24,length(wfreqs));
    bphase_amp_p_late(blockid,:,:) = nan(24,length(wfreqs));

    for ch = 1:24
        ch
        for ww = 1:length(wfreqs)
            [rho,p] = circ_corrcl(lfp_samps_alpha_phase(:,ch),ampgram(:,ww,ch));
            phase_amp_corr(blockid,ch,ww) = rho;
            phase_amp_p(blockid,ch,ww) = p;
            [rho,p] = circ_corrcl(lfp_samps_alpha_phase(in_beg==1,ch),ampgram(in_beg==1,ww,ch));
            phase_amp_corr_beg(blockid,ch,ww) = rho;
            phase_amp_p_beg(blockid,ch,ww) = p;
            [rho,p] = circ_corrcl(lfp_samps_alpha_phase(in_late==1,ch),ampgram(in_late==1,ww,ch));
            phase_amp_corr_late(blockid,ch,ww) = rho;
            phase_amp_p_late(blockid,ch,ww) = p;
            
%             [rho,p] = circ_corrcl(lfp_samps_alpha_phase(:,alpha_ch),ampgram(:,ww,ch));
%             phase_amp_corr(blockid,ch,ww) = rho;
%             phase_amp_p(blockid,ch,ww) = p;
            [rho,p] = circ_corrcl(lfp_samps_alpha_phase(in_beg==1,alpha_ch),ampgram(in_beg==1,ww,ch));
            bphase_amp_corr_beg(blockid,ch,ww) = rho;
            bphase_amp_p_beg(blockid,ch,ww) = p;
            [rho,p] = circ_corrcl(lfp_samps_alpha_phase(in_late==1,alpha_ch),ampgram(in_late==1,ww,ch));
            bphase_amp_corr_late(blockid,ch,ww) = rho;
            bphase_amp_p_late(blockid,ch,ww) = p;
            
            [rho,p] = circ_corrcl(lfp_samps_gamma_phase(:,ch),ampgram(:,ww,ch));
            gphase_amp_corr(blockid,ch,ww) = rho;
            gphase_amp_p(blockid,ch,ww) = p;
            [rho,p] = circ_corrcl(lfp_samps_gamma_phase(in_beg==1,ch),ampgram(in_beg==1,ww,ch));
            gphase_amp_corr_beg(blockid,ch,ww) = rho;
            gphase_amp_p_beg(blockid,ch,ww) = p;
            [rho,p] = circ_corrcl(lfp_samps_gamma_phase(in_late==1,ch),ampgram(in_late==1,ww,ch));
            gphase_amp_corr_late(blockid,ch,ww) = rho;
            gphase_amp_p_late(blockid,ch,ww) = p;
        end
    end
    
    
    
    
    %%
%     clear phase_dep*
%     for p = 1:n_phase_bins
%         cur_phase_set = find(lfp_samps_alpha_phase >= phase_bins(p) & lfp_samps_alpha_phase < phase_bins(p+1));
%         cur_phase_set_beg = cur_phase_set(in_beg(cur_phase_set) == 1);
%         cur_phase_set_mid = cur_phase_set(in_mid(cur_phase_set) == 1);
%         cur_phase_set_late = cur_phase_set(in_late(cur_phase_set) == 1);
%         
%         phase_dep_rates(p,:) = mean(spike_rates(:,cur_phase_set),2);
%         phase_dep_rates_beg(p,:) = mean(spike_rates(:,cur_phase_set_beg),2);
%         phase_dep_rates_mid(p,:) = mean(spike_rates(:,cur_phase_set_mid),2);
%         phase_dep_rates_late(p,:) = mean(spike_rates(:,cur_phase_set_late),2);
%        
%         phase_dep_pow(p,:,:) = mean(ampgram(cur_phase_set,:,:));
%         phase_dep_pow_beg(p,:,:) = mean(ampgram(cur_phase_set_beg,:,:));
%         phase_dep_pow_mid(p,:,:) = mean(ampgram(cur_phase_set_mid,:,:));
%         phase_dep_pow_late(p,:,:) = mean(ampgram(cur_phase_set_late,:,:));
%         
%         cur_gphase_set = find(lfp_samps_gamma_phase >= phase_bins(p) & lfp_samps_gamma_phase < phase_bins(p+1));
%         cur_gphase_set_beg = cur_gphase_set(in_beg(cur_gphase_set) == 1);
%         cur_gphase_set_mid = cur_gphase_set(in_mid(cur_gphase_set) == 1);
%         cur_gphase_set_late = cur_gphase_set(in_late(cur_gphase_set) == 1);
%         gphase_dep_rates(p,:) = mean(spike_rates(:,cur_gphase_set),2);
%         gphase_dep_rates_beg(p,:) = mean(spike_rates(:,cur_gphase_set_beg),2);
%         gphase_dep_rates_mid(p,:) = mean(spike_rates(:,cur_gphase_set_mid),2);
%         gphase_dep_rates_late(p,:) = mean(spike_rates(:,cur_gphase_set_late),2);
%         gphase_dep_pow(p,:,:) = mean(ampgram(cur_gphase_set,:,:));
%         gphase_dep_pow_beg(p,:,:) = mean(ampgram(cur_gphase_set_beg,:,:));
%         gphase_dep_pow_mid(p,:,:) = mean(ampgram(cur_gphase_set_mid,:,:));
%         gphase_dep_pow_late(p,:,:) = mean(ampgram(cur_gphase_set_late,:,:));
%    end
%     all_phase_dep_rates(blockid,:,:) = phase_dep_rates;
%     all_phase_dep_rates_beg(blockid,:,:) = phase_dep_rates_beg;
%     all_phase_dep_rates_mid(blockid,:,:) = phase_dep_rates_mid;
%     all_phase_dep_rates_late(blockid,:,:) = phase_dep_rates_late;
%     
%     all_gphase_dep_rates(blockid,:,:) = gphase_dep_rates;
%     all_gphase_dep_rates_beg(blockid,:,:) = gphase_dep_rates_beg;
%     all_gphase_dep_rates_mid(blockid,:,:) = gphase_dep_rates_mid;
%     all_gphase_dep_rates_late(blockid,:,:) = gphase_dep_rates_late;
% 
%     all_phase_dep_pow(blockid,:,:,:) = phase_dep_pow;
%     all_phase_dep_pow_beg(blockid,:,:,:) = phase_dep_pow_beg;
%     all_phase_dep_pow_mid(blockid,:,:,:) = phase_dep_pow_mid;
%     all_phase_dep_pow_late(blockid,:,:,:) = phase_dep_pow_late;
%  
%     all_gphase_dep_pow(blockid,:,:,:) = gphase_dep_pow;
%     all_gphase_dep_pow_beg(blockid,:,:,:) = gphase_dep_pow_beg;
%     all_gphase_dep_pow_mid(blockid,:,:,:) = gphase_dep_pow_mid;
%     all_gphase_dep_pow_late(blockid,:,:,:) = gphase_dep_pow_late;
end

%%
% avg_phase_dep_rates_late = squeeze(mean(all_phase_dep_rates_late));
% avg_gphase_dep_rates_late = squeeze(mean(all_gphase_dep_rates_late));
% avg_phase_dep_rates_beg = squeeze(mean(all_phase_dep_rates_beg));
% avg_gphase_dep_rates_beg = squeeze(mean(all_gphase_dep_rates_beg));
% 
% avg_phase_dep_pow_late = squeeze(mean(all_phase_dep_pow_late));
% avg_gphase_dep_pow_late = squeeze(mean(all_gphase_dep_pow_late));
% avg_phase_dep_pow_beg = squeeze(mean(all_phase_dep_pow_beg));
% avg_gphase_dep_pow_beg = squeeze(mean(all_gphase_dep_pow_beg));

%%
save alpha_phase_analysis_v2 phase_* bphase* gphase_* su_* mu_*

%%
% phase_dep_pow_mod_beg = squeeze(max(avg_phase_dep_pow_beg)) - squeeze(min(avg_phase_dep_pow_beg));
% phase_dep_pow_mod_late = squeeze(max(avg_phase_dep_pow_late)) - squeeze(min(avg_phase_dep_pow_late));
% 
% gphase_dep_pow_mod_beg = squeeze(max(avg_gphase_dep_pow_beg)) - squeeze(min(avg_gphase_dep_pow_beg));
% gphase_dep_pow_mod_late = squeeze(max(avg_gphase_dep_pow_late)) - squeeze(min(avg_gphase_dep_pow_late));
% 
% 
% 
% %%
% norm_phase_dep_rates_late = bsxfun(@minus,avg_phase_dep_rates_late,mean(avg_phase_dep_rates_late));
% norm_phase_dep_rates_beg = bsxfun(@minus,avg_phase_dep_rates_beg,mean(avg_phase_dep_rates_beg));
% norm_gphase_dep_rates_late = bsxfun(@minus,avg_gphase_dep_rates_late,mean(avg_gphase_dep_rates_late));
% norm_gphase_dep_rates_beg = bsxfun(@minus,avg_gphase_dep_rates_beg,mean(avg_gphase_dep_rates_beg));

%%
avg_su_alpha_lock = squeeze(mean(su_alpha_lock));
avg_su_alpha_lock_beg = squeeze(mean(su_alpha_lock_beg));
avg_su_alpha_lock_late = squeeze(mean(su_alpha_lock_late));
avg_mu_alpha_lock = squeeze(mean(mu_alpha_lock));
avg_mu_alpha_lock_beg = squeeze(mean(mu_alpha_lock_beg));
avg_mu_alpha_lock_late = squeeze(mean(mu_alpha_lock_late));
avg_su_gamma_lock = squeeze(mean(su_gamma_lock));
avg_su_gamma_lock_beg = squeeze(mean(su_gamma_lock_beg));
avg_su_gamma_lock_late = squeeze(mean(su_gamma_lock_late));
avg_mu_gamma_lock = squeeze(mean(mu_gamma_lock));
avg_mu_gamma_lock_beg = squeeze(mean(mu_gamma_lock_beg));
avg_mu_gamma_lock_late = squeeze(mean(mu_gamma_lock_late));

su_probes = [Blocks{1}.suprobes];
mu_probes = [Blocks{1}.muprobes];
all_probes = [su_probes mu_probes];
[~,probe_ord] = sort(all_probes);

all_gamma_lock = [avg_su_gamma_lock; avg_mu_gamma_lock];
all_gamma_lock_beg = [avg_su_gamma_lock_beg; avg_mu_gamma_lock_beg];
all_gamma_lock_late = [avg_su_gamma_lock_late; avg_mu_gamma_lock_late];

all_alpha_lock = [avg_su_alpha_lock; avg_mu_alpha_lock];
all_alpha_lock_beg = [avg_su_alpha_lock_beg; avg_mu_alpha_lock_beg];
all_alpha_lock_late = [avg_su_alpha_lock_late; avg_mu_alpha_lock_late];

avg_phase_amp_corr = squeeze(mean(phase_amp_corr));
avg_phase_amp_corr_beg = squeeze(mean(phase_amp_corr_beg));
avg_phase_amp_corr_late = squeeze(mean(phase_amp_corr_late));
avg_bphase_amp_corr_beg = squeeze(mean(bphase_amp_corr_beg));
avg_bphase_amp_corr_late = squeeze(mean(bphase_amp_corr_late));

avg_gphase_amp_corr = squeeze(mean(gphase_amp_corr));
avg_gphase_amp_corr_beg = squeeze(mean(gphase_amp_corr_beg));
avg_gphase_amp_corr_late = squeeze(mean(gphase_amp_corr_late));
