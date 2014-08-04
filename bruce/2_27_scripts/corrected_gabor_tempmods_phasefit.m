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

cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_data_v3

cd /Users/James/Data/bruce/2_27_12/stimrecon
load sing_eye_pgabortrack_d2
load sing_eye_pgabortrack_fin_est_d2
cd ~/Data/bruce/2_27_12/ExptA

SDIM = length(xpatch_inds);
cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_image_patches_corrected

NT = length(used_fixs);
X_resh = reshape(X,NT,SDIM^2);
% dt = .01;
dt = .0025;

cellids = [1 2 3 4 5 6 7 8 10];

spk_cnts = [spk_cnts(:,cellids)];

n_used_cells = size(spk_cnts,2);

%%
fix_stimout = zeros(n_used_cells,length(used_fixs));
all_model_time_axis = [];
all_model_blockids = [];
all_model_fixids = [];
spikes_binned = [];
time_since_fix = [];
for blockid = 1:3;
    fprintf('Block %d of %d\n',blockid,3);
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
        for c = 1:10
            if c <= 10
                temp = histc(Blocks{blockid}.spktimes{c},cur_tedges);
            else
                temp = histc(Blocks{blockid}.mutimes{c-10},cur_tedges);
            end
            temp_binned(c,:) = temp(1:end-1);
        end
        spikes_binned = [spikes_binned; temp_binned'];
        
        time_since_fix = [time_since_fix cur_tcents-start_T];
        
    end
end

%%
Fs = 1000;
dsf = 2;
Fsd = Fs/dsf;
niqf = Fs/2;
[b,a] = butter(2,[1 200]/niqf);
[b_gam,a_gam] = butter(2,[25 50]/niqf);
[b_the,a_the] = butter(2,[6 12]/niqf);

% scales = [11];
% % scales = logspace(log10(4),log10(80),10);
% wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
% nwfreqs = length(wfreqs);

all_lfp_amp_g = [];
all_lfp_phase_g = [];
all_lfp_amp_t = [];
all_lfp_phase_t = [];
all_used_inds = [];
for blockid = 1:3;
    fprintf('Block %d of %d\n',blockid,3);
    
    %%
    cd /Users/James/Data/bruce/2_27_12
    load(sprintf('lemM232A.5%d.lfp.mat',blockid));
    %get start times of each LFP trial
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
        lfp_time = [lfp_time linspace(lfp_trial_start(i),lfp_trial_start(i)+lfp_dur(i),size(LFP.Trials(i).LFP,1))];
        lfp_samps = [lfp_samps; LFP.Trials(i).LFP];
    end
    
    lfp_samps_g = filtfilt(b_gam,a_gam,lfp_samps);
    lfp_sampsd_g = zscore(downsample(lfp_samps_g,dsf));
    lfp_samps_t = filtfilt(b_the,a_the,lfp_samps);
    lfp_sampsd_t = zscore(downsample(lfp_samps_t,dsf));
    clear lfp_samps
    
    %     lfp_samps = filtfilt(b,a,lfp_samps);
    %     lfp_samps_hf = filtfilt(b_hf,a_hf,lfp_samps);
    %     lfp_sampsd = zscore(downsample(lfp_samps,dsf));
    %     lfp_samps_hf = abs(hilbert(lfp_samps_hf));
    %     lfp_sampsd_hf = downsample(lfp_samps_hf,dsf);
    %     lfp_sampsd_hf = zscore(lfp_sampsd_hf);
    %     clear lfp_samps lfp_samps_hf
    
    lfp_timed = downsample(lfp_time,dsf)';
    
    cur_all_model = find(all_model_blockids==blockid);
    
    %     interp_ampgram = zeros(24,length(cur_all_model),length(wfreqs));
    %     interp_phasegram = zeros(24,length(cur_all_model),length(wfreqs));
    %     for cc = 1:24
    %         fprintf('Channel %d of %d\n',cc,24);
    %         temp = cwt(lfp_sampsd(:,cc),scales,'cmor1-1');
    %         temp = temp';
    %         tampgram = abs(temp);
    %         tphasegram = angle(temp);
    %         interp_ampgram(cc,:,:) = interp1(lfp_timed,tampgram,all_model_time_axis(cur_all_model));
    %         interp_phasegram(cc,:,:) = interp1(lfp_timed,tphasegram,all_model_time_axis(cur_all_model));
    %     end
    
    %     interp_hf = interp1(lfp_timed,lfp_sampsd_hf,all_model_time_axis(cur_all_model));
    interp_lfp_g = interp1(lfp_timed,lfp_sampsd_g,all_model_time_axis(cur_all_model));
    interp_lfp_t = interp1(lfp_timed,lfp_sampsd_t,all_model_time_axis(cur_all_model));
    used_tpoints = find(~isnan(interp_lfp_g(:,1)));
    
    h_tr = hilbert(interp_lfp_g(used_tpoints,:));
    interp_phase = angle(h_tr);
%     interp_amp = abs(h_tr);
    all_lfp_phase_g = [all_lfp_phase_g; interp_phase];
%     all_lfp_amp_g = [all_lfp_amp_g; interp_amp];

    h_tr = hilbert(interp_lfp_t(used_tpoints,:));
    interp_phase = angle(h_tr);
%     interp_amp = abs(h_tr);
    all_lfp_phase_t = [all_lfp_phase_t; interp_phase];
%     all_lfp_amp_t = [all_lfp_amp_t; interp_amp];

    %     used_tpoints = find(~isnan(interp_ampgram(1,:,1)));
    %     all_lfp_phase = cat(2,all_lfp_phase, interp_phasegram(:,used_tpoints,:));
    %     all_lfp_amp = cat(2,all_lfp_amp, interp_ampgram(:,used_tpoints,:));
    %     all_lfp_hf = [all_lfp_hf; interp_hf(used_tpoints,:)];
    all_used_inds = [all_used_inds; cur_all_model(used_tpoints)'];
    
end

%% Compute TBR time-since fix onset
max_tsf = 0.7; nbins = 30;
used_tsf = time_since_fix;
used_tsf(used_tsf > max_tsf) = [];
tax = prctile(used_tsf,linspace(100/nbins,100-100/nbins,nbins));

% tax = logspace(log10(dt),log10(max_tsf),nbins);
Tmat = tbrep(time_since_fix,tax);

%%
xv_frac = 0;
rperm = randperm(length(used_fixs));
xv_fixs = rperm(1:round(xv_frac*length(used_fixs)));
% xv_inds = find(ismember(all_model_fixids,xv_fixs));
% tr_inds = find(~ismember(all_model_fixids,xv_fixs));

xv_inds = all_used_inds(ismember(all_model_fixids(all_used_inds),xv_fixs));
tr_inds = all_used_inds(~ismember(all_model_fixids(all_used_inds),xv_fixs));
tr_inds_lfp = find(~ismember(all_model_fixids(all_used_inds),xv_fixs));
xv_inds_lfp = find(ismember(all_model_fixids(all_used_inds),xv_fixs));

% xv_fixs = find(all_stim_filtered(used_fixs)==0);
% xv_inds = find(ismember(all_model_fixids,xv_fixs));
% tr_inds = find(~ismember(all_model_fixids,xv_fixs));
%% COMPUTE SACCADE TRIG AVG FIRING RATES
sac_trg_rate = zeros(n_used_cells,length(tax));
sac_trg_cnt = zeros(n_used_cells,length(tax));
sac_trg_mphase_g = zeros(n_used_cells,length(tax));
sac_trg_kphase_g = zeros(n_used_cells,length(tax));
sac_trg_mphase_t = zeros(n_used_cells,length(tax));
sac_trg_kphase_t = zeros(n_used_cells,length(tax));
taxb = [0 tax];
for i = 1:length(tax)
    curset = find(time_since_fix(all_used_inds) >= taxb(i) & time_since_fix(all_used_inds) < taxb(i+1));
    curset2 = all_used_inds(curset);
    for c = 1:n_used_cells
        cur_probe = Blocks{1}.suprobes(cellids(c));
        sac_trg_rate(c,i) = mean(spikes_binned(curset2,c));
        sac_trg_cnt(c,i) = sum(spikes_binned(curset2,c));
%         for ww = 1:length(wfreqs)
%             sac_trg_mphase(c,ww,i) = circ_mean(all_lfp_phase(cur_probe,curset,ww)',spikes_binned(curset2,c));
%             sac_trg_kphase(c,ww,i) = circ_kappa(all_lfp_phase(cur_probe,curset,ww)',spikes_binned(curset2,c));
            sac_trg_mphase_g(c,i) = circ_mean(all_lfp_phase_g(curset,cur_probe),spikes_binned(curset2,c));
            sac_trg_kphase_g(c,i) = circ_kappa(all_lfp_phase_g(curset,cur_probe),spikes_binned(curset2,c));
            sac_trg_mphase_t(c,i) = circ_mean(all_lfp_phase_t(curset,cur_probe),spikes_binned(curset2,c));
            sac_trg_kphase_t(c,i) = circ_kappa(all_lfp_phase_t(curset,cur_probe),spikes_binned(curset2,c));
%         end
    end
    n_occ(i) = length(curset);
end
avg_rates = mean(spikes_binned)/dt;
sac_trg_rate = sac_trg_rate/dt;

%%
clear mod_*
sac_hist = zeros(n_used_cells,length(tax));
stim_kern = cell(n_used_cells,1);
slope2_lambda = 150;
l1_pen = 0;
% stim_mod_predrate = zeros(9,length(time_since_fix));
% sac_hist_predrate = zeros(9,length(time_since_fix));
stim_mod_predrate = zeros(9,length(tr_inds));
sac_hist_predrate = zeros(9,length(tr_inds));
clear net_mod* phase_mod* amp_mod* phasenet_mod*cl sac_hist* stim_mod*
for c = 1:n_used_cells
    % c = 6;
    fprintf('Cell %d of %d\n',c,9);
    
    un_spk_cnts = length(unique(spikes_binned(tr_inds,c))) - 1;
    spkbs = [];
    for i = 1:un_spk_cnts
        curset = find(spikes_binned(tr_inds,cellids(c)) == i);
        spkbs = [spkbs; repmat(curset,i,1)];
    end
    spkbs = sort(spkbs);
    Robs = spikes_binned(tr_inds,cellids(c));
    Robs_xv = spikes_binned(xv_inds,cellids(c));
    
    %%
%     disp('Fitting sac hist model')
%     %first fit sachist model
%     Xmat = Tmat;
%     W0 = zeros(size(Xmat,2),1);
%     NT = length(tr_inds);
%     
%     %null model
%     avg_rate = repmat(mean(Robs),length(Robs),1);
%     null_LL(c) = sum(-Robs.*log(avg_rate)+avg_rate)/sum(Robs);
%     avg_rate_xv = repmat(mean(Robs),length(Robs_xv),1);
%     null_xvLL(c) = sum(-Robs_xv.*log(avg_rate_xv)+avg_rate_xv)/sum(Robs_xv);
%     
%     %initialize parameters
%     silent = 1;
%     lamrange = [];
%     lamrange2 = [slope2_lambda 1 length(tax)];
%     nlxrange = [tax];
%     Pcon = [];
%     Pmono = [];
%     hold_const = [];
%     NLtype = 0;
%     llist = [];
%     
%     [fitp,grad] = GLMsolve_jmm( Xmat(tr_inds,:), spkbs, W0, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const, NLtype, 1e-6,nlxrange);
%     sac_hist_k(c,:) = fitp.k(1:end-1);
%     sac_hist_const(c) = fitp.k(end);
%     sac_hist_LP(c) = fitp.LP;
%     
%     cur_genfun = Tmat(tr_inds,:)*sac_hist_k(c,:)' + sac_hist_const(c);
%     cur_predrate = log(1+exp(cur_genfun));
%     sac_hist_LL(c) = sum(-Robs.*log(cur_predrate) + cur_predrate)/sum(Robs);
% %     cur_genfun_xv = Tmat(xv_inds,:)*sac_hist_k(c,:)' + sac_hist_const(c);
% %     cur_predrate_xv = log(1+exp(cur_genfun_xv));
% %     sac_hist_xvLL(c) = sum(-Robs_xv.*log(cur_predrate_xv) + cur_predrate_xv)/sum(Robs_xv);
% %     cur_genfun = Tmat*sac_hist_k(c,:)' + sac_hist_const(c);
%     sac_hist_predrate(c,:) = log(1+exp(cur_genfun));
    
    %%
%     disp('Fitting total stim model')
%     X_resh = reshape(X,size(X,1),SDIM^2);
%     cur_mask1 = get_pgabor_mask(gabor_params_fin(c,1:6),0,[SDIM SDIM]);
%     cur_mask2 = get_pgabor_mask(gabor_params_fin(c,1:6),pi/2,[SDIM SDIM]);
%     mask1_out = X_resh*cur_mask1(:);
%     mask2_out = X_resh*cur_mask2(:);
%     lin_out = gabor_params_fin(c,8)*mask1_out(all_model_fixids) + gabor_params_fin(c,9)*mask2_out(all_model_fixids);
%     energy_out = gabor_params_fin(c,7)*sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
%     total_out = lin_out + energy_out;
%     total_out = zscore(total_out);
%     
%     lamrange2 = [slope2_lambda 1 length(tax);
%         slope2_lambda length(tax)+1 length(tax)*2];
%     nlxrange = [tax tax];
%     Xmat = Tmat;
%     Xmat = [Xmat bsxfun(@times,Tmat,total_out)];
%     
%     W0 = zeros(size(Xmat,2),1);
%     llist = [];
%     [fitp,grad] = GLMsolve_jmm( Xmat(tr_inds,:), spkbs, W0, 1, lamrange, lamrange2,Pcon, Pmono, llist, [], NLtype,1e-5,nlxrange);
%     stim_mod_const(c) = fitp.k(end);
%     stim_mod_sachist(c,:) = fitp.k(1:length(tax));
%     stim_mod_tot(c,:) = fitp.k(length(tax)+1:length(tax)*2);
%     
%     cur_genfun = Xmat(tr_inds,:)*fitp.k(1:end-1) + stim_mod_const(c);
%     cur_predrate = log(1+exp(cur_genfun));
%     stim_mod_LL(c) = sum(-Robs.*log(cur_predrate) + cur_predrate)/sum(Robs);
% %     cur_genfun_xv = Xmat(xv_inds,:)*fitp.k(1:end-1) + stim_mod_const(c);
% %     cur_predrate_xv = log(1+exp(cur_genfun_xv));
% %     stim_mod_predrate_xv(c,:) = cur_predrate_xv;
% %     stim_mod_xvLL(c) = sum(-Robs_xv.*log(cur_predrate_xv) + cur_predrate_xv)/sum(Robs_xv);
% %     cur_genfun = Xmat*fitp.k(1:end-1) + stim_mod_const(c);
%     
%     stim_mod_predrate(c,:) = log(1+exp(cur_genfun));
    
    %%
%     disp('Fitting network stim model')
%     X_resh = reshape(X,size(X,1),SDIM^2);
%     cur_mask1 = get_pgabor_mask(gabor_params_fin(c,1:6),0,[SDIM SDIM]);
%     cur_mask2 = get_pgabor_mask(gabor_params_fin(c,1:6),pi/2,[SDIM SDIM]);
%     mask1_out = X_resh*cur_mask1(:);
%     mask2_out = X_resh*cur_mask2(:);
%     lin_out = gabor_params_fin(c,8)*mask1_out(all_model_fixids) + gabor_params_fin(c,9)*mask2_out(all_model_fixids);
%     energy_out = gabor_params_fin(c,7)*sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
%     total_out = lin_out + energy_out;
%     total_out = zscore(total_out);
%     
%     %     other_units = setdiff(1:24,cellids(c));
%     %     norm_spikes_binned = zscore(spikes_binned);
%     %     network_rate = mean(norm_spikes_binned(:,other_units),2);
%     cur_probe = Blocks{1}.suprobes(cellids(c));
%     if cur_probe == 24
%         neighbs = 23;
%     else
%         neighbs = [cur_probe-1 cur_probe+1];
%     end
%     network_rate = zeros(size(total_out));
%     network_rate(all_used_inds) = mean(all_lfp_hf(:,neighbs),2);
%     
%     lamrange2 = [slope2_lambda 1 length(tax);
%         slope2_lambda length(tax)+1 length(tax)*2;
%         slope2_lambda length(tax)*2+1 length(tax)*3];
%     nlxrange = [tax tax tax];
%     Xmat = Tmat;
%     Xmat = [Xmat bsxfun(@times,Tmat,total_out)];
%     Xmat = [Xmat bsxfun(@times,Tmat,network_rate)];
%     
%     W0 = zeros(size(Xmat,2),1);
%     llist = [];
%     [fitp,grad] = GLMsolve_jmm( Xmat(tr_inds,:), spkbs, W0, 1, lamrange, lamrange2,Pcon, Pmono, llist, [], NLtype,1e-5,nlxrange);
%     net_mod_const(c) = fitp.k(end);
%     net_mod_sachist(c,:) = fitp.k(1:length(tax));
%     net_mod_stim(c,:) = fitp.k(length(tax)+1:length(tax)*2);
%     net_mod_net(c,:) = fitp.k(length(tax)*2+1:length(tax)*3);
%     
%     cur_genfun = Xmat(tr_inds,:)*fitp.k(1:end-1) + net_mod_const(c);
%     cur_predrate = log(1+exp(cur_genfun));
%     net_mod_LL(c) = sum(-Robs.*log(cur_predrate) + cur_predrate)/sum(Robs);
%     cur_genfun_xv = Xmat(xv_inds,:)*fitp.k(1:end-1) + net_mod_const(c);
%     cur_predrate_xv = log(1+exp(cur_genfun_xv));
%     net_mod_predrate_xv(c,:) = cur_predrate_xv;
%     net_mod_xvLL(c) = sum(-Robs_xv.*log(cur_predrate_xv) + cur_predrate_xv)/sum(Robs_xv);
%     cur_genfun = Xmat*fitp.k(1:end-1) + net_mod_const(c);
%     
%     net_mod_predrate(c,:) = log(1+exp(cur_genfun));
    
    %%
    %     disp('Fitting amp stim model')
    %     X_resh = reshape(X,size(X,1),SDIM^2);
    %     cur_mask1 = get_pgabor_mask(gabor_params_fin(c,1:6),0,[SDIM SDIM]);
    %     cur_mask2 = get_pgabor_mask(gabor_params_fin(c,1:6),pi/2,[SDIM SDIM]);
    %     mask1_out = X_resh*cur_mask1(:);
    %     mask2_out = X_resh*cur_mask2(:);
    %     lin_out = gabor_params_fin(c,8)*mask1_out(all_model_fixids) + gabor_params_fin(c,9)*mask2_out(all_model_fixids);
    %     energy_out = gabor_params_fin(c,7)*sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
    %     total_out = lin_out + energy_out;
    %     total_out = zscore(total_out);
    %
    %     cur_probe = Blocks{1}.suprobes(cellids(c));
    %     used_freqs = [1 2 3 4 5];
    %     cur_lfp_amp = squeeze(all_lfp_amp(cur_probe,tr_inds_lfp,used_freqs));
    %     cur_lfp_amp = zscore(cur_lfp_amp);
    %
    %     silent = 0;
    %     lamrange = [];
    %     Pcon = [];
    %     Pmono = [];
    %     hold_const = [];
    %     NLtype = 0;
    %     llist = [];
    %     lamrange2 = [slope2_lambda 1 length(tax);
    %         slope2_lambda length(tax)+1 length(tax)*2];
    %     nlxrange = [tax tax];
    %     Xmat = Tmat(tr_inds,:);
    %     Xmat = [Xmat bsxfun(@times,Tmat(tr_inds,:),total_out(tr_inds))];
    %     for i = 1:length(used_freqs)
    %         Xmat = [Xmat bsxfun(@times,Tmat(tr_inds,:),cur_lfp_amp(:,i))];
    %         nlxrange = [nlxrange tax];
    %         lamrange2 = [lamrange2; slope2_lambda length(tax)*(i+1)+1 length(tax)*(i+2)];
    %     end
    %
    %     W0 = zeros(size(Xmat,2),1);
    %     NT = length(tr_inds);
    %
    %     llist = [];
    %     [fitp,grad] = GLMsolve_jmm( Xmat, spkbs, W0, silent, lamrange, lamrange2,Pcon, Pmono, llist, [], NLtype,1e-5,nlxrange);
    %     amp_mod_const(c) = fitp.k(end);
    %     amp_mod_sachist(c,:) = fitp.k(1:length(tax));
    %     amp_mod_stim(c,:) = fitp.k(length(tax)+1:length(tax)*2);
    %     amp_mod_amp(c,:,:) = reshape(fitp.k(length(tax)*2+1:end-1),length(tax),length(used_freqs));
    %
    %     cur_genfun = Xmat*fitp.k(1:end-1) + phase_mod_const(c);
    %     cur_predrate = log(1+exp(cur_genfun));
    %     amp_mod_LL(c) = sum(-Robs.*log(cur_predrate) + cur_predrate)/sum(Robs);
    
    %%
    disp('Fitting phase stim model')
    X_resh = reshape(X,size(X,1),SDIM^2);
    cur_mask1 = get_pgabor_mask(gabor_params_fin(c,1:6),0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(gabor_params_fin(c,1:6),pi/2,[SDIM SDIM]);
    mask1_out = X_resh*cur_mask1(:);
    mask2_out = X_resh*cur_mask2(:);
    lin_out = gabor_params_fin(c,8)*mask1_out(all_model_fixids) + gabor_params_fin(c,9)*mask2_out(all_model_fixids);
    energy_out = gabor_params_fin(c,7)*sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
    total_out = lin_out + energy_out;
    total_out = zscore(total_out);
    
    cur_probe = Blocks{1}.suprobes(cellids(c));
%     used_freqs = [5];
%     cur_phase_range = zeros(length(used_freqs),length(tr_inds));
%     for ww = 1:length(used_freqs)
%         
%         %         cur_phase = all_lfp_phase(cur_probe,tr_inds_lfp,used_freqs(ww));
%         cur_phase = all_lfp_phase(tr_inds_lfp,cur_probe);
%         %         cur_amp = all_lfp_amp(cur_probe,tr_inds_lfp,used_freqs(ww));
%         %         cur_avg_phases = zeros(size(cur_phase));
%         %         for i = 1:length(tax)
%         %             curset = find(time_since_fix(tr_inds) >= taxb(i) & time_since_fix(tr_inds) < taxb(i+1));
%         % %             cur_avg_phases(curset) = sac_trg_mphase(c,used_freqs(ww),i);
%         %             cur_avg_phases(curset) = sac_trg_mphase(c,i);
%         %         end
%         cur_avg_phases = circ_mean(cur_phase,Robs);
%         cur_phase_range(ww,:) = cos(cur_phase - cur_avg_phases);
%     end
        
cur_phase = all_lfp_phase_g(tr_inds_lfp,cur_probe);
cur_avg_phases_g = circ_mean(cur_phase,Robs);
cur_phase_range_g = cos(cur_phase - cur_avg_phases_g);
% cur_phase = all_lfp_phase_t(tr_inds_lfp,cur_probe);
% cur_avg_phases_t = circ_mean(cur_phase,Robs);
% cur_phase_range_t = cos(cur_phase - cur_avg_phases_t);
    
    silent = 1;
    lamrange = [];
    Pcon = [];
    Pmono = [];
    hold_const = [];
    NLtype = 0;
    llist = [];
    lamrange2 = [slope2_lambda 1 length(tax);
        slope2_lambda length(tax)+1 length(tax)*2];
    nlxrange = [tax tax];
    Xmat = Tmat(tr_inds,:);
    Xmat = [Xmat bsxfun(@times,Tmat(tr_inds,:),total_out(tr_inds))];
%     for i = 1:length(used_freqs)
%         Xmat = [Xmat bsxfun(@times,Tmat(tr_inds,:),cur_phase_range(i,:)')];
%         nlxrange = [nlxrange tax];
%         lamrange2 = [lamrange2; slope2_lambda length(tax)*(i+1)+1 length(tax)*(i+2)];
%     end
Xmat = [Xmat bsxfun(@times,Tmat(tr_inds,:),cur_phase_range_g)];
% Xmat = [Xmat bsxfun(@times,Tmat(tr_inds,:),cur_phase_range_t)];
nlxrange = [nlxrange tax];
% nlxrange = [nlxrange tax];
lamrange2 = [lamrange2; slope2_lambda length(tax)*2+1 length(tax)*3];
% lamrange2 = [lamrange2; slope2_lambda length(tax)*3+1 length(tax)*4];

    W0 = zeros(size(Xmat,2),1);
    NT = length(tr_inds);
    
    llist = [];
    [fitp,grad] = GLMsolve_jmm( Xmat, spkbs, W0, silent, lamrange, lamrange2,Pcon, Pmono, llist, [], NLtype,1e-5,nlxrange);
    phase_mod_const(c) = fitp.k(end);
    phase_mod_sachist(c,:) = fitp.k(1:length(tax));
    phase_mod_stim(c,:) = fitp.k(length(tax)+1:length(tax)*2);
    phase_mod_phase_g(c,:) = fitp.k(length(tax)*2+1:length(tax)*3);
%     phase_mod_phase_t(c,:) = fitp.k(length(tax)*3+1:length(tax)*4);
    
    cur_genfun = Xmat*fitp.k(1:end-1) + phase_mod_const(c);
    cur_predrate = log(1+exp(cur_genfun));
    phase_mod_LL(c) = sum(-Robs.*log(cur_predrate) + cur_predrate)/sum(Robs);
    
%     cur_phase_range = zeros(length(used_freqs),length(xv_inds));
%     for ww = 1:length(used_freqs)
%         cur_phase = all_lfp_phase(cur_probe,xv_inds,used_freqs(ww));
%         cur_avg_phases = zeros(size(cur_phase));
%         for i = 1:length(tax)
%             curset = find(time_since_fix(xv_inds) >= taxb(i) & time_since_fix(xv_inds) < taxb(i+1));
%             cur_avg_phases(curset) = sac_trg_mphase(c,used_freqs(ww),i);
%         end
%         cur_phase_range(ww,:) = cos(cur_phase - cur_avg_phases);
%     end
%     Xmat = Tmat(xv_inds,:);
%     Xmat = [Xmat bsxfun(@times,Tmat(xv_inds,:),total_out(xv_inds))];
%     for i = 1:length(used_freqs)
%         Xmat = [Xmat bsxfun(@times,Tmat(xv_inds,:),cur_phase_range(i,:)')];
%     end
%     cur_genfun_xv = Xmat*fitp.k(1:end-1) + phase_mod_const(c);
%     cur_predrate_xv = log(1+exp(cur_genfun_xv));
%     phase_mod_predrate_xv(c,:) = cur_predrate_xv;
%     phase_mod_xvLL(c) = sum(-Robs_xv.*log(cur_predrate_xv) + cur_predrate_xv)/sum(Robs_xv);
    
    %%
%     disp('Fitting phase net stim model')
%     X_resh = reshape(X,size(X,1),SDIM^2);
%     cur_mask1 = get_pgabor_mask(gabor_params_fin(c,1:6),0,[SDIM SDIM]);
%     cur_mask2 = get_pgabor_mask(gabor_params_fin(c,1:6),pi/2,[SDIM SDIM]);
%     mask1_out = X_resh*cur_mask1(:);
%     mask2_out = X_resh*cur_mask2(:);
%     lin_out = gabor_params_fin(c,8)*mask1_out(all_model_fixids) + gabor_params_fin(c,9)*mask2_out(all_model_fixids);
%     energy_out = gabor_params_fin(c,7)*sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
%     total_out = lin_out + energy_out;
%     total_out = zscore(total_out);
%     
%     cur_probe = Blocks{1}.suprobes(cellids(c));
%     used_freqs = [5];
%     cur_phase_range = zeros(length(used_freqs),length(tr_inds));
%     for ww = 1:length(used_freqs)
%         
%         cur_phase = all_lfp_phase(cur_probe,tr_inds_lfp,used_freqs(ww));
%         cur_amp = all_lfp_amp(cur_probe,tr_inds_lfp,used_freqs(ww));
%         cur_avg_phases = zeros(size(cur_phase));
%         for i = 1:length(tax)
%             curset = find(time_since_fix(tr_inds) >= taxb(i) & time_since_fix(tr_inds) < taxb(i+1));
%             cur_avg_phases(curset) = sac_trg_mphase(c,used_freqs(ww),i);
%         end
%         cur_phase_range(ww,:) = cos(cur_phase - cur_avg_phases);
%     end
%     
%     %     other_units = setdiff(1:24,cellids(c));
%     %     norm_spikes_binned = zscore(spikes_binned);
%     %     network_rate = mean(norm_spikes_binned(:,other_units),2);
%     cur_probe = Blocks{1}.suprobes(cellids(c));
%     if cur_probe == 24
%         neighbs = 23;
%     else
%         neighbs = [cur_probe-1 cur_probe+1];
%     end
%     network_rate = zeros(size(total_out));
%     network_rate(all_used_inds) = mean(all_lfp_hf(:,neighbs),2);
%     
%     silent = 1;
%     lamrange = [];
%     Pcon = [];
%     Pmono = [];
%     hold_const = [];
%     NLtype = 0;
%     llist = [];
%     lamrange2 = [slope2_lambda 1 length(tax);
%         slope2_lambda length(tax)+1 length(tax)*2];
%     nlxrange = [tax tax];
%     Xmat = Tmat(tr_inds,:);
%     Xmat = [Xmat bsxfun(@times,Tmat(tr_inds,:),total_out(tr_inds))];
%     for i = 1:length(used_freqs)
%         Xmat = [Xmat bsxfun(@times,Tmat(tr_inds,:),cur_phase_range(i,:)')];
%         nlxrange = [nlxrange tax];
%         lamrange2 = [lamrange2; slope2_lambda length(tax)*(i+1)+1 length(tax)*(i+2)];
%     end
%     Xmat = [Xmat bsxfun(@times,Tmat(tr_inds,:),network_rate(tr_inds))];
%     nlxrange = [nlxrange tax];
%     lamrange2 = [lamrange2; slope2_lambda length(tax)*(length(used_freqs)+2)+1 length(tax)*(length(used_freqs)+3)];
%     
%     W0 = zeros(size(Xmat,2),1);
%     NT = length(tr_inds);
%     
%     llist = [];
%     [fitp,grad] = GLMsolve_jmm( Xmat, spkbs, W0, silent, lamrange, lamrange2,Pcon, Pmono, llist, [], NLtype,1e-5,nlxrange);
%     phasenet_mod_const(c) = fitp.k(end);
%     phasenet_mod_sachist(c,:) = fitp.k(1:length(tax));
%     phasenet_mod_stim(c,:) = fitp.k(length(tax)+1:length(tax)*2);
%     %         phase_mod_phase(c,:,:) = reshape(fitp.k(length(tax)*2+1:end-1),length(tax),length(used_freqs));
%     phasenet_mod_phase(c,:) = fitp.k(length(tax)*2+1:length(tax)*3);
%     phasenet_mod_net(c,:) = fitp.k(length(tax)*3+1:end-1);
%     
%     cur_genfun = Xmat*fitp.k(1:end-1) + phasenet_mod_const(c);
%     cur_predrate = log(1+exp(cur_genfun));
%     phasenet_mod_LL(c) = sum(-Robs.*log(cur_predrate) + cur_predrate)/sum(Robs);
%     
%     cur_phase_range = zeros(length(used_freqs),length(xv_inds));
%     for ww = 1:length(used_freqs)
%         cur_phase = all_lfp_phase(cur_probe,xv_inds,used_freqs(ww));
%         cur_avg_phases = zeros(size(cur_phase));
%         for i = 1:length(tax)
%             curset = find(time_since_fix(xv_inds) >= taxb(i) & time_since_fix(xv_inds) < taxb(i+1));
%             cur_avg_phases(curset) = sac_trg_mphase(c,used_freqs(ww),i);
%         end
%         cur_phase_range(ww,:) = cos(cur_phase - cur_avg_phases);
%     end
%     Xmat = Tmat(xv_inds,:);
%     Xmat = [Xmat bsxfun(@times,Tmat(xv_inds,:),total_out(xv_inds))];
%     for i = 1:length(used_freqs)
%         Xmat = [Xmat bsxfun(@times,Tmat(xv_inds,:),cur_phase_range(i,:)')];
%     end
%     Xmat = [Xmat bsxfun(@times,Tmat(xv_inds,:),network_rate(xv_inds))];
%     cur_genfun_xv = Xmat*fitp.k(1:end-1) + phasenet_mod_const(c);
%     cur_predrate_xv = log(1+exp(cur_genfun_xv));
%     phasenet_mod_predrate_xv(c,:) = cur_predrate_xv;
%     phasenet_mod_xvLL(c) = sum(-Robs_xv.*log(cur_predrate_xv) + cur_predrate_xv)/sum(Robs_xv);
    
end

%%
for c = 1:9
    
    a = corrcoef(sac_hist_predrate(c,xv_inds),spikes_binned(xv_inds,c));
    sac_hist_corr(c) = a(2,1);
    a = corrcoef(stim_mod_predrate_xv(c,:),spikes_binned(xv_inds,c));
    stim_hist_corr(c) = a(2,1);
    a = corrcoef(net_mod_predrate_xv(c,:),spikes_binned(xv_inds,c));
    net_corr(c) = a(2,1);
    a = corrcoef(phase_mod_predrate_xv(c,:),spikes_binned(xv_inds,c));
    phase_corr(c) = a(2,1);
    a = corrcoef(phasenet_mod_predrate_xv(c,:),spikes_binned(xv_inds,c));
    phasenet_corr(c) = a(2,1);
    
    sm_rate = smooth(spikes_binned(xv_inds,c),20);
    a = corrcoef(sm_rate,spikes_binned(xv_inds,c));
    expected_corr(c) = a(2,1);
end

%%
% for c = 1:9;
%     fprintf('Cell %d of %d\n',c,9);
% un_spk_cnts = length(unique(spikes_binned(tr_inds,c))) - 1;
% spkbs = [];
% for i = 1:un_spk_cnts
%     curset = find(spikes_binned(tr_inds,c) == i);
%     spkbs = [spkbs; repmat(curset,i,1)];
% end
% spkbs = sort(spkbs);
% Robs = spikes_binned(tr_inds,c);
% Robs_xv = spikes_binned(xv_inds,c);
%
% init_params = gabor_params_fin(c,:);
% cur_gain = Tmat*stim_mod_tot(c,:)';
% cur_offset = Tmat*stim_mod_sachist(c,:)';
% LB = [-5 -5 0 6 2 0.5 -Inf -Inf -Inf -Inf];
% UB = [10 10 pi SDIM(1)/2 SDIM(1)/6 2 Inf Inf Inf Inf];
% hold_const = [1 1 1 1 1 1 1 1 1 0];
% [gabor_params,LL] = fit_temporal_gabor_params(init_params,X_resh,cur_gain(tr_inds),cur_offset(tr_inds),...
%     Robs,all_model_fixids(tr_inds),[SDIM SDIM],hold_const,LB,UB);
%
% % hold_const = [1 1 1 1 1 1 0 0 0 0];
% % [gabor_params,LL] = fit_temporal_gabor_params(gabor_params,X_resh,cur_gain(tr_inds),cur_offset(tr_inds),...
% %     Robs,all_model_fixids(tr_inds),[SDIM SDIM],hold_const,LB,UB);
%
% hold_const = [1 1 0 0 0 0 0 0 0 0];
% [new_gabor_params(c,:),LL] = fit_temporal_gabor_params(gabor_params,X_resh,cur_gain(tr_inds),cur_offset(tr_inds),...
%     Robs,all_model_fixids(tr_inds),[SDIM SDIM],hold_const,LB,UB);
%
% cur_mask1 = get_pgabor_mask(new_gabor_params(c,1:6),0,[SDIM SDIM]);
% cur_mask2 = get_pgabor_mask(new_gabor_params(c,1:6),pi/2,[SDIM SDIM]);
% mask1_out = X_resh*cur_mask1(:);
% mask2_out = X_resh*cur_mask2(:);
% lin_out = new_gabor_params(c,8)*mask1_out(all_model_fixids) + new_gabor_params(c,9)*mask2_out(all_model_fixids);
% energy_out = new_gabor_params(c,7)*sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
% total_out = lin_out + energy_out;
% total_out = zscore(total_out);
% lamrange2 = [slope2_lambda 1 length(tax);
%     slope2_lambda length(tax)+1 length(tax)*2];
% nlxrange = [tax tax];
% Xmat = Tmat;
% Xmat = [Xmat bsxfun(@times,Tmat,total_out)];
% llist = [l1_pen (length(tax)+1):size(Xmat,2)];
% [fitp,grad] = GLMsolve_jmm( Xmat(tr_inds,:), spkbs, W0, 1, lamrange, lamrange2,Pcon, Pmono, llist, [], NLtype,1e-5,nlxrange);
% rstim_mod_const(c) = fitp.k(end);
% rstim_mod_sachist(c,:) = fitp.k(1:length(tax));
% rstim_mod_tot(c,:) = fitp.k(length(tax)+1:length(tax)*2);
%
% cur_genfun = Xmat(tr_inds,:)*fitp.k(1:end-1) + rstim_mod_const(c);
% cur_predrate = log(1+exp(cur_genfun));
% rstim_mod_LL(c) = sum(-Robs.*log(cur_predrate) + cur_predrate)/sum(Robs);
% cur_genfun_xv = Xmat(xv_inds,:)*fitp.k(1:end-1) + stim_mod_const(c);
% cur_predrate_xv = log(1+exp(cur_genfun_xv));
% rstim_mod_xvLL(c) = sum(-Robs_xv.*log(cur_predrate_xv) + cur_predrate_xv)/sum(Robs_xv);
% cur_genfun = Xmat*fitp.k(1:end-1) + rstim_mod_const(c);
% rstim_mod_predrate(c,:) = log(1+exp(cur_genfun));
% end
%%
stim_pred_diff = zeros(n_used_cells,length(tax));
obs_diff = zeros(n_used_cells,length(tax));
taxb = [0 tax];
for i = 1:length(tax)
    curset = find(time_since_fix >= taxb(i) & time_since_fix < taxb(i+1));
    for c = 1:n_used_cells
        %         stim_pred_diff(c,i) =mean(abs(stim_mod_predrate(c,curset)-sac_hist_predrate(c,curset)));
        stim_pred_diff(c,i) =sqrt(mean((stim_mod_predrate(c,curset)-sac_hist_predrate(c,curset)).^2));
        %         obs_diff(c,i) = mean(abs(spikes_binned(curset,c)'-sac_hist_predrate(c,curset)));
        obs_diff(c,i) = sqrt(mean((spikes_binned(curset,c)'-sac_hist_predrate(c,curset)).^2));
    end
    n_occ(i) = length(curset);
end
stim_pred_diff = stim_pred_diff/dt;
obs_diff = obs_diff/dt;

%%
cd /Users/James/James_scripts/bruce/modelfits/
save gabor_tempmodfits_phasemods stim_mod* sac_* null* tax stim_mod_predrate all_model_* fano* time_since_fix phase_mod*

%%
clf
cur_cell = 6;
subplot(3,1,1)
plot(tax,sac_trg_rate(cur_cell,:))
xlim([0 0.5])
subplot(3,1,2)
plot(tax,sac_hist_k(cur_cell,:),'r')
hold on
plot(tax,stim_mod_sachist(cur_cell,:),'k')
title('Constant')
xlim([0 0.5])
subplot(3,1,3)
plot(tax,stim_mod_tot(cur_cell,:),'k')
hold on
% plot(tax,sep_mod_lin(cur_cell,:),'b')
% plot(tax,sep_mod_en(cur_cell,:),'r')
hold on
xlim([0 0.5])

%%
used = 1:n_used_cells;
sac_LL_imp = (null_LL(used) - sac_hist_LL(used))/log(2);
sac_xvLL_imp = (null_xvLL(used) - sac_hist_xvLL(used))/log(2);

stim_LL_imp = (null_LL(used) - stim_mod_LL(used))/log(2);
stim_xvLL_imp = (null_xvLL(used) - stim_mod_xvLL(used))/log(2);

% sep_LL_imp = null_LL(used) - sep_mod_LL(used);
% sep_xvLL_imp = null_xvLL(used) - sep_mod_xvLL(used);

%%
% taxold = tax;
% max_tsf = 0.75; nbins = 25;
% used_tsf = time_since_fix;
% used_tsf(used_tsf > max_tsf) = [];
% tax = prctile(used_tsf,linspace(100/nbins,100-100/nbins,nbins));
stim_mod_predrate_hz = stim_mod_predrate/dt;
taxb = [0 tax];
bin_widths = taxb(2:end)-taxb(1:end-1);
rel_spike_binned = nan(n_used_cells,length(used_fixs),length(tax));
combined_spike_binned = nan(n_used_cells,length(used_fixs),length(tax));
sac_spike_binned = nan(n_used_cells,length(used_fixs),length(tax));
for blockid = 1:4;
    fprintf('Block %d of %d\n',blockid,4);
    cur_set = find(blockids(used_fixs) == blockid);
    cur_model_set = find(all_model_blockids==blockid);
    for i = 1:length(cur_set)
        start_T = all_fix_start_times(used_fixs(cur_set(i)));
        end_T = all_fix_stop_times(used_fixs(cur_set(i)));
        cur_tedges = start_T + taxb;
        model_t_interp = round(interp1(all_model_time_axis(cur_model_set),1:length(cur_model_set),start_T+tax));
        cur_tcents = 0.5*cur_tedges(1:end-1)+0.5*cur_tedges(2:end);
        cur_bad = find(cur_tedges > end_T);
        cur_bad2 = find(cur_tcents > end_T | isnan(model_t_interp));
        for c = 1:n_used_cells
            if c <= 9
                temp = histc(Blocks{blockid}.spktimes{cellids(c)},cur_tedges);
            else
                temp = histc(Blocks{blockid}.mutimes{muaids(c-9)},cur_tedges);
            end
            temp(cur_bad) = nan;
            %             rel_spike_binned(c,cur_set(i),:) = temp(1:end-1);
            rel_spike_binned(c,cur_set(i),:) = temp(1:end-1);
            
            model_t_interp(cur_bad2) = 1;
            combined_pred_cnts = stim_mod_predrate(c,cur_model_set(model_t_interp));
            combined_pred_cnts(cur_bad2) = nan;
            combined_spike_binned(c,cur_set(i),:) = combined_pred_cnts.*bin_widths;
            %             combined_spike_binned(c,cur_set(i),:) = combined_pred_cnts;
            
            %             sac_pred_cnts = sac_pred_r(cur_model_set(model_t_interp),c);
            %             sac_pred_cnts(cur_bad2) = nan;
            %             sac_spike_binned(c,cur_set(i),:) = sac_pred_cnts.*bin_widths';
        end
    end
end
% combined_residual = rel_spike_binned - combined_spike_binned;
% combined_resvar = squeeze(nanvar(combined_residual,[],2));
% orivar = squeeze(nanvar(rel_spike_binned,[],2));
% combined_mean = squeeze(nanmean(combined_spike_binned,2));
% ori_mean = squeeze(nanmean(rel_spike_binned,2));
% combined_fano = combined_resvar./ori_mean;
% ori_fano = orivar./ori_mean;

cur_corr = zeros(n_used_cells,length(tax));
for i = 1:length(tax)
    uset = find(~isnan(combined_spike_binned(1,:,i)));
    for c = 1:n_used_cells
        [a,b] = corrcoef(rel_spike_binned(c,uset,i),combined_spike_binned(c,uset,i));
        cur_corr(c,i) = a(2,1);
        cur_p(c,i) = b(2,1);
    end
end

%%
for cur_cell = 1:n_used_cells;
    % cur_cell = 8
    % clf
    [~,peakloc] = max(sac_trg_rate(cur_cell,:));
    
    subplot(3,2,1)
    plot(tax,sac_trg_rate(cur_cell,:))
    axis tight
    xlim([0 0.5])
    yl = ylim();
    line(tax([peakloc peakloc]),yl,'color','k','linestyle','--')
    ylabel('Avreage rate (Hz)','fontsize',14)
    title('Fixation triggered average rate','fontsize',14)
    
    subplot(3,2,3)
    plot(tax,sac_hist_k(cur_cell,:)+sac_hist_const(cur_cell),'r')
    hold on
    plot(tax,stim_mod_sachist(cur_cell,:)+stim_mod_const(cur_cell),'k')
    title('Offset','fontsize',14)
    axis tight
    xlim([0 0.5])
    ylabel('Offset term','fontsize',14)
    legend('Without Stim','With Stim')
    yl = ylim();
    line(tax([peakloc peakloc]),yl,'color','k','linestyle','--')
    
    subplot(3,2,5)
    plot(tax,stim_mod_tot(cur_cell,:),'k')
    title('Stimulus tuning gain','fontsize',14)
    axis tight
    xlim([0 0.5])
    ylabel('Stimulus filter gain (1/z)','fontsiz',14)
    xlabel('Time since fixation onset (s)','fontsize',14)
    yl = ylim();
    line(tax([peakloc peakloc]),yl,'color','k','linestyle','--')
    
    subplot(3,2,2)
    % plot(tax,obs_diff(cur_cell,:))
    hold on
    plot(tax,stim_pred_diff(cur_cell,:),'b')
    % legend('Observed','Model')
    axis tight
    xlim([0 0.5])
    ylabel('RMS rate difference (Hz)','fontsize',14)
    title('Stim-Sac Model Difference','fontsize',14)
    yl = ylim();
    line(tax([peakloc peakloc]),yl,'color','k','linestyle','--')
    
    subplot(3,2,4)
    plot(tax,fano_fac(cur_cell,:),'ko-')
    hold on
    plot(tax,smooth(fano_fac(cur_cell,:),3),'r')
    % plot(tax,combined_fano(cur_cell,:),'r')
    % plot(tax,ori_fano(cur_cell,:),'k')
    axis tight
    xlim([0 0.5])
    ylabel('Fano Factor','fontsize',14)
    title('Fano Factor','fontsize',14)
    yl = ylim();
    line(tax([peakloc peakloc]),yl,'color','k','linestyle','--')
    
    subplot(3,2,6)
    plot(tax,cur_corr(cur_cell,:).^2,'ko-')
    hold on
    plot(tax,smooth(cur_corr(cur_cell,:).^2,3),'r')
    title('Stimulus prediction R2','fontsize',14)
    hold on
    % plot(tax,smooth(cur_corr(cur_cell,:).^2,3,'moving'),'r')
    axis tight
    xlim([0 0.5])
    xlabel('Time since fixation onset (s)','fontsize',14)
    ylabel('Spike count R^2','fontsize',14)
    yl = ylim();
    line(tax([peakloc peakloc]),yl,'color','k','linestyle','--')
    
    
    fillPage(gcf,'margins',[0 0 0 0],'papersize',[12 14]);
    if cur_cell < 10
        print(sprintf('Cell%d_modfit',cellids(cur_cell)),'-dpdf','-painters');
    else
        print(sprintf('MUA%d_modfit',muaids(cur_cell-9)),'-dpdf','-painters');
    end
    close
end

%%
cur_cell = 6;
xl = [0 0.4];
n_stims = 10;
stim_vals = linspace(-3,3,n_stims);
clf
clear cur_rate
cur_tax = linspace(0,0.5,200);
test_tax = tbrep(cur_tax,tax);
for i = 1:n_stims
    cur_gen = stim_mod_sachist(cur_cell,:)*test_tax';
    cur_gen = cur_gen + stim_mod_tot(cur_cell,:)*test_tax'*stim_vals(i);
    cur_gen = cur_gen + stim_mod_const(cur_cell);
    cur_rate(i,:) = log(1+exp(cur_gen));
end
norm_crate = bsxfun(@rdivide,cur_rate,max(cur_rate,[],2));
cur_rate = cur_rate/dt;

for i = 1:n_stims
    subplot(2,1,1)
    plot(cur_tax,norm_crate(i,:),'color',cmap(i,:))
    hold on
    subplot(2,1,2)
    plot(cur_tax,cur_rate(i,:),'color',cmap(i,:))
    hold on
end
% subplot(2,1,1)
% plot(tax,sac_trg_rate(cur_cell,:)/max(sac_trg_rate(cur_cell,:)),'color','k','linewidth',2)

subplot(2,1,1)
xlim(xl);
xlabel('Time since fixation onset (s)','fontsize',14)
ylabel('Relative predicted rate','fontsize',14)
subplot(2,1,2)
xlim(xl);
xlabel('Time since fixation onset (s)','fontsize',14)
ylabel('Predicted rate (Hz)','fontsize',14)