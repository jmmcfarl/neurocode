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
[bb,aa] = butter(2,120/niqf,'low');
scales = logspace(log10(5),log10(60),30);
scales = [scales 70 80 90];
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

% new_dt = .0025;
new_dt = .004;

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
new_lfps = [];
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
    clear lfp_sampsd
    for ll = 1:24
        fprintf('Processing channel %d of %d\n',ll,24);
        lfp_sampsd(:,ll) = decimate(lfp_samps(:,ll),dsf);
        temp = cwt(lfp_sampsd(:,ll),scales,'cmor1-1');
        cur_phasegram(:,:,ll) = angle(temp)';
        cur_ampgram(:,:,ll) = abs(temp)';
    end
    
    cur_set = find(new_model_blockids==blockid);
    %     interp_phasegram = interp1(lfp_timed,cur_phasegram,new_model_time_axis(cur_set));
    unwr_phasegram = unwrap(cur_phasegram);
    interp_phasegram = interp1(lfp_timed,unwr_phasegram,new_model_time_axis(cur_set));
    interp_phasegram = mod(interp_phasegram+pi,2*pi)-pi;
    
    interp_lfps = interp1(lfp_timed,lfp_sampsd,new_model_time_axis(cur_set));
    
    interp_ampgram = interp1(lfp_timed,cur_ampgram,new_model_time_axis(cur_set));
    
    new_lfps = cat(1,new_lfps,interp_lfps);
    new_phasegram = cat(1,new_phasegram,interp_phasegram);
    new_ampgram = cat(1,new_ampgram,interp_ampgram);
    %     all_t = [all_t; lfp_timed(:)];
end

new_ampgram_n = bsxfun(@rdivide,new_ampgram,nanstd(new_ampgram));

%%

cur_used_fixs = unique(new_model_fixids);
n_fixs = length(cur_used_fixs);

xv_frac = 0.2;
n_xv_trials = round(n_fixs*xv_frac);
xv_set = randperm(n_fixs);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(new_model_fixids,xv_set));
tr_inds = find(~ismember(new_model_fixids,xv_set))';


trial_start_inds = 1 + find(diff(new_model_fixids) ~= 0);
trial_stop_inds = find(diff(new_model_fixids) ~= 0);
trial_start_inds = [1 trial_start_inds];
trial_stop_inds = [trial_stop_inds length(new_model_fixids)];
%%
cur_dt = 0.008;
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
bad_inds = find(time_since_fix > max(tent_centers) | is_bad' == 1);
xv_inds(ismember(xv_inds,bad_inds)) = [];
tr_inds(ismember(tr_inds,bad_inds)) = [];

trial_inds = zeros(size(new_model_time_axis));
trial_inds(trial_start_inds) = 1;
trial_Tmat = zeros(length(new_model_time_axis),ntents);
for i = 1:ntents
    trial_Tmat(:,i) = conv(trial_inds,tbmat(i,:),'same');
end

beg_dur = round(0.15/new_dt);
late_indicator = zeros(size(new_model_time_axis));
early_indicator = zeros(size(new_model_time_axis));
for i = 1:length(trial_start_inds)
    cur_inds = (beg_dur+trial_start_inds(i)):trial_stop_inds(i);
    late_indicator(cur_inds) = 1;
    cur_inds = trial_start_inds(i):(beg_dur+trial_start_inds(i));
    early_indicator(cur_inds) = 1;
end

NT = length(new_model_time_axis);

%%
close all
NL_type = 1; %exp

reg_params.dl1_ind = 0;
reg_params.dl1_dep = 0;
reg_params.dl2_ind = 800;
reg_params.dl2_dep = 800;
reg_params.dl2_freq_ind = 10000;
reg_params.dl2_freq_dep = 10000;
reg_params.dl2_time_ind = 200;
reg_params.dl2_time_dep = 200;
reg_params.l1_ind = 0;
reg_params.l1_dep = 0;
reg_params.l1t_ind = 0;
reg_params.l1t_dep = 0;
reg_params.is_phase = 0;

% reg_params2.dl1_ind = 0;
% reg_params2.dl1_dep = 0;
% reg_params2.dl2_ind = 10000;
% reg_params2.dl2_dep = 10000;
% reg_params2.dl2_freq_ind = 1000;
% reg_params2.dl2_freq_dep = 1000;
% reg_params2.dl2_time_ind = 200;
% reg_params2.dl2_time_dep = 200;
% reg_params2.l1_ind = 0;
% reg_params2.l1_dep = 0;
% reg_params2.l1t_ind = 0;
% reg_params2.l1t_dep = 0;
% reg_params2.is_phase = 0;

silent = 1;
NT = length(new_model_time_axis);
% used_inds = 1:NT;

nbins = 30;
pax = linspace(-pi,pi,nbins+1);
doub_pax = [pax(1:end-1) 2*pi+pax(1:end-1)];

phase_set = [reshape(cos(new_phasegram),NT,length(wfreqs)*24) reshape(sin(new_phasegram),NT,length(wfreqs)*24)];
ampphase_set = [reshape(new_ampgram_n,NT,length(wfreqs)*24).*reshape(cos(new_phasegram),NT,length(wfreqs)*24) ...
    reshape(new_ampgram_n,NT,length(wfreqs)*24).*reshape(sin(new_phasegram),NT,length(wfreqs)*24)];
% phase_elec_set = [repmat(1:length(use_lfps),1,length(wfreqs)) repmat(1:length(use_lfps),1,length(wfreqs))];
% phase_elec_set = phase_elec_set(:);
phase_elec_set = ones(length(wfreqs),1)*(1:24);
phase_elec_set = [phase_elec_set(:); phase_elec_set(:)];
phase_freq_set = wfreqs'*(ones(1,24));
phase_freq_set = [phase_freq_set(:); phase_freq_set(:)];

%%
use_elecs = 1:24;
use_pset = find(ismember(phase_elec_set,use_elecs));
Xmat = [phase_set(:,use_pset)];
Xmat = bsxfun(@times,Xmat,late_indicator');
Xmat = [Xmat trial_Tmat];

aXmat = [ampphase_set(:,use_pset)];
aXmat = bsxfun(@times,aXmat,late_indicator');
aXmat = [aXmat trial_Tmat];

%%
for cc = [3]
% cc=16;    
    fprintf('Cell %d of %d\n',cc,24);
    
    Robs = spikes_binned(tr_inds,cc);
    tr_spkbns = convert_to_spikebins(Robs);
    
    X_resh = reshape(X,size(X,1),SDIM^2);
    cur_mask1 = get_pgabor_mask(gabor_params_fin(cc,1:6),0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(gabor_params_fin(cc,1:6),pi/2,[SDIM SDIM]);
    mask1_out = X_resh*cur_mask1(:);
    mask2_out = X_resh*cur_mask2(:);
    
    energy_out = gabor_params_fin(cc,7)*sqrt(mask1_out(new_model_fixids).^2 + mask2_out(new_model_fixids).^2);
    lin_out = gabor_params_fin(cc,8)*mask1_out(new_model_fixids) + gabor_params_fin(cc,9)*mask1_out(new_model_fixids);
    
    total_out = lin_out + energy_out;
        total_out = zscore(total_out);

%     stim_params = [0,0,ntents];
%     klen = size(trial_Tmat(tr_inds,:),2);
%     K0 = zeros(klen+1,1);
%     [fitp_to_ns] = fit_GLM_phase_model(trial_Tmat(tr_inds,:), Robs, K0, 1, stim_params, reg_params,1,NL_type);
%     to_ns_filt(cc,:) = fitp_to_ns.k((1):(ntents));
%     to_ns_avg(cc,:) = exp(fitp_to_ns.k(1:ntents)+fitp_to.k(end))';
% 
%     Tmat = [trial_Tmat(tr_inds,:) bsxfun(@times,trial_Tmat(tr_inds,:),total_out(tr_inds))];
%     stim_params = [0,0,ntents];
%     klen = size(Tmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_to] = fit_GLM_phase_model(Tmat, Robs, K0, 1, stim_params, reg_params,1,NL_type);
%     to_ind_filt(cc,:) = fitp_to.k((1:ntents));
%     to_dep_filt(cc,:) = fitp_to.k((ntents+1):end-1);

%     nearest_probe = all_probes(cc);
%     cur_phase = squeeze(new_phasegram(uset,:,nearest_probe));
%     Pmat = nan(length(uset),length(wfreqs)*nbins);
%     for ww = 1:length(wfreqs)
%         cur_tb = tbrep(cur_phase(:,ww),pax);
%         cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
%         cur_tb(:,end) = [];
%         cur_set = ((ww-1)*nbins+1):ww*nbins;
%         Pmat(:,cur_set) = cur_tb;
%     end
%     
%    Pmat = [bsxfun(@times,Pmat,late_indicator(uset)') trial_Tmat(uset,:)];
%     stim_params = [nbins,length(wfreqs),ntents];
%     klen = size(Pmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_dphase] = fit_GLM_phase_model(Pmat, Robs, K0,silent, stim_params, reg_params2,1,NL_type);
%     phase_pfilt(cc,:) = fitp_dphase.k(1:nbins*length(wfreqs));
%     phase_tfilt(cc,:) = fitp_dphase.k(nbins*length(wfreqs)+1:end-1);
    
    cur_Xmat = [Xmat bsxfun(@times,trial_Tmat,total_out)];
    stim_params = [length(wfreqs),length(use_elecs), ntents];
%     %     stim_params = [length(wfreqs),length(use_elecs)];
%     klen = size(cur_Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_fphase] = fit_GLM_phase_model(cur_Xmat(tr_inds,:), Robs, K0,silent, stim_params, reg_params,1,NL_type);
%     sinphase_cfilt(cc,:) = fitp_fphase.k(1:length(use_elecs)*length(wfreqs));
%     sinphase_sfilt(cc,:) = fitp_fphase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
%     sinphase_ind_tfilt(cc,:) = fitp_fphase.k((length(use_elecs)*length(wfreqs)*2+1):(length(use_elecs)*length(wfreqs)*2+ntents));
%     sinphase_dep_tfilt(cc,:) = fitp_fphase.k((length(use_elecs)*length(wfreqs)*2+ntents+1):end-1);
%     phase_ampkern(cc,:) = sqrt(sinphase_cfilt(cc,:).^2 + sinphase_sfilt(cc,:).^2);
    
%      klen = size(cur_Xmat,2);
%     K0 = zeros(klen+1,1);
%    cur_aXmat = [aXmat bsxfun(@times,trial_Tmat, total_out)];
%     [fitp_aphase] = fit_GLM_phase_model(cur_aXmat(tr_inds,:), Robs, K0,silent, stim_params, reg_params,1,NL_type);
%     ampphase_cfilt(cc,:) = fitp_aphase.k(1:length(use_elecs)*length(wfreqs));
%     ampphase_sfilt(cc,:) = fitp_aphase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
%     amphase_ind_tfilt(cc,:) = fitp_aphase.k((length(use_elecs)*length(wfreqs)*2+1):(length(use_elecs)*length(wfreqs)*2+ntents));
%     ampphase_dep_tfilt(cc,:) = fitp_aphase.k((length(use_elecs)*length(wfreqs)*2+ntents+1):end-1);
%     ampphase_ampkern(cc,:) = sqrt(ampphase_cfilt(cc,:).^2 + ampphase_sfilt(cc,:).^2);

    cur_use_elecs = all_probes(cc);
    cur_use_pset = find(ismember(phase_elec_set,cur_use_elecs));
    
%     se_aXmat = [ampphase_set(:,cur_use_pset)];
% %     se_aXmat = bsxfun(@times,se_aXmat,late_indicator');
%     se_aXmat = [se_aXmat trial_Tmat bsxfun(@times,trial_Tmat,total_out)];
%     klen = size(se_aXmat,2);
%     K0 = zeros(klen+1,1);
%     cur_stim_params = [length(wfreqs),length(cur_use_elecs), ntents];
%     [fitp_seaphase] = fit_GLM_phase_model(se_aXmat(tr_inds,:), Robs, K0,silent, cur_stim_params, reg_params,1,NL_type);
%     seampphase_cfilt(cc,:) = fitp_seaphase.k(1:length(wfreqs));
%     seampphase_sfilt(cc,:) = fitp_seaphase.k((length(wfreqs)+1):length(wfreqs)*2);
%     seampphase_ind_tfilt(cc,:) = fitp_seaphase.k((length(wfreqs)*2+1):(length(wfreqs)*2+ntents));
%     seampphase_dep_tfilt(cc,:) = fitp_seaphase.k((length(wfreqs)*2+ntents+1):end-1);
%     seampphase_ampkern(cc,:) = sqrt(seampphase_cfilt(cc,:).^2 + seampphase_sfilt(cc,:).^2);

    se_Xmat = [phase_set(:,cur_use_pset)];
%     se_Xmat = bsxfun(@times,se_Xmat,late_indicator');
    se_Xmat = [se_Xmat trial_Tmat bsxfun(@times,trial_Tmat,total_out)];
    klen = size(se_aXmat,2);
    K0 = zeros(klen+1,1);
    cur_stim_params = [length(wfreqs),length(cur_use_elecs), ntents];
    [fitp_sephase] = fit_GLM_phase_model(se_Xmat(tr_inds,:), Robs, K0,silent, cur_stim_params, reg_params,1,NL_type);
    sephase_cfilt(cc,:) = fitp_sephase.k(1:length(wfreqs));
    sephase_sfilt(cc,:) = fitp_sephase.k((length(wfreqs)+1):length(wfreqs)*2);
    sephase_ind_tfilt(cc,:) = fitp_sephase.k((length(wfreqs)*2+1):(length(wfreqs)*2+ntents));
    sephase_dep_tfilt(cc,:) = fitp_sephase.k((length(wfreqs)*2+ntents+1):end-1);
    sephase_ampkern(cc,:) = sqrt(sephase_cfilt(cc,:).^2 + sephase_sfilt(cc,:).^2);

    
%     ese_aXmat = [ampphase_set(:,cur_use_pset)];
%     ese_aXmat = bsxfun(@times,ese_aXmat,early_indicator');
%     ese_aXmat = [ese_aXmat trial_Tmat bsxfun(@times,trial_Tmat,total_out)];
%     klen = size(ese_aXmat,2);
%     K0 = zeros(klen+1,1);
%     cur_stim_params = [length(wfreqs),length(cur_use_elecs), ntents];
%     [fitp_eseaphase] = fit_GLM_phase_model(ese_aXmat(tr_inds,:), Robs, K0,silent, cur_stim_params, reg_params,1,NL_type);
%     eseampphase_cfilt(cc,:) = fitp_eseaphase.k(1:length(wfreqs));
%     eseampphase_sfilt(cc,:) = fitp_eseaphase.k((length(wfreqs)+1):length(wfreqs)*2);
%     eseamphase_ind_tfilt(cc,:) = fitp_eseaphase.k((length(wfreqs)*2+1):(length(wfreqs)*2+ntents));
%     eseampphase_dep_tfilt(cc,:) = fitp_eseaphase.k((length(wfreqs)*2+ntents+1):end-1);
%     eseampphase_ampkern(cc,:) = sqrt(eseampphase_cfilt(cc,:).^2 + eseampphase_sfilt(cc,:).^2);
%      early_pred_rate = se_aXmat*fitp_eseaphase.k(1:end-1)+fitp_eseaphase.k(end);
   
    
    xv_Robs = spikes_binned(xv_inds,cc);
    
%     Tmat = [trial_Tmat(xv_inds,:)];
%     xv_pred_rate = Tmat*fitp_to_ns.k(1:end-1) + fitp_to_ns.k(end);
%     if NL_type == 0
%         xv_pred_rate = log(1+exp(xv_pred_rate));
%     else
%         xv_pred_rate = exp(xv_pred_rate);
%     end
%     xv_to_ns_LL(cc) = -sum(xv_Robs.*log(xv_pred_rate)-xv_pred_rate)/sum(xv_Robs);
% 
%         Tmat = [Tmat bsxfun(@times,trial_Tmat(xv_inds,:),total_out(xv_inds))];
%     xv_pred_rate = Tmat*fitp_to.k(1:end-1) + fitp_to.k(end);
%     if NL_type == 0
%         xv_pred_rate = log(1+exp(xv_pred_rate));
%     else
%         xv_pred_rate = exp(xv_pred_rate);
%     end
%     xv_to_LL(cc) = -sum(xv_Robs.*log(xv_pred_rate)-xv_pred_rate)/sum(xv_Robs);
% 
%         xv_pred_rate = cur_Xmat(xv_inds,:)*fitp_fphase.k(1:end-1) + fitp_fphase.k(end);
%         if NL_type == 0
%             xv_pred_rate = log(1+exp(xv_pred_rate));
%         else
%             xv_pred_rate = exp(xv_pred_rate);
%         end
%         xv_sinphase_LL(cc) = -sum(xv_Robs.*log(xv_pred_rate)-xv_pred_rate)/sum(xv_Robs);
        
%         xv_pred_rate = cur_aXmat(xv_inds,:)*fitp_aphase.k(1:end-1) + fitp_aphase.k(end);
%         if NL_type == 0
%             xv_pred_rate = log(1+exp(xv_pred_rate));
%         else
%             xv_pred_rate = exp(xv_pred_rate);
%         end
%         xv_ampphase_LL(cc) = -sum(xv_Robs.*log(xv_pred_rate)-xv_pred_rate)/sum(xv_Robs);

%         xv_pred_rate = se_aXmat(xv_inds,:)*fitp_seaphase.k(1:end-1) + fitp_seaphase.k(end);
%         if NL_type == 0
%             xv_pred_rate = log(1+exp(xv_pred_rate));
%         else
%             xv_pred_rate = exp(xv_pred_rate);
%         end
%         xv_seampphase_LL(cc) = -sum(xv_Robs.*log(xv_pred_rate)-xv_pred_rate)/sum(xv_Robs);
% 

%         xv_pred_rate = se_Xmat(xv_inds,:)*fitp_sephase.k(1:end-1) + fitp_sephase.k(end);
%         if NL_type == 0
%             xv_pred_rate = log(1+exp(xv_pred_rate));
%         else
%             xv_pred_rate = exp(xv_pred_rate);
%         end
%         xv_sephase_LL(cc) = -sum(xv_Robs.*log(xv_pred_rate)-xv_pred_rate)/sum(xv_Robs);

%                 xv_epred_rate = ese_aXmat(xv_inds,:)*fitp_eseaphase.k(1:end-1) + fitp_eseaphase.k(end);
%         if NL_type == 0
%             xv_epred_rate = log(1+exp(xv_epred_rate));
%         else
%             xv_epred_rate = exp(xv_epred_rate);
%         end
%         xv_eseampphase_LL(cc) = -sum(xv_Robs.*log(xv_epred_rate)-xv_epred_rate)/sum(xv_Robs);

        avg_rate = mean(Robs);
        null_pred = avg_rate*ones(size(xv_Robs));
        xv_null_LL(cc) = -sum(xv_Robs.*log(null_pred)-null_pred)/sum(xv_Robs);

end

%%
cd ~/Data/bruce/2_27_12
save free_viewing_phase_models_late_v5 *filt *kern wfreqs ntents tent_centers nbins pax xv*
%%
close all
load ./free_viewing_phase_models_late_v4
to_ns_imp = (xv_null_LL - xv_to_ns_LL)/log(2);
to_imp = (xv_null_LL - xv_to_LL)/log(2);
seampphase_imp = (xv_null_LL - xv_seampphase_LL)/log(2);
phase_imp = (xv_null_LL - xv_sinphase_LL)/log(2);
ampphase_imp = (xv_null_LL - xv_ampphase_LL)/log(2);
ucells = 1:9;
ucells = 10:18;

boxplot([to_ns_imp(ucells); to_imp(ucells); seampphase_imp(ucells); phase_imp(ucells); ampphase_imp(ucells)]');
ylim([0 0.7]);

load ./free_viewing_phase_models_all_v4
to_ns_imp = (xv_null_LL - xv_to_ns_LL)/log(2);
to_imp = (xv_null_LL - xv_to_LL)/log(2);
seampphase_imp = (xv_null_LL - xv_seampphase_LL)/log(2);
phase_imp = (xv_null_LL - xv_sinphase_LL)/log(2);
ampphase_imp = (xv_null_LL - xv_ampphase_LL)/log(2);
figure
boxplot([to_ns_imp(ucells); to_imp(ucells); seampphase_imp(ucells); phase_imp(ucells); ampphase_imp(ucells)]');
ylim([0 0.7]);

%%
close all
seampphase_ampkern = sqrt(seampphase_cfilt.^2+seampphase_sfilt.^2);
seampphase_phasekern = -atan2(seampphase_cfilt,seampphase_sfilt)+pi/2;
sephase_ampkern = sqrt(sephase_cfilt.^2+sephase_sfilt.^2);
sephase_phasekern = -atan2(sephase_cfilt,sephase_sfilt)+pi/2;

phase_ax = linspace(-pi,pi,50);
for unit = [3 8 9 16]
    unit
% unit = 4;
phase_mod = nan(length(wfreqs),length(phase_ax));
ampphase_mod = nan(length(wfreqs),length(phase_ax));
for i = 1:length(wfreqs)
   phase_mod(i,:) = cos(phase_ax -  sephase_phasekern(unit,i))*sephase_ampkern(unit,i);
   ampphase_mod(i,:) = cos(phase_ax -  seampphase_phasekern(unit,i))*seampphase_ampkern(unit,i);
end
% subplot(2,2,1)
% pcolor(phase_ax,wfreqs,phase_mod);shading interp; set(gca,'yscale','log')
subplot(2,1,1)
pcolor(phase_ax,wfreqs,phase_mod);shading interp; set(gca,'yscale','log')
subplot(2,1,2)
temp = reshape(ampphase_ampkern(unit,:),length(wfreqs),24)';
% temp2 = reshape(phase_ampkern(unit,:),length(wfreqs),24)';
pcolor(wfreqs,1:24,reshape(ampphase_ampkern(unit,:),length(wfreqs),24)');shading interp
 set(gca,'xscale','log')
 
pause
close all
end
%%
close all
ucells = 1:9;
load ./free_viewing_phase_models_late_v4
shadedErrorBar(wfreqs,nanmean(seampphase_ampkern(ucells,:)),nanstd(seampphase_ampkern(ucells,:))/sqrt(length(ucells)),{'color','b'});
hold on
shadedErrorBar(wfreqs,nanmean(sephase_ampkern(ucells,:)),nanstd(sephase_ampkern(ucells,:))/sqrt(length(ucells)),{'color','r'});
xlim(wfreqs([end 1]))


%%
% close all
% ucells = 1:9;
% load ./free_viewing_phase_models_late_v4
% late_xv_imp = xv_null_LL - xv_seampphase_LL;
% late_ampkern = seampphase_ampkern;
% shadedErrorBar(wfreqs,nanmean(seampphase_ampkern(ucells,:)),nanstd(seampphase_ampkern(ucells,:))/sqrt(length(ucells)),{'color','b'});
% hold on
% load ./free_viewing_phase_models_early_v4
% early_xv_imp = xv_null_LL - xv_seampphase_LL;
% shadedErrorBar(wfreqs,nanmean(seampphase_ampkern(ucells,:)),nanstd(seampphase_ampkern(ucells,:))/sqrt(length(ucells)),{'color','r'});
% early_ampkern = seampphase_ampkern;
% xlim(wfreqs([end 1]))
%%
for cc = 1:n_used_cells
    pcolor(wfreqs,1:24,reshape(phase_ampkern(cc,:),length(wfreqs),24)')
shading flat
% all_probes(cc)
pause
clf
end

%%
for cc = 1:n_used_cells
% % full_phase_filt = [sinphase_cfilt(cc,:) sinphase_sfilt(cc,:)];
% % full_phasemod_out(cc,:) = phase_set(:,use_pset)*full_phase_filt';
% % full_phasemod_out(cc,:) = full_phasemod_out(cc,:) - nanmean(full_phasemod_out(cc,:));

full_ampphase_filt = [ampphase_cfilt(cc,:) ampphase_sfilt(cc,:)];
% curset = bsxfun(@times,ampphase_set(:,use_pset),late_indicator');
curset =ampphase_set(:,use_pset);
full_ampphasemod_out(cc,:) = curset*full_ampphase_filt';
full_ampphasemod_out(cc,:) = full_ampphasemod_out(cc,:) - nanmean(full_ampphasemod_out(cc,:));
end

%%
all_use_late_inds = all_use_inds(late_indicator(all_use_inds) == 1);
for cur_freq = 1:length(wfreqs)
    cur_freq
    cur_fset = find(phase_freq_set==wfreqs(cur_freq));
    for cc = 1:n_used_cells
        full_ampphase_filt = [ampphase_cfilt(cc,cur_fset(1:24)) ampphase_sfilt(cc,cur_fset(1:24))];
        curset =ampphase_set(:,cur_fset);
        fdep_ampphasemod_out(cc,:) = curset*full_ampphase_filt';
    end
    fcmat(cur_freq,:,:) = corr(fdep_ampphasemod_out(:,all_use_late_inds)');
    
end
%%
[~,ord] = sort(all_probes);
for i = 1:length(trial_start_inds)
    cur_set = trial_start_inds(i):trial_stop_inds(i);
%     subplot(2,1,1)
%     imagesc(new_model_time_axis(cur_set)-new_model_time_axis(cur_set(1)),1:length(all_probes),full_phasemod_out(ord,cur_set));
%     caxis([-2 1.5])
%     xlim([0 0.5])
%     subplot(2,1,2)
    imagesc(new_model_time_axis(cur_set)-new_model_time_axis(cur_set(1)),1:length(all_probes),full_ampphasemod_out(ord,cur_set));
    caxis([-2 1.5])
    xlim([0 0.5])
    
    pause
    clf
end

%%
lags = 0:round(Fsd*0.6);
counter = zeros(length(lags),1);
phase_avg = nan(length(trial_start_inds),n_used_cells,length(lags));
amp_avg = nan(length(trial_start_inds),n_used_cells,length(lags));
for i = 1:length(trial_start_inds)
    cur_set = trial_start_inds(i):trial_stop_inds(i);
    cur_set(length(lags)+1:end) = [];
    cl = length(cur_set);
    phase_avg(i,:,1:cl) = full_phasemod_out(:,cur_set);
    amp_avg(i,:,1:cl) = full_ampphasemod_out(:,cur_set);
    counter(1:cl) = counter(1:cl) + 1;
end
amp_avg = squeeze(nanmean(amp_avg));
phase_avg = squeeze(nanmean(phase_avg));

amp_avg_all = nan(24,length(lags));
phase_avg_all = nan(24,length(lags));
   
amp_avg_all(all_probes,:) = amp_avg;
phase_avg_all(all_probes,:) = phase_avg;
    
save phaseavg_outs amp_avg_all phase_avg_all
    %%
    cd /home/james/Data/bruce/2_27_12
    load ./free_viewing_phase_models_late_v3
    fv_ampkern = sqrt(sinphase_cfilt.^2 + sinphase_sfilt.^2);
    
    cd ~/Data/bruce/2_27_12/M232
load ./random_bar_phase_models_v2
rb_ampkern = sqrt(sinphase_cfilt.^2 + sinphase_sfilt.^2);

for cc = 1:24
    cur = find(all_probes==cc);
    if ~isempty(cur)
        subplot(2,1,1)
            pcolor(wfreqs,1:24,reshape(fv_ampkern(cur,:),length(wfreqs),24)')
shading flat
    end
    subplot(2,1,2)
            pcolor(wfreqs,1:24,reshape(rb_ampkern(cc,:),length(wfreqs),24)')
    shading flat
    cc
    pause
    clf
end
    