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

%%
xv_frac = 0.2;
n_xv_trials = round(n_fixs*xv_frac);
xv_set = randperm(n_fixs);
xv_set(n_xv_trials+1:end) = [];
xv_inds = uset(ismember(new_model_fixids(uset),xv_set));
tr_inds = uset(~ismember(uset,xv_inds))';

%%
close all
NL_type = 1; %exp

reg_params.dl1_ind = 1000;
reg_params.dl1_dep = 5000;
reg_params.dl2_ind = 3000;
reg_params.dl2_dep = 60000;
reg_params.dl2_freq_ind = 50;
reg_params.dl2_freq_dep = 1500;
reg_params.dl2_time_ind = 100;
reg_params.dl2_time_dep = 600;
reg_params.l1_ind = 0;
reg_params.l1_dep = 0;
reg_params.l1t_ind = 0;
reg_params.l1t_dep = 0;
reg_params.is_phase = 1;

reg_params2.dl1_ind = 0;
reg_params2.dl1_dep = 0;
reg_params2.dl2_ind = 25000;
reg_params2.dl2_dep = 25000;
reg_params2.dl2_freq_ind = 15000;
reg_params2.dl2_freq_dep = 15000;
reg_params2.dl2_time_ind = 100;
reg_params2.dl2_time_dep = 600;
reg_params2.l1_ind = 0;
reg_params2.l1_dep = 0;
reg_params2.l1t_ind = 0;
reg_params2.l1t_dep = 0;
reg_params2.is_phase = 0;

silent = 1;
nbins = 30;
pax = linspace(-pi,pi,nbins+1);
doub_pax = [pax(1:end-1) 2*pi+pax(1:end-1)];

new_phase_set = [reshape(cos(new_phasegram),length(trial_inds),length(wfreqs)*24) reshape(sin(new_phasegram),length(trial_inds),length(wfreqs)*24)];
phase_elec_set = [repmat(1:24,1,length(wfreqs)) repmat(1:24,1,length(wfreqs))];
phase_elec_set = phase_elec_set(:);

for cc = 1:n_used_cells
    fprintf('Fitting GEM: Cell %d of %d\n',cc,n_used_cells);
    
    Robs = spikes_binned(tr_inds,cc);
    spkbs = convert_to_spikebins(Robs);
    
    
    X_resh = reshape(X,size(X,1),SDIM^2);
    cur_mask1 = get_pgabor_mask(gabor_params_fin(cc,1:6),0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(gabor_params_fin(cc,1:6),pi/2,[SDIM SDIM]);
    mask1_out = X_resh*cur_mask1(:);
    mask2_out = X_resh*cur_mask2(:);
    energy_out = gabor_params_fin(cc,7)*sqrt(mask1_out(new_model_fixids).^2 + mask2_out(new_model_fixids).^2);
    lin_out = gabor_params_fin(cc,8)*mask1_out(new_model_fixids) + gabor_params_fin(cc,9)*mask1_out(new_model_fixids);
    
    total_out = lin_out + energy_out;
    total_out = zscore(total_out);
    
    stim_params = [0,0,ntents];
    Tmat = [trial_Tmat(tr_inds,:) bsxfun(@times,trial_Tmat(tr_inds,:),total_out(tr_inds))];
    Xmat = [Tmat];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_to] = fit_GLM_phase_model(Xmat, Robs, K0, 1, stim_params, reg_params,1,NL_type);
    stim_ind_to_filt(cc,:) = fitp_to.k((1):(ntents));
    stim_dep_to_filt(cc,:) = fitp_to.k((ntents+1):end-1);
    
    stim_params = [0,0,ntents];
    Tmat = [trial_Tmat(tr_inds,:)];
    Xmat = [Tmat];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_to_ns] = fit_GLM_phase_model(Xmat, Robs, K0, 1, stim_params, reg_params,1,NL_type);
    stim_ind_to_ns_filt(cc,:) = fitp_to_ns.k((1):(ntents));
    
    %     nearest_probe = cc;
    if all_probes(cc) < 24
        nearest_probe = all_probes(cc)+1;
    else
        nearest_probe = all_probes(cc)-1;
    end
    cur_alpha_phase = squeeze(new_phasegram(tr_inds,:,nearest_probe));
    Pmat = nan(length(tr_inds),length(wfreqs)*nbins);
    for ww = 1:length(wfreqs)
        cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        cur_set = ((ww-1)*nbins+1):ww*nbins;
        Pmat(:,cur_set) = cur_tb;
        %         Pmat = [Pmat cur_tb];
    end
    
    tr_late_indicator = logical(late_indicator(tr_inds))';
    aPmat = bsxfun(@times,Pmat,tr_late_indicator);
    
    Tmat = [trial_Tmat(tr_inds,:) bsxfun(@times,trial_Tmat(tr_inds,:),total_out(tr_inds))];
    Xmat = [aPmat Tmat];
    stim_params = [nbins,length(wfreqs),ntents,1];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_phase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params,1,NL_type);
    stim_ind_phase_pfilt(cc,:) = fitp_phase.k(1:nbins*length(wfreqs));
    stim_ind_phase_tfilt(cc,:) = fitp_phase.k((nbins*length(wfreqs)+1):(nbins*length(wfreqs)+ntents));
    stim_dep_phase_tfilt(cc,:) = fitp_phase.k((nbins*length(wfreqs)+ntents+1):end-1);
    
    Xmat = [aPmat];
    stim_params = [nbins,length(wfreqs)];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_po] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
    stim_ind_phase_po(cc,:) = fitp_po.k(1:nbins*length(wfreqs));
    
    %     resh_kern = reshape(stim_ind_phase_pfilt(cc,:),nbins,length(wfreqs));
    %     ms_kern = bsxfun(@minus,resh_kern,mean(resh_kern));
    %     ms_kern = bsxfun(@rdivide,ms_kern,std(ms_kern));
    %     ms_kern = ms_kern(:);
    %
    %     phase_outs = nan(length(tr_inds),length(wfreqs),24);
    %     for ll = 1:24
    %         fprintf('Channel %d of %d\n',ll,24);
    %         if ll ~= all_probes(cc)
    %             cur_alpha_phase = squeeze(new_phasegram(tr_inds,:,ll));
    %             for ww = 1:length(wfreqs)
    %                 cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
    %                 cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
    %                 cur_tb(:,end) = [];
    %                 cur_set = ((ww-1)*nbins+1):ww*nbins;
    %                 phase_outs(:,ww,ll) = cur_tb*ms_kern(cur_set);
    %             end
    %         end
    %     end
    %     phase_outs(:,:,all_probes(cc)) = [];
    %     phase_outs = reshape(phase_outs,length(tr_inds),length(wfreqs)*23);
    %     phase_outs = bsxfun(@times,phase_outs,tr_late_indicator);
    %
    %     Xmat = [phase_outs Tmat];
    %     stim_params = [length(wfreqs),23,ntents,1];
    %     klen = size(Xmat,2);
    %     K0 = zeros(klen+1,1);
    %     [fitp_fphase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params2,1,NL_type);
    %     stim_ind_fphase_pfilt(cc,:) = fitp_fphase.k(1:23*length(wfreqs));
    %     stim_ind_fphase_tfilt(cc,:) = fitp_fphase.k((23*length(wfreqs)+1):(23*length(wfreqs)+ntents));
    %     stim_dep_fphase_tfilt(cc,:) = fitp_fphase.k((23*length(wfreqs)+ntents+1):end-1);
    
    use_elecs = 1:24;
    use_elecs(cc) = [];
    use_set = find(ismember(phase_elec_set,use_elecs));
    Xmat = [new_phase_set(tr_inds,use_set) Tmat];
    stim_params = [length(wfreqs),length(use_elecs),ntents,2];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_fphase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params2,1,NL_type);
    %     stim_ind_fphase_pfilt(cc,:) = fitp_fphase.k(1:23*length(wfreqs));
    %     stim_ind_fphase_tfilt(cc,:) = fitp_fphase.k((23*length(wfreqs)+1):(23*length(wfreqs)+ntents));
    %     stim_dep_fphase_tfilt(cc,:) = fitp_fphase.k((23*length(wfreqs)+ntents+1):end-1);
    stim_ind_sinphase_cfilt(cc,:) = fitp_fphase.k(1:length(use_elecs)*length(wfreqs));
    stim_ind_sinphase_sfilt(cc,:) = fitp_fphase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
    stim_ind_sinphase_tfilt(cc,:) = fitp_fphase.k((length(use_elecs)*length(wfreqs)*2+1):(length(use_elecs)*length(wfreqs)*2+ntents));
    stim_dep_sinphase_tfilt(cc,:) = fitp_fphase.k((length(use_elecs)*length(wfreqs)*2+ntents+1):end-1);
    
    
    
    xv_Robs = spikes_binned(xv_inds,cc);
    
    cur_alpha_phase = squeeze(new_phasegram(xv_inds,:,nearest_probe));
    Pmat = nan(length(xv_inds),length(wfreqs)*nbins);
    for ww = 1:length(wfreqs)
        cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        cur_set = ((ww-1)*nbins+1):ww*nbins;
        Pmat(:,cur_set) = cur_tb;
        %         Pmat = [Pmat cur_tb];
    end
    cur_late_indicator = logical(late_indicator(xv_inds))';
    
    %     phase_outs = nan(length(xv_inds),length(wfreqs),24);
    %     for ll = 1:24
    %         fprintf('Channel %d of %d\n',ll,24);
    %         if ll ~= all_probes(cc)
    %             cur_alpha_phase = squeeze(new_phasegram(xv_inds,:,ll));
    %             for ww = 1:length(wfreqs)
    %                 cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
    %                 cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
    %                 cur_tb(:,end) = [];
    %                 cur_set = ((ww-1)*nbins+1):ww*nbins;
    %                 phase_outs(:,ww,ll) = cur_tb*ms_kern(cur_set);
    %             end
    %         end
    %     end
    %     phase_outs(:,:,all_probes(cc)) = [];
    %     phase_outs = reshape(phase_outs,length(xv_inds),length(wfreqs)*23);
    %     phase_outs = bsxfun(@times,phase_outs,cur_late_indicator);
    
    
    Tmat = [trial_Tmat(xv_inds,:) bsxfun(@times,trial_Tmat(xv_inds,:),total_out(xv_inds))];
    Xmat = [Tmat];
    xv_to_pred_rate = Xmat*fitp_to.k(1:end-1) + fitp_to.k(end);
    if NL_type == 0
        xv_to_pred_rate = log(1+exp(xv_to_pred_rate));
    else
        xv_to_pred_rate = exp(xv_to_pred_rate);
    end
    
    aPmat = bsxfun(@times,Pmat,cur_late_indicator);
    Xmat = [aPmat Tmat];
    xv_phase_pred_rate = Xmat*fitp_phase.k(1:end-1) + fitp_phase.k(end);
    if NL_type == 0
        xv_phase_pred_rate = log(1+exp(xv_phase_pred_rate));
    else
        xv_phase_pred_rate = exp(xv_phase_pred_rate);
    end
    
    %     Xmat = [phase_outs Tmat];
    Xmat = [new_phase_set(xv_inds,use_set) Tmat];
    xv_fphase_pred_rate = Xmat*fitp_fphase.k(1:end-1) + fitp_fphase.k(end);
    if NL_type == 0
        xv_fphase_pred_rate = log(1+exp(xv_fphase_pred_rate));
    else
        xv_fphase_pred_rate = exp(xv_fphase_pred_rate);
    end
    
    Tmat = [trial_Tmat(xv_inds,:)];
    Xmat = [Tmat];
    xv_to_ns_pred_rate = Xmat*fitp_to_ns.k(1:end-1) + fitp_to_ns.k(end);
    if NL_type == 0
        xv_to_ns_pred_rate = log(1+exp(xv_to_ns_pred_rate));
    else
        xv_to_ns_pred_rate = exp(xv_to_ns_pred_rate);
    end
    
    xv_phase_LL(cc) = -sum(xv_Robs.*log(xv_phase_pred_rate) - xv_phase_pred_rate)/sum(xv_Robs);
    xv_phase_LL_late(cc) = -sum(xv_Robs(cur_late_indicator).*log(xv_phase_pred_rate(cur_late_indicator)) - ...
        xv_phase_pred_rate(cur_late_indicator))/sum(xv_Robs(cur_late_indicator));
    
    xv_fphase_LL(cc) = -sum(xv_Robs.*log(xv_fphase_pred_rate) - xv_fphase_pred_rate)/sum(xv_Robs);
    xv_fphase_LL_late(cc) = -sum(xv_Robs(cur_late_indicator).*log(xv_fphase_pred_rate(cur_late_indicator)) - ...
        xv_fphase_pred_rate(cur_late_indicator))/sum(xv_Robs(cur_late_indicator));
    
    xv_to_LL(cc) = -sum(xv_Robs.*log(xv_to_pred_rate) - xv_to_pred_rate)/sum(xv_Robs);
    xv_to_LL_late(cc) = -sum(xv_Robs(cur_late_indicator).*log(xv_to_pred_rate(cur_late_indicator)) - ...
        xv_to_pred_rate(cur_late_indicator))/sum(xv_Robs(cur_late_indicator));
    
    xv_to_ns_LL(cc) = -sum(xv_Robs.*log(xv_to_ns_pred_rate) - xv_to_ns_pred_rate)/sum(xv_Robs);
    xv_to_ns_LL_late(cc) = -sum(xv_Robs(cur_late_indicator).*log(xv_to_ns_pred_rate(cur_late_indicator)) - ...
        xv_to_ns_pred_rate(cur_late_indicator))/sum(xv_Robs(cur_late_indicator));
    
    avg_rate = mean(Robs);
    null_pred = avg_rate*ones(size(xv_Robs));
    xv_null_LL(cc) = -sum(xv_Robs.*log(null_pred)-null_pred)/sum(xv_Robs);
    avg_rate = mean(Robs(tr_late_indicator));
    null_pred = avg_rate*ones(size(xv_Robs));
    xv_null_LL_late(cc) = -sum(xv_Robs(cur_late_indicator).*log(null_pred(cur_late_indicator)) -...
        null_pred(cur_late_indicator))/sum(xv_Robs(cur_late_indicator));
    
end

%%
save free_viewing_phase_models stim* wfreqs pax nbins xv* all_probes

%%
cd /home/james/Data/bruce/2_27_12
load ./free_viewing_phase_models
fv_ampkern = sqrt(stim_ind_sinphase_cfilt.^2+stim_ind_sinphase_sfilt.^2);
fv_wfreqs = wfreqs;
fv_nbins = nbins;

cd /home/james/Data/bruce/2_27_12/M232/
load ./bar_phase_models_ftime_v2
stim_ind_ampkern = sqrt(stim_ind_sinphase_cfilt.^2+stim_ind_sinphase_sfilt.^2);

for cc = 1:length(all_probes)
    subplot(2,1,1)
    pcolor(fv_wfreqs,1:23,reshape(fv_ampkern(cc,:),length(fv_wfreqs),23)');shading flat
    
    subplot(2,1,2)
    for dd = 1:24
        fprintf('%d vs %d\n',all_probes(cc),dd);
        pcolor(wfreqs,1:23,reshape(stim_ind_ampkern(dd,:),length(wfreqs),23)');shading flat
        
        pause
    end
    
    
end


%%
close all
stim_ind_ampkern = sqrt(stim_ind_sinphase_cfilt.^2+stim_ind_sinphase_sfilt.^2);
stim_ind_phasekern = atan2(stim_ind_sinphase_sfilt,stim_ind_sinphase_cfilt);
for cc = 1:24
    all_probes(cc)
    temp = stim_ind_phase_pfilt(cc,:);
    temp = temp-mean(temp);
    %     temp2 = stim_ind_dphase_pfilt(cc,:);
    %     temp2 = temp2-mean(temp2);
    subplot(3,1,1)
    pcolor(pax(1:end-1),wfreqs,reshape(temp,nbins,length(wfreqs))');shading flat
    %         set(gca,'yscale','log')
    %     caxis([-0.05 0.05])
    %     subplot(2,2,2)
    %     pcolor(pax(1:end-1),wfreqs,reshape(temp2,nbins,length(wfreqs))');shading flat
    %         set(gca,'yscale','log')
    %     caxis([-0.05 0.05])
    %     colorbar
    subplot(3,1,2)
    pcolor(wfreqs,1:23,reshape(stim_ind_ampkern(cc,:),length(wfreqs),23)');shading flat;
    subplot(3,1,3)
    pcolor(wfreqs,1:23,reshape(stim_ind_phasekern(cc,:),length(wfreqs),23)');shading flat;
    %     set(gca,'xscale','log')
    %     colorbar
    %     caxis([0 0.005])
    %     subplot(2,2,4)
    %     plot(tent_centers*new_dt,stim_ind_phase_tfilt(cc,:))
    %     hold on
    %     plot(tent_centers*new_dt,stim_ind_sac_kern(cc,:),'r')
    %     axis tight
    pause
    clf
end


