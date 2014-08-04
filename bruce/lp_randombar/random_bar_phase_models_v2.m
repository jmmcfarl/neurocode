clear all
close all

cd ~/Data/bruce/2_27_12/M232
% load ./random_bar_eyedata.mat
load ./random_bar_eyedata_ftime.mat

load ./lemM232Expts.mat
%%
rf_cent = [4.73 -2.2];
axis_or = 50*pi/180; %in radians

Fs = 1000;
dsf = 3;
Fsd = Fs/dsf;
niqf = Fs/2;
[bb,aa] = butter(2,120/niqf,'low');
scales = logspace(log10(5),log10(60),30);
scales = [scales 70 80 90];
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

new_dt = .0025;
dsfrac = 4;
%% COMPUTE PERPINDICULAR EYE POSITION
for ee = 1:length(bar_expts)
    
    avg_eyepos = 0.5*all_eye_vals_cal{ee}(:,1:2) + 0.5*all_eye_vals_cal{ee}(:,3:4);
    
    fix_r = sqrt(sum(avg_eyepos.^2,2));
    fix_angle = atan2(avg_eyepos(:,2),avg_eyepos(:,1));
    angle_diff = fix_angle - axis_or;
    fix_perp{ee} = fix_r.*sin(angle_diff);
    
end

%%
load ./CellList.mat
good_sus = find(all(CellList(bar_expts,:,1) > 0));
% use_sus = find(all(CellList(bar_expts,:,1) > 0));
use_sus = 1:24;
% used_expts{1} = 1:6;
% used_expts{2} = 1:5;
% used_expts{3} = 1:6;
% used_expts{4} =
%%
load ./un_bar_pos.mat
un_bar_pos(1:2) = [];

flen = 20;
full_Xmat = [];
full_spkbinned = [];
full_exptvec = [];
full_exptvec_new = [];
full_trialvec = [];
full_taxis = [];
full_taxis_new = [];
full_old_t_inds = [];
full_bar_pos = [];
for ee = 1:length(bar_expts);
    fprintf('Processing expt %d of %d\n',ee,length(bar_expts));
    fname = sprintf('Expt%dClusterTimes.mat',bar_expts(ee));
    load(fname);
    
    trial_durs{ee} = [Expts{bar_expts(ee)}.Trials(:).dur];
    n_trials(ee) = length(Expts{bar_expts(ee)}.Trials);
    all_t_axis = [];
    all_used_inds = [];
    all_t_axis_new = [];
    all_old_t_inds = [];
    all_used_inds_new = [];
    all_expt_vec = [];
    all_trial_vec = [];
    all_bar_Op = [];
    all_bar_e0 = [];
    all_bar_Xmat = [];
    all_binned_spikes = [];
    all_used_inds = [];
    for tt = 1:n_trials(ee)
        
        cur_t_axis = [Expts{bar_expts(ee)}.Trials(tt).Start]/1e4;
        cur_bar_Op = [Expts{bar_expts(ee)}.Trials(tt).Op];
        %         cur_bar_e0 = [Expts{bar_expts(ee)}.Trials(tt).e0];
        
        cur_t_edges = [cur_t_axis; Expts{bar_expts(ee)}.Trials(tt).End(end)/1e4];
        cur_t_edges_new = cur_t_edges(1):new_dt:cur_t_edges(end);
        cur_t_axis_new = cur_t_edges_new(1:end-1);
        
        cur_binned_spks = nan(length(cur_t_axis_new),length(use_sus));
        for cc = 1:length(use_sus)
            cur_hist = histc(Clusters{use_sus(cc)}.times,cur_t_edges_new);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
        end
        
        %         cur_interp_eye_perp = interp1(all_eye_ts{ee},fix_perp{ee},cur_t_axis);
        %         cur_bar_Op_cor = cur_bar_Op - cur_interp_eye_perp;
        
        cur_bar_mat = zeros(length(cur_t_axis),length(un_bar_pos));
        for bb = 1:length(un_bar_pos)
            cur_set = find(cur_bar_Op==un_bar_pos(bb));
            cur_bar_mat(cur_set,bb) = 1;
            %             cur_bar_mat(cur_set,bb) = cur_bar_e0(cur_set);
        end
        bar_Xmat = makeStimRows(cur_bar_mat,flen);
        
        cur_used_inds = ones(length(cur_t_axis),1);
        cur_used_inds(1:flen) = 0;
        cur_used_inds_new = ones(length(cur_t_axis_new),1);
        cur_used_inds_new(1:(flen*dsfrac)) = 0;
        
        cur_used_inds_new(end-round(0.1/new_dt):end) = 0; %reduce phase estimation edge artifacts
        
        all_used_inds_new = [all_used_inds_new; cur_used_inds_new(:)];
        all_t_axis = [all_t_axis; cur_t_axis(:)];
        all_used_inds = [all_used_inds; cur_used_inds(:)];
        all_t_axis_new = [all_t_axis_new; cur_t_axis_new(:)];
        %         all_old_t_inds = [all_old_t_inds; old_t_inds(:)];
        all_binned_spikes = [all_binned_spikes; cur_binned_spks];
        all_expt_vec = [all_expt_vec; ones(length(cur_t_axis_new),1)*bar_expts(ee)];
        all_trial_vec = [all_trial_vec; ones(length(cur_t_axis),1)*tt];
        all_bar_Op = [all_bar_Op; cur_bar_Op(:)];
        %         all_bar_e0 = [all_bar_e0; cur_bar_e0(:)];
        all_bar_Xmat = [all_bar_Xmat; bar_Xmat];
    end
    
    cur_blink_start_times = all_eye_ts{ee}(all_blink_startinds{ee});
    cur_blink_stop_times = all_eye_ts{ee}(all_blink_stopinds{ee});
    cur_blink_start_inds = round(interp1(all_t_axis,1:length(all_t_axis),cur_blink_start_times));
    cur_blink_stop_inds = round(interp1(all_t_axis,1:length(all_t_axis),cur_blink_stop_times));
    cur_poss_blinks = find(~isnan(cur_blink_start_inds) & ~isnan(cur_blink_stop_inds));
    
    all_blink_inds = zeros(size(all_t_axis));
    for i = 1:length(cur_poss_blinks)
        all_blink_inds(cur_blink_start_inds(cur_poss_blinks(i)):cur_blink_stop_inds(cur_poss_blinks(i))) = 1;
    end
    
    cur_blink_start_inds = round(interp1(all_t_axis_new,1:length(all_t_axis_new),cur_blink_start_times));
    cur_blink_stop_inds = round(interp1(all_t_axis_new,1:length(all_t_axis_new),cur_blink_stop_times));
    cur_poss_blinks = find(~isnan(cur_blink_start_inds) & ~isnan(cur_blink_stop_inds));
    all_blink_inds_new = zeros(size(all_t_axis_new));
    for i = 1:length(cur_poss_blinks)
        all_blink_inds_new(cur_blink_start_inds(cur_poss_blinks(i)):cur_blink_stop_inds(cur_poss_blinks(i))) = 1;
    end
    
    all_used_inds(all_blink_inds == 1) = 0;
    all_used_inds = logical(all_used_inds);
    all_used_inds_new(all_blink_inds_new == 1) = 0;
    all_used_inds_new = logical(all_used_inds_new);
    
    full_Xmat = [full_Xmat; all_bar_Xmat(all_used_inds,:)];
    full_spkbinned = [full_spkbinned; all_binned_spikes(all_used_inds_new,:)];
    full_exptvec = [full_exptvec; ones(sum(all_used_inds),1)*ee];
    full_exptvec_new = [full_exptvec_new; ones(sum(all_used_inds_new),1)*ee];
    full_taxis = [full_taxis; all_t_axis(all_used_inds)];
    full_taxis_new = [full_taxis_new; all_t_axis_new(all_used_inds_new)];
    %     full_old_t_inds = [full_old_t_inds; all_old_t_inds(all_used_inds_new)];
    full_trialvec = [full_trialvec; all_trial_vec(all_used_inds)];
    full_bar_pos = [full_bar_pos; all_bar_Op(all_used_inds)];
end

%%
full_lfps = [];
full_phasegrams = [];
for ee = 1:length(bar_expts);
    fprintf('Expt %d of %d\n',ee,length(bar_expts));
    fname = sprintf('lemM232A.%d.lfp.mat',bar_expts(ee));
    load(fname);
    
    Fs = 1/LFP.Header.CRsamplerate;
    
    n_trials(ee) = length(LFP.Trials);
    %     lfp_trial_starts = [LFP.Trials(:).Start]/1e4;
    lfp_trial_starts = [LFP.Trials(:).ftime]/1e4;
    lfp_trial_ends = [LFP.Trials(:).End]/1e4;
    expt_lfp_t_axis = [];
    expt_lfps = [];
    for tt = 1:n_trials(ee)
        
        cur_npts = size(LFP.Trials(tt).LFP,1);
        cur_t_end(tt) = lfp_trial_starts(tt)+(cur_npts-1)/Fs;
        cur_t_axis = lfp_trial_starts(tt):1/Fs:cur_t_end(tt);
        
        if ~isempty(expt_lfp_t_axis)
            cur_sp = find(cur_t_axis > max(expt_lfp_t_axis),1,'first');
        else
            cur_sp = 1;
        end
        cur_t_axis = cur_t_axis(cur_sp:end);
        cur_LFP = [LFP.Trials(tt).LFP];
        cur_LFP = filtfilt(bb,aa,cur_LFP);
        
        cur_LFP = cur_LFP(cur_sp:end,:);
        
        
        expt_lfp_t_axis = [expt_lfp_t_axis; cur_t_axis(:)];
        expt_lfps = [expt_lfps; cur_LFP];
    end
    
    seg_start_inds = [1; 1+find(diff(expt_lfp_t_axis) > 2/Fs)];
    seg_stop_inds = [find(diff(expt_lfp_t_axis) > 2/Fs); length(expt_lfp_t_axis)];
    res_t_axis = [];
    res_phasegram = [];
    res_lfps = [];
    for jj = 1:length(seg_start_inds)
        cur_res_inds = seg_start_inds(jj):seg_stop_inds(jj);
        cur_res_t = expt_lfp_t_axis(seg_start_inds(jj)):dsf/Fs:expt_lfp_t_axis(seg_stop_inds(jj));
        res_t_axis = [res_t_axis; cur_res_t(:)];
        
        interp_lfps = interp1(expt_lfp_t_axis(cur_res_inds),expt_lfps(cur_res_inds,:),cur_res_t);
        
        cur_phasegram = nan(length(cur_res_t),length(wfreqs),24);
        for ll = 1:24
            temp = cwt(interp_lfps(:,ll),scales,'cmor1-1');
            cur_phasegram(:,:,ll) = angle(temp)';
        end
        res_phasegram = cat(1,res_phasegram,cur_phasegram);
        res_lfps = [res_lfps; interp_lfps];
    end
    
    cur_set = find(full_exptvec_new==ee);
    interp_lfps = interp1(res_t_axis,res_lfps,full_taxis_new(cur_set));
    unwr_phasegram = unwrap(res_phasegram);
    interp_phasegrams = interp1(res_t_axis,unwr_phasegram,full_taxis_new(cur_set));
    interp_phasegrams = mod(interp_phasegrams+pi,2*pi)-pi;
    %     interp_phasegrams = interp1(res_t_axis,res_phasegram,full_taxis_new(cur_set));
    
    full_lfps = [full_lfps; interp_lfps];
    full_phasegrams = cat(1,full_phasegrams, interp_phasegrams);
end

%% create mapping to old taxis (with stimulus stored)
full_old_inds = [];
ind_cnt = 0;
for ee = 1:length(bar_expts)
    cur_set = find(full_exptvec==ee);
    cur_set_new = find(full_exptvec_new==ee);
    cur_inds = round(interp1(full_taxis(cur_set),1:length(cur_set),full_taxis_new(cur_set_new)));
    full_old_inds = [full_old_inds; cur_inds + ind_cnt];
    ind_cnt = ind_cnt + length(cur_set);
end

%%
[c,ia,ic] = unique([full_exptvec full_trialvec],'rows');
n_trials = length(ia);

xv_frac = 0.2;
n_xv_trials = round(n_trials*xv_frac);
xv_set = randperm(n_trials);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(ic,xv_set));
tr_inds = find(~ismember(ic,xv_inds))';

xv_inds_new = find(ismember(full_old_inds,xv_inds));
tr_inds_new = find(ismember(full_old_inds,tr_inds));

stim_Xmat = full_Xmat;
%%

load ./random_bar_unit_models_ftime.mat
%%
close all
NL_type = 1; %exp

reg_params.dl1_ind = 0;
reg_params.dl1_dep = 0;
reg_params.dl2_ind = 20000;
reg_params.dl2_dep = 20000;
reg_params.dl2_freq_ind = 20000;
reg_params.dl2_freq_dep = 20000;
reg_params.dl2_time_ind = 0;
reg_params.dl2_time_dep = 0;
reg_params.l1_ind = 0;
reg_params.l1_dep = 0;
reg_params.l1t_ind = 0;
reg_params.l1t_dep = 0;
reg_params.is_phase = 0;

reg_params2.dl1_ind = 3000;
reg_params2.dl1_dep = 0;
reg_params2.dl2_ind = 30000;
reg_params2.dl2_dep = 0;
reg_params2.dl2_freq_ind = 300;
reg_params2.dl2_freq_dep = 0;
reg_params2.dl2_time_ind = 0;
reg_params2.dl2_time_dep = 0;
reg_params2.l1_ind = 0;
reg_params2.l1_dep = 0;
reg_params2.is_phase = 1;

silent = 1;

nbins = 30;
pax = linspace(-pi,pi,nbins+1);
doub_pax = [pax(1:end-1) 2*pi+pax(1:end-1)];

new_phase_set = [reshape(cos(full_phasegrams),length(full_taxis_new),length(wfreqs)*24) reshape(sin(full_phasegrams),length(full_taxis_new),length(wfreqs)*24)];
% phase_elec_set = [repmat(1:24,1,length(wfreqs)) repmat(1:24,1,length(wfreqs))];
phase_elec_set = ones(length(wfreqs),1)*(1:24);
phase_elec_set = [phase_elec_set(:); phase_elec_set(:)];

use_elecs = 1:24;
use_set = find(ismember(phase_elec_set,use_elecs));
Xmat = [new_phase_set(tr_inds_new,use_set)];
Xmatxv = [new_phase_set(xv_inds_new,use_set)];

%%
for cc = 1:length(use_sus)
    fprintf('Cell %d of %d\n',cc,length(use_sus));
    
    Robs = full_spkbinned(tr_inds_new,cc);
    tr_spkbns = convert_to_spikebins(Robs);
    
    glm_kern = get_k_mat(glm_fit(cc));
    stim_out = stim_Xmat*glm_kern;
    stim_out_interp = stim_out(full_old_inds(tr_inds_new));
    
    klen = size(stim_out_interp,2);
    K0 = zeros(klen+1,1);
    [fitp_so,grad] = GLMsolve_jmm(stim_out_interp, tr_spkbns, K0, 1, [], [],[], [], [], [], NL_type);
    
%     %     if cc < 22
%     %         nearest_probe = cc+1;
%     %     else
%     %         nearest_probe = cc-1;
%     %     end
%     nearest_probe = cc;
%     cur_alpha_phase = squeeze(full_phasegrams(tr_inds_new,:,nearest_probe));
%     Pmat = nan(length(tr_inds_new),length(wfreqs)*nbins);
%     for ww = 1:length(wfreqs)
%         cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
%         cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
%         cur_tb(:,end) = [];
%         cur_set = ((ww-1)*nbins+1):ww*nbins;
%         Pmat(:,cur_set) = cur_tb;
%     end
    
%     Xmat = [Pmat];
%     stim_params = [nbins,length(wfreqs)];
%     klen = size(Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_dphase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params2,0,NL_type);
%     phase_pfilt(cc,:) = fitp_dphase.k(1:nbins*length(wfreqs));
    
%     %     use_elecs = cc;
%     use_elecs = 1:24;
%     %     use_elecs(cc) = [];
%     use_set = find(ismember(phase_elec_set,use_elecs));
    stim_params = [length(wfreqs),length(use_elecs)];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_fphase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
    sinphase_cfilt(cc,:) = fitp_fphase.k(1:length(use_elecs)*length(wfreqs));
    sinphase_sfilt(cc,:) = fitp_fphase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
    ampkern(cc,:) = sqrt(sinphase_cfilt(cc,:).^2 + sinphase_sfilt(cc,:).^2);

%     nXmat = [Xmat stim_out_interp];
%         klen = size(nXmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_nfphase] = fit_GLM_phase_model(nXmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
%     stsinphase_cfilt(cc,:) = fitp_nfphase.k(1:length(use_elecs)*length(wfreqs));
%     stsinphase_sfilt(cc,:) = fitp_nfphase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
%     stampkern(cc,:) = sqrt(stsinphase_cfilt(cc,:).^2 + stsinphase_sfilt(cc,:).^2);

%         for ww = 1:length(wfreqs)
%             use_phase = squeeze(full_phasegrams(tr_inds_new,ww,cc));
%             Xmat = [cos(use_phase) sin(use_phase)];
%             b(ww,:) = glmfit(Xmat,Robs,'poisson');
%         end
%         ind_amp(cc,:) = sqrt(b(:,2).^2 + b(:,3).^2);
    
    
    xv_Robs = full_spkbinned(xv_inds_new,cc);
    stim_out_interp = stim_out(full_old_inds(xv_inds_new));
    
    xv_so_pred_rate = stim_out_interp*fitp_so.k(1:end-1) + fitp_so.k(end);
    if NL_type == 1
        xv_so_pred_rate = exp(xv_so_pred_rate);
    else
        xv_so_pred_rate = log(1+exp(xv_so_pred_rate));
    end
    xv_so_LL(cc) = -sum(xv_Robs.*log(xv_so_pred_rate)-xv_so_pred_rate)/sum(xv_Robs);
    
%     for ww = 1:length(wfreqs)
%         use_phase = squeeze(full_phasegrams(xv_inds_new,ww,cc));
%         Xmat = [cos(use_phase) sin(use_phase)];
%         cur_out = Xmat*b(ww,2:end)'+b(ww,1);
%         cur_out = exp(cur_out);
%         ind_xvLL(cc,ww) = -sum(xv_Robs.*log(cur_out)-cur_out)/sum(xv_Robs);
%     end
%     ind_amp(cc,:) = sqrt(b(:,2).^2 + b(:,3).^2);
    
%         cur_alpha_phase = squeeze(full_phasegrams(xv_inds_new,:,nearest_probe));
%         Pmat = nan(length(xv_inds_new),length(wfreqs)*nbins);
%         for ww = 1:length(wfreqs)
%             cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
%             cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
%             cur_tb(:,end) = [];
%             cur_set = ((ww-1)*nbins+1):ww*nbins;
%             Pmat(:,cur_set) = cur_tb;
%         end
%     
%         Xmat = [Pmat];
%         xv_phase_pred_rate = Xmat*fitp_dphase.k(1:end-1) + fitp_dphase.k(end);
%         if NL_type == 0
%         xv_phase_pred_rate = log(1+exp(xv_phase_pred_rate));
%         else
%         xv_phase_pred_rate = exp(xv_phase_pred_rate);
%         end
%         xv_phase_LL(cc) = -sum(xv_Robs.*log(xv_phase_pred_rate)-xv_phase_pred_rate)/sum(xv_Robs);
    
    
%     Xmat = [new_phase_set(xv_inds_new,use_set)];
    xv_fphase_pred_rate = Xmatxv*fitp_fphase.k(1:end-1) + fitp_fphase.k(end);
    if NL_type == 0
        xv_fphase_pred_rate = log(1+exp(xv_fphase_pred_rate));
    else
        xv_fphase_pred_rate = exp(xv_fphase_pred_rate);
    end
    xv_fphase_LL(cc) = -sum(xv_Robs.*log(xv_fphase_pred_rate)-xv_fphase_pred_rate)/sum(xv_Robs);

%     nXmatxv = [Xmatxv stim_out_interp];
%         xv_fphase_pred_rate = nXmatxv*fitp_nfphase.k(1:end-1) + fitp_nfphase.k(end);
%     if NL_type == 0
%         xv_fphase_pred_rate = log(1+exp(xv_fphase_pred_rate));
%     else
%         xv_fphase_pred_rate = exp(xv_fphase_pred_rate);
%     end
%     xv_nfphase_LL(cc) = -sum(xv_Robs.*log(xv_fphase_pred_rate)-xv_fphase_pred_rate)/sum(xv_Robs);
%
    %
    avg_rate = mean(Robs);
    null_pred = avg_rate*ones(size(xv_Robs));
    xv_null_LL(cc) = -sum(xv_Robs.*log(null_pred)-null_pred)/sum(xv_Robs);
    
end

%%
cd ~/Data/bruce/2_27_12/M232
save random_bar_phase_models_v2 xv_* *filt wfreqs pax nbins
%%
% n_prc_bins = 10;
% prc_bin_edges = linspace(0,1,n_prc_bins+1);
% prc_bin_edges = 1+floor(prc_bin_edges*length(tr_inds_new));
% prc_bin_edges(prc_bin_edges > length(tr_inds_new)) = length(tr_inds_new);
% 
% n_phase_bins = 8;
% phase_bin_edges = linspace(-pi,pi,n_phase_bins+1);
% 
% clear np_phase_stim
% cc=19;
% % for cc = 1:24
% fprintf('Cell %d of %d\n',cc,24);
% Robs = full_spkbinned(tr_inds_new,cc);
% tr_spkbns = convert_to_spikebins(Robs);
% 
% glm_kern = get_k_mat(glm_fit(cc));
% stim_out = stim_Xmat*glm_kern;
% stim_out_interp = stim_out(full_old_inds(tr_inds_new));
% 
% phase_rate = zeros(length(wfreqs),n_phase_bins);
% phase_lock = zeros(length(wfreqs),1);
% for use_freq = 1:length(wfreqs)
%     use_phase = squeeze(full_phasegrams(tr_inds_new,use_freq,cc));
%     for pp = 1:n_phase_bins
%         phase_set = find(use_phase >= phase_bin_edges(pp) & use_phase < phase_bin_edges(pp+1));
%         phase_rate(use_freq,pp) = mean(Robs(phase_set));
%     end
%     phase_lock(use_freq) = circ_kappa(use_phase(tr_spkbns));
% end
% 
% 
% np_stim = zeros(n_prc_bins,1);
% for i = 1:n_prc_bins
%     cur_set = sort_ord(prc_bin_edges(i):prc_bin_edges(i+1));
%     np_stim(i) = mean(Robs(cur_set));
% end
% for use_freq = 1:length(wfreqs)
%     use_phase = squeeze(full_phasegrams(tr_inds_new,use_freq,cc));
%     [~,sort_ord] = sort(stim_out_interp);
%     np_phase_stim(use_freq,:,:) = zeros(n_prc_bins,n_phase_bins);
%     for i = 1:n_prc_bins
%         cur_set = sort_ord(prc_bin_edges(i):prc_bin_edges(i+1));
%         cur_phase = use_phase(cur_set);
%         for pp = 1:n_phase_bins
%             phase_set = find(cur_phase >= phase_bin_edges(pp) & cur_phase < phase_bin_edges(pp+1));
%             np_phase_stim(use_freq,i,pp) = mean(Robs(cur_set(phase_set)));
%         end
%     end
%     
% end

% figure(1)
% figure(2)
% figure(3)
% pcolor(pax(1:end-1),wfreqs,reshape(phase_pfilt(cc,:),nbins,length(wfreqs))');shading flat
% for i = 1:33
%     figure(1)
%     plot(wfreqs,ind_amp(cc,:));hold on
%     plot(wfreqs(i),ind_amp(cc,i),'ro')
%     figure(2)
%     cur_ps = squeeze(np_phase_stim(i,:,:));
%     norm_ps = bsxfun(@rdivide,cur_ps,mean(cur_ps,2));
%     norm_ps2 = bsxfun(@rdivide,cur_ps,mean(cur_ps));
%     pred_out = np_stim*phase_rate(i,:);
%     subplot(2,2,1)
%     imagesc(cur_ps)
%     subplot(2,2,2)
%     imagesc(norm_ps)
%     subplot(2,2,3)
%     imagesc(pred_out)
%     subplot(2,2,4)
%     imagesc(norm_ps2)
%     
%     
%     pause
%     figure(1);clf;figure(2);clf;
% end

%%
close all
stim_ind_ampkern = sqrt(sinphase_cfilt.^2 + sinphase_sfilt.^2);
stim_ind_phasekern = -atan2(sinphase_cfilt,sinphase_sfilt)+pi/2;
for cc = 1:24
    cc
    temp = phase_pfilt(cc,:);
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
% %     colorbar
    subplot(3,1,2)
    pcolor(wfreqs,1:24,reshape(stim_ind_ampkern(cc,:),length(wfreqs),24)');shading flat
    
    subplot(3,1,3)
    pcolor(wfreqs,1:24,reshape(stim_ind_phasekern(cc,:),length(wfreqs),24)');shading flat
    %     set(gca,'xscale','log')
%     colorbar
%     caxis([0 0.005])
%     subplot(2,2,4)
%     plot(tent_centers*new_dt,stim_ind_dphase_tfilt(cc,:))
%     hold on
%     plot(tent_centers*new_dt,stim_ind_sac_kern(cc,:),'r')
% %     cc
%     fprintf('Phase imp: %.3f\n',xv_phase_imp(cc));
%     fprintf('Depth imp: %.3f\n',xv_fphase_imp(cc));

    pause
    clf
end
%%
close all
stim_ind_ampkern = sqrt(sinphase_cfilt.^2 + sinphase_sfilt.^2);
stim_dep_ampkern = sqrt(stsinphase_cfilt.^2 + stsinphase_sfilt.^2);
for cc = 1:24
    cc
    subplot(2,1,1)
    pcolor(wfreqs,1:24,reshape(stim_ind_ampkern(cc,:),length(wfreqs),24)');shading flat
    
    subplot(2,1,2)
    pcolor(wfreqs,1:24,reshape(stim_dep_ampkern(cc,:),length(wfreqs),24)');shading flat

    pause
    clf
end