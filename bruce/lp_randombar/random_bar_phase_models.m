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
all_fix_inds = [];
all_first_inds = [];
all_second_inds = [];
all_small_inds = [];
for ee = 1:length(bar_expts)
    expt_trial_starts = [Expts{bar_expts(ee)}.Trials(:).TrialStart]/1e4;
    expt_trial_end = [Expts{bar_expts(ee)}.Trials(:).TrueEnd]/1e4;
    
    sac_start_times = all_eye_ts{ee}(all_sac_startinds{ee});
    sac_stop_times = all_eye_ts{ee}(all_sac_stopinds{ee});
    sac_amps = all_sac_peakamps{ee};
    
    cur_sac_set = find(all_expt_num==ee);
    if length(cur_sac_set) ~= length(sac_start_times)
        error('saccade mismatch')
    end
    intrial_sac_set = find(~isnan(all_trial_num(cur_sac_set)));
    infirst_sac_set = find(ismember(cur_sac_set,all_isfirst));
    insecond_sac_set = find(ismember(cur_sac_set,all_issecond));
    micro_sac_set = find(sac_amps < 50);
    
%     fix_times = sac_stop_times(intrial_sac_set);
%     first_fix_times = sac_stop_times(infirst_sac_set);
%     second_fix_times = sac_stop_times(insecond_sac_set);
    fix_times = sac_start_times(intrial_sac_set);
    first_fix_times = sac_start_times(infirst_sac_set);
    second_fix_times = sac_start_times(insecond_sac_set);
    small_fix_times = sac_start_times(micro_sac_set);
    
    cur_e_set = find(full_exptvec_new == ee);
    
    fix_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),fix_times));
    first_fix_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),first_fix_times));
    second_fix_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),second_fix_times));
    small_fix_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),small_fix_times));

    fix_inds(isnan(fix_inds)) = [];
    first_fix_inds(isnan(first_fix_inds)) = [];
    second_fix_inds(isnan(second_fix_inds)) = [];
    small_fix_inds(isnan(small_fix_inds)) = [];
    
    all_fix_inds = [all_fix_inds; cur_e_set(fix_inds(:))];
    all_first_inds = [all_first_inds; cur_e_set(first_fix_inds(:))];
    all_second_inds = [all_second_inds; cur_e_set(second_fix_inds(:))];
    all_small_inds = [all_small_inds; cur_e_set(small_fix_inds(:))];
end

used_sac_inds = sort([all_first_inds; all_second_inds]);
% used_sac_inds = sort(all_fix_inds);
% small_set = setdiff(all_fix_inds,used_sac_inds);

%%
cur_dt = 0.01;
% flen_t = 0.5;
tent_centers = [0:cur_dt:0.8];
tent_centers = round(tent_centers/new_dt);
tbmat = construct_tent_bases(tent_centers,1);
[ntents,tblen] = size(tbmat);

shift = round(0.3/new_dt);
tbmat = [zeros(ntents,tblen-2*shift-1) tbmat];
tent_centers = tent_centers-shift;

trial_inds = zeros(size(full_taxis_new));
trial_inds(used_sac_inds) = 1;
trial_Tmat = zeros(length(full_taxis_new),ntents);
for i = 1:ntents
    trial_Tmat(:,i) = conv(trial_inds,tbmat(i,:),'same');
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

reg_params.dl1_ind = 2000;
reg_params.dl1_dep = 0;
reg_params.dl2_ind = 20000;
reg_params.dl2_dep = 0;
reg_params.dl2_freq_ind = 200;
reg_params.dl2_freq_dep = 0;
reg_params.dl2_time_ind = 0;
reg_params.dl2_time_dep = 0;
reg_params.l1_ind = 0;
reg_params.l1_dep = 0;
reg_params.is_phase = 1;

reg_params_s.dl1_ind = 2000;
reg_params_s.dl1_dep = 0;
reg_params_s.dl2_ind = 20000;
reg_params_s.dl2_dep = 0;
reg_params_s.dl2_freq_ind = 200;
reg_params_s.dl2_freq_dep = 0;
reg_params_s.dl2_time_ind = 400;
reg_params_s.dl2_time_dep = 0;
reg_params_s.l1_ind = 0;
reg_params_s.l1_dep = 0;
reg_params_s.l1t_ind = 0;
reg_params_s.l1t_dep = 0;
reg_params_s.is_phase = 1;
% 
reg_params2.dl1_ind = 0;
reg_params2.dl1_dep = 0;
reg_params2.dl2_ind = 40000;
reg_params2.dl2_dep = 40000;
reg_params2.dl2_freq_ind = 20000;
reg_params2.dl2_freq_dep = 20000;
reg_params2.dl2_time_ind = 0;
reg_params2.dl2_time_dep = 0;
reg_params2.l1_ind = 0;
reg_params2.l1_dep = 0;
reg_params2.l1t_ind = 0;
reg_params2.l1t_dep = 0;
reg_params2.is_phase = 0;

silent = 1;

nbins = 30;
pax = linspace(-pi,pi,nbins+1);
doub_pax = [pax(1:end-1) 2*pi+pax(1:end-1)];

new_phase_set = [reshape(cos(full_phasegrams),length(full_taxis_new),length(wfreqs)*24) reshape(sin(full_phasegrams),length(full_taxis_new),length(wfreqs)*24)];
phase_elec_set = [repmat(1:24,1,length(wfreqs)) repmat(1:24,1,length(wfreqs))];
phase_elec_set = phase_elec_set(:);
%
for cc = 1:length(use_sus)
    fprintf('Cell %d of %d\n',cc,length(use_sus));
    
    Robs = full_spkbinned(tr_inds_new,cc);
    tr_spkbns = convert_to_spikebins(Robs);
    
    
    glm_kern = get_k_mat(glm_fit(cc));
    stim_out = stim_Xmat*glm_kern;
    stim_out_interp = stim_out(full_old_inds(tr_inds_new));
    
    Xmat = [stim_out_interp];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_so,grad] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], [],[], [], [], [], NL_type);

    Xmat = [trial_Tmat(tr_inds_new,:) stim_out_interp];
    klen = size(Xmat,2);
    lamrange2 = [400 1 ntents 0];
    K0 = zeros(klen+1,1);
    [fitp_sac,grad] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], lamrange2,[], [], [], [], NL_type);
    stim_ind_sac_kern(cc,:) = fitp_sac.k(1:ntents);
    
    if cc < 24
        nearest_probe = cc+1;
    else
        nearest_probe = cc-1;
    end
% nearest_probe = cc;
    cur_alpha_phase = squeeze(full_phasegrams(tr_inds_new,:,nearest_probe));
    Pmat = nan(length(tr_inds_new),length(wfreqs)*nbins);
    for ww = 1:length(wfreqs)
        cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        cur_set = ((ww-1)*nbins+1):ww*nbins;
        Pmat(:,cur_set) = cur_tb;
    end
    
    Xmat = [Pmat trial_Tmat(tr_inds_new,:) stim_out_interp];
    stim_params = [nbins,length(wfreqs),ntents,1];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_dphase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params_s,1,NL_type);
    stim_ind_dphase_pfilt(cc,:) = fitp_dphase.k(1:nbins*length(wfreqs));
    stim_ind_dphase_tfilt(cc,:) = fitp_dphase.k((nbins*length(wfreqs)+1):(nbins*length(wfreqs)+ntents));
    
    Xmat = [Pmat stim_out_interp];
    stim_params = [nbins,length(wfreqs)];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_phase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
    stim_ind_phase_pfilt(cc,:) = fitp_phase.k(1:nbins*length(wfreqs));
    stim_phase_g(cc) = fitp_phase.k(end-1);
    
    use_elecs = 1:24;
    use_elecs(cc) = [];
    use_set = find(ismember(phase_elec_set,use_elecs)); 
    Xmat = [new_phase_set(tr_inds_new,use_set) stim_out_interp];
    stim_params = [length(wfreqs),length(use_elecs)];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_fphase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params2,0,NL_type);
    stim_ind_sinphase_cfilt(cc,:) = fitp_fphase.k(1:length(use_elecs)*length(wfreqs));
    stim_ind_sinphase_sfilt(cc,:) = fitp_fphase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
    stim_fphase_g(cc) = fitp_fphase.k(end-1);
    
    
    xv_Robs = full_spkbinned(xv_inds_new,cc);
    stim_out_interp = stim_out(full_old_inds(xv_inds_new));

    Xmat = [stim_out_interp];
    xv_so_pred_rate = Xmat*fitp_so.k(1:end-1) + fitp_so.k(end);
    if NL_type == 1
    xv_so_pred_rate = exp(xv_so_pred_rate);
    else
        xv_so_pred_rate = log(1+exp(xv_so_pred_rate));
    end
    xv_so_LL(cc) = -sum(xv_Robs.*log(xv_so_pred_rate)-xv_so_pred_rate)/sum(xv_Robs);
%     
    Xmat = [trial_Tmat(xv_inds_new,:) stim_out_interp];
    xv_sac_pred_rate = Xmat*fitp_sac.k(1:end-1) + fitp_sac.k(end);
    if NL_type == 1
        xv_sac_pred_rate = exp(xv_sac_pred_rate);
    else
        xv_sac_pred_rate = log(1+exp(xv_sac_pred_rate));
    end
    xv_sac_LL(cc) = -sum(xv_Robs.*log(xv_sac_pred_rate)-xv_sac_pred_rate)/sum(xv_Robs);

    cur_alpha_phase = squeeze(full_phasegrams(xv_inds_new,:,nearest_probe));
    Pmat = nan(length(xv_inds_new),length(wfreqs)*nbins);
    for ww = 1:length(wfreqs)
        cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        cur_set = ((ww-1)*nbins+1):ww*nbins;
        Pmat(:,cur_set) = cur_tb;
    end

    Xmat = [Pmat stim_out_interp];
    xv_phase_pred_rate = Xmat*fitp_phase.k(1:end-1) + fitp_phase.k(end);
    if NL_type == 0
    xv_phase_pred_rate = log(1+exp(xv_phase_pred_rate));
    else
    xv_phase_pred_rate = exp(xv_phase_pred_rate);
    end
    xv_phase_LL(cc) = -sum(xv_Robs.*log(xv_phase_pred_rate)-xv_phase_pred_rate)/sum(xv_Robs);
    
        Xmat = [Pmat trial_Tmat(xv_inds_new,:) stim_out_interp];
    xv_dphase_pred_rate = Xmat*fitp_dphase.k(1:end-1) + fitp_dphase.k(end);
    if NL_type == 0
    xv_dphase_pred_rate = log(1+exp(xv_dphase_pred_rate));
    else
    xv_dphase_pred_rate = exp(xv_dphase_pred_rate);
    end
    xv_dphase_LL(cc) = -sum(xv_Robs.*log(xv_dphase_pred_rate)-xv_dphase_pred_rate)/sum(xv_Robs);
    
    Xmat = [new_phase_set(xv_inds_new,use_set) stim_out_interp];
    xv_fphase_pred_rate = Xmat*fitp_fphase.k(1:end-1) + fitp_fphase.k(end);
    if NL_type == 0
    xv_fphase_pred_rate = log(1+exp(xv_fphase_pred_rate));
    else
    xv_fphase_pred_rate = exp(xv_fphase_pred_rate);
    end
    xv_fphase_LL(cc) = -sum(xv_Robs.*log(xv_fphase_pred_rate)-xv_fphase_pred_rate)/sum(xv_Robs);
%     
%     
    avg_rate = mean(Robs);
    null_pred = avg_rate*ones(size(xv_Robs));
    xv_null_LL(cc) = -sum(xv_Robs.*log(null_pred)-null_pred)/sum(xv_Robs);

end
%%
xv_so_imp = xv_null_LL - xv_so_LL;
xv_phase_imp = xv_so_LL - xv_phase_LL;
xv_fphase_imp = xv_so_LL - xv_fphase_LL;
xv_sacp_imp = xv_phase_LL - xv_dphase_LL;
xv_sac_imp = xv_so_LL - xv_sac_LL;

boxplot([xv_so_imp(:) xv_phase_imp(:) xv_fphase_imp(:) xv_sac_imp(:)],{'stim','phase','depth-phase','sac'})
%%
save bar_phase_models_ftime_v2 stim_* tent_* new_dt xv_* wfreqs pax nbins 
%%
close all
stim_ind_ampkern = sqrt(stim_ind_sinphase_cfilt.^2+stim_ind_sinphase_sfilt.^2);
stim_ind_phasekern = -atan2(stim_ind_sinphase_cfilt,stim_ind_sinphase_sfilt)+pi/2;
for cc = 1:24
    cc
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
% %     colorbar
    subplot(3,1,2)
    pcolor(wfreqs,1:23,reshape(stim_ind_ampkern(cc,:),length(wfreqs),23)');shading flat
    subplot(3,1,3)
    pcolor(wfreqs,1:23,reshape(stim_ind_phasekern(cc,:),length(wfreqs),23)');shading flat
    %     set(gca,'xscale','log')
%     colorbar
%     caxis([0 0.005])
%     subplot(2,2,4)
%     plot(tent_centers*new_dt,stim_ind_dphase_tfilt(cc,:))
%     hold on
%     plot(tent_centers*new_dt,stim_ind_sac_kern(cc,:),'r')
%     cc
    fprintf('Phase imp: %.3f\n',xv_phase_imp(cc));
    fprintf('Depth imp: %.3f\n',xv_fphase_imp(cc));

    pause
    clf
end

%%
%%
close all
stim_ind_ampkern = sqrt(stim_ind_sinphase_cfilt.^2+stim_ind_sinphase_sfilt.^2);
stim_ind_phasekern = atan2(stim_ind_sinphase_sfilt,stim_ind_sinphase_cfilt);
for cc = 1:24
    temp = stim_ind_phase_pfilt(cc,:);
    temp = temp-mean(temp);
    subplot(2,4,[1 2 5 6])
    pcolor(180/pi*pax(1:end-1),wfreqs,reshape(temp,nbins,length(wfreqs))');shading flat
        set(gca,'yscale','log')
        xlabel('Phase deg','fontsize',16)
        ylabel('Frequency (Hz)','fontsize',16)
        title('Phase kernel (Non-parametric')
%     caxis([-0.03 0.03])
    subplot(2,4,3)
    pcolor(wfreqs,1:24,reshape((stim_ind_ampkern(cc,:)),length(wfreqs),24)');shading flat
%     caxis([0 0.01])
        set(gca,'xscale','log')
        xlabel('Frequency (Hz)','fontsize',16)
        ylabel('Channel','fontsize',16)
        title('Amplitude')
%     colorbar
    subplot(2,4,4)
    pcolor(wfreqs,1:24,reshape(stim_ind_phasekern(cc,:),length(wfreqs),24)');shading flat
        set(gca,'xscale','log')
        xlabel('Frequency (Hz)','fontsize',16)
        ylabel('Channel','fontsize',16)
        title('Preferred phase')
    subplot(2,4,7)
    pcolor(wfreqs,1:24,reshape(stim_ind_sinphase_cfilt(cc,:),length(wfreqs),24)');shading flat
%      caxis([-0.0075 0.0075])
        set(gca,'xscale','log')
        xlabel('Frequency (Hz)','fontsize',16)
        ylabel('Channel','fontsize',16)
        title('Cosine kernel')
   subplot(2,4,8)
    pcolor(wfreqs,1:24,reshape(stim_ind_sinphase_sfilt(cc,:),length(wfreqs),24)');shading flat
%      caxis([-0.0075 0.0075])
        set(gca,'xscale','log')
        xlabel('Frequency (Hz)','fontsize',16)
        ylabel('Channel','fontsize',16)
        title('Sin kernel')
        
%         fillPage(gcf,'papersize',[20 10])
%         fname = sprintf('Unit_%d_phasemod2',cc);
%         print('-dpdf',fname);
% close
     cc
     pause
    clf
end


%%
close all
for cc = 1:24

  % Phase kernel
  phs = reshape(stim_ind_phasekern(cc,:),length(wfreqs),24)';
  phs = fliplr(flipud(phs));

  % Amplitude kernel
  amp = reshape(stim_ind_ampkern(cc,:),length(wfreqs),24)';
  amp = amp/max(amp(:));
  amp = fliplr(flipud(amp));

  combmap = zeros(24,length(wfreqs),3);

  % Phase between red and blue
  combmap(:,:,1) = amp.* (0.5+0.5*cos(phs));
  combmap(:,:,3) = amp.* (0.5-0.5*cos(phs));

  image(combmap)

  cc
  pause
  clf
end
%%

for cc = 1:24

  % Phase filter
  frqphs = reshape(stim_ind_phase_pfilt(cc,:),nbins,length(wfreqs))';
  % Make 2 cycles
  frqphs = [frqphs frqphs];
  frqphs = frqphs-min(frqphs(:));

  % Try stretching to same time scale by upsampling
  dt = 0.001; 
  Tmax = 0.27*2;
  NT = ceil(Tmax/dt);
  NF = length(wfreqs);

  stretcher = zeros(NF,NT);

  for f = 1:length(wfreqs)
    tres = 1/wfreqs(NF+1-f) / (nbins);
    bins = floor((0:NT-1)*(dt/tres))+1;
    bins = bins(bins <= (2*nbins));
    stretcher(NF+1-f,1:length(bins)) = frqphs(NF+1-f,bins);
  end

  imagesc(stretcher)
  cc
  pause
  clf
end

%%
close all
NL_type = 1; %exp 

reg_params.dl1_ind = 40;
reg_params.dl1_dep = 0;
reg_params.dl2_ind = 1000;
reg_params.dl2_dep = 0;
reg_params.dl2_freq_ind = 0;
reg_params.dl2_freq_dep = 0;
reg_params.dl2_time_ind = 0;
reg_params.dl2_time_dep = 0;
reg_params.l1_ind = 0;
reg_params.l1_dep = 0;
reg_params.is_phase = 1;

use_freqs = 20;

silent = 1;

nbins = 30;
pax = linspace(-pi,pi,nbins+1);
doub_pax = [pax(1:end-1) 2*pi+pax(1:end-1)];

new_phase_set = [reshape(cos(full_phasegrams(:,use_freqs,:)),length(full_taxis_new),length(use_freqs)*24) reshape(sin(full_phasegrams(:,use_freqs,:)),length(full_taxis_new),length(use_freqs)*24)];
phase_elec_set = [repmat(1:24,1,length(use_freqs)) repmat(1:24,1,length(use_freqs))];
phase_elec_set = phase_elec_set(:);
%
for cc = 1:length(use_sus)
    fprintf('Cell %d of %d\n',cc,length(use_sus));
    
    Robs = full_spkbinned(tr_inds_new,cc);
    tr_spkbns = convert_to_spikebins(Robs);
    
    
    glm_kern = get_k_mat(glm_fit(cc));
    stim_out = stim_Xmat*glm_kern;
    stim_out_interp = stim_out(full_old_inds(tr_inds_new));
    
    
    if cc < 24
        nearest_probe = cc+1;
    else
        nearest_probe = cc-1;
    end
    cur_alpha_phase = squeeze(full_phasegrams(tr_inds_new,use_freqs,nearest_probe));
    Pmat = nan(length(tr_inds_new),length(use_freqs)*nbins);
    for ww = 1:length(use_freqs)
        cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        cur_set = ((ww-1)*nbins+1):ww*nbins;
        Pmat(:,cur_set) = cur_tb;
    end
       
    Xmat = [Pmat];
    stim_params = [nbins,length(use_freqs)];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_phase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
    phase_k(cc,:) = fitp_phase.k(1:nbins*length(use_freqs));
    
    use_elecs = nearest_probe;
    use_set = find(ismember(phase_elec_set,use_elecs)); 
    Xmat = [new_phase_set(tr_inds_new,use_set)];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_fphase,grad] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], [],[], [], [], [], NL_type);
    fphase_kern(cc,:) = fitp_fphase.k(1:end-1);
    
        xv_Robs = full_spkbinned(xv_inds_new,cc);
    
    cur_alpha_phase = squeeze(full_phasegrams(xv_inds_new,use_freqs,nearest_probe));
    Pmat = nan(length(xv_inds_new),length(use_freqs)*nbins);
    for ww = 1:length(use_freqs)
        cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        cur_set = ((ww-1)*nbins+1):ww*nbins;
        Pmat(:,cur_set) = cur_tb;
    end

    Xmat = [Pmat];
    xv_phase_pred_rate = Xmat*fitp_phase.k(1:end-1) + fitp_phase.k(end);
    if NL_type == 0
    xv_phase_pred_rate = log(1+exp(xv_phase_pred_rate));
    else
    xv_phase_pred_rate = exp(xv_phase_pred_rate);
    end
    xv_phase_LL(cc) = -sum(xv_Robs.*log(xv_phase_pred_rate)-xv_phase_pred_rate)/sum(xv_Robs);
        
    Xmat = [new_phase_set(xv_inds_new,use_set)];
    xv_fphase_pred_rate = Xmat*fitp_fphase.k(1:end-1) + fitp_fphase.k(end);
    if NL_type == 0
    xv_fphase_pred_rate = log(1+exp(xv_fphase_pred_rate));
    else
    xv_fphase_pred_rate = exp(xv_fphase_pred_rate);
    end
    xv_fphase_LL(cc) = -sum(xv_Robs.*log(xv_fphase_pred_rate)-xv_fphase_pred_rate)/sum(xv_Robs);
    avg_rate = mean(Robs);
    null_pred = avg_rate*ones(size(xv_Robs));
    xv_null_LL(cc) = -sum(xv_Robs.*log(null_pred)-null_pred)/sum(xv_Robs);

end
