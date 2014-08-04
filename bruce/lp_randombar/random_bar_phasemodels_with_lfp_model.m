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
dsf = 5;
Fsd = Fs/dsf;
niqf = Fs/2;
[bb,aa] = butter(2,[1 80]/niqf);

new_dt = .005;
dsfrac = 2;

%%
load ./CellList.mat
good_sus = find(all(CellList(bar_expts,:,1) > 0));
use_sus = 1:24;
%%
load ./un_bar_pos.mat
un_bar_pos(1:2) = [];

flen = 40;
full_Xmat = [];
full_spkbinned = [];
full_exptvec = [];
full_exptvec_new = [];
full_trialvec_new = [];
full_taxis = [];
full_taxis_new = [];
full_old_t_inds = [];
full_bar_pos = [];
full_used_inds = [];
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
    all_bar_Xmat = [];
    all_binned_spikes = [];
    all_used_inds = [];
    for tt = 1:n_trials(ee)
        
        cur_t_axis = [Expts{bar_expts(ee)}.Trials(tt).Start]/1e4;
        cur_bar_Op = [Expts{bar_expts(ee)}.Trials(tt).Op];
        
        cur_bar_Op_new = repmat(cur_bar_Op,1,dsfrac)';
        cur_bar_Op_new = cur_bar_Op_new(:);
        cur_bar_Op_new = cur_bar_Op_new(2:end-1);
        
        cur_t_edges = [cur_t_axis; Expts{bar_expts(ee)}.Trials(tt).End(end)/1e4];
        cur_t_edges_new = cur_t_edges(1):new_dt:cur_t_edges(end);
        cur_t_axis_new = cur_t_edges_new(1:end-1);
        
        cur_binned_spks = nan(length(cur_t_axis_new),length(use_sus));
        for cc = 1:length(use_sus)
            cur_hist = histc(Clusters{use_sus(cc)}.times,cur_t_edges_new);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
        end
        
        cur_bar_mat = zeros(length(cur_t_axis_new),length(un_bar_pos));
        for b = 1:length(un_bar_pos)
            cur_set = find(cur_bar_Op_new==un_bar_pos(b));
            cur_bar_mat(cur_set,b) = 1;
        end
        bar_Xmat = makeStimRows(cur_bar_mat,flen);
        
        %         cur_used_inds = ones(length(cur_t_axis),1);
        %         cur_used_inds(1:flen) = 0;
        cur_used_inds_new = ones(length(cur_t_axis_new),1);
        cur_used_inds_new(1:flen) = 0;
        
        
        cur_used_inds_new(end-round(0.2/new_dt):end) = 0; %reduce phase estimation edge artifacts
        
        all_used_inds_new = [all_used_inds_new; cur_used_inds_new(:)];
        all_t_axis = [all_t_axis; cur_t_axis(:)];
        %         all_used_inds = [all_used_inds; cur_used_inds(:)];
        all_t_axis_new = [all_t_axis_new; cur_t_axis_new(:)];
        all_binned_spikes = [all_binned_spikes; cur_binned_spks];
        all_expt_vec = [all_expt_vec; ones(length(cur_t_axis_new),1)*bar_expts(ee)];
        all_trial_vec = [all_trial_vec; ones(length(cur_t_axis_new),1)*tt];
        all_bar_Op = [all_bar_Op; cur_bar_Op_new(:)];
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
    
    %     all_used_inds(all_blink_inds == 1) = 0;
    %     all_used_inds = logical(all_used_inds);
    all_used_inds_new(all_blink_inds_new == 1) = 0;
    all_used_inds_new = logical(all_used_inds_new);
    
    full_used_inds = [full_used_inds; all_used_inds_new];
    full_Xmat = [full_Xmat; all_bar_Xmat];
    full_spkbinned = [full_spkbinned; all_binned_spikes];
    full_exptvec = [full_exptvec; ones(size(all_used_inds_new))*ee];
    full_exptvec_new = [full_exptvec_new; ones(size(all_used_inds_new))*ee];
%     full_taxis = [full_taxis; all_t_axis];
    full_taxis_new = [full_taxis_new; all_t_axis_new];
    %     full_old_t_inds = [full_old_t_inds; all_old_t_inds(all_used_inds_new)];
    full_trialvec_new = [full_trialvec_new; all_trial_vec];
    full_bar_pos = [full_bar_pos; all_bar_Op];
end

%%
full_lfps = [];
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
    res_lfps = [];
    for jj = 1:length(seg_start_inds)
        cur_res_inds = seg_start_inds(jj):seg_stop_inds(jj);
        cur_res_t = expt_lfp_t_axis(seg_start_inds(jj)):dsf/Fs:expt_lfp_t_axis(seg_stop_inds(jj));
        res_t_axis = [res_t_axis; cur_res_t(:)];
        
        interp_lfps = interp1(expt_lfp_t_axis(cur_res_inds),expt_lfps(cur_res_inds,:),cur_res_t);
        res_lfps = [res_lfps; interp_lfps];
    end
    
    cur_set = find(full_exptvec_new==ee);
    interp_lfps = interp1(res_t_axis,res_lfps,full_taxis_new(cur_set));
    
    full_lfps = [full_lfps; interp_lfps];
end
% full_lfps = zscore(full_lfps);

%%
lin_X = zeros(length(full_taxis_new),length(bar_expts)-1);
for i = 1:length(bar_expts)-1
    cur_set = find(full_exptvec_new==i);
    lin_X(cur_set,i) = 1;
end

[c,ia,ic] = unique([full_exptvec_new full_trialvec_new],'rows');
n_trials = length(ia);

rp = randperm(n_trials);
rand_trial_vec = rp(ic);
[~,ind_shuff] = sort(rand_trial_vec);

xv_frac = 0.2;
n_xv_trials = round(n_trials*xv_frac);
xv_tset = randperm(n_trials);
xv_tset(n_xv_trials+1:end) = [];
tr_tset = find(~ismember(1:n_trials,xv_tset));
xv_inds = find(ismember(ic,xv_tset));
tr_inds = find(ismember(ic,tr_tset));

tr_inds(full_used_inds(tr_inds) ~= 1) = [];
xv_inds(full_used_inds(tr_inds) ~= 1) = [];
full_inds = sort([tr_inds; xv_inds]);
%%
n_bar_pos = length(un_bar_pos);
stim_params.spatial_dims = 1;
stim_params.sdim = n_bar_pos;
stim_params.flen = flen;
klen = n_bar_pos*flen;
clear defmod

stim_pred = full_lfps;
lfp_resid = full_lfps;
for cc = 1:24
    fprintf('Cell %d of %d\n',cc,24);
    Yobs = full_lfps(:,cc);
    defmod.lambda_L1x = 0;
    defmod.lambda_d2T = 100;
    defmod.lambda_d2X = 2;
    defmod.lambda_d2XT = 10;
    
    kern_types{1} = 'lin';
    init_kerns = 0.01*randn(klen,1);
    glm = createGNM(init_kerns,1,kern_types,defmod,stim_params,'gauss');
    lfp_glm_fit(cc) = fitGNM_filters(glm,full_Xmat(tr_inds,:),Yobs(tr_inds),'none',[],1e-4,1e-6,0);
    
    
    if xv_frac > 0
        [glm_sse, ~, ~, prate] = getLL_GNM(lfp_glm_fit(cc),full_Xmat(xv_inds,:),Yobs(xv_inds),'none');
        yobs_var = var(Yobs(xv_inds));
        xv_glm_r2(cc) = 1-glm_sse/yobs_var;
    end
    
    [glm_sse, ~, ~, stim_pred(full_inds,cc)] = getLL_GNM(lfp_glm_fit(cc),full_Xmat(full_inds,:),Yobs(full_inds),'none');
    lfp_resid(full_inds,cc) = Yobs(full_inds) - stim_pred(full_inds,cc);
end

%%
cd /home/james/Data/bruce/2_27_12/M232
save lfp_models *_inds lfp_resid lfp_glm_fit stim_pred
%%
params.tapers = [4 7];
params.Fs = Fsd;
win = 2;
params.fpass = [2 70];
params.trialave = 1;
for cc = 1:24
    Yobs = full_lfps(full_inds,cc);
    [Cmn(cc,:),~,~,~,~,f] = coherencysegc(Yobs,stim_pred(full_inds,cc),win,params);
end

%%
scales = logspace(log10(4),log10(100),30);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

all_phasegrams = nan(length(full_taxis_new),length(wfreqs),24);
all_ampgrams = nan(length(full_taxis_new),length(wfreqs),24);
all_rphasegrams = nan(length(full_taxis_new),length(wfreqs),24);
all_rampgrams = nan(length(full_taxis_new),length(wfreqs),24);
for ii = 1:n_trials
    ii
    cur_trial_inds = find(ic == ii);
    for cc = 1:24
        temp = cwt(stim_pred(cur_trial_inds,cc),scales,'cmor1-1');
        all_phasegrams(cur_trial_inds,:,cc) = angle(temp)';
        all_ampgrams(cur_trial_inds,:,cc) = abs(temp)';
        temp = cwt(lfp_resid(cur_trial_inds,cc),scales,'cmor1-1');
        all_rphasegrams(cur_trial_inds,:,cc) = angle(temp)';
        all_rampgrams(cur_trial_inds,:,cc) = abs(temp)';
    end
    
end

all_ampgrams = bsxfun(@rdivide,all_ampgrams,std(all_ampgrams));
all_rampgrams = bsxfun(@rdivide,all_rampgrams,std(all_rampgrams));

%%
cd /home/james/Data/bruce/2_27_12/M232
load ./random_bar_phasemods_stimmods

%%
fNT = length(full_taxis_new);
new_phase_set = [reshape(all_ampgrams,fNT,length(wfreqs)*24).*reshape(cos(all_phasegrams),fNT,length(wfreqs)*24) ...
    reshape(all_ampgrams,fNT,length(wfreqs)*24).*reshape(sin(all_phasegrams),fNT,length(wfreqs)*24)];
new_rphase_set = [reshape(all_rampgrams,fNT,length(wfreqs)*24).*reshape(cos(all_rphasegrams),fNT,length(wfreqs)*24) ...
    reshape(all_rampgrams,fNT,length(wfreqs)*24).*reshape(sin(all_rphasegrams),fNT,length(wfreqs)*24)];


phase_elec_set = ones(length(wfreqs),1)*(1:24);
phase_elec_set = [phase_elec_set(:); phase_elec_set(:)];

use_elecs = 1:24;
use_set = find(ismember(phase_elec_set,use_elecs));
%%
NL_type = 0; %exp

reg_params.dl2_freq = 20; %20000
reg_params.dl2_ch =  400; %20000
reg_params.dl2_time = 100; %500
reg_params.dl_time = 20;
reg_params.dl_freq = 5; %100
reg_params.dl_ch =  200; %100
reg_params.d2_phase = 25;
reg_params.d2_time = 1;

silent = 1;

%%
for cc = 1:24
    fprintf('Cell %d of %d\n',cc,24);
    
    Robs = full_spkbinned(tr_inds,cc);
    tr_spkbns = convert_to_spikebins(Robs);
    
    glm_kern = get_k_mat(glm_fit(cc));
    glm_mat = reshape(glm_kern,30,n_bar_pos);
    glm_mat = [zeros(flen-30,n_bar_pos); glm_mat];
    glm_kern = glm_mat(:);
    
    stim_out = full_Xmat*glm_kern;
    stim_out_interp = stim_out(tr_inds);
    
    cur_Xmat = [lin_X(tr_inds,:)];
    klen = size(cur_Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_null,grad] = GLMsolve_jmm(cur_Xmat, tr_spkbns, K0, 1, [], [],[], [], [], [], NL_type);
    
    cur_Xmat = [stim_out_interp lin_X(tr_inds,:)];
    klen = size(cur_Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_so,grad] = GLMsolve_jmm(cur_Xmat, tr_spkbns, K0, 1, [], [],[], [], [], [], NL_type);
    
    use_elecs = 1:24;
    %     if mod(cc,2)==0
    %         use_elecs = 2:2:24;
    %     else
    %         use_elecs = 1:2:24;
    %     end
    use_set = find(ismember(phase_elec_set,use_elecs));
        
    cur_Xmat = [new_phase_set(tr_inds,use_set) stim_out_interp lin_X(tr_inds,:)];
    stim_params = [length(wfreqs),length(use_elecs)];
    klen = size(cur_Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_post(cc)] = fit_GLM_sinphase_model(cur_Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
    post_sinphase_cfilt(cc,:) = fitp_post(cc).k(1:length(use_elecs)*length(wfreqs));
    post_sinphase_sfilt(cc,:) = fitp_post(cc).k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
    post_ampkern(cc,:) = sqrt(post_sinphase_cfilt(cc,:).^2 + post_sinphase_sfilt(cc,:).^2);
    
    cur_Xmat = [new_rphase_set(tr_inds,use_set) stim_out_interp lin_X(tr_inds,:)];
    stim_params = [length(wfreqs),length(use_elecs)];
    klen = size(cur_Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_rpost(cc)] = fit_GLM_sinphase_model(cur_Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
    rpost_sinphase_cfilt(cc,:) = fitp_rpost(cc).k(1:length(use_elecs)*length(wfreqs));
    rpost_sinphase_sfilt(cc,:) = fitp_rpost(cc).k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
    rpost_ampkern(cc,:) = sqrt(rpost_sinphase_cfilt(cc,:).^2 + rpost_sinphase_sfilt(cc,:).^2);
    
    cur_Xmat = [new_phase_set(tr_inds,use_set) lin_X(tr_inds,:)];
    stim_params = [length(wfreqs),length(use_elecs)];
    klen = size(cur_Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_po(cc)] = fit_GLM_sinphase_model(cur_Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
    po_sinphase_cfilt(cc,:) = fitp_po(cc).k(1:length(use_elecs)*length(wfreqs));
    po_sinphase_sfilt(cc,:) = fitp_po(cc).k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
    po_ampkern(cc,:) = sqrt(po_sinphase_cfilt(cc,:).^2 + po_sinphase_sfilt(cc,:).^2);
    
    cur_Xmat = [new_rphase_set(tr_inds,use_set) lin_X(tr_inds,:)];
    stim_params = [length(wfreqs),length(use_elecs)];
    klen = size(cur_Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_rpo(cc)] = fit_GLM_sinphase_model(cur_Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
    rpo_sinphase_cfilt(cc,:) = fitp_rpo(cc).k(1:length(use_elecs)*length(wfreqs));
    rpo_sinphase_sfilt(cc,:) = fitp_rpo(cc).k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
    rpo_ampkern(cc,:) = sqrt(rpo_sinphase_cfilt(cc,:).^2 + rpo_sinphase_sfilt(cc,:).^2);
      
    xv_Robs = full_spkbinned(xv_inds,cc);
    stim_out_interp = stim_out(xv_inds);
    
    cur_Xmat = [stim_out_interp lin_X(xv_inds,:)];
    xv_so_pred_rate = cur_Xmat*fitp_so.k(1:end-1) + fitp_so.k(end);
    if NL_type == 1
        xv_so_pred_rate = exp(xv_so_pred_rate);
    else
        xv_so_pred_rate = log(1+exp(xv_so_pred_rate));
    end
    xv_so_LL(cc) = -sum(xv_Robs.*log(xv_so_pred_rate)-xv_so_pred_rate)/sum(xv_Robs);
    
    cur_Xmat = [new_phase_set(xv_inds,use_set) stim_out_interp lin_X(xv_inds,:)];
    xv_fphase_pred_rate = cur_Xmat*fitp_post(cc).k(1:end-1) + fitp_post(cc).k(end);
    if NL_type == 0
        xv_fphase_pred_rate = log(1+exp(xv_fphase_pred_rate));
    else
        xv_fphase_pred_rate = exp(xv_fphase_pred_rate);
    end
    xv_post_LL(cc) = -sum(xv_Robs.*log(xv_fphase_pred_rate)-xv_fphase_pred_rate)/sum(xv_Robs);
    
    cur_Xmat = [new_rphase_set(xv_inds,use_set) stim_out_interp lin_X(xv_inds,:)];
    xv_fphase_pred_rate = cur_Xmat*fitp_rpost(cc).k(1:end-1) + fitp_rpost(cc).k(end);
    if NL_type == 0
        xv_fphase_pred_rate = log(1+exp(xv_fphase_pred_rate));
    else
        xv_fphase_pred_rate = exp(xv_fphase_pred_rate);
    end
    xv_rpost_LL(cc) = -sum(xv_Robs.*log(xv_fphase_pred_rate)-xv_fphase_pred_rate)/sum(xv_Robs);
    
    cur_Xmat = [new_phase_set(xv_inds,use_set) lin_X(xv_inds,:)];
    xv_fphase_pred_rate = cur_Xmat*fitp_po(cc).k(1:end-1) + fitp_po(cc).k(end);
    if NL_type == 0
        xv_fphase_pred_rate = log(1+exp(xv_fphase_pred_rate));
    else
        xv_fphase_pred_rate = exp(xv_fphase_pred_rate);
    end
    xv_po_LL(cc) = -sum(xv_Robs.*log(xv_fphase_pred_rate)-xv_fphase_pred_rate)/sum(xv_Robs);
    
    cur_Xmat = [new_rphase_set(xv_inds,use_set) lin_X(xv_inds,:)];
    xv_fphase_pred_rate = cur_Xmat*fitp_rpo(cc).k(1:end-1) + fitp_rpo(cc).k(end);
    if NL_type == 0
        xv_fphase_pred_rate = log(1+exp(xv_fphase_pred_rate));
    else
        xv_fphase_pred_rate = exp(xv_fphase_pred_rate);
    end
    xv_rpo_LL(cc) = -sum(xv_Robs.*log(xv_fphase_pred_rate)-xv_fphase_pred_rate)/sum(xv_Robs);

    cur_Xmat = [lin_X(xv_inds,:)];
    xv_null_pred_rate = cur_Xmat*fitp_null.k(1:end-1) + fitp_null.k(end);
    if NL_type == 0
        xv_null_pred_rate = log(1+exp(xv_null_pred_rate));
    else
        xv_null_pred_rate = exp(xv_null_pred_rate);
    end
    xv_null_LL(cc) = -sum(xv_Robs.*log(xv_null_pred_rate)-xv_null_pred_rate)/sum(xv_Robs);
    
end

%%
save pred_lfp_phasemods xv_* wfreqs fit* *kern *filt xv_inds tr_inds
%%
stim_imp = (xv_null_LL - xv_so_LL)/log(2);
phase_imp = (xv_null_LL - xv_po_LL)/log(2);
rphase_imp = (xv_null_LL - xv_rpo_LL)/log(2);
stimphase_imp = (xv_null_LL - xv_post_LL)/log(2);
rstimphase_imp = (xv_null_LL - xv_rpost_LL)/log(2);

sus = good_sus;
mus = setdiff(1:24,sus);
mus(mus==6) = [];
figure
subplot(2,1,1)
boxplot([stim_imp(sus)' phase_imp(sus)' rphase_imp(sus)' stimphase_imp(sus)' rstimphase_imp(sus)'],...
    {'stim','stim-dep LFP','stim-ind LFP','stim + stim-dep LFP','stim + stim-ind LFP'});
ylim([0 1.1])
ylabel('LL improvement (bits/spk)','fontsize',16)
title('SUs')
subplot(2,1,2)
boxplot([stim_imp(mus)' phase_imp(mus)' rphase_imp(mus)' stimphase_imp(mus)' rstimphase_imp(mus)'],...
    {'stim','stim-dep LFP','stim-ind LFP','stim + stim-dep LFP','stim + stim-ind LFP'});
ylim([0 1.1])
ylabel('LL improvement (bits/spk)','fontsize',16)
title('MUs')

%%
close all
for cc = 1:24
    
    subplot(2,1,1)
    pcolor(wfreqs,1:25,[reshape(po_ampkern(cc,:)',length(wfreqs),24)'; zeros(1,length(wfreqs))]);shading interp; set(gca,'xscale','log');
    colorbar
     subplot(2,1,2)
    pcolor(wfreqs,1:25,[reshape(rpo_ampkern(cc,:)',length(wfreqs),24)'; zeros(1,length(wfreqs))]);shading interp; set(gca,'xscale','log');
    colorbar
   pause
    clf
end

%%
for cc = 1:24
    w_kern(cc,:,:) = reshape(po_ampkern(cc,:)',length(wfreqs),24);
    w_kern2(cc,:,:) = reshape(rpo_ampkern(cc,:)',length(wfreqs),24);
end
%%
full_eye_speed = [];
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
    big_sac_set = [infirst_sac_set; insecond_sac_set];
    micro_sac_set = find(sac_amps' < 60 & ~isnan(all_trial_num(cur_sac_set)));
    micro_sac_set(ismember(micro_sac_set,big_sac_set)) = [];

    sac_times = sac_start_times(intrial_sac_set);
    first_sac_times = sac_start_times(infirst_sac_set);
    second_sac_times = sac_start_times(insecond_sac_set);
    micro_sac_times = sac_start_times(micro_sac_set);
    
    cur_e_set = find(full_exptvec_new == ee);
    
        interp_eye_speed = interp1(all_eye_ts{ee},all_eye_speed{ee},full_taxis_new(cur_e_set));
    full_eye_speed = [full_eye_speed; interp_eye_speed];
        
    sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),sac_times));
    sac_inds(isnan(sac_inds)) = [];
    micro_sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),micro_sac_times));
    micro_sac_inds(isnan(micro_sac_inds)) = [];
    first_sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),first_sac_times));
    first_sac_inds(isnan(first_sac_inds)) = [];
    second_sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),second_sac_times));
    second_sac_inds(isnan(second_sac_inds)) = [];
    
%     all_sac_inds = [all_sac_inds; cur_e_set(sac_inds(:))];
%     all_microsac_inds = [all_microsac_inds; cur_e_set(micro_sac_inds(:))];
%     all_firstsac_inds = [all_firstsac_inds; cur_e_set(first_sac_inds(:))];
%     all_secondsac_inds = [all_secondsac_inds; cur_e_set(second_sac_inds(:))];
end

% all_sac_inds = [];
% all_microsac_inds = [];
% all_firstsac_inds = [];
% all_secondsac_inds = [];
% for ee = 1:length(bar_expts)
%     expt_trial_starts = [Expts{bar_expts(ee)}.Trials(:).TrialStart]/1e4;
%     expt_trial_end = [Expts{bar_expts(ee)}.Trials(:).TrueEnd]/1e4;
%     
%     sac_start_times = all_eye_ts{ee}(all_sac_startinds{ee});
%     sac_stop_times = all_eye_ts{ee}(all_sac_stopinds{ee});
%     sac_amps = all_sac_peakamps{ee};
%     
%     cur_sac_set = find(all_expt_num==ee);
%     if length(cur_sac_set) ~= length(sac_start_times)
%         error('saccade mismatch')
%     end
%     intrial_sac_set = find(~isnan(all_trial_num(cur_sac_set)));
% 
%     infirst_sac_set = find(ismember(cur_sac_set,all_isfirst));
%     insecond_sac_set = find(ismember(cur_sac_set,all_issecond));
%     big_sac_set = [infirst_sac_set; insecond_sac_set];
%     micro_sac_set = find(sac_amps' < 60 & ~isnan(all_trial_num(cur_sac_set)));
%     micro_sac_set(ismember(micro_sac_set,big_sac_set)) = [];
% 
%     sac_times = sac_start_times(intrial_sac_set);
%     first_sac_times = sac_start_times(infirst_sac_set);
%     second_sac_times = sac_start_times(insecond_sac_set);
%     micro_sac_times = sac_start_times(micro_sac_set);
%     
%     cur_e_set = find(full_exptvec_new == ee);
%     
%     sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),sac_times));
%     sac_inds(isnan(sac_inds)) = [];
%     micro_sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),micro_sac_times));
%     micro_sac_inds(isnan(micro_sac_inds)) = [];
%     first_sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),first_sac_times));
%     first_sac_inds(isnan(first_sac_inds)) = [];
%     second_sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),second_sac_times));
%     second_sac_inds(isnan(second_sac_inds)) = [];
%     
%     all_sac_inds = [all_sac_inds; cur_e_set(sac_inds(:))];
%     all_microsac_inds = [all_microsac_inds; cur_e_set(micro_sac_inds(:))];
%     all_firstsac_inds = [all_firstsac_inds; cur_e_set(first_sac_inds(:))];
%     all_secondsac_inds = [all_secondsac_inds; cur_e_set(second_sac_inds(:))];
% end
% all_bigsac_inds = sort([all_firstsac_inds; all_secondsac_inds]);
% 
% %%
% cur_dt = 0.01;
% tent_centers = [0:cur_dt:0.7];
% tent_centers = round(tent_centers/new_dt);
% tbmat = construct_tent_bases(tent_centers,1);
% [ntents,tblen] = size(tbmat);
% 
% shift = round(0.2/new_dt);
% tbmat = [zeros(ntents,tblen-2*shift-1) tbmat];
% tent_centers = tent_centers-shift;
% 
% trial_sac_inds = zeros(size(full_taxis_new));
% trial_sac_inds(all_bigsac_inds) = 1;
% trial_sac_mat = zeros(length(full_taxis_new),ntents);
% for i = 1:ntents
%     trial_sac_mat(:,i) = conv(trial_sac_inds,tbmat(i,:),'same');
% end
% 
% trial_msac_inds = zeros(size(full_taxis_new));
% trial_msac_inds(all_microsac_inds) = 1;
% trial_msac_mat = zeros(length(full_taxis_new),ntents);
% for i = 1:ntents
%     trial_msac_mat(:,i) = conv(trial_msac_inds,tbmat(i,:),'same');
% end
% 
% %%
% for cc = 1:length(use_sus)
%     fprintf('Cell %d of %d\n',cc,length(use_sus));
%     
%     Robs = full_spkbinned(tr_inds,cc);
%     tr_spkbns = convert_to_spikebins(Robs);
%    
%     glm_kern = get_k_mat(glm_fit(cc));
%     glm_mat = reshape(glm_kern,30,n_bar_pos);
%     glm_mat = [zeros(flen-30,n_bar_pos); glm_mat];
%     glm_kern = glm_mat(:);
%     
%     stim_out = full_Xmat*glm_kern;
%     stim_out_interp = stim_out(tr_inds);
% 
%     cur_Xmat = [trial_sac_mat(tr_inds,:) trial_msac_mat(tr_inds,:) stim_out_interp lin_X(tr_inds,:)];
%     lamrange = [reg_params.dl_time 1 ntents 0; reg_params.dl_time ntents + 1 2*ntents 0];
%     lamrange2 = [reg_params.dl2_time 1 ntents 0; reg_params.dl2_time ntents+1 2*ntents 0];
%     klen = size(cur_Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_sac,grad] = GLMsolve_jmm(cur_Xmat, tr_spkbns, K0, silent, lamrange, lamrange2,[], [], [], [], NL_type);
%     so_sac_kern(cc,:) = fitp_sac.k(1:ntents);
%     so_msac_kern(cc,:) = fitp_sac.k(ntents+1:2*ntents);
% 
%     
%     xv_Robs = full_spkbinned(xv_inds,cc);
%     stim_out_interp = stim_out(xv_inds);
%     
%     cur_Xmat = [trial_sac_mat(xv_inds,:) trial_msac_mat(xv_inds,:) stim_out_interp lin_X(xv_inds,:)];
%     xv_so_pred_rate = cur_Xmat*fitp_sac.k(1:end-1) + fitp_sac.k(end);
%     if NL_type == 1
%         xv_so_pred_rate = exp(xv_so_pred_rate);
%     else
%         xv_so_pred_rate = log(1+exp(xv_so_pred_rate));
%     end
%     xv_sac_LL(cc) = -sum(xv_Robs.*log(xv_so_pred_rate)-xv_so_pred_rate)/sum(xv_Robs);
%     
% end
% 
%%
n_bar_pos = length(un_bar_pos);
use_bar_pos = full_bar_pos; 
use_bar_mat = zeros(length(full_taxis_new),n_bar_pos);
for b = 1:n_bar_pos
    cur_set = find(use_bar_pos == un_bar_pos(b));
    use_bar_mat(cur_set,b) = 1;
end
use_bar_pos(use_bar_pos < -10) = nan;

use_t_axis = (1:length(full_taxis_new))*new_dt;
%%
examp_cell = 14;

for i = 20:n_trials
    i
    cur_set = find(ic == i);
    cur_set(full_used_inds(cur_set) ~= 1) = [];
    
    subplot(4,1,1)
    imagesc(use_t_axis(cur_set),un_bar_pos,use_bar_mat(cur_set,:)');colormap(gray);
    xlim(use_t_axis(cur_set([1 end])))
    xlabel('Time (s)','fontsize',16)
    ylabel('Bar position (deg)','fontsize',16)
     set(gca,'fontsize',14)
    
    subplot(4,1,2)
    plot(use_t_axis(cur_set),full_lfps(cur_set,examp_cell),'k')
    xlim(use_t_axis(cur_set([1 end])))
    xlabel('Time (s)','fontsize',16)
    ylabel('LFP amplitude (z)','fontsize',16)
     set(gca,'fontsize',14)
    box off
   
    subplot(4,1,3)
    plot(use_t_axis(cur_set),full_lfps(cur_set,examp_cell),'k--')
    hold on
    plot(use_t_axis(cur_set),stim_pred(cur_set,examp_cell),'b')
     plot(use_t_axis(cur_set),lfp_resid(cur_set,examp_cell)-3,'r')
     xlim(use_t_axis(cur_set([1 end])))
    xlabel('Time (s)','fontsize',16)
    ylabel('LFP amplitude (z)','fontsize',16);
    set(gca,'fontsize',14)
    box off

    subplot(4,1,4)
        plot(use_t_axis(cur_set),full_eye_speed(cur_set),'k')
    xlim(use_t_axis(cur_set([1 end])))
     xlabel('Time (s)','fontsize',16)
    ylabel('Eye speed (deg/s)','fontsize',16)
    set(gca,'fontsize',14)
    box off

    pause
    clf
end