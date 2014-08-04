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
dsf = 4;
Fsd = Fs/dsf;
niqf = Fs/2;
[bb,aa] = butter(2,120/niqf,'low');
scales = logspace(log10(6),log10(60),25);
scales = [scales 70 80];
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

new_dt = .005;
dsfrac = 2;
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
use_sus = 1:24;
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
full_ampgrams = [];
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
    
    %compute CSD
    vars.Fs = Fs;
    vars.BrainBound = 1;
    vars.ChanSep = 0.05;
     vars.diam = 0.5; %0.5
    CSD = PettersenCSD(expt_lfps','spline',vars)';
    
    
    seg_start_inds = [1; 1+find(diff(expt_lfp_t_axis) > 2/Fs)];
    seg_stop_inds = [find(diff(expt_lfp_t_axis) > 2/Fs); length(expt_lfp_t_axis)];
    res_t_axis = [];
    res_phasegram = [];
    res_ampgram = [];
    res_lfps = [];
    for jj = 1:length(seg_start_inds)
        cur_res_inds = seg_start_inds(jj):seg_stop_inds(jj);
        cur_res_t = expt_lfp_t_axis(seg_start_inds(jj)):dsf/Fs:expt_lfp_t_axis(seg_stop_inds(jj));
        res_t_axis = [res_t_axis; cur_res_t(:)];
        
%         interp_lfps = interp1(expt_lfp_t_axis(cur_res_inds),expt_lfps(cur_res_inds,:),cur_res_t);
        interp_lfps = interp1(expt_lfp_t_axis(cur_res_inds),CSD(cur_res_inds,:),cur_res_t);
        
        cur_phasegram = nan(length(cur_res_t),length(wfreqs),24);
        cur_ampgram = nan(length(cur_res_t),length(wfreqs),24);
        for ll = 1:24
            temp = cwt(interp_lfps(:,ll),scales,'cmor1-1');
            cur_phasegram(:,:,ll) = angle(temp)';
            cur_ampgram(:,:,ll) = abs(temp)';
        end
        res_phasegram = cat(1,res_phasegram,cur_phasegram);
        res_ampgram = cat(1,res_ampgram,cur_ampgram);
        res_lfps = [res_lfps; interp_lfps];
    end
    
    cur_set = find(full_exptvec_new==ee);
    interp_lfps = interp1(res_t_axis,res_lfps,full_taxis_new(cur_set));
    unwr_phasegram = unwrap(res_phasegram);
    interp_phasegrams = interp1(res_t_axis,unwr_phasegram,full_taxis_new(cur_set));
    interp_phasegrams = mod(interp_phasegrams+pi,2*pi)-pi;
    %     interp_phasegrams = interp1(res_t_axis,res_phasegram,full_taxis_new(cur_set));
    interp_ampgrams = interp1(res_t_axis,res_ampgram,full_taxis_new(cur_set));
    
    full_lfps = [full_lfps; interp_lfps];
    full_phasegrams = cat(1,full_phasegrams, interp_phasegrams);
    full_ampgrams = cat(1,full_ampgrams, interp_ampgrams);
end

%%
full_ampgrams = bsxfun(@rdivide,full_ampgrams,std(full_ampgrams));

%%
all_sac_inds = [];
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

    sac_times = sac_start_times(intrial_sac_set);

    cur_e_set = find(full_exptvec_new == ee);
    
    sac_inds = round(interp1(full_taxis_new(cur_e_set),1:length(cur_e_set),sac_times));

    sac_inds(isnan(sac_inds)) = [];
    
    all_sac_inds = [all_sac_inds; cur_e_set(sac_inds(:))];
end

%%
cur_dt = 0.01;
% flen_t = 0.5;
tent_centers = [0:cur_dt:0.7];
tent_centers = round(tent_centers/new_dt);
tbmat = construct_tent_bases(tent_centers,1);
[ntents,tblen] = size(tbmat);

shift = round(0.2/new_dt);
tbmat = [zeros(ntents,tblen-2*shift-1) tbmat];
tent_centers = tent_centers-shift;

trial_inds = zeros(size(full_taxis_new));
trial_inds(all_sac_inds) = 1;
trial_Tmat = zeros(length(full_taxis_new),ntents);
for i = 1:ntents
    trial_Tmat(:,i) = conv(trial_inds,tbmat(i,:),'same');
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

load ./random_bar_unit_models_ftime.mat glm_fit

%%
[c,ia,ic] = unique([full_exptvec full_trialvec],'rows');
n_trials = length(ia);

xv_frac = 0.2;
n_xv_trials = round(n_trials*xv_frac);
xv_tset = randperm(n_trials);
xv_tset(n_xv_trials+1:end) = [];
tr_tset = find(~ismember(1:n_trials,xv_tset));
xv_inds = find(ismember(ic,xv_tset));
tr_inds = find(ismember(ic,tr_tset))';

xv_inds_new = find(ismember(full_old_inds,xv_inds));
tr_inds_new = find(ismember(full_old_inds,tr_inds));

stim_Xmat = full_Xmat;
%%
fNT = length(full_taxis_new);
% new_phase_set = [reshape(cos(full_phasegrams),length(full_taxis_new),length(wfreqs)*24) reshape(sin(full_phasegrams),length(full_taxis_new),length(wfreqs)*24)];
new_phase_set = [reshape(full_ampgrams,fNT,length(wfreqs)*24).*reshape(cos(full_phasegrams),fNT,length(wfreqs)*24) ...
    reshape(full_ampgrams,fNT,length(wfreqs)*24).*reshape(sin(full_phasegrams),fNT,length(wfreqs)*24)];
phase_elec_set = ones(length(wfreqs),1)*(1:24);
phase_elec_set = [phase_elec_set(:); phase_elec_set(:)];

use_elecs = 1:24;
use_set = find(ismember(phase_elec_set,use_elecs));
Xmat = [new_phase_set(tr_inds_new,use_set)];
Xmatxv = [new_phase_set(xv_inds_new,use_set)];

tr_sac_inds = find(trial_inds(tr_inds_new)==1);
%%
NL_type = 0; %exp

% reg_params.dl2_freq = 5000;
% reg_params.dl2_ch = 5000;
reg_params.dl2_freq = 5000; %20000
reg_params.dl2_ch =  10000; %20000
reg_params.dl2_time = 300; %500
reg_params.d2_phase = 10;
reg_params.d2_time = 10;

silent = 1;

%%
for cc = 1:length(use_sus)
fprintf('Cell %d of %d\n',cc,length(use_sus));
    
    Robs = full_spkbinned(tr_inds_new,cc);
    tr_spkbns = convert_to_spikebins(Robs);
    
    glm_kern = get_k_mat(glm_fit(cc));
    stim_out = stim_Xmat*glm_kern;
    stim_out_interp = stim_out(full_old_inds(tr_inds_new));
    
    cur_Xmat = [stim_out_interp];
    klen = size(cur_Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_so,grad] = GLMsolve_jmm(cur_Xmat, tr_spkbns, K0, 1, [], [],[], [], [], [], NL_type);

    cur_Xmat = [trial_Tmat(tr_inds_new,:) stim_out_interp];
    lamrange2 = [300 1 ntents 0];
    klen = size(cur_Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_sac,grad] = GLMsolve_jmm(cur_Xmat, tr_spkbns, K0, silent, [], lamrange2,[], [], [], [], NL_type);
    so_sac_kern(cc,:) = fitp_sac.k(1:ntents);
    predrate_out = cur_Xmat*fitp_sac.k(1:end-1) + fitp_sac.k(end);
    [st_avg_sacmod(cc,:),cur_lags] = get_event_trig_avg(exp(predrate_out),tr_sac_inds,round(0.2/new_dt),round(0.5/new_dt));
    
    
    cur_Xmat = [Xmat trial_Tmat(tr_inds_new,:) stim_out_interp];
    use_elecs = 1:24;
    use_set = find(ismember(phase_elec_set,use_elecs));
    stim_params = [length(wfreqs),length(use_elecs) ntents];
    klen = size(cur_Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_fphase] = fit_GLM_sinphase_model(cur_Xmat, Robs, K0,silent, stim_params, reg_params,1,NL_type);
    sinphase_cfilt(cc,:) = fitp_fphase.k(1:length(use_elecs)*length(wfreqs));
    sinphase_sfilt(cc,:) = fitp_fphase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
    ampkern(cc,:) = sqrt(sinphase_cfilt(cc,:).^2 + sinphase_sfilt(cc,:).^2);
    phase_sac_kern(cc,:) = fitp_fphase.k(length(use_elecs)*length(wfreqs)*2+1:length(use_elecs)*length(wfreqs)*2+ntents);

    phasemod_out = Xmat*fitp_fphase.k(1:2*length(use_elecs)*length(wfreqs));
    predrate_out = cur_Xmat*fitp_fphase.k(1:end-1) + fitp_fphase.k(end);
    [st_avg_phasemod(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,tr_sac_inds,round(0.2/new_dt),round(0.5/new_dt));
    [st_avg_allmod(cc,:),cur_lags] = get_event_trig_avg(exp(predrate_out),tr_sac_inds,round(0.2/new_dt),round(0.5/new_dt));
    
    xv_Robs = full_spkbinned(xv_inds_new,cc);
    stim_out_interp = stim_out(full_old_inds(xv_inds_new));
    
    xv_so_pred_rate = stim_out_interp*fitp_so.k(1:end-1) + fitp_so.k(end);
    if NL_type == 1
        xv_so_pred_rate = exp(xv_so_pred_rate);
    else
        xv_so_pred_rate = log(1+exp(xv_so_pred_rate));
    end
    xv_so_LL(cc) = -sum(xv_Robs.*log(xv_so_pred_rate)-xv_so_pred_rate)/sum(xv_Robs);
    
    cur_Xmat = [trial_Tmat(xv_inds_new,:) stim_out_interp];
    xv_so_pred_rate = cur_Xmat*fitp_sac.k(1:end-1) + fitp_so.k(end);
    if NL_type == 1
        xv_so_pred_rate = exp(xv_so_pred_rate);
    else
        xv_so_pred_rate = log(1+exp(xv_so_pred_rate));
    end
    xv_sac_LL(cc) = -sum(xv_Robs.*log(xv_so_pred_rate)-xv_so_pred_rate)/sum(xv_Robs);

    cur_Xmat = [Xmatxv trial_Tmat(xv_inds_new,:) stim_out_interp];
    xv_fphase_pred_rate = cur_Xmat*fitp_fphase.k(1:end-1) + fitp_fphase.k(end);
    if NL_type == 0
        xv_fphase_pred_rate = log(1+exp(xv_fphase_pred_rate));
    else
        xv_fphase_pred_rate = exp(xv_fphase_pred_rate);
    end
    xv_fphase_LL(cc) = -sum(xv_Robs.*log(xv_fphase_pred_rate)-xv_fphase_pred_rate)/sum(xv_Robs);

    avg_rate = mean(Robs);
    null_pred = avg_rate*ones(size(xv_Robs));
    xv_null_LL(cc) = -sum(xv_Robs.*log(null_pred)-null_pred)/sum(xv_Robs);
    
    fprintf('Stim only: %.4f\n',(xv_null_LL(cc) - xv_so_LL(cc))/log(2));
    fprintf('Sac + Stim: %.4f\n',(xv_null_LL(cc) - xv_sac_LL(cc))/log(2));
    fprintf('Stim + phase: %.4f\n',(xv_null_LL(cc) - xv_fphase_LL(cc))/log(2));
    
end

%%
cd ~/Data/bruce/2_27_12/M232
save random_bar_phase_models_csd xv_* *filt wfreqs ampkern reg_params

%%
close all
for cc = 1:24
    subplot(3,2,1);hold on
    plot(tent_centers*new_dt,so_sac_kern(cc,:),'r')
    plot(cur_lags*new_dt,st_avg_phasemod(cc,:),'b')
    plot(tent_centers*new_dt,phase_sac_kern(cc,:),'k')
    subplot(3,2,3)
    pcolor(wfreqs,1:24,reshape(ampkern(cc,:)',length(wfreqs),24)');shading interp; set(gca,'xscale','log');
    ca = max([max(abs(sinphase_cfilt(cc,:))) max(abs(sinphase_sfilt(cc,:)))]);
    subplot(3,2,5)
    pcolor(wfreqs,1:24,reshape(sinphase_cfilt(cc,:)',length(wfreqs),24)');shading interp; set(gca,'xscale','log');
    caxis(0.8*[-ca ca])
    subplot(3,2,6)
    pcolor(wfreqs,1:24,reshape(sinphase_sfilt(cc,:)',length(wfreqs),24)');shading interp; set(gca,'xscale','log');
    caxis(0.8*[-ca ca])
    cc
    pause
    clf
end