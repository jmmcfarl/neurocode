clear all
close all

cd ~/Data/bruce/2_27_12/M232
load ./random_bar_eyedata_ftime.mat

load ./lemM232Expts.mat

load ./CellList.mat
good_sus = find(all(CellList(bar_expts,:,1) > 0));
use_sus = 1:24;

%%
Fs = 1000; %assume sample rate exactly 1 kHz for LFP!
dsf = 3;
Fsd = Fs/dsf;
niqf = Fsd/2;

maxlag = round(Fsd*0.5); %time range around stimulus to compute triggered average
maxlag_spk = round(Fsd*0.5); %time range around stimulus to compute triggered average

%bandpass filter for LFP
bb_lcf = 1;
[b_bb,a_bb] = butter(2,[bb_lcf]/niqf,'high');

all_lfps = [];
all_lfp_t_axis = [];
all_expt = [];
all_binned_spks = [];
for ee = 1:length(bar_expts);
    fprintf('Expt %d of %d\n',ee,length(bar_expts));
    fname = sprintf('lemM232A.%d.lfp.mat',bar_expts(ee));
    load(fname);
    fname = sprintf('Expt%dClusterTimes.mat',bar_expts(ee));
    load(fname);
    
    %     Fs = 1/LFP.Header.CRsamplerate;
    
    n_trials(ee) = length(LFP.Trials);
    %     lfp_trial_starts = [LFP.Trials(:).Start]/1e4;
    lfp_trial_starts = [LFP.Trials(:).ftime]/1e4;
    lfp_trial_ends = [LFP.Trials(:).End]/1e4;
    
    expt_lfp_t_axis = [];
    expt_lfps = [];
    expt_binned_spks = [];
    cur_npts = nan(n_trials(ee),1);
    cur_t_end = nan(n_trials(ee),1);
    for tt = 1:n_trials(ee)
        
        cur_npts(tt) = size(LFP.Trials(tt).LFP,1);
        cur_t_end(tt) = lfp_trial_starts(tt)+(cur_npts(tt)-1)/Fs;
        cur_t_axis = lfp_trial_starts(tt):1/Fs:cur_t_end(tt);
        
        if ~isempty(expt_lfp_t_axis)
            cur_sp = find(cur_t_axis > max(expt_lfp_t_axis),1,'first');
        else
            cur_sp = 1;
        end
        cur_t_axis = cur_t_axis(cur_sp:end);
        cur_LFP = [LFP.Trials(tt).LFP];
        cur_LFP = cur_LFP(cur_sp:end,:);
        cur_LFP_d = [];
        for ch = 1:24
            cur_LFP_d = [cur_LFP_d; decimate(cur_LFP(:,ch),dsf)'];
        end
        cur_LFP_d = filtfilt(b_bb,a_bb,cur_LFP_d);
        
        cur_t_axis_d = downsample(cur_t_axis,dsf);
        expt_lfp_t_axis = [expt_lfp_t_axis; cur_t_axis_d(:)];
        expt_lfps = [expt_lfps; cur_LFP_d'];
        
        cur_binned_spks = nan(length(cur_t_axis_d),length(use_sus));
        for cc = 1:length(use_sus)
            cur_spikes = Clusters{use_sus(cc)}.times;
            cur_spikes = cur_spikes(cur_spikes >= cur_t_axis_d(1) & cur_spikes < cur_t_axis_d(end));
            cur_binned_spks(:,cc) = hist(cur_spikes,cur_t_axis_d);
        end
        expt_binned_spks = [expt_binned_spks; cur_binned_spks];
        
    end
    all_lfps = [all_lfps; expt_lfps];
    all_lfp_t_axis = [all_lfp_t_axis; expt_lfp_t_axis];
    all_expt = [all_expt; ee*ones(length(expt_lfp_t_axis),1)];
    all_binned_spks = [all_binned_spks; expt_binned_spks];
end

%%
load ./un_bar_pos.mat
un_bar_pos(1:2) = [];
flen = 20;

all_trial_start = [];
all_trial_stop = [];
all_trial_expt = [];
all_stim_starts = [];
all_stim_start_times = [];
all_stim_trials = [];
all_bar_locs = [];
all_stim_class =[];
all_rel_start_times = [];
all_bar_Xmat = [];
all_stim_expts = [];
for ee = 1:length(bar_expts)
    expt_trial_starts = [Expts{bar_expts(ee)}.Trials(:).TrialStart]/1e4;
    expt_trial_end = [Expts{bar_expts(ee)}.Trials(:).TrueEnd]/1e4;
    
    %     buff = round(0.2*100);
    buff = 1;
    cur_expt_start_times = [];
    cur_rel_times = [];
    cur_expt_bar_locs = [];
    cur_bar_Xmat = [];
    cur_stim_trial_num = [];
    for t = 1:length(expt_trial_starts)
        cur_expt_start_times = [cur_expt_start_times; Expts{bar_expts(ee)}.Trials(t).Start];
        cur_rel_times = [cur_rel_times; Expts{bar_expts(ee)}.Trials(t).Start/1e4 - expt_trial_starts(t)];
        temp_bar_locs = Expts{bar_expts(ee)}.Trials(t).Op;
        cur_expt_bar_locs = [cur_expt_bar_locs; temp_bar_locs];
        cur_stim_trial_num = [cur_stim_trial_num; ones(length(temp_bar_locs),1)*t];
        
        cur_bar_mat = zeros(length(temp_bar_locs),length(un_bar_pos));
        for bb = 1:length(un_bar_pos)
            cur_set = find(temp_bar_locs==un_bar_pos(bb));
            cur_bar_mat(cur_set,bb) = 1;
        end
        cur_bar_Xmat = [cur_bar_Xmat; makeStimRows(cur_bar_mat,flen)];
    end
    cur_expt_start_times = cur_expt_start_times/1e4;
    
    %     sac_start_times = all_eye_ts{ee}(all_sac_startinds{ee});
    %     sac_stop_times = all_eye_ts{ee}(all_sac_stopinds{ee});
    %
    %     cur_sac_set = find(all_expt_num==ee);
    %     if length(cur_sac_set) ~= length(sac_start_times)
    %         error('saccade mismatch')
    %     end
    %     intrial_sac_set = find(~isnan(all_trial_num(cur_sac_set)));
    %     infirst_sac_set = find(ismember(cur_sac_set,all_isfirst));
    %     insecond_sac_set = find(ismember(cur_sac_set,all_issecond));
    
    cur_lfp_set = find(all_expt == ee);
    
    stim_start_inds = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),cur_expt_start_times));
    rid_stims = find(stim_start_inds < maxlag | stim_start_inds > length(cur_lfp_set) - maxlag);
    rid_stims = [rid_stims 1+find(diff(cur_expt_start_times) == 0)];
    stim_start_inds(rid_stims) = []; cur_expt_bar_locs(rid_stims) = []; 
    cur_rel_times(rid_stims) = []; cur_expt_start_times(rid_stims) = [];
    cur_stim_trial_num(rid_stims) = [];
    cur_bar_Xmat(rid_stims,:) = [];
    
    stim_class = zeros(size(cur_expt_start_times));
    stim_class(cur_rel_times < 0.8) = 1;
    stim_class(cur_rel_times > 0.8 & cur_rel_times < 1.5) = 2;
    
    all_stim_class = [all_stim_class(:); stim_class];
    all_stim_starts = [all_stim_starts(:); cur_lfp_set(stim_start_inds)];
    all_stim_start_times = [all_stim_start_times; cur_expt_start_times];
    all_stim_trials = [all_stim_trials; cur_stim_trial_num];
    all_stim_expts = [all_stim_expts(:); ones(length(stim_start_inds),1)*ee];
    all_bar_locs = [all_bar_locs(:); cur_expt_bar_locs];
    all_trial_start = [all_trial_start expt_trial_starts];
    all_trial_stop = [all_trial_stop expt_trial_end];
    all_trial_expt = [all_trial_expt ones(1,length(expt_trial_starts))*ee];
    all_rel_start_times = [all_rel_start_times; cur_rel_times];
    all_bar_Xmat = [all_bar_Xmat; cur_bar_Xmat];
end

%%
spk_sm = round(Fsd*0.01);
sm_binned_spks = zeros(size(all_binned_spks));
for i = 1:length(use_sus)
    sm_binned_spks(:,i) = jmm_smooth_1d_cor(all_binned_spks(:,i),spk_sm);
    
end

%%
load ./random_bar_unit_models_ftime.mat
all_pred_rate = [];
for ee = 1:length(bar_expts)
    cur_set_lfp = find(all_expt==ee);
    cur_set = find(all_stim_expts==ee);
    pred_rate_set = nan(length(cur_set_lfp),length(use_sus));
    for i = 1:length(use_sus)
        [~, ~, ~, prate] = getLL_GNM(glm_fit(i),all_bar_Xmat(cur_set,:),[],'none');
        interp_rate = interp1(all_stim_start_times(cur_set),prate,all_lfp_t_axis(cur_set_lfp));
        pred_rate_set(:,i) = interp_rate;
    end
    all_pred_rate = [all_pred_rate; pred_rate_set];
end
%%
lags = (-maxlag:maxlag);
lags_spk = (-maxlag_spk:maxlag_spk);

st_bounds = [0.3 1.7];

[un_bar_pos,ia,ic] = unique(all_bar_locs);

stim_trig_avgs = nan(length(un_bar_pos),length(lags),24);
stim_trig_avg_spks = nan(length(un_bar_pos),length(lags_spk),length(use_sus));
stim_trig_avg_pred = nan(length(un_bar_pos),length(lags_spk),length(use_sus));
for i = 1:length(un_bar_pos)
    fprintf('Bar pos %d of %d\n',i,length(un_bar_pos));
    cur_barset = find(all_bar_locs==un_bar_pos(i) & all_rel_start_times >= st_bounds(1) & all_rel_start_times <= st_bounds(2));
    for j = 1:length(lags)
        stim_trig_avgs(i,j,:) =  mean(all_lfps(all_stim_starts(cur_barset)+lags(j),:));
    end
    for j = 1:length(lags_spk)
        stim_trig_avg_spks(i,j,:) =  mean(sm_binned_spks(all_stim_starts(cur_barset)+lags_spk(j),:));
        stim_trig_avg_pred(i,j,:) =  nanmean(all_pred_rate(all_stim_starts(cur_barset)+lags_spk(j),:));
    end
end
ms_stim_trig_avgs = bsxfun(@minus,stim_trig_avgs,mean(stim_trig_avgs));

%%
clear lin_pred
for cur_ch = 1:24;
    lin_pred(cur_ch,:) = zeros(1,length(all_lfp_t_axis));
    for i = 1:length(un_bar_pos)
        cur_lin_filt = squeeze(ms_stim_trig_avgs(i,:,cur_ch));
        %                    cur_lin_filt = cur_lin_filt - mean(cur_lin_filt);
        cur_barset = find(all_bar_locs==un_bar_pos(i));
        null_sig = zeros(size(all_lfp_t_axis));
        null_sig(all_stim_starts(cur_barset)) = 1;
        cur_pred = conv(null_sig,cur_lin_filt,'same');
        lin_pred(cur_ch,:) = lin_pred(cur_ch,:) + cur_pred';
    end
end
lin_pred = lin_pred';
lin_res = all_lfps-lin_pred;

%% COMPUTE OBSERVED AND PREDICTED PHASE AND AMPLITUDE SPECTRA
% scales = logspace(log10(5.5),log10(80),25);
% % scales = [scales];
% wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
scales = logspace(log10(5),log10(60),30);
scales = [scales 70 80 90];
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

all_res_phasegram = [];
all_res_ampgram = [];
all_pred_phasegram = [];
all_pred_ampgram = [];
for ee = 1:length(bar_expts)
    fprintf('Computing CWT Expt %d of %d\n',ee,length(bar_expts));
    cur_lfp_set = find(all_expt==ee);
    
    cur_phasegram = nan(length(cur_lfp_set),length(wfreqs),24);
    cur_ampgram = nan(length(cur_lfp_set),length(wfreqs),24);
    pred_phasegram = nan(length(cur_lfp_set),length(wfreqs),24);
    pred_ampgram = nan(length(cur_lfp_set),length(wfreqs),24);
    for ll = 1:24
        temp = cwt(lin_res(cur_lfp_set,ll),scales,'cmor1-1');
        cur_phasegram(:,:,ll) = angle(temp)';
        cur_ampgram(:,:,ll) = abs(temp)';
        temp = cwt(lin_pred(cur_lfp_set,ll),scales,'cmor1-1');
        pred_phasegram(:,:,ll) = angle(temp)';
        pred_ampgram(:,:,ll) = abs(temp)';
    end
    
    all_res_phasegram = cat(1,all_res_phasegram,cur_phasegram);
    all_res_ampgram = cat(1,all_res_ampgram,cur_ampgram);
    all_pred_phasegram = cat(1,all_pred_phasegram,pred_phasegram);
    all_pred_ampgram = cat(1,all_pred_ampgram,pred_ampgram);
    
end
all_res_n_ampgram = bsxfun(@rdivide,all_res_ampgram,std(all_res_ampgram));
all_pred_n_ampgram = bsxfun(@rdivide,all_pred_ampgram,std(all_pred_ampgram));

%% COMPUTE USED INDS
buffer = round(0.2*Fsd);
use_inds = [];
interp_inds = [];
for ee = 1:length(bar_expts)
   
    cur_lfp_set = find(all_expt==ee);
    cur_trial_set = find(all_trial_expt==ee);
    trial_start_inds = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),all_trial_start(cur_trial_set)));
    trial_stop_inds = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),all_trial_stop(cur_trial_set)));
    bad = find(isnan(trial_start_inds) | isnan(trial_stop_inds));
    trial_start_inds(bad) = []; trial_stop_inds(bad) = [];
    cur_use_inds = [];
    for i = 1:length(trial_start_inds)
        cur_use_inds = [cur_use_inds (trial_start_inds(i)+buffer):(trial_stop_inds(i)-buffer)];
    end
    use_inds = [use_inds; cur_lfp_set(cur_use_inds(:))];
    
    cur_set = find(all_stim_expts==ee);
    cur_int = round(interp1(all_stim_start_times(cur_set),1:length(cur_set),all_lfp_t_axis(cur_lfp_set(cur_use_inds))));
    %     cur_int = round(interp1(all_lfp_t_axis(cur_lfp_set),1:length(cur_lfp_set),all_stim_start_times(cur_set)));
    uu = find(~isnan(cur_int));
    temp = nan(size(cur_int));
    temp(uu) = cur_set(cur_int(uu));
    interp_inds = [interp_inds; temp(:)];
end

use_inds(isnan(interp_inds)) = [];
interp_inds(isnan(interp_inds)) = [];

[c,ia,ic] = unique([all_stim_trials all_stim_expts],'rows');
n_trials = length(ia);

xv_frac = 0.2;
n_xv_trials = round(n_trials*xv_frac);
xv_set = randperm(n_trials);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(ic,xv_set));
tr_inds = find(~ismember(ic,xv_set))';

xv_inds_new = find(ismember(interp_inds,xv_inds));
tr_inds_new = find(ismember(interp_inds,tr_inds));

%%
NL_type = 1;
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


res_phase_set = [reshape(cos(all_res_phasegram(use_inds,:,:)),length(use_inds),length(wfreqs)*24) reshape(sin(all_res_phasegram(use_inds,:,:)),length(use_inds),length(wfreqs)*24)];
pred_phase_set = [reshape(cos(all_pred_phasegram(use_inds,:,:)),length(use_inds),length(wfreqs)*24) reshape(sin(all_pred_phasegram(use_inds,:,:)),length(use_inds),length(wfreqs)*24)];

res_ampphase_set = [reshape(all_res_n_ampgram(use_inds,:,:),length(use_inds),length(wfreqs)*24).*reshape(cos(all_res_phasegram(use_inds,:,:)),length(use_inds),length(wfreqs)*24) ...
    reshape(all_res_n_ampgram(use_inds,:,:),length(use_inds),length(wfreqs)*24).*reshape(sin(all_res_phasegram(use_inds,:,:)),length(use_inds),length(wfreqs)*24)];
pred_ampphase_set = [reshape(all_pred_n_ampgram(use_inds,:,:),length(use_inds),length(wfreqs)*24).*reshape(cos(all_pred_phasegram(use_inds,:,:)),length(use_inds),length(wfreqs)*24) ...
    reshape(all_pred_n_ampgram(use_inds,:,:),length(use_inds),length(wfreqs)*24).*reshape(sin(all_pred_phasegram(use_inds,:,:)),length(use_inds),length(wfreqs)*24)];

% phase_elec_set = [repmat(1:24,1,length(wfreqs)) repmat(1:24,1,length(wfreqs))];
% phase_elec_set = phase_elec_set(:);
phase_elec_set = ones(length(wfreqs),1)*(1:24);
phase_elec_set = [phase_elec_set(:); phase_elec_set(:)];

%%
silent = 1;
for cc = 1:24
    fprintf('Cell %d of %d\n',cc,24);
    Robs = all_binned_spks(use_inds(tr_inds_new),cc);
    tr_spkbns = convert_to_spikebins(Robs);
   
    glm_kern = get_k_mat(glm_fit(cc));
    stim_out = all_bar_Xmat*glm_kern;
    stim_out_interp = stim_out(interp_inds(tr_inds_new));
    cur_mean = mean(stim_out_interp);
    cur_std = std(stim_out_interp);
    stim_out_interp = (stim_out_interp - cur_mean)/cur_std;
    
    Xmat = [stim_out_interp];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_so,grad] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], [],[], [], [], [], NL_type);

    use_elecs = 1:24;
%     use_elecs(cc) = [];
    use_set = find(ismember(phase_elec_set,use_elecs));
        
    Xmat = [res_ampphase_set(tr_inds_new,use_set) stim_out_interp];
    stim_params = [length(wfreqs),length(use_elecs)];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_raphase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
    res_ampphase_cfilt(cc,:) = fitp_raphase.k(1:length(use_elecs)*length(wfreqs));
    res_ampphase_sfilt(cc,:) = fitp_raphase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
    res_ampphase_g(cc) = fitp_raphase.k(end-1);

        Xmat = [pred_ampphase_set(tr_inds_new,use_set) stim_out_interp];
    stim_params = [length(wfreqs),length(use_elecs)];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_paphase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
    pred_ampphase_cfilt(cc,:) = fitp_paphase.k(1:length(use_elecs)*length(wfreqs));
    pred_ampphase_sfilt(cc,:) = fitp_paphase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
    pred_ampphase_g(cc) = fitp_paphase.k(end-1);

%     Xmat = [bsxfun(@times,res_ampphase_set(tr_inds_new,use_set),stim_out_interp) phasemod_out stim_out_interp];
%     stim_params = [length(wfreqs),length(use_elecs)];
%     klen = size(Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_raphase_f] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
%     res_ampphasef_cfilt(cc,:) = fitp_raphase_f.k(1:length(use_elecs)*length(wfreqs));
%     res_ampphasef_sfilt(cc,:) = fitp_raphase_f.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);

%     Xmat = [bsxfun(@times,res_ampphase_set(tr_inds_new,use_set),stim_out_interp) stim_out_interp];
%     stim_params = [length(wfreqs),length(use_elecs)];
%     klen = size(Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_raphase_m] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
%     res_ampphasem_cfilt(cc,:) = fitp_raphase_m.k(1:length(use_elecs)*length(wfreqs));
%     res_ampphasem_sfilt(cc,:) = fitp_raphase_m.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
    
    %%
    xv_Robs = all_binned_spks(use_inds(xv_inds_new),cc);
    stim_out_interp = stim_out(interp_inds(xv_inds_new));
    stim_out_interp = (stim_out_interp - cur_mean)/cur_std;

    Xmat = [stim_out_interp];
    xv_so_pred_rate = Xmat*fitp_so.k(1:end-1) + fitp_so.k(end);
    if NL_type == 1
        xv_so_pred_rate = exp(xv_so_pred_rate);
    else
        xv_so_pred_rate = log(1+exp(xv_so_pred_rate));
    end
    xv_so_LL(cc) = -sum(xv_Robs.*log(xv_so_pred_rate)-xv_so_pred_rate)/sum(xv_Robs);
           
    
    Xmat = [res_ampphase_set(xv_inds_new,use_set) stim_out_interp];
    xv_raphase_pred_rate = Xmat*fitp_raphase.k(1:end-1) + fitp_raphase.k(end);
    if NL_type == 0
        xv_raphase_pred_rate = log(1+exp(xv_raphase_pred_rate));
    else
        xv_raphase_pred_rate = exp(xv_raphase_pred_rate);
    end
    xv_raphase_LL(cc) = -sum(xv_Robs.*log(xv_raphase_pred_rate)-xv_raphase_pred_rate)/sum(xv_Robs);

        
    Xmat = [pred_ampphase_set(xv_inds_new,use_set) stim_out_interp];
    xv_raphase_pred_rate = Xmat*fitp_paphase.k(1:end-1) + fitp_paphase.k(end);
    if NL_type == 0
        xv_raphase_pred_rate = log(1+exp(xv_raphase_pred_rate));
    else
        xv_raphase_pred_rate = exp(xv_raphase_pred_rate);
    end
    xv_paphase_LL(cc) = -sum(xv_Robs.*log(xv_raphase_pred_rate)-xv_raphase_pred_rate)/sum(xv_Robs);

%     phasemod_out = res_ampphase_set(xv_inds_new,use_set)*phase_kern';
%     Xmat = [bsxfun(@times,res_ampphase_set(xv_inds_new,use_set),stim_out_interp) phasemod_out stim_out_interp];
%     xv_raphase_pred_rate = Xmat*fitp_raphase_f.k(1:end-1) + fitp_raphase_f.k(end);
%     if NL_type == 0
%         xv_raphase_pred_rate = log(1+exp(xv_raphase_pred_rate));
%     else
%         xv_raphase_pred_rate = exp(xv_raphase_pred_rate);
%     end
%     xv_raphasef_LL(cc) = -sum(xv_Robs.*log(xv_raphase_pred_rate)-xv_raphase_pred_rate)/sum(xv_Robs);
   
%     phasemod_out = res_ampphase_set(xv_inds_new,use_set)*phase_kern';
%     Xmat = [bsxfun(@times,res_ampphase_set(xv_inds_new,use_set),stim_out_interp) stim_out_interp];
%     xv_raphase_pred_rate = Xmat*fitp_raphase_m.k(1:end-1) + fitp_raphase_m.k(end);
%     if NL_type == 0
%         xv_raphase_pred_rate = log(1+exp(xv_raphase_pred_rate));
%     else
%         xv_raphase_pred_rate = exp(xv_raphase_pred_rate);
%     end
%     xv_raphasem_LL(cc) = -sum(xv_Robs.*log(xv_raphase_pred_rate)-xv_raphase_pred_rate)/sum(xv_Robs);

    avg_rate = mean(Robs);
    null_pred = avg_rate*ones(size(xv_Robs));
    xv_null_LL(cc) = -sum(xv_Robs.*log(null_pred)-null_pred)/sum(xv_Robs);


end

%%
cd ~/Data/bruce/2_27_12/M232
save random_bar_spec_mods *filt xv_*
%%
use_c = 1:24;
xv_so_imp = xv_null_LL(use_c) - xv_so_LL(use_c);
xv_raphase_imp = xv_so_LL(use_c) - xv_raphase_LL(use_c);
xv_paphase_imp = xv_so_LL(use_c) - xv_paphase_LL(use_c);
% xv_raphasem_imp = xv_so_LL(use_c) - xv_raphasem_LL(use_c);
% xv_raphasef_imp = xv_so_LL(use_c) - xv_raphasef_LL(use_c);

%%
res_ampphase_ampkern = sqrt(res_ampphase_cfilt.^2+res_ampphase_sfilt.^2);
res_ampphase_phasekern = -atan2(res_ampphase_cfilt,res_ampphase_sfilt)+pi/2;
pred_ampphase_ampkern = sqrt(pred_ampphase_cfilt.^2+pred_ampphase_sfilt.^2);
pred_ampphase_phasekern = -atan2(pred_ampphase_cfilt,pred_ampphase_sfilt)+pi/2;

close all
for cc = 1:24
    subplot(2,2,1)
    pcolor(wfreqs,1:24,reshape(res_ampphase_ampkern(cc,:),length(wfreqs),24)');shading flat
%     set(gca,'xscale','log')
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Channel','fontsize',16)
    colorbar
%     title('Cosine kernel')
    subplot(2,2,3)
    pcolor(wfreqs,1:24,reshape(pred_ampphase_ampkern(cc,:),length(wfreqs),24)');shading flat
%     set(gca,'xscale','log')
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Channel','fontsize',16)
%     title('Sin kernel')
    colorbar

cc
    pause
    clf
end

%%
cd ~/Data/bruce/2_27_12/M232
load random_bar_spec_mods
res_ampphase_ampkern = sqrt(res_ampphase_cfilt.^2+res_ampphase_sfilt.^2);
pred_ampphase_ampkern = sqrt(pred_ampphase_cfilt.^2+pred_ampphase_sfilt.^2);
cd ~/Data/bruce/2_27_12
load ./free_viewing_phase_models_late_v3
fv_ampkern = sqrt(sinphase_cfilt.^2 + sinphase_sfilt.^2);
cd ~/Data/bruce/2_27_12/
load Blocks.mat
cellids = [1 2 3 4 5 6 7 8 10];
muaids = [1 2 4 5 6 11 12 13 14];
su_probes = [Blocks{1}.suprobes(cellids)];
mu_probes = [Blocks{1}.muprobes(muaids)];
all_probes = [su_probes mu_probes];

cd ~/Data/bruce/2_27_12/M232
load ./random_bar_phase_models_v2
rb_ampkern = sqrt(sinphase_cfilt.^2 + sinphase_sfilt.^2);

close all
for cc = 1:24
    subplot(2,2,1)
    cm = max(abs(res_ampphase_ampkern(cc,:)));
    pcolor(wfreqs,1:24,reshape(res_ampphase_ampkern(cc,:),length(wfreqs),24)');shading flat
%     set(gca,'xscale','log')
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Channel','fontsize',16)
    caxis([0 cm]*0.8)
    colorbar
%     title('Cosine kernel')
    subplot(2,2,3)
    pcolor(wfreqs,1:24,reshape(pred_ampphase_ampkern(cc,:),length(wfreqs),24)');shading flat
%     set(gca,'xscale','log')
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Channel','fontsize',16)
%     title('Sin kernel')
    caxis([0 cm]*0.8)
    colorbar

    subplot(2,2,2)
    pcolor(wfreqs,1:24,reshape(rb_ampkern(cc,:),length(wfreqs),24)');shading flat
%     set(gca,'xscale','log')
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Channel','fontsize',16)
        caxis([0 cm]*0.8)
colorbar

    %     title('Cosine kernel')
    cur_cell = find(all_probes==cc);
    if ~isempty(cur_cell)
    subplot(2,2,4)
pcolor(wfreqs,1:24,reshape(fv_ampkern(cur_cell,:),length(wfreqs),24)');shading flat
%     set(gca,'xscale','log')
        caxis([0 cm]*0.8)
    xlabel('Frequency (Hz)','fontsize',16)
    ylabel('Channel','fontsize',16)
%     title('Sin kernel')
    colorbar
    end

cc
    pause
    clf
end
%%
reg_params.dl1_ind = 0;
reg_params.dl1_dep = 0;
reg_params.dl2_ind = 500;
reg_params.dl2_dep = 500;
reg_params.dl2_freq_ind = 0;
reg_params.dl2_freq_dep = 0;
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


nbins = 30;
pax = linspace(-pi,pi,nbins+1);
doub_pax = [pax(1:end-1) 2*pi+pax(1:end-1)];

for cc = 1:24
        fprintf('Cell %d of %d\n',cc,24);
    Robs = all_binned_spks(use_inds(tr_inds_new),cc);
    tr_spkbns = convert_to_spikebins(Robs);

    use_elecs =cc;
    use_set = find(ismember(phase_elec_set,use_elecs));

    Xmat = [res_phase_set(tr_inds_new,use_set)];
    stim_params = [length(wfreqs),length(use_elecs)];
    klen = size(Xmat,2);
    K0 = randn(klen+1,1);
    [fitp_raphase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
    cfilt(cc,:) = fitp_raphase.k(1:length(use_elecs)*length(wfreqs));
    sfilt(cc,:) = fitp_raphase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
    
    
    cur_phase = squeeze(all_res_phasegram(use_inds(tr_inds_new),:,cc));
    Pmat = nan(length(tr_inds_new),length(wfreqs)*nbins);
    for ww = 1:length(wfreqs)
        cur_tb = tbrep(cur_phase(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        cur_set = ((ww-1)*nbins+1):ww*nbins;
        Pmat(:,cur_set) = cur_tb;
    end
    stim_params = [nbins,length(wfreqs)];
    klen = size(Pmat,2);
    K0 = zeros(klen+1,1);
    [fitp_phase] = fit_GLM_phase_model(Pmat, Robs, K0,silent, stim_params, reg_params2,0,NL_type);
    phase_pfilt(cc,:) = fitp_phase.k(1:nbins*length(wfreqs));

%     full_b = glmfit(Xmat,Robs,'poisson');
%     full_amp(cc,:) = sqrt(full_b(2:26).^2 + full_b(27:end).^2);
% 
%     phase_lock(cc,:) = zeros(length(wfreqs),1);
%     for use_freq = 1:length(wfreqs)
%         use_phase = squeeze(all_res_phasegram(tr_inds_new,use_freq,cc));
%         phase_lock(cc,use_freq) = circ_kappa(use_phase(tr_spkbns));
%     end

    for ww = 1:length(wfreqs)
        use_phase = squeeze(all_res_phasegram(use_inds(tr_inds_new),ww,cc));
        Xmat = [cos(use_phase) sin(use_phase)];
        b(ww,:) = glmfit(Xmat,Robs,'poisson');
    end
    ind_amp(cc,:) = sqrt(b(:,2).^2 + b(:,3).^2);
    
    xv_Robs = all_binned_spks(use_inds(xv_inds_new),cc);
    
    Xmat = [res_phase_set(xv_inds_new,use_set)];
    xv_phase_pred_rate = Xmat*fitp_raphase.k(1:end-1) + fitp_raphase.k(end);
    if NL_type == 0
        xv_phase_pred_rate = log(1+exp(xv_phase_pred_rate));
    else
        xv_phase_pred_rate = exp(xv_phase_pred_rate);
    end
    cur_xvLL(cc) = -sum(xv_Robs.*log(xv_phase_pred_rate)-xv_phase_pred_rate)/sum(xv_Robs);

    cur_phase = squeeze(all_res_phasegram(use_inds(xv_inds_new),:,cc));
    Pmat = nan(length(xv_inds_new),length(wfreqs)*nbins);
    for ww = 1:length(wfreqs)
        cur_tb = tbrep(cur_phase(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        cur_set = ((ww-1)*nbins+1):ww*nbins;
        Pmat(:,cur_set) = cur_tb;
    end
    xv_phase_pred_rate = Pmat*fitp_phase.k(1:end-1) + fitp_phase.k(end);
    if NL_type == 0
        xv_phase_pred_rate = log(1+exp(xv_phase_pred_rate));
    else
        xv_phase_pred_rate = exp(xv_phase_pred_rate);
    end
    cur_pxvLL(cc) = -sum(xv_Robs.*log(xv_phase_pred_rate)-xv_phase_pred_rate)/sum(xv_Robs);

    for ww = 1:length(wfreqs)
        use_phase = squeeze(all_res_phasegram(use_inds(xv_inds_new),ww,cc));
        Xmat = [cos(use_phase) sin(use_phase)];
        cur_out = Xmat*b(ww,2:end)'+b(ww,1);
        cur_out = exp(cur_out);
        ind_xvLL(cc,ww) = -sum(xv_Robs.*log(cur_out)-cur_out)/sum(xv_Robs);
    end
    ind_amp(cc,:) = sqrt(b(:,2).^2 + b(:,3).^2);
    
end

ampkern = sqrt(cfilt.^2+sfilt.^2);
%%
n_prc_bins = 10;
prc_bin_edges = linspace(0,1,n_prc_bins+1);
prc_bin_edges = 1+floor(prc_bin_edges*length(tr_inds_new));
prc_bin_edges(prc_bin_edges > length(tr_inds_new)) = length(tr_inds_new);

n_phase_bins = 8;
phase_bin_edges = linspace(-pi,pi,n_phase_bins+1);

cc=6;
% for cc = 1:24
fprintf('Cell %d of %d\n',cc,24);
Robs = all_binned_spks(use_inds(tr_inds_new),cc);
tr_spkbns = convert_to_spikebins(Robs);

glm_kern = get_k_mat(glm_fit(cc));
stim_out = all_bar_Xmat*glm_kern;
stim_out_interp = stim_out(interp_inds(tr_inds_new));
cur_mean = mean(stim_out_interp);
cur_std = std(stim_out_interp);
stim_out_interp = (stim_out_interp - cur_mean)/cur_std;

phase_rate = zeros(length(wfreqs),n_phase_bins);
phase_lock = zeros(length(wfreqs),1);
for use_freq = 1:length(wfreqs)
   use_phase = squeeze(all_res_phasegram(tr_inds_new,use_freq,cc));
   for pp = 1:n_phase_bins
        phase_set = find(use_phase >= phase_bin_edges(pp) & use_phase < phase_bin_edges(pp+1));
        phase_rate(use_freq,pp) = mean(Robs(phase_set));
   end
   phase_lock(use_freq) = circ_kappa(use_phase(tr_spkbns));
end







use_freq = 25;
use_phase = squeeze(all_res_phasegram(tr_inds_new,use_freq,cc));
[~,sort_ord] = sort(stim_out_interp);
np_phase_stim = zeros(n_prc_bins,n_phase_bins);
for i = 1:n_prc_bins
   cur_set = sort_ord(prc_bin_edges(i):prc_bin_edges(i+1)); 
   cur_phase = use_phase(cur_set);
   for pp = 1:n_phase_bins
      phase_set = find(cur_phase >= phase_bin_edges(pp) & cur_phase < phase_bin_edges(pp+1));
      np_phase_stim(i,pp) = mean(Robs(cur_set(phase_set)));
   end
end

Xmat = [stim_out_interp];
klen = size(Xmat,2);
K0 = zeros(klen+1,1);
[fitp_so,grad] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], [],[], [], [], [], NL_type);

use_elecs = 1:24;
%     use_elecs(cc) = [];
use_set = find(ismember(phase_elec_set,use_elecs));

%     Xmat = [res_ampphase_set(tr_inds_new,use_set) stim_out_interp];
%     stim_params = [length(wfreqs),length(use_elecs)];
%     klen = size(Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_raphase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
%     res_ampphase_cfilt(cc,:) = fitp_raphase.k(1:length(use_elecs)*length(wfreqs));
%     res_ampphase_sfilt(cc,:) = fitp_raphase.k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
%     res_ampphase_g(cc) = fitp_raphase.k(end-1);

phase_kern = [res_ampphase_cfilt(cc,:) res_ampphase_sfilt(cc,:)];
phasemod_out = res_ampphase_set(tr_inds_new,use_set)*phase_kern';


% end