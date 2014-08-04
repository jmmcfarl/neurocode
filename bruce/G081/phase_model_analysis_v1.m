clear all
% close all
cd ~/Data/bruce/G081/
load jbeG081Expts.mat
load ./CellList.mat
good_sus = find(all(CellList(:,:,1) > 0));

stim_fs = 100; %in Hz
Fs = 1e3;
dsf = 6;Fsd = Fs/dsf;
niqf = Fs/2;
[filt_b,filt_a] = butter(2,[1 60]/niqf);
use_lfps = [1:3:96];
use_lfps = sort(unique([use_lfps(:); good_sus(:)]));
use_units = 1:96;

% use_lfps = [1 63];
backlag = round(0*Fsd);
forwardlag = round(4*Fsd);
lags = (-backlag:forwardlag);

% scales = logspace(log10(0.2),log10(4.5),25);
scales = logspace(log10(0.4),log10(7.5),35);
scales = scales*60/dsf;
% scales = scales*2;
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

beg_buffer = round(Fsd*0.2); %don't use the first X data after start of a trial.
end_buffer = round(Fsd*0.2);

%% PARSE TRIAL DATA STRUCTURES
for i = 1:length(Expts)
    if strcmp(Expts{i}.Header.expname,'grating.OpXseRC') | strcmp(Expts{i}.Header.expname,'grating.OpRC')
        is_bar_expt(i) = 1;
    else
        is_bar_expt(i) = 0;
    end
    
    if strcmp(Expts{i}.Stimvals.Bs,'image')
        expt_image_back(i) = 1;
    else
        expt_image_back(i) = 0;
    end
    
    expt_sim_sacs(i) = Expts{i}.Stimvals.ijump;
    expt_bar_ori(i) = Expts{i}.Stimvals.or;
    
end
expt_bar_ori(expt_bar_ori == -45) = 135;

load ./all_un_bar_pos
n_bar_pos = size(all_un_bar_pos,1);


%% DETERMINE SET OF EXPERIMENTS TO USE
flen = 12;
bar_oris = [0];
un_bar_pos = all_un_bar_pos(:,1);

fprintf('Analyzing %d ori expts\n',bar_oris);

%     cur_expt_set = find(is_bar_expt == 1 & expt_image_back == 1 & expt_bar_ori == bar_oris(EE));

% %expts with X deg bars and gray back (sim sacs)
% cur_expt_set = find(is_bar_expt==1 & expt_image_back == 0 & expt_bar_ori == bar_oris);

%expts with X deg bars and any back (including sim sacs)
cur_expt_set = find(is_bar_expt==1 & expt_bar_ori == bar_oris);

cur_expt_set(cur_expt_set<=8) = []; %this expt had continuous bar luminance
cur_expt_set(cur_expt_set == 13) = []; %problem with LFP data?
% cur_expt_set(cur_expt_set >= 46 & cur_expt_set <= 51) = []; %problem with image background
cur_expt_set(ismember(cur_expt_set,[46 48 49 51])) = []; %problem with image background
cur_expt_set(cur_expt_set > 60) = []; %no rect

%% COMPUTE TRIAL DATA
all_stim_times = [];
all_phase = [];
all_Op = [];
all_bar_mat = [];

all_t_axis = [];
all_rel_stimes = [];
all_rel_etimes = [];
all_binned_spks = [];
all_used_inds = [];
all_exptvec = [];
all_trialvec = [];
all_trial_start_times = [];
all_trial_end_times = [];
for ee = 1:length(cur_expt_set);
    fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
    cur_expt = cur_expt_set(ee);
    fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
    load(fname);
    
    n_trials = length(Expts{cur_expt}.Trials);
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    
    all_trial_start_times = [all_trial_start_times; trial_start_times(:)];
    all_trial_end_times = [all_trial_end_times; trial_end_times(:)];
    for tt = 1:n_trials
        cur_stim_times = Expts{cur_expt}.Trials(tt).Start/1e4;
        cur_Op = Expts{cur_expt}.Trials(tt).Op;
        cur_phase = Expts{cur_expt}.Trials(tt).ph;
        
        cur_t_edges = [(Expts{cur_expt}.Trials(tt).Start(1)/1e4):1/Fsd:(Expts{cur_expt}.Trials(tt).End(end)/1e4)];
        cur_t_cents = cur_t_edges(1:end-1) + 2/Fsd;
        cur_binned_spks = nan(length(cur_t_cents),length(use_units));
        for cc = 1:length(use_units)
            cur_hist = histc(Clusters{use_units(cc)}.times,cur_t_edges);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
        end
        
        cur_used_inds = ones(length(cur_t_cents),1);
        cur_used_inds(1:beg_buffer) = 0;
        cur_used_inds(end-end_buffer+1:end) = 0;
        
        all_stim_times = [all_stim_times; cur_stim_times];
        all_t_axis = [all_t_axis; cur_t_cents(:)];
        all_rel_stimes = [all_rel_stimes; cur_t_cents(:)- trial_start_times(tt)];
        all_rel_etimes = [all_rel_etimes; trial_end_times(tt) - cur_t_cents(:)];
        all_Op = [all_Op; cur_Op];
        all_phase = [all_phase; cur_phase'];
        all_used_inds = [all_used_inds; cur_used_inds];
        %         all_bar_mat = [all_bar_mat; bar_Xmat];
        all_binned_spks = [all_binned_spks; cur_binned_spks];
        all_exptvec = [all_exptvec; ones(length(cur_t_cents),1)*ee];
        all_trialvec = [all_trialvec; ones(length(cur_t_cents),1)*tt];
    end
    
end

cur_bar_mat = zeros(length(all_stim_times),n_bar_pos);
for bb = 1:n_bar_pos
    cur_set = find(all_Op==un_bar_pos(bb));
    pset = all_phase(cur_set) == 0;
    nset = all_phase(cur_set) == pi;
    cur_bar_mat(cur_set(pset),bb) = 1;
    cur_bar_mat(cur_set(nset),bb) = -1;
    %                 cur_bar_mat(cur_set,bb) = 1;
end

all_bar_mat = makeStimRows(cur_bar_mat,flen);

%%
all_interp_phasegrams = [];
all_interp_ampgrams = [];
all_interp_Vmat = [];
for ee = 1:length(cur_expt_set);
    
    %%
    fprintf('Analyzing block %d of %d\n',ee,length(cur_expt_set));
    
    Vmat = [];
    phasegrams = [];
    ampgrams = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d of %d\n',ll,length(use_lfps));
        filename = sprintf('Expt%d_p%d_despkLFP.mat',cur_expt_set(ee),use_lfps(ll));
        load(filename);
        
        V_down = filtfilt(filt_b,filt_a,despiked_lfp);
        V_down = downsample(V_down,dsf);
        Vmat = [Vmat V_down];
        
        temp = cwt(V_down,scales,'cmor1-1');
        phasegrams = cat(3,phasegrams,angle(temp)');
        ampgrams = cat(3,ampgrams,abs(temp)');
    end
    lfp_t_axis = downsample(t_axd,dsf);
    
    cur_expt_inds = find(all_exptvec == ee);
    
    unwr_phasegram = unwrap(phasegrams);
    interp_phasegrams = interp1(lfp_t_axis,unwr_phasegram,all_t_axis(cur_expt_inds));
    interp_phasegrams = mod(interp_phasegrams+pi,2*pi)-pi;
    interp_ampgrams = interp1(lfp_t_axis,ampgrams,all_t_axis(cur_expt_inds));
    
    interp_Vmat = interp1(lfp_t_axis,Vmat,all_t_axis(cur_expt_inds));
    
    all_interp_Vmat = cat(1,all_interp_Vmat,interp_Vmat);
    all_interp_phasegrams = cat(1,all_interp_phasegrams,interp_phasegrams);
    all_interp_ampgrams = cat(1,all_interp_ampgrams,interp_ampgrams);
end

%%
all_interp_ampgrams = bsxfun(@rdivide,all_interp_ampgrams,std(all_interp_ampgrams));

%%
trial_set = unique(all_trialvec);
n_trials = length(trial_set);

xv_frac = 0.25;
n_xv_trials = round(n_trials*xv_frac);
xv_set = randperm(n_trials);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(all_trialvec,xv_set));
tr_inds = find(~ismember(all_trialvec,xv_set))';
tr_inds(all_used_inds(tr_inds) == 0) = [];
xv_inds(all_used_inds(xv_inds) == 0) = [];
%%
NT = length(all_t_axis);
new_phase_set = [reshape(all_interp_ampgrams,NT,length(wfreqs)*length(use_lfps)).*reshape(cos(all_interp_phasegrams),NT,length(wfreqs)*length(use_lfps)) ...
    reshape(all_interp_ampgrams,NT,length(wfreqs)*length(use_lfps)).*reshape(sin(all_interp_phasegrams),NT,length(wfreqs)*length(use_lfps))];
% new_phase_set = [reshape(cos(all_interp_phasegrams),NT,length(wfreqs)*length(use_lfps)) ...
%     reshape(sin(all_interp_phasegrams),NT,length(wfreqs)*length(use_lfps))];

% phase_elec_set = ones(length(wfreqs),1)*(use_lfps);
% phase_elec_set = [phase_elec_set(:); phase_elec_set(:)];

[phase_freq_set,phase_elec_set] = meshgrid(wfreqs,use_lfps);
phase_elec_set = phase_elec_set';
phase_freq_set = phase_freq_set';
phase_elec_set = [phase_elec_set(:); phase_elec_set(:)];
phase_freq_set = [phase_freq_set(:); phase_freq_set(:)];

%%
NL_type = 0; %exp

reg_params.dl2_freq = 300; %20000
reg_params.dl2_ch =  00; %20000
reg_params.dl2_time = 0; %500
reg_params.dl_freq = 100; %20000
reg_params.dl_ch =  0; %20000
reg_params.dl_time = 0; %500
reg_params.d2_phase = 100;
reg_params.d2_time = 0;

reg_params2.dl2_freq = 300; %20000
reg_params2.dl2_ch =  0; %20000
reg_params2.dl_freq = 100; %20000
reg_params2.dl_ch =  0; %20000
reg_params2.dl2_time = 0; %500
reg_params2.d2_phase = 20000;
reg_params2.d2_time = 0;
reg_params2.dl_time = 0;

silent = 1;
%%
for cc = 1:length(use_lfps)
%%
fprintf('Cell %d\n',cc);
    
    Robs = all_binned_spks(tr_inds,use_lfps(cc));
    tr_spkbns = convert_to_spikebins(Robs);
    
    use_elecs = use_lfps;
    use_set = find(ismember(phase_elec_set,use_elecs));
    cur_Xmat = [new_phase_set(tr_inds,use_set)];
    stim_params = [length(wfreqs),length(use_elecs)];
    klen = size(cur_Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_po(cc)] = fit_GLM_sinphase_model(cur_Xmat, Robs, K0,silent, stim_params, reg_params2,0,NL_type);
    po_sinphase_cfilt(cc,:) = fitp_po(cc).k(1:length(use_elecs)*length(wfreqs));
    po_sinphase_sfilt(cc,:) = fitp_po(cc).k((length(use_elecs)*length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
    po_ampkern(cc,:) = sqrt(po_sinphase_cfilt(cc,:).^2 + po_sinphase_sfilt(cc,:).^2);
    
    use_elecs2 = use_lfps(cc);
    use_set2 = find(ismember(phase_elec_set,use_elecs2));
    cur_Xmat = [new_phase_set(tr_inds,use_set2)];
    stim_params = [length(wfreqs),length(use_elecs2)];
    klen = size(cur_Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_spo(cc)] = fit_GLM_sinphase_model(cur_Xmat, Robs, K0,silent, stim_params, reg_params,0,NL_type);
    spo_sinphase_cfilt(cc,:) = fitp_spo(cc).k(1:length(use_elecs2)*length(wfreqs));
    spo_sinphase_sfilt(cc,:) = fitp_spo(cc).k((length(use_elecs2)*length(wfreqs)+1):length(use_elecs2)*length(wfreqs)*2);
    spo_ampkern(cc,:) = sqrt(spo_sinphase_cfilt(cc,:).^2 + spo_sinphase_sfilt(cc,:).^2);
    
    Robs_xv = all_binned_spks(xv_inds,use_lfps(cc));
    
    cur_Xmat = [new_phase_set(xv_inds,use_set)];
    xv_po_pred_rate = cur_Xmat*fitp_po(cc).k(1:end-1) + fitp_po(cc).k(end);
    if NL_type == 0
        xv_po_pred_rate = log(1+exp(xv_po_pred_rate));
    else
        xv_po_pred_rate = exp(xv_po_pred_rate);
    end
    xvLL_po(cc) = -sum(Robs_xv.*log(xv_po_pred_rate)-xv_po_pred_rate)/sum(Robs_xv);

    cur_Xmat = [new_phase_set(xv_inds,use_set2)];
    xv_spo_pred_rate = cur_Xmat*fitp_spo(cc).k(1:end-1) + fitp_spo(cc).k(end);
    if NL_type == 0
        xv_spo_pred_rate = log(1+exp(xv_spo_pred_rate));
    else
        xv_spo_pred_rate = exp(xv_spo_pred_rate);
    end
    xvLL_spo(cc) = -sum(Robs_xv.*log(xv_spo_pred_rate)-xv_spo_pred_rate)/sum(Robs_xv);

    avg_rate = mean(Robs);
    null_prate = avg_rate*ones(size(Robs_xv));
    xvLL_null(cc) = -sum(Robs_xv.*log(null_prate)-null_prate)/sum(Robs_xv);
    
    xv_imp(cc) = (xvLL_null(cc) - xvLL_po(cc))/log(2);
    xv_imps(cc) = (xvLL_null(cc) - xvLL_spo(cc))/log(2);
    fprintf('xv imp: %.4f  xv imps: %.4f\n',xv_imp(cc),xv_imps(cc));
    
end
phasekern = 180/pi*(atan2(po_sinphase_cfilt,po_sinphase_sfilt)+pi);

%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)
load ./jbeG081.em.mat
all_eye_vals = [];
all_eye_speed = [];
all_eye_ts = [];
for ee = 1:length(cur_expt_set);
    cur_set = find(all_exptvec==ee);
    [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data_new(Expt,all_t_axis(cur_set([1 end])));

    eye_dt = median(diff(eye_ts_interp));
    eye_fs = 1/eye_dt;
    lEyeXY = eye_vals_interp(:,1:2);
    rEyeXY = eye_vals_interp(:,3:4);
    clear sm_avg_eyepos eye_vel
    sm_avg_eyepos(:,1) = smooth(lEyeXY(:,1),3);
    sm_avg_eyepos(:,2) = smooth(lEyeXY(:,2),3);
    eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/eye_dt;
    eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/eye_dt;

    eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);

    all_eye_vals = [all_eye_vals; lEyeXY rEyeXY];
    all_eye_speed = [all_eye_speed; eye_speed];
    all_eye_ts = [all_eye_ts; eye_ts_interp'];
end
back_pts = 1 + find(diff(all_eye_ts) <= 0);
double_samples = [];
for i = 1:length(back_pts)
    next_forward = find(all_eye_ts > all_eye_ts(back_pts(i)-1),1,'first');
    double_samples = [double_samples back_pts(i):next_forward];
end
all_eye_ts(double_samples) = [];
all_eye_speed(double_samples) = [];
all_eye_vals(double_samples,:) = [];
eye_fs = 1/median(diff(all_eye_ts));

interp_eye_vals = interp1(all_eye_ts,all_eye_vals,all_t_axis);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_t_axis);

orth_eye_pos = interp_eye_vals(:,2);
par_eye_pos = interp_eye_vals(:,1);
bad_pts = find(abs(orth_eye_pos) > 1); %cheap way of detecting blinks because the right eye signal was unreliable even for that at times

all_eye_intrial = zeros(size(all_eye_ts));
for i = 1:length(all_trial_start_times)
    curset = find(all_eye_ts >= all_trial_start_times(i) & all_eye_ts <= all_trial_end_times(i));
    all_eye_intrial(curset) = 1;
end

%compute times saccades
sac_thresh = 10;
min_isi = 0.05;
orig_saccade_inds = 1 + find(all_eye_speed(1:end-1) < sac_thresh & all_eye_speed(2:end) > sac_thresh);
orig_saccade_inds = unique(orig_saccade_inds);
isis = [Inf; diff(orig_saccade_inds)]/eye_fs;
double_sacs = find(isis < min_isi);
orig_saccade_inds(double_sacs) = [];
orig_saccade_inds(all_eye_intrial(orig_saccade_inds) == 0) = [];
interp_saccade_inds = round(interp1(all_t_axis,1:length(all_t_axis),all_eye_ts(orig_saccade_inds)));

peri_thresh = 3;
sac_starts = nan(size(interp_saccade_inds));
sac_stops = nan(size(interp_saccade_inds));
for i = 1:length(interp_saccade_inds)
    cur_up = find(interp_eye_speed(1:interp_saccade_inds(i)) < peri_thresh,1,'last');
    cur_down = find(interp_eye_speed(interp_saccade_inds(i):end) < peri_thresh,1,'first');
    if ~isempty(cur_up) & ~isempty(cur_down)
        sac_starts(i) = cur_up;
        sac_stops(i) = interp_saccade_inds(i) + cur_down;
    end
end
bad = find(isnan(sac_starts));
interp_saccade_inds(bad) = []; sac_starts(bad) = []; sac_stops(bad) = [];

sac_starts = sac_starts - round(Fsd*0.01);
sac_stops = sac_stops + round(Fsd*0.01);
pre_pos = par_eye_pos(sac_starts);
post_pos = par_eye_pos(sac_stops);

sac_dist = abs(post_pos-pre_pos);
if bar_oris == 0
    outgoing_sacs = interp_saccade_inds(pre_pos > -0.9 & post_pos < -0.9);
    returning_sacs = interp_saccade_inds(pre_pos < -1 & post_pos > -1);
elseif bar_oris == 135
    outgoing_sacs = interp_saccade_inds(pre_pos > -1 & post_pos < -1);
    returning_sacs = interp_saccade_inds(pre_pos < -1.1 & post_pos > -1.2);
elseif bar_oris == 45
    outgoing_sacs = interp_saccade_inds(pre_pos > -1 & post_pos < -1);
    returning_sacs = interp_saccade_inds(pre_pos < -1.1 & post_pos > -1);
elseif bar_oris == 90
    outgoing_sacs = interp_saccade_inds(pre_pos > -1 & post_pos < -1);
    returning_sacs = interp_saccade_inds(pre_pos < -1 & post_pos > -1);
end
msacs = interp_saccade_inds(sac_dist' < 1); msacs(ismember(msacs,outgoing_sacs)) = []; msacs(ismember(msacs,returning_sacs)) = [];

use_out = find(all_rel_stimes(outgoing_sacs) > 0.8 & all_rel_stimes(outgoing_sacs) < 1.3);
use_ret = find(all_rel_stimes(returning_sacs) > 1.4 & all_rel_stimes(returning_sacs) < 1.9);

outgoing_sacs = outgoing_sacs(use_out);
returning_sacs = returning_sacs(use_ret);


trial_sac_inds = zeros(size(all_t_axis));
trial_sac_inds(interp_saccade_inds) = 1;
tr_sac_inds = find(trial_sac_inds(tr_inds)==1);
xv_sac_inds = find(trial_sac_inds(xv_inds)==1);

trial_msac_inds = zeros(size(all_t_axis));
trial_msac_inds(msacs) = 1;
tr_msac_inds = find(trial_msac_inds(tr_inds)==1);
xv_msac_inds = find(trial_msac_inds(xv_inds)==1);

trial_fsac_inds = zeros(size(all_t_axis));
trial_fsac_inds(outgoing_sacs) = 1;
tr_fsac_inds = find(trial_fsac_inds(tr_inds)==1);
xv_fsac_inds = find(trial_fsac_inds(xv_inds)==1);

trial_ssac_inds = zeros(size(all_t_axis));
trial_ssac_inds(returning_sacs) = 1;
tr_ssac_inds = find(trial_ssac_inds(tr_inds)==1);
xv_ssac_inds = find(trial_ssac_inds(xv_inds)==1);

%%
clear *sac_trig_avg*
all_binned_spks_norm = bsxfun(@rdivide,all_binned_spks,mean(all_binned_spks(tr_inds,:)));
sm_sig = round(Fsd*0.005);
for cc = 1:96
    all_binned_spks_norm(:,cc) = jmm_smooth_1d_cor(all_binned_spks_norm(:,cc),sm_sig);
end
for cc = 1:length(use_lfps)
    [sac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(tr_inds,use_lfps(cc)),tr_sac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [sac_trig_avg_V(cc,:),cur_lags] = get_event_trig_avg(all_interp_Vmat(tr_inds,cc),tr_sac_inds,round(0.4*Fsd),round(0.5*Fsd));

    [msac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(tr_inds,use_lfps(cc)),tr_msac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [msac_trig_avg_V(cc,:),cur_lags] = get_event_trig_avg(all_interp_Vmat(tr_inds,cc),tr_msac_inds,round(0.4*Fsd),round(0.5*Fsd));

    [fsac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(tr_inds,use_lfps(cc)),tr_fsac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [fsac_trig_avg_V(cc,:),cur_lags] = get_event_trig_avg(all_interp_Vmat(tr_inds,cc),tr_fsac_inds,round(0.4*Fsd),round(0.5*Fsd));

    [ssac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(tr_inds,use_lfps(cc)),tr_ssac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [ssac_trig_avg_V(cc,:),cur_lags] = get_event_trig_avg(all_interp_Vmat(tr_inds,cc),tr_ssac_inds,round(0.4*Fsd),round(0.5*Fsd));


    [xvsac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(xv_inds,use_lfps(cc)),xv_sac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [xvsac_trig_avg_V(cc,:),cur_lags] = get_event_trig_avg(all_interp_Vmat(xv_inds,cc),xv_sac_inds,round(0.4*Fsd),round(0.5*Fsd));

    [xvmsac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(xv_inds,use_lfps(cc)),xv_msac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [xvmsac_trig_avg_V(cc,:),cur_lags] = get_event_trig_avg(all_interp_Vmat(xv_inds,cc),xv_msac_inds,round(0.4*Fsd),round(0.5*Fsd));

    [xvfsac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(xv_inds,use_lfps(cc)),xv_fsac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [xvfsac_trig_avg_V(cc,:),cur_lags] = get_event_trig_avg(all_interp_Vmat(xv_inds,cc),tr_fsac_inds,round(0.4*Fsd),round(0.5*Fsd));

    [xvssac_trig_avg_rate(cc,:),cur_lags] = get_event_trig_avg(all_binned_spks_norm(xv_inds,use_lfps(cc)),xv_ssac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [xvssac_trig_avg_V(cc,:),cur_lags] = get_event_trig_avg(all_interp_Vmat(xv_inds,cc),xv_ssac_inds,round(0.4*Fsd),round(0.5*Fsd));
end

%%
% ufreq_set = find(wfreqs <= 60);
ufreq_set = find(wfreqs <= 60 & wfreqs > 4);
% ufreq_set2 = find(wfreqs <= 10);
% ufreq_set3 = find(wfreqs <= 80);
clear st_*
for cc = 1:length(use_lfps)
    cc
%     use_elecs = use_lfps(cc);
        use_elecs = use_lfps;
    use_set = find(ismember(phase_elec_set,use_elecs));
%     cur_Xmat = [new_phase_set(tr_inds,use_set)];
    cur_kern = fitp_po(cc).k(1:2*length(use_elecs)*length(wfreqs));
    
    use_freqs = find(ismember(phase_freq_set(use_set),wfreqs(ufreq_set)));
     cur_Xmat = [new_phase_set(tr_inds,use_set)];
   phasemod_out = exp(cur_Xmat(:,use_freqs)*cur_kern(use_freqs));
    [st_avg_phasemod_po(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,tr_sac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [st_avg_mphasemod_po(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,tr_msac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [st_avg_fphasemod_po(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,tr_fsac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [st_avg_sphasemod_po(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,tr_ssac_inds,round(0.4*Fsd),round(0.5*Fsd));
      cur_Xmat = [new_phase_set(xv_inds,use_set)];
   phasemod_out = exp(cur_Xmat(:,use_freqs)*cur_kern(use_freqs));
   [xvst_avg_phasemod_po(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,xv_sac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [xvst_avg_mphasemod_po(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,xv_msac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [xvst_avg_fphasemod_po(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,xv_fsac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [xvst_avg_sphasemod_po(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,xv_ssac_inds,round(0.4*Fsd),round(0.5*Fsd));

    use_elecs = use_lfps(cc);
    use_set = find(ismember(phase_elec_set,use_elecs));
    cur_kern = fitp_spo(cc).k(1:2*length(use_elecs)*length(wfreqs));
    
    use_freqs = find(ismember(phase_freq_set(use_set),wfreqs(ufreq_set)));
     cur_Xmat = [new_phase_set(tr_inds,use_set)];
   phasemod_out = exp(cur_Xmat(:,use_freqs)*cur_kern(use_freqs));
    [st_avg_phasemod_pos(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,tr_sac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [st_avg_mphasemod_pos(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,tr_msac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [st_avg_fphasemod_pos(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,tr_fsac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [st_avg_sphasemod_pos(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,tr_ssac_inds,round(0.4*Fsd),round(0.5*Fsd));
      cur_Xmat = [new_phase_set(xv_inds,use_set)];
   phasemod_out = exp(cur_Xmat(:,use_freqs)*cur_kern(use_freqs));
   [xvst_avg_phasemod_pos(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,xv_sac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [xvst_avg_mphasemod_pos(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,xv_msac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [xvst_avg_fphasemod_pos(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,xv_fsac_inds,round(0.4*Fsd),round(0.5*Fsd));
    [xvst_avg_sphasemod_pos(cc,:),cur_lags] = get_event_trig_avg(phasemod_out,xv_ssac_inds,round(0.4*Fsd),round(0.5*Fsd));
end

%%
% close all
% sm_sig = round(0.0025*Fsd);
% sm_sig2 = round(0.01*Fsd);
% sm_sig3 = round(0.02*Fsd);
% for cc = 1:length(use_lfps)
%     cc
%     subplot(3,1,1)
%     plot(cur_lags/Fsd,jmm_smooth_1d_cor(sac_trig_avg_rate(cc,:),sm_sig))
% hold on
%     plot(cur_lags/Fsd,jmm_smooth_1d_cor(sac_trig_avg_rate(cc,:),sm_sig2),'k','linewidth',2)
%     plot(cur_lags/Fsd,jmm_smooth_1d_cor(sac_trig_avg_rate(cc,:),sm_sig3),'r')
% xlim([-0.2 0.5]); grid on
%    subplot(3,1,2)
%     plot(cur_lags/Fsd,st_avg_phasemod_po(cc,:))
%     hold on
%     plot(cur_lags/Fsd,st_avg_phasemod_po2(cc,:),'r')
%     plot(cur_lags/Fsd,st_avg_phasemod_po3(cc,:),'k','linewidth',2)
%     xlim([-0.2 0.5]); grid on
%    subplot(3,1,3)
%     plot(cur_lags/Fsd,sac_trig_avg_V(cc,:))
% xlim([-0.2 0.5]); grid on
%     pause
%     clf
% end
% 
cd ~/Data/bruce/7_15_12/G034/
load ./gabor_tracking_varmeans
gabor_params_used = gabor_params_f{2};
load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;

cd ~/Data/bruce/G081

%fit smoothed retinotopic surface
interp_x = nan(96,1);
interp_y = nan(96,1);
tempinds = zeros(10,10);
for i = 1:96
    tempx(Y_pos(i),X_pos(i)) = gabor_params_used(i,1);
    tempy(Y_pos(i),X_pos(i)) = gabor_params_used(i,2);
    tempinds(Y_pos(i),X_pos(i)) = i;
end
weights = ones(10,10);
weights(tempx==0) = 0;
xpos_interp = smoothn(tempx,weights,'robust');
ypos_interp = smoothn(tempy,weights,'robust');
used_inds = find(weights == 1);
tempinds = tempinds(used_inds);
interp_x(tempinds) = xpos_interp(used_inds);
interp_y(tempinds) = ypos_interp(used_inds);
xi = linspace(min(interp_x),max(interp_x),50);
yi = linspace(min(interp_y),max(interp_y),50);
[Xi,Yi] = meshgrid(xi,yi);

id_mat = nan(10,10);
for i = 1:10
    for j = 1:10
        cur = find(X_pos==j&Y_pos==i);
        if ~isempty(cur)
            id_mat(i,j) = cur;
        end
    end
end
use_ids = find(~isnan(id_mat));

% [~,ord] = sort(Y_pos);
% [~,ord] = sort(X_pos);
% close all
% 
el_fdist = sqrt(interp_x.^2 + interp_y.^2);
el_fdist = el_fdist(use_lfps);
[temp,ord] = sort(el_fdist);
%%
 xv_imp = (xvLL_null - xvLL_po)/log(2);
 avg_rates = mean(all_binned_spks(tr_inds,use_lfps));
close all
sm_sig = ceil(0.0025*Fsd);
xl = [-0.15 0.4];
    f1 = figure;
for cc = 1:length(use_lfps)
    use_lfps(ord(cc))
%     f1 = figure;
    subplot(3,4,1);hold on
    plot(cur_lags/Fsd,sac_trig_avg_rate(ord(cc),:),'r')
    plot(cur_lags/Fsd,xvsac_trig_avg_rate(ord(cc),:),'b')
 xlim(xl); grid on
    subplot(3,4,2);hold on
    plot(cur_lags/Fsd,msac_trig_avg_rate(ord(cc),:),'r')
    plot(cur_lags/Fsd,xvmsac_trig_avg_rate(ord(cc),:),'b')
 xlim(xl); grid on
   subplot(3,4,3);hold on
    plot(cur_lags/Fsd,fsac_trig_avg_rate(ord(cc),:),'r')
     plot(cur_lags/Fsd,xvfsac_trig_avg_rate(ord(cc),:),'b')
xlim(xl); grid on
    subplot(3,4,4);hold on
    plot(cur_lags/Fsd,ssac_trig_avg_rate(ord(cc),:),'r')
     plot(cur_lags/Fsd,xvssac_trig_avg_rate(ord(cc),:),'b')
xlim(xl); grid on
    subplot(3,4,5); hold on
%     plot(cur_lags/Fsd,st_avg_sphasemod_po(cc,:),'k')
    plot(cur_lags/Fsd,(st_avg_phasemod_po(ord(cc),:)),'r')
    plot(cur_lags/Fsd,(xvst_avg_phasemod_po(ord(cc),:)),'b','linewidth',2)
%     plot(cur_lags/Fsd,exp(xvst_avg_phasemod_pos(ord(cc),:)),'k')
 xlim(xl); grid on
   subplot(3,4,6); hold on
%     plot(cur_lags/Fsd,st_avg_mphasemod_po(cc,:))
    plot(cur_lags/Fsd,(st_avg_mphasemod_po(ord(cc),:)),'r')
    plot(cur_lags/Fsd,(xvst_avg_mphasemod_po(ord(cc),:)),'b','linewidth',2)
%     plot(cur_lags/Fsd,exp(xvst_avg_mphasemod_pos(ord(cc),:)),'k')
    hold on
 xlim(xl); grid on
   subplot(3,4,7);hold on
%     plot(cur_lags/Fsd,st_avg_fphasemod_po(cc,:),'r')
    plot(cur_lags/Fsd,(st_avg_fphasemod_po(ord(cc),:)),'r')
     plot(cur_lags/Fsd,(xvst_avg_fphasemod_po(ord(cc),:)),'b','linewidth',2)
%      plot(cur_lags/Fsd,exp(xvst_avg_fphasemod_pos(ord(cc),:)),'k')
 xlim(xl); grid on
   subplot(3,4,8); hold on
%     plot(cur_lags/Fsd,st_avg_sphasemod_po(cc,:),'k')
    plot(cur_lags/Fsd,(st_avg_sphasemod_po(ord(cc),:)),'r')
    plot(cur_lags/Fsd,(xvst_avg_sphasemod_po(ord(cc),:)),'b','linewidth',2)
%     plot(cur_lags/Fsd,exp(xvst_avg_sphasemod_pos(ord(cc),:)),'k')
 xlim(xl); grid on
    subplot(3,4,9); hold on
    plot(cur_lags/Fsd,sac_trig_avg_V(ord(cc),:),'r')
    plot(cur_lags/Fsd,xvsac_trig_avg_V(ord(cc),:),'b')
 xlim(xl); grid on
  subplot(3,4,10); hold on
    plot(cur_lags/Fsd,msac_trig_avg_V(ord(cc),:),'r')
    plot(cur_lags/Fsd,xvmsac_trig_avg_V(ord(cc),:))
 xlim(xl); grid on
   subplot(3,4,11); hold on
    plot(cur_lags/Fsd,fsac_trig_avg_V(ord(cc),:),'r')
    plot(cur_lags/Fsd,xvfsac_trig_avg_V(ord(cc),:),'b')
 xlim(xl); grid on
   subplot(3,4,12); hold on
    plot(cur_lags/Fsd,ssac_trig_avg_V(ord(cc),:),'r')
    plot(cur_lags/Fsd,xvssac_trig_avg_V(ord(cc),:),'b')
 xlim(xl); grid on
    
% f2 = figure;
%         subplot(2,1,1)
%     plot(wfreqs,po_ampkern(ord(cc),:))
%     hold on
%     xlim(wfreqs([end 1]))
%     subplot(2,1,2)
%     plot(wfreqs,phasekern(ord(cc),:))
%     hold on
%     ylim([0 360])
%     xlim(wfreqs([end 1]))

    fprintf('phase imp %.4f\n',xv_imp(ord(cc)));
    
    set(f1,'Position',[1 200 2000 1000]);shg
%     set(f2,'Position',[1100 400 500 500]);shg
pause;clf
%     close all
end

%%
close all
for cc = 1:length(use_lfps)
    subplot(2,1,1)
    plot(wfreqs,po_ampkern(ord(cc),:))
    hold on
    xlim(wfreqs([end 1]))
    subplot(2,1,2)
    plot(wfreqs,phasekern(ord(cc),:))
    hold on
    ylim([0 360])
    xlim(wfreqs([end 1]))
    pause
    clf
end