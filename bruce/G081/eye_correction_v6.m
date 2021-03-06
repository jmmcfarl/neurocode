clear all
% close all
cd ~/Data/bruce/G081/
load jbeG081Expts.mat
addpath('~/James_scripts/bruce/G081/');
%%
stim_fs = 100; %in Hz
use_sus = 1:96;

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
beg_buffer = round(stim_fs*0.15); %don't use the first X data after start of a trial.
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
all_rel_stimes = [];
all_rel_etimes = [];
all_phase = [];
all_Op = [];
all_bar_mat = [];
all_used_inds = [];
all_binned_spks = [];
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
        
        cur_t_edges = [cur_stim_times; Expts{cur_expt}.Trials(tt).End(end)/1e4];
        cur_binned_spks = nan(length(cur_stim_times),length(use_sus));
        for cc = 1:length(use_sus)
            cur_hist = histc(Clusters{use_sus(cc)}.times,cur_t_edges);
            cur_binned_spks(:,cc) = cur_hist(1:end-1);
        end
        
        cur_used_inds = ones(length(cur_stim_times),1);
        cur_used_inds(1:flen) = 0;
        cur_used_inds(1:beg_buffer) = 0;
        
        all_stim_times = [all_stim_times; cur_stim_times];
        all_rel_stimes = [all_rel_stimes; cur_stim_times- trial_start_times(tt)];
        all_rel_etimes = [all_rel_etimes; trial_end_times(tt) - cur_stim_times];
        all_Op = [all_Op; cur_Op];
        all_phase = [all_phase; cur_phase'];
        all_used_inds = [all_used_inds; cur_used_inds];
        %         all_bar_mat = [all_bar_mat; bar_Xmat];
        all_binned_spks = [all_binned_spks; cur_binned_spks];
        all_exptvec = [all_exptvec; ones(size(cur_stim_times))*ee];
        all_trialvec = [all_trialvec; ones(size(cur_stim_times))*tt];
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

%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)
load ./jbeG081.em.mat
all_eye_vals = [];
all_eye_speed = [];
all_eye_ts = [];
for ee = 1:length(cur_expt_set);
    cur_set = find(all_exptvec==ee);
    [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data_new(Expt,all_stim_times(cur_set([1 end])));
    
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

interp_eye_vals = interp1(all_eye_ts,all_eye_vals,all_stim_times);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_stim_times);

orth_eye_pos = interp_eye_vals(:,2);
par_eye_pos = interp_eye_vals(:,1);
% orth_eye_pos = interp_eye_vals(:,1);
% orth_eye_pos = interp_eye_vals(:,1)*cos(pi/4) + interp_eye_vals(:,2)*sin(pi/4);
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
interp_saccade_inds = round(interp1(all_stim_times,1:length(all_stim_times),all_eye_ts(orig_saccade_inds)));

% peri_thresh = 3;
% sac_starts = nan(size(interp_saccade_inds));
% sac_stops = nan(size(interp_saccade_inds));
% for i = 1:length(interp_saccade_inds)
%     cur_up = find(interp_eye_speed(1:interp_saccade_inds(i)) < peri_thresh,1,'last');
%     cur_down = find(interp_eye_speed(interp_saccade_inds(i):end) < peri_thresh,1,'first');
%     if ~isempty(cur_up) & ~isempty(cur_down)
%         sac_starts(i) = cur_up;
%         sac_stops(i) = interp_saccade_inds(i) + cur_down;
%     end
% end
% bad = find(isnan(sac_starts));
% interp_saccade_inds(bad) = []; interp_saccade_inds(bad) = []; interp_saccade_inds(bad) = [];
% 
% sac_starts = sac_starts - 1;
% sac_stops = sac_stops + 1;
% bad = find(sac_starts <= 0 | sac_starts > length(orth_eye_pos));
% interp_saccade_inds(bad) = []; sac_starts(bad) = []; sac_stops(bad) = [];
% pre_pos = orth_eye_pos(sac_starts);
% post_pos = orth_eye_pos(sac_stops);
% sac_amp = post_pos-pre_pos;

%% SET UP USEABLE TRIALS
[c,ia,ic] = unique([all_exptvec all_trialvec],'rows');
n_trials = length(ia);
bad_trials = unique(ic(bad_pts)); %trials with putative blinks

tr_set = 1:n_trials;
tr_set = find(ismember(1:n_trials,tr_set));

tr_inds = find(ismember(ic,tr_set));

tr_inds(all_used_inds(tr_inds) == 0) = [];

Xmat = all_bar_mat;
%% double sample Xmat
bar_dx = 0.08;
up_fac = 2;
bar_dxd = bar_dx/up_fac;
SDIM = n_bar_pos;
NT = length(tr_inds);
resh_X = reshape(0.5*Xmat(tr_inds,:)',[flen SDIM NT]);
if up_fac > 1
    fprintf('Up-sampling X-matrix\n');
    dresh_X = [];
    for i = 1:SDIM-1
        dresh_X = cat(2,dresh_X,resh_X(:,i,:));
        %for interpolating bar positions
        dresh_X = cat(2,dresh_X,0.5*resh_X(:,i,:) + 0.5*resh_X(:,i+1,:));
    end
    dresh_X = cat(2,dresh_X,resh_X(:,end,:));
    
    up_SDIM = up_fac*SDIM-1;
else
    up_SDIM = SDIM;
    dresh_X = resh_X;
end
X_z = reshape(dresh_X,flen*up_SDIM,NT)';

%% LOAD INITIAL MODELS
load ./CellList.mat
good_sus = find(all(CellList(:,:,1) > 0));

% load ./eye_track_initial_models_0deg2.mat
load ./eye_track_initial_models_90deg.mat
% load ./eye_track_initial_models_0deg_gray.mat
%% DETERMINE WHICH CELLS HAVE USABLE MODELS
NUNITS = 96;
xvLL_imp = (null_xvLL - init_xvLL)/log(2);
xvLL_thresh = 0.001;

useable_cells = find(xvLL_imp >= xvLL_thresh);

%% DEFINE POINTS TO ALLOW RAPID CHANGES (TRIAL STARTS AFTER SACCADES)
trial_start_inds = [1; find(diff(all_trialvec(tr_inds)) ~= 0) + 1];

sac_thresh = 10;
temp = zeros(length(all_stim_times),1);
temp(interp_saccade_inds) = 1;
temp = temp(tr_inds);
saccade_inds = find(temp==1);

%push the effects of saccades forward in time
sac_shift = round(0.05*stim_fs);

use_prior = logical(zeros(size(tr_inds)));
use_prior(trial_start_inds) = 1;
for i = 1:length(saccade_inds)
    next_trial = trial_start_inds(find(trial_start_inds > saccade_inds(i),1,'first'));
    %mark a point (forward in time) as a saccade if it occurs within the
    %same trial as the saccade
    if next_trial > saccade_inds(i) + sac_shift
        use_prior(saccade_inds(i) + sac_shift) = 1;
    end
end
chunk_boundaries = find(diff(use_prior) > 0);
chunk_start_inds = [1; chunk_boundaries+1];
chunk_stop_inds = [chunk_boundaries; NT];
chunk_durs = (chunk_stop_inds - chunk_start_inds)/stim_fs;
chunk_ids = nan(NT,1);
for i = 1:length(chunk_durs)
    chunk_ids(chunk_start_inds(i):chunk_stop_inds(i)) = i;
end

rsac_inds = zeros(size(tr_inds));
rsac_inds(saccade_inds) = 1;

sim_sac_inds1 = find(all_rel_stimes(tr_inds(1:end-1)) < 0.7 & all_rel_stimes(tr_inds(2:end)) >= 0.7 & expt_sim_sacs(cur_expt_set(all_exptvec(tr_inds(1:end-1))))' > 0);
sim_sac_inds2 = find(all_rel_stimes(tr_inds(1:end-1)) < 1.4 & all_rel_stimes(tr_inds(2:end)) >= 1.4 & expt_sim_sacs(cur_expt_set(all_exptvec(tr_inds(1:end-1))))' > 0);
sim_sac_inds = unique([sim_sac_inds1; sim_sac_inds2]);
ssac_inds = zeros(size(tr_inds));
ssac_inds(sim_sac_inds) = 1;

%% CREATE EVENT PREDICTORS FOR REAL AND SIM SACCADES (POOLING MICRO AND MACRO SACS)
dt = 1/stim_fs;
tent_centers = [0:dt:0.4];
tent_centers = round(tent_centers/dt);
tbmat = construct_tent_bases(tent_centers,1);
[ntents,tblen] = size(tbmat);

shift = round(0.05/dt); %look this far before measured saccade times for modulation.
tbmat = [zeros(ntents,tblen-2*shift-1) tbmat];
tent_centers = tent_centers-shift;

%saccade indicator matrix
rsac_Tmat = zeros(NT,ntents);
for i = 1:ntents
    rsac_Tmat(:,i) = conv(rsac_inds,tbmat(i,:),'same');
end

%indicator variable for experiment number
un_expts = unique(all_exptvec);
expt_Tmat = zeros(NT,length(un_expts));
for i = 1:length(un_expts)
    temp_set = find(all_exptvec(tr_inds)==i);
    expt_Tmat(temp_set,i) = 1;
end

%% FIT EVENT KERNELS (SACCADES AND EXPERIMENT NUMBER)
for cc = 1:96
    fprintf('Fitting event kernels for cell %d of %d\n',cc,96);
    
    Robs = all_binned_spks(tr_inds,cc);
    tr_spkbns = convert_to_spikebins(Robs);
    
    lamrange2 = [1000 1 ntents 0]; %smoothness and saccade kernels
    llist = []; %l1 prior
    %prediction of stimulus model
    [nll, pnll, lpen, prate, g] = getLL_GNM(upquad_fit(cc),X_z,tr_spkbns,'none');
    Xmat = [rsac_Tmat expt_Tmat g];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    K0(end-1) = 1;
    constraints = [klen klen+1]; %don't refit these coefficients
    [fitp_sac_kerns(cc)] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], lamrange2,[], [], llist, constraints, 0);
    init_rsac_kern(cc,:) = fitp_sac_kerns(cc).k(1:ntents);
    init_expt_kern(cc,:) = fitp_sac_kerns(cc).k(ntents+1:end-2);
    pre_offset(cc) = upquad_fit(cc).spk_theta;
end

%% INITIALIZE TRANSITION PRIORS FOR HMM
NT = length(tr_inds);
sp_dx = 0.08/up_fac;
max_shift = 8*up_fac;
dshift = 1;
shifts = -max_shift:dshift:max_shift;
n_shifts = length(shifts);
zero_frame = find(shifts==0);

%overall prior on shifts
eps_prior_sigma = 0.15; %0.2
leps_prior = -(shifts*sp_dx).^2/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
lA_tflip = repmat(leps_prior,n_shifts,1);

cdist = squareform(pdist(shifts'*sp_dx));
deps_sigma = 0.015; %.015
lA = -cdist.^2/(2*deps_sigma^2);
lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize

%% ESTIMATE LL for each shift in each stimulus frame
tr_set = useable_cells;
n_tr_chs = length(tr_set);

Robs = all_binned_spks(tr_inds,tr_set);

up_klen = up_SDIM*flen;

%create filter banks for all units
lfilt_bank = zeros(length(use_sus),up_klen);
qfilt1_bank = zeros(length(use_sus),up_klen);
qfilt2_bank = zeros(length(use_sus),up_klen);
for cc = 1:length(use_sus)
    cur_k = get_k_mat(upquad_fit(cc));
    lfilt_bank(cc,:) = cur_k(:,1);
    qfilt1_bank(cc,:) = cur_k(:,2);
    if use_mod(cc) == 3
        qfilt2_bank(cc,:) = cur_k(:,3);
    end
end
lfilt_bank = reshape(lfilt_bank',[flen,up_SDIM,96]);
qfilt1_bank = reshape(qfilt1_bank',[flen,up_SDIM,96]);
qfilt2_bank = reshape(qfilt2_bank',[flen,up_SDIM,96]);
lfilt_bank = lfilt_bank(:,:,tr_set);
qfilt1_bank = qfilt1_bank(:,:,tr_set);
qfilt2_bank = qfilt2_bank(:,:,tr_set);

%initialize filter shift matrices
shifted_lfilt_bank = nan(klen,n_tr_chs);
shifted_qfilt1_bank = nan(klen,n_tr_chs);
shifted_qfilt2_bank = nan(klen,n_tr_chs);

%indicator predictions
rsac_out = rsac_Tmat*init_rsac_kern(tr_set,:)';
expt_out = expt_Tmat*init_expt_kern(tr_set,:)';

%precompute LL at all shifts for all units
LLs = nan(NT,length(tr_set),n_shifts);
for xx = 1:length(shifts)
    fprintf('Shift %d of %d\n',xx,n_shifts);
    d2 = dist_shift2d(lfilt_bank,shifts(xx),2,0);
    shifted_lfilt_bank = reshape(d2,up_klen,n_tr_chs);
    d2 = dist_shift2d(qfilt1_bank,shifts(xx),2,0);
    shifted_qfilt1_bank = reshape(d2,up_klen,n_tr_chs);
    d2 = dist_shift2d(qfilt2_bank,shifts(xx),2,0);
    shifted_qfilt2_bank = reshape(d2,up_klen,n_tr_chs);
    
    %outputs of stimulus models at current X-matrix shift
    lin_out = X_z*shifted_lfilt_bank;
    quad_out = (X_z*shifted_qfilt1_bank).^2 + (X_z*shifted_qfilt2_bank).^2;
    gfun = lin_out + quad_out;
    gfun = bsxfun(@plus,gfun,pre_offset(tr_set));
        
    %add contributions from event kernels
    gfun = gfun + rsac_out + expt_out;
    
    too_large = gfun > 50;
    pred_rate = log(1+exp(gfun));
    pred_rate(too_large) = gfun(too_large);
    
    pred_rate(pred_rate < 1e-20) = 1e-20;
    
    LLs(:,:,xx) = Robs.*log(pred_rate) - pred_rate;
end

%% DO INITIAL EYE INFERENCE AND MODEL RE-FITTING FOR XV SET
resh_X = reshape(X_z',[flen up_SDIM NT]);
resh_X_sh = zeros(size(resh_X));

poss_xv_set = find(xvLL_imp > 0.005);
for XV = 1:length(poss_xv_set)
    xv_cell = poss_xv_set(XV);
    
    cur_tr_inds = find(tr_set~=xv_cell);
    
    fprintf('Using XV cell %d,  %d of %d\n',xv_cell,XV,length(poss_xv_set));
    frame_LLs = squeeze(sum(LLs(:,cur_tr_inds,:),2));
    
    %%
    lB = frame_LLs;
    lalpha=zeros(NT,n_shifts);
    lbeta = zeros(NT,n_shifts);
    lscale=zeros(NT,1); %initialize rescaling parameters
    %compute rescaled forward messages
    lalpha(1,:) = leps_prior + lB(1,:);
    lscale(1)=logsumexp(lalpha(1,:));
    lalpha(1,:) = lalpha(1,:) - lscale(1);
    for t=2:NT
        %         if mod(t,1000)==0
        %             fprintf('%d of %d\n',t,NT);
        %         end
        if use_prior(t)
            cur_lA = lA_tflip;
        else
            cur_lA = lA;
        end
        lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + lB(t,:);
        lscale(t) = logsumexp(lalpha(t,:));
        lalpha(t,:)= lalpha(t,:) - lscale(t);
    end
    
    %compute rescaled backward messages
    lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
    for t=NT-1:-1:1
        %         if mod(t,1000)==0
        %             fprintf('%d of %d\n',t,NT);
        %         end
        if use_prior(t+1)
            cur_lA = lA_tflip;
        else
            cur_lA = lA;
        end
        lf1 = lbeta(t+1,:) + lB(t+1,:);
        lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
    end
    lgamma= lalpha + lbeta;
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    
    gamma = exp(lgamma);
    post_mean(:,XV) = sum(bsxfun(@times,gamma,shifts),2);
    post_std(:,XV) = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - post_mean(:,XV).^2);
    
    %% RECONSTRUCT MAP STIMULUS (OLD METHOD)
    [max_post,max_loc] = max(lgamma,[],2);
    
    shift_cor(:,XV) = shifts(max_loc);
    % shift_cor(:,XV) = post_mean;
    for ii = 1:NT
        d2 = dist_shift2d(resh_X(:,:,ii), -shift_cor(ii,XV),2,0);
        resh_X_sh(:,:,ii) = d2;
    end
    X_sh = reshape(resh_X_sh,flen*up_SDIM,NT)';
    
    %% REFIT MODELS
    upklen = flen*up_SDIM;
    stim_params.spatial_dims = 1;
    stim_params.sdim = up_SDIM;
    stim_params.flen = flen;
    tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,xv_cell));
    
    quad_fit_sh(xv_cell) = fitGNM_filters(upquad_fit(xv_cell),X_sh,tr_spkbns,'none',[],1e-4,1e-6,1);
    
    ref_LL(1,xv_cell) = quad_fit_sh(xv_cell).LL;
    
    fprintf('Original: %.4f  New: %.4f\n',null_LL(xv_cell)-up_LL(xv_cell),null_LL(xv_cell)-ref_LL(1,xv_cell));
    
    %%
    [nll, pnll, lpen, prate, g] = getLL_GNM(quad_fit_sh(xv_cell),X_sh,tr_spkbns,'none');
    Xmat = [rsac_Tmat expt_Tmat g];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    K0(end-1) = 1;
    [fitp_sac_kerns(xv_cell)] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], lamrange2,[], [], [], [klen klen+1], 0);
    rev_rsac_kern(xv_cell,:) = fitp_sac_kerns(xv_cell).k(1:ntents);
    rev_expt_kern(xv_cell,:) = fitp_sac_kerns(xv_cell).k(ntents+1:end-2);
    
end


%% NOW FOR ALL NON XV CELLS (EYE-TRACKING USING ALL CELLS).
for jjj = 1
    cur_tr_inds = 1:length(tr_set);
    non_xv_cells = setdiff(tr_set,poss_xv_set);
    
    fprintf('Using all cells\n');
    frame_LLs = squeeze(sum(LLs(:,cur_tr_inds,:),2));
    
    lB = frame_LLs;
    lalpha=zeros(NT,n_shifts);
    lbeta = zeros(NT,n_shifts);
    lscale=zeros(NT,1); %initialize rescaling parameters
    %compute rescaled forward messages
    lalpha(1,:) = leps_prior + lB(1,:);
    lscale(1)=logsumexp(lalpha(1,:));
    lalpha(1,:) = lalpha(1,:) - lscale(1);
    for t=2:NT
        if use_prior(t)
            cur_lA = lA_tflip;
        else
            cur_lA = lA;
        end
        lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + lB(t,:);
        lscale(t) = logsumexp(lalpha(t,:));
        lalpha(t,:)= lalpha(t,:) - lscale(t);
    end
    
    %compute rescaled backward messages
    lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
    for t=NT-1:-1:1
        if use_prior(t+1)
            cur_lA = lA_tflip;
        else
            cur_lA = lA;
        end
        lf1 = lbeta(t+1,:) + lB(t+1,:);
        lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
    end
    lgamma= lalpha + lbeta;
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    
    gamma = exp(lgamma);
    fpost_mean = sum(bsxfun(@times,gamma,shifts),2);
    fpost_std = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - fpost_mean.^2);
    
    [max_post,max_loc] = max(lgamma,[],2);
    
    fshift_cor = shifts(max_loc);
    % shift_cor(:,XV) = post_mean;
    for ii = 1:NT
        d2 = dist_shift2d(resh_X(:,:,ii), -fshift_cor(ii),2,0);
        resh_X_sh(:,:,ii) = d2;
    end
    X_sh = reshape(resh_X_sh,flen*up_SDIM,NT)';
    
    for cc = 1:length(non_xv_cells)
        fprintf('Fitting non xv-cell %d of %d\n',cc,length(non_xv_cells));
        
        tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,non_xv_cells(cc)));
        
        quad_fit_sh(non_xv_cells(cc)) = fitGNM_filters(upquad_fit(non_xv_cells(cc)),X_sh,tr_spkbns,'none',[],1e-4,1e-6,1);
        
        ref_LL(1,non_xv_cells(cc)) = quad_fit_sh(non_xv_cells(cc)).LL;
        
        fprintf('Original: %.4f  New: %.4f\n',null_LL(non_xv_cells(cc))-up_LL(non_xv_cells(cc)),null_LL(non_xv_cells(cc))-ref_LL(1,non_xv_cells(cc)));
        
        [nll, pnll, lpen, prate, g] = getLL_GNM(quad_fit_sh(non_xv_cells(cc)),X_sh,tr_spkbns,'none');
        Xmat = [rsac_Tmat expt_Tmat g];
        klen = size(Xmat,2);
        K0 = zeros(klen+1,1);
        K0(end-1) = 1;
        [fitp_sac_kerns(non_xv_cells(cc))] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], lamrange2,[], [], [], [klen klen+1], 0);
        rev_rsac_kern(non_xv_cells(cc),:) = fitp_sac_kerns(non_xv_cells(cc)).k(1:ntents);
        rev_expt_kern(non_xv_cells(cc),:) = fitp_sac_kerns(non_xv_cells(cc)).k(ntents+1:end-2);
    end
end

%% PLOT IMPROVEMENT
old_imp = (null_LL- up_LL)/log(2);
new_imp = (null_LL(poss_xv_set) - ref_LL(1,poss_xv_set))/log(2);
new_imp_nonxv = (null_LL(non_xv_cells) - ref_LL(1,non_xv_cells))/log(2);
mdl = LinearModel.fit(old_imp(poss_xv_set)',new_imp');
xx = linspace(0,0.2,100);
[ypred,pred_errs] = predict(mdl,xx');
figure;hold on
plot(xx,ypred,'k','linewidth',2);plot(xx,pred_errs,'r--')
plot(old_imp(poss_xv_set),new_imp,'o')
plot(old_imp(non_xv_cells),new_imp_nonxv,'ko')
line([0 0.2],[0 0.2])
xlim([0 0.25]); ylim([0 0.25])

beta = mdl.Coefficients.Estimate;
fprintf('%.3f-fold improvement on xval\n',beta(2));
mdl = LinearModel.fit(old_imp(non_xv_cells)',new_imp_nonxv');
beta = mdl.Coefficients.Estimate;
fprintf('%.3f-fold improvement on others\n',beta(2));

%% INITIALIZE FOR ITERATIVE EYE_CORRECTIONS
it_shift_cor{1} = shift_cor;
it_quad_fit{1} = quad_fit_sh;
it_rsac_kern{1} = rev_rsac_kern;
it_expt_kern{1} = rev_expt_kern;
it_post_mean{1} = post_mean;
it_post_std{1} = post_std;
it_fpost_mean{1} = fpost_mean;
it_fpost_std{1} = fpost_std;
it_fshift_cor{1} = fshift_cor;
%% NOW ITERATE

n_ITER = 15;
for nn = 2:n_ITER+1
    
    fprintf('Eye correction iteration %d of %d\n',nn-1,n_ITER);
    
    lfilt_bank = zeros(length(tr_set),up_klen);
    qfilt1_bank = zeros(length(tr_set),up_klen);
    qfilt2_bank = zeros(length(tr_set),up_klen);
    for cc = 1:length(tr_set)
        cur_k = get_k_mat(it_quad_fit{nn-1}(tr_set(cc)));
        lfilt_bank(cc,:) = cur_k(:,1);
        qfilt1_bank(cc,:) = cur_k(:,2);
        if use_mod(tr_set(cc)) == 3
            qfilt2_bank(cc,:) = cur_k(:,3);
        end
    end
    lfilt_bank = reshape(lfilt_bank',[flen,up_SDIM,length(tr_set)]);
    qfilt1_bank = reshape(qfilt1_bank',[flen,up_SDIM,length(tr_set)]);
    qfilt2_bank = reshape(qfilt2_bank',[flen,up_SDIM,length(tr_set)]);
    
    shifted_lfilt_bank = nan(klen,n_tr_chs);
    shifted_qfilt1_bank = nan(klen,n_tr_chs);
    shifted_qfilt2_bank = nan(klen,n_tr_chs);
    
    rsac_out = rsac_Tmat*init_rsac_kern(tr_set,:)';
    expt_out = expt_Tmat*init_expt_kern(tr_set,:)';
    
    LLs = nan(NT,length(tr_set),n_shifts);
    for xx = 1:length(shifts)
        fprintf('Shift %d of %d\n',xx,n_shifts);
        d2 = dist_shift2d(lfilt_bank,shifts(xx),2,0);
        shifted_lfilt_bank = reshape(d2,up_klen,n_tr_chs);
        d2 = dist_shift2d(qfilt1_bank,shifts(xx),2,0);
        shifted_qfilt1_bank = reshape(d2,up_klen,n_tr_chs);
        d2 = dist_shift2d(qfilt2_bank,shifts(xx),2,0);
        shifted_qfilt2_bank = reshape(d2,up_klen,n_tr_chs);
        
        lin_out = X_z*shifted_lfilt_bank;
        quad_out = (X_z*shifted_qfilt1_bank).^2 + (X_z*shifted_qfilt2_bank).^2;
        gfun = lin_out + quad_out;
        gfun = bsxfun(@plus,gfun,pre_offset(tr_set));
        
        %     %weight stimulus term
        %     gfun = bsxfun(@times,gfun,init_gweight(tr_set));
        
        %add contributions from event kernels
        gfun = gfun + rsac_out + expt_out;
        
        too_large = gfun > 50;
        pred_rate = log(1+exp(gfun));
        pred_rate(too_large) = gfun(too_large);
        
        pred_rate(pred_rate < 1e-20) = 1e-20;
        
        LLs(:,:,xx) = Robs.*log(pred_rate) - pred_rate;
    end
    
    for XV = 1:length(poss_xv_set)
        xv_cell = poss_xv_set(XV);
        
        cur_tr_inds = find(tr_set~=xv_cell);
        
        fprintf('Using XV cell %d,  %d of %d\n',xv_cell,XV,length(poss_xv_set));
        frame_LLs = squeeze(sum(LLs(:,cur_tr_inds,:),2));
        
        %%
        lB = frame_LLs;
        lalpha=zeros(NT,n_shifts);
        lbeta = zeros(NT,n_shifts);
        lscale=zeros(NT,1); %initialize rescaling parameters
        %compute rescaled forward messages
        lalpha(1,:) = leps_prior + lB(1,:);
        lscale(1)=logsumexp(lalpha(1,:));
        lalpha(1,:) = lalpha(1,:) - lscale(1);
        for t=2:NT
            if use_prior(t)
                cur_lA = lA_tflip;
            else
                cur_lA = lA;
            end
            lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + lB(t,:);
            lscale(t) = logsumexp(lalpha(t,:));
            lalpha(t,:)= lalpha(t,:) - lscale(t);
        end
        
        %compute rescaled backward messages
        lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
        for t=NT-1:-1:1
            if use_prior(t)
                cur_lA = lA_tflip;
            else
                cur_lA = lA;
            end
            lf1 = lbeta(t+1,:) + lB(t+1,:);
            lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
        end
        lgamma= lalpha + lbeta;
        lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
        
        gamma = exp(lgamma);
        post_mean(:,XV) = sum(bsxfun(@times,gamma,shifts),2);
        post_std(:,XV) = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - post_mean(:,XV).^2);
        
        %% RECONSTRUCT MAP STIMULUS (OLD METHOD)
        [max_post,max_loc] = max(lgamma,[],2);
        shift_cor(:,XV) = shifts(max_loc);
        
        resh_X_sh = zeros(size(resh_X));
        for ii = 1:NT
            d2 = dist_shift2d(resh_X(:,:,ii), -shift_cor(ii,XV),2,0);
            resh_X_sh(:,:,ii) = d2;
        end
        X_sh = reshape(resh_X_sh,flen*up_SDIM,NT)';
        
        %% REFIT MODELS
        tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,xv_cell));
        
        it_quad_fit{nn}(xv_cell) = fitGNM_filters(it_quad_fit{nn-1}(xv_cell),X_sh,tr_spkbns,'none',[],1e-4,1e-6,1);
        
        ref_LL(nn,xv_cell) = it_quad_fit{nn}(xv_cell).LL;
        
        fprintf('Original: %.4f  First: %.4f  Previous: %.4f   Current: %.4f\n',...
            null_LL(xv_cell)-up_LL(xv_cell),null_LL(xv_cell)-it_quad_fit{1}(xv_cell).LL,...
            null_LL(xv_cell)-it_quad_fit{nn-1}(xv_cell).LL,...
            null_LL(xv_cell)-it_quad_fit{nn}(xv_cell).LL);
        
        %%
        [nll, pnll, lpen, prate, g] = getLL_GNM(it_quad_fit{nn}(xv_cell),X_sh,tr_spkbns,'none');
        Xmat = [rsac_Tmat expt_Tmat g];
        klen = size(Xmat,2);
        K0 = zeros(klen+1,1);
        K0(end-1) = 1;
        [fitp_sac_kerns(xv_cell)] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], lamrange2,[], [], [], [klen klen+1], 0);
        it_rsac_kern{nn}(xv_cell,:) = fitp_sac_kerns(xv_cell).k(1:ntents);
        it_expt_kern{nn}(xv_cell,:) = fitp_sac_kerns(xv_cell).k(ntents+1:end-2);
        
    end
    
    %% NOW FOR ALL NON XV CELLS
    for jjj = 1
        cur_tr_inds = 1:length(tr_set);
        non_xv_cells = setdiff(tr_set,poss_xv_set);
        
        fprintf('Using all cells\n');
        frame_LLs = squeeze(sum(LLs(:,cur_tr_inds,:),2));
        
        lB = frame_LLs;
        lalpha=zeros(NT,n_shifts);
        lbeta = zeros(NT,n_shifts);
        lscale=zeros(NT,1); %initialize rescaling parameters
        %compute rescaled forward messages
        lalpha(1,:) = leps_prior + lB(1,:);
        lscale(1)=logsumexp(lalpha(1,:));
        lalpha(1,:) = lalpha(1,:) - lscale(1);
        for t=2:NT
            if use_prior(t)
                cur_lA = lA_tflip;
            else
                cur_lA = lA;
            end
            lalpha(t,:) = logmulexp(lalpha(t-1,:),cur_lA) + lB(t,:);
            lscale(t) = logsumexp(lalpha(t,:));
            lalpha(t,:)= lalpha(t,:) - lscale(t);
        end
        
        %compute rescaled backward messages
        lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
        for t=NT-1:-1:1
            if use_prior(t+1)
                cur_lA = lA_tflip;
            else
                cur_lA = lA;
            end
            lf1 = lbeta(t+1,:) + lB(t+1,:);
            lbeta(t,:) = logmulexp(lf1,cur_lA') - lscale(t);
        end
        lgamma= lalpha + lbeta;
        lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
        
        gamma = exp(lgamma);
        it_fpost_mean{nn} = sum(bsxfun(@times,gamma,shifts),2);
        it_fpost_std{nn} = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - it_fpost_mean{nn}.^2);
        
        [max_post,max_loc] = max(lgamma,[],2);
        
        it_fshift_cor{nn} = shifts(max_loc);
        % shift_cor(:,XV) = post_mean;
        for ii = 1:NT
            d2 = dist_shift2d(resh_X(:,:,ii), -it_fshift_cor{nn}(ii),2,0);
            resh_X_sh(:,:,ii) = d2;
        end
        X_sh = reshape(resh_X_sh,flen*up_SDIM,NT)';
        
        for cc = 1:length(non_xv_cells)
            fprintf('Fitting non xv-cell %d of %d\n',cc,length(non_xv_cells));
            
            tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,non_xv_cells(cc)));
            
            it_quad_fit{nn}(non_xv_cells(cc)) = fitGNM_filters(it_quad_fit{nn-1}(non_xv_cells(cc)),X_sh,tr_spkbns,'none',[],1e-4,1e-6,1);
            
            ref_LL(nn,non_xv_cells(cc)) = it_quad_fit{nn}(non_xv_cells(cc)).LL;
            
            fprintf('Original: %.4f  First: %.4f  Previous: %.4f   Current: %.4f\n',...
                null_LL(non_xv_cells(cc))-up_LL(non_xv_cells(cc)),null_LL(non_xv_cells(cc))-it_quad_fit{1}(non_xv_cells(cc)).LL,...
                null_LL(non_xv_cells(cc))-it_quad_fit{nn-1}(non_xv_cells(cc)).LL,...
                null_LL(non_xv_cells(cc))-it_quad_fit{nn}(non_xv_cells(cc)).LL);
            
            [nll, pnll, lpen, prate, g] = getLL_GNM(it_quad_fit{nn}(non_xv_cells(cc)),X_sh,tr_spkbns,'none');
            Xmat = [rsac_Tmat expt_Tmat g];
            klen = size(Xmat,2);
            K0 = zeros(klen+1,1);
            K0(end-1) = 1;
            [fitp_sac_kerns(non_xv_cells(cc))] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], lamrange2,[], [], [], [klen klen+1], 0);
            it_rsac_kern{nn}(non_xv_cells(cc),:) = fitp_sac_kerns(non_xv_cells(cc)).k(1:ntents);
            it_expt_kern{nn}(non_xv_cells(cc),:) = fitp_sac_kerns(non_xv_cells(cc)).k(ntents+1:end-2);
        end
    end
    
    %%
    it_post_mean{nn} = post_mean;
    it_post_std{nn} = post_std;
    it_shift_cor{nn} = shift_cor;
    
    %% PLOT PROGRESS
    old_imp = (null_LL(poss_xv_set)- up_LL(poss_xv_set))/log(2);
    old_imp_nonxv = (null_LL(non_xv_cells) - up_LL(non_xv_cells))/log(2);
    new_imp = (null_LL(poss_xv_set) - ref_LL(end,poss_xv_set))/log(2);
    new_imp_nonxv = (null_LL(non_xv_cells) - ref_LL(end,non_xv_cells))/log(2);
    mdl = LinearModel.fit(old_imp',new_imp');
    beta = mdl.Coefficients.Estimate;
    fprintf('%.3f-fold overall improvement on xval\n',beta(2));
    mdl2 = LinearModel.fit(old_imp_nonxv',new_imp_nonxv');
    beta = mdl2.Coefficients.Estimate;
    fprintf('%.3f-fold overall improvement on others\n',beta(2));

%         xx = linspace(0,0.5,100);
%     [ypred,pred_errs] = predict(mdl,xx');
%     figure;
%     subplot(2,1,1)
%     hold on
%     plot(xx,ypred,'k','linewidth',2);plot(xx,pred_errs,'r--')
%     plot(old_imp,new_imp,'o')
%     plot(old_imp_nonxv,new_imp_nonxv,'k*')
%     line([0 0.2],[0 0.2])
%     xlim([0 0.5]); ylim([0 0.5])
    
    
    prev_imp = (null_LL(poss_xv_set) - ref_LL(end-1,poss_xv_set))/log(2);
    prev_imp_nonxv = (null_LL(non_xv_cells) - ref_LL(end-1,non_xv_cells))/log(2);
    mdl = LinearModel.fit(prev_imp',new_imp');
    beta = mdl.Coefficients.Estimate;
    fprintf('%.3f-fold previous improvement on xval\n',beta(2));
    mdl2 = LinearModel.fit(prev_imp_nonxv',new_imp_nonxv');
    beta = mdl2.Coefficients.Estimate;
    fprintf('%.3f-fold previous improvement on others\n',beta(2));

%     xx = linspace(0,0.2,100);
%     [ypred,pred_errs] = predict(mdl,xx');
%     subplot(2,1,2)
%     hold on
%     plot(xx,ypred,'k','linewidth',2);plot(xx,pred_errs,'r--')
%     plot(prev_imp,new_imp,'o')
%     plot(prev_imp_nonxv,new_imp_nonxv,'k*')
%     line([0 0.2],[0 0.2])
%     xlim([0 0.25]); ylim([0 0.25])

    
    %% SAVE PROGRESS
    eye_times = all_stim_times(tr_inds);
    save full_eye_correct_v6_90deg *_LL it* tr_set non_xv_cells tr_inds eye_times
    
end
%%
eye_times = all_stim_times(tr_inds);
save full_eye_correct_v6_90deg *_LL it* tr_set non_xv_cells tr_inds eye_times

%%
pre_improve = (null_LL-old_LL)/log(2);
post_improve = (null_LL - ref_LL(end,:))/log(2);

figure
% plot(pre_improve,post_improve,'o','markersize',8)
line([0 0.5],[0 0.5],'color','k')
hold on
plot(pre_improve(tr_set),post_improve(tr_set),'r*')
plot(pre_improve(xv_set),post_improve(xv_set),'k*')
pp = polyfit(pre_improve(tr_set),post_improve(tr_set),1);
xx = linspace(0.001,0.5,100);
plot(xx,polyval(pp,xx),'r')
pp = polyfit(pre_improve(xv_set),post_improve(xv_set),1);
plot(xx,polyval(pp,xx),'k')
xlim([0 0.6])
ylim([0 0.6])
% set(gca,'xscale','log')
% set(gca,'yscale','log')
xlabel('Pre-correction LL improvement (bits/spk)','fontsize',16)
ylabel('Post-correction LL improvement (bits/spk)','fontsize',16)

xv_lim = xv_set; xv_lim(xv_lim==13) = [];
clear cur_slope
for i = 1:n_ITER+1
    cur_post_improve = (null_LL - ref_LL(i,:))/log(2);
    pp = polyfit(pre_improve,cur_post_improve,1);
    cur_slope(i) = pp(1);
    pp = polyfit(pre_improve(xv_set),cur_post_improve(xv_set),1);
    cur_slope_xv(i) = pp(1);
    pp = polyfit(pre_improve(xv_lim),cur_post_improve(xv_lim),1);
    cur_slope_xvl(i) = pp(1);
    pp = polyfit(pre_improve(tr_set),cur_post_improve(tr_set),1);
    cur_slope_tr(i) = pp(1);
end
%%

%%
dt = 0.01;
up_bar_pos = (-floor(up_SDIM/2):floor(up_SDIM/2))*0.04;
close all
for cc = tr_set
%     if ismember(cc,tr_set)
%         fprintf('TR Cell %d of %d\n',cc,96);
%     else
%         fprintf('XV Cell %d of %d\n',cc,96);
%     end
    k = get_k_mat(upquad_fit(cc));
    subplot(3,2,1)
    %     imagesc(up_bar_pos,-(0:flen)*dt,[zeros(flen,zpads) reshape(k(:,1),flen,up_SDIM) zeros(flen,zpads)]);
    imagesc(up_bar_pos,-(0:flen)*dt,reshape(k(:,1),flen,up_SDIM));
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    title('Pre Linear')
    xlabel('Relative bar position (deg)')
    ylabel('Time lag (s)')
    subplot(3,2,3)
    %     imagesc(up_bar_pos,-(0:flen)*dt,[zeros(flen,zpads) reshape(k(:,2),flen,up_SDIM) zeros(flen,zpads)])
    imagesc(up_bar_pos,-(0:flen)*dt,reshape(k(:,2),flen,up_SDIM))
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    title('Pre Quad')
    xlabel('Relative bar position (deg)')
    ylabel('Time lag (s)')
    if use_mod(cc) == 3
    subplot(3,2,5)
    %     imagesc(up_bar_pos,-(0:flen)*dt,[zeros(flen,zpads) reshape(k(:,2),flen,up_SDIM) zeros(flen,zpads)])
    imagesc(up_bar_pos,-(0:flen)*dt,reshape(k(:,3),flen,up_SDIM))
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    title('Pre Quad')
    xlabel('Relative bar position (deg)')
    ylabel('Time lag (s)')
    end
    
    k = get_k_mat(it_quad_fit{end}(cc));
    subplot(3,2,2)
    imagesc(up_bar_pos,-(0:flen)*dt,reshape(k(:,1),flen,up_SDIM));
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    title('Post Linear')
    xlabel('Relative bar position (deg)')
    ylabel('Time lag (s)')
    subplot(3,2,1)
    caxis([-ca ca]);colorbar
    
    subplot(3,2,4)
    imagesc(up_bar_pos,-(0:flen)*dt,reshape(k(:,2),flen,up_SDIM))
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    title('Post Quad')
    xlabel('Relative bar position (deg)')
    ylabel('Time lag (s)')
    subplot(3,2,3)
    caxis([-ca ca]);colorbar
    
    
    if use_mod(cc) == 3
    subplot(3,2,6)
    imagesc(up_bar_pos,-(0:flen)*dt,reshape(k(:,3),flen,up_SDIM))
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    title('Post Quad')
    xlabel('Relative bar position (deg)')
    ylabel('Time lag (s)')
    subplot(3,2,5)
    caxis([-ca ca]);colorbar

    caxis([-ca ca]);colorbar
    end
        orig_imp = (null_LL(cc) - up_LL(cc))/log(2);
    new_imp = (null_LL(cc) - ref_LL(end,cc))/log(2);
    fprintf('Orig: %.4f  New: %.4f\n',orig_imp,new_imp);

    pause
    clf
end

%%
% for cc = use_sus
%     if ismember(cc,tr_set)
%         temp_str = 'TR';
%         fprintf('Fitting training cell %d of %d\n',cc,length(use_sus));
%     elseif ismember(cc,xv_set)
%         temp_str = 'XV';
%         fprintf('Fitting xv cell %d of %d\n',cc,length(use_sus));
%     else
%         temp_str = 'CR';
%         fprintf('Fitting crap cell %d of %d\n',cc,length(use_sus));
%     end
%     tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,(cc)));
%     big_quad(cc) = add_gnm_filter(it_quad_fit{nn}(cc),0.01*randn(upklen,1),-1,'quad');
%     big_quad(cc) = add_gnm_filter(big_quad(cc),0.01*randn(upklen,1),1,'quad');
%
%     big_quad(cc) = fitGNM_filters(big_quad(cc),X_sh,tr_spkbns,'none',[],1e-4,1e-6,1);
%     big_quad_LL(cc) = big_quad(cc).LL;
%
%     fprintf('Original: %.3f  Previous: %.3f   Current: %.3f\n',it_quad_fit{1}(cc).LL,it_quad_fit{nn}(cc).LL,big_quad_LL(cc));
%
% end
%
%% %DECONVOLVE ESTIMATE OF TRUE EYE POSITION
inferred_ev1 = zeros(size(all_stim_times));
inferred_ev1(tr_inds) = it_post_mean{1};
inferred_ev2 = zeros(size(all_stim_times));
inferred_ev2(tr_inds) = it_post_mean{end};
inferred_eye_cors1 = conv(inferred_ev1,avg_tfilt_kern,'same');
inferred_eye_cors1 = inferred_eye_cors1(tr_inds);
inferred_eye_cors2 = conv(inferred_ev2,avg_tfilt_kern,'same');
inferred_eye_cors2 = inferred_eye_cors2(tr_inds);

%align to saccades
for i = length(saccade_inds):-1:1
    cur_inds = saccade_inds(i):(saccade_inds(i)+sac_shift);
    inferred_eye_cors1(cur_inds) = it_post_mean{1}(cur_inds(end));
    inferred_eye_cors2(cur_inds) = it_post_mean{end}(cur_inds(end));
    cur_back_inds = (saccade_inds(i)-(flen-sac_shift)):(saccade_inds(i)-1);
    inferred_eye_cors1(cur_back_inds) = it_post_mean{1}(saccade_inds(i));
    inferred_eye_cors2(cur_back_inds) = it_post_mean{end}(saccade_inds(i));
end
for i = length(trial_start_inds):-1:2
    cur_back_inds = (trial_start_inds(i)-flen):(trial_start_inds(i)-1);
    cur_back_inds(cur_back_inds <= 0) = [];
    inferred_eye_cors1(cur_back_inds) = it_post_mean{1}(trial_start_inds(i)-1);
    inferred_eye_cors2(cur_back_inds) = it_post_mean{end}(trial_start_inds(i)-1);
end

%%
close all
figure; hold on
shadedErrorBar(all_stim_times(tr_inds),it_fpost_mean{end}*.04,it_fpost_std{end}*.04,{'r'})
% shadedErrorBar(all_stim_times(tr_inds),it_post_mean{end}*.04,it_inferred_std{4}*.04,{'b'},0)
% shadedErrorBar(all_stim_times(tr_inds),inferred_eye_cors1*.04,it_post_std{1}*.04,{'r'})
% shadedErrorBar(all_stim_times(tr_inds),inferred_eye_cors2*.04,it_post_std{end}*.04,{'b'},0)
% plot(all_stim_times(tr_inds),it_fshift_cor{end}*0.04,'k','linewidth',2)

% plot(all_stim_times(tr_inds),rsac_inds,'r.-')
% plot(all_stim_times(tr_inds),shift_cor*0.04,'.-')
% plot(all_stim_times(tr_inds),interp_eye_vals(tr_inds,2),'k.-');
% plot(all_stim_times(tr_inds),interp_eye_vals(tr_inds,4),'g.-');
% plot(all_stim_times(tr_inds),0.5*interp_eye_vals(tr_inds,4)+0.5*interp_eye_vals(tr_inds,2),'c.-');

plot(all_stim_times(tr_inds),interp_eye_vals(tr_inds,2),'k.-');
plot(all_stim_times(tr_inds),interp_eye_vals(tr_inds,4),'g.-');
plot(all_stim_times(tr_inds),0.5*interp_eye_vals(tr_inds,2)+0.5*interp_eye_vals(tr_inds,4),'c.-');

plot(all_stim_times(interp_saccade_inds),0.1*ones(size(interp_saccade_inds)),'k*','markersize',16)
ylim([-0.75 0.75])
yl = ylim();
for i = 1:length(all_trial_start_times)
    line(all_trial_start_times([i i]),yl,'color',[0.2 0.9 0.2],'linewidth',2)
    line(all_trial_end_times([i i]),yl,'color','r','linewidth',2)
end

for i =500:length(all_trial_start_times)
    xlim([all_trial_start_times(i)-0.1 all_trial_end_times(i)+0.1])
    pause
end
%%
figure
plot(all_stim_times(tr_inds),it_fpost_mean{end}*0.04)
hold on
plot(all_stim_times(tr_inds),interp_eye_vals(tr_inds,2),'k.-');
plot(all_stim_times(tr_inds),interp_eye_vals(tr_inds,4),'g.-');
plot(all_stim_times(interp_saccade_inds),0.1*ones(size(interp_saccade_inds)),'k*','markersize',16)
% plot(all_stim_times(tr_inds(saccade_inds)),0.1*ones(size(saccade_inds)),'g*','markersize',16)

for i = 1:length(saccade_inds)
    pre_inds = (saccade_inds(i) - 6):(saccade_inds(i) - 2);
    post_inds = (saccade_inds(i) + 2):(saccade_inds(i) + 6);
    pre_obspos(i,:) = mean(interp_eye_vals(tr_inds(pre_inds),:));
    post_obspos(i,:) = mean(interp_eye_vals(tr_inds(post_inds),:));
    
% plot(all_stim_times(tr_inds(pre_inds)),interp_eye_vals(tr_inds(pre_inds),2),'ro');
% plot(all_stim_times(tr_inds(pre_inds)),interp_eye_vals(tr_inds(pre_inds),4),'ro');
% plot(all_stim_times(tr_inds(post_inds)),interp_eye_vals(tr_inds(post_inds),2),'go');
% plot(all_stim_times(tr_inds(post_inds)),interp_eye_vals(tr_inds(post_inds),4),'go');
% xlim(all_stim_times(tr_inds(saccade_inds(i))) + [-1 1])
% pause

pre_inds = (saccade_inds(i) -1):(saccade_inds(i) + 3);
    post_inds = (saccade_inds(i) + 7):(saccade_inds(i) + 11);
    pre_infpos(i) = mean(it_fpost_mean{end}(pre_inds,:))*.04;
    post_infpos(i) = mean(it_fpost_mean{end}(post_inds,:))*.04;   
    peri_inds = (saccade_inds(i) + 4):(saccade_inds(i) + 6);
    post_infunc(i) = mean(it_fpost_std{end}(peri_inds));

end
amp_obs = post_obspos - pre_obspos;
amp_inf = post_infpos - pre_infpos;

%%
figure
imagesc(1:NT,shifts*0.04,gamma')
hold on
plot(inferred_eye_cors*0.04,'r-','linewidth',2)
plot(interp_eye_vals(tr_inds,2),'w-','linewidth',2)

%%
pre_improve = (null_LL-old_LL)/log(2);
pre_improve2 = (null_LL-nocor_LL)/log(2);
post_improve = (null_LL - ref_LL(end,:))/log(2);

figure
% plot(pre_improve,post_improve,'o','markersize',8)
line([0 0.5],[0 0.5],'color','k')
hold on
plot(pre_improve,post_improve,'ro')
plot(pre_improve2,post_improve,'o')
pp = polyfit(pre_improve,post_improve,1);
xx = linspace(0.001,0.5,100);
plot(xx,polyval(pp,xx),'g')
xlim([0 0.5])
ylim([0 0.5])
% set(gca,'xscale','log')
% set(gca,'yscale','log')
xlabel('Pre-correction LL improvement (bits/spk)','fontsize',16)
ylabel('Post-correction LL improvement (bits/spk)','fontsize',16)

%% RECOMPUTE BASE LLs
for cc = 1:length(tr_set)
    fprintf('Fitting cell %d of %d\n',cc,length(tr_set));
    
    tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,tr_set(cc)));
    
    new_base_mod(tr_set(cc)) = fitGNM_filters(it_quad_fit{end}(tr_set(cc)),X_z,tr_spkbns,'none',[],1e-4,1e-6,1);
    
    new_base_LL(tr_set(cc)) = new_base_mod(tr_set(cc)).LL;
    
end
%%
old_imp = (null_LL(poss_xv_set)- new_base_LL(poss_xv_set))/log(2);
old_imp_nonxv = (null_LL(non_xv_cells) - new_base_LL(non_xv_cells))/log(2);
new_imp = (null_LL(poss_xv_set) - ref_LL(end,poss_xv_set))/log(2);
new_imp_nonxv = (null_LL(non_xv_cells) - ref_LL(end,non_xv_cells))/log(2);
mdl = LinearModel.fit(old_imp',new_imp');
beta = mdl.Coefficients.Estimate;
fprintf('%.3f-fold overall improvement on xval\n',beta(2));
mdl2 = LinearModel.fit(old_imp_nonxv',new_imp_nonxv');
beta = mdl2.Coefficients.Estimate;
fprintf('%.3f-fold overall improvement on others\n',beta(2));

xx = linspace(0,0.5,100);
[ypred,pred_errs] = predict(mdl,xx');
figure;
hold on
plot(xx,ypred,'k','linewidth',2);plot(xx,pred_errs,'r--')
plot(old_imp,new_imp,'o')
plot(old_imp_nonxv,new_imp_nonxv,'k*')
line([0 0.2],[0 0.2])
xlim([0 1]); ylim([0 1])

%%
% save nocorr_mod_135 nocor_LL test_nocor_mod
save full_eye_correct_v6_0deg *_LL it* tr_set non_xv_cells tr_inds eye_times new_*
