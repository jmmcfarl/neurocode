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

%% USE ONLY GRAY BACKGROUND DATA
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

sac_thresh = 10;
min_isi = 0.05;
orig_saccade_inds = 1 + find(all_eye_speed(1:end-1) < sac_thresh & all_eye_speed(2:end) > sac_thresh);
orig_saccade_inds = unique(orig_saccade_inds);
isis = [Inf; diff(orig_saccade_inds)]/eye_fs;
double_sacs = find(isis < min_isi);
orig_saccade_inds(double_sacs) = [];
orig_saccade_inds(all_eye_intrial(orig_saccade_inds) == 0) = [];
interp_saccade_inds = round(interp1(all_stim_times,1:length(all_stim_times),all_eye_ts(orig_saccade_inds)));


%% IF YOU WANT TO DO CROSS-VALIDATION ON INITIAL MODEL FITS, OTHERWISE,
[c,ia,ic] = unique([all_exptvec all_trialvec],'rows');
n_trials = length(ia);
bad_trials = unique(ic(bad_pts)); %trials with putative blinks

xv_frac = 0.2;
xv_set = randperm(n_trials);
xv_set(round(n_trials*xv_frac)+1:end) = [];

tr_set = setdiff(1:n_trials,xv_set);

tr_inds = find(ismember(ic,tr_set));
xv_inds = find(ismember(ic,xv_set));
tr_inds(all_used_inds(tr_inds) == 0) = [];
xv_inds(all_used_inds(xv_inds) == 0) = [];

%% double sample Xmat
bar_dx = 0.08;
up_fac = 2;
bar_dxd = bar_dx/up_fac;
SDIM = n_bar_pos;
full_NT = length(all_stim_times);
NT = length(tr_inds);
xvNT = length(xv_inds);
resh_X = reshape(all_bar_mat',[flen SDIM full_NT]);

fprintf('Up-sampling X-matrix\n');
dresh_X = [resh_X(:,1,:)];
for i = 1:SDIM-1
    dresh_X = cat(2,dresh_X,resh_X(:,i,:));
    %for interpolating bar positions
    dresh_X = cat(2,dresh_X,resh_X(:,i,:) + resh_X(:,i+1,:));
    %adding zeros for absent bar positions
    %         dresh_X = cat(2,dresh_X,zeros(flen,1,NT));
end
dresh_X = cat(2,dresh_X,resh_X(:,end,:));
dresh_X = cat(2,dresh_X,resh_X(:,end,:));

up_SDIM = up_fac*SDIM-1;

%%
up_xax = 1:up_SDIM;
tent_centers = [2:2:up_SDIM-1];
tbmat = construct_tent_bases_v2(tent_centers,up_xax);
tbmat(1,1) = 1;
tbmat(end,end) = 1;
ntents = length(tent_centers);

X_tb = nan(flen,full_NT,ntents);
for i = 1:flen
   X_tb(i,:,:) = squeeze(dresh_X(i,:,:))'*tbmat'; 
end
X_tb = permute(X_tb,[1 3 2]);

X_tb = reshape(X_tb,flen*ntents,full_NT)';

%%
SDIM = ntents;
klen = flen*SDIM;
stim_params.spatial_dims = 1;
stim_params.sdim = SDIM;
stim_params.flen = flen;

for cc = 1:96
    fprintf('Fitting cell %d of %d\n',cc,96);
    tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,cc));
    
    clear defmod
    defmod.lambda_L2x = 0;
    defmod.lambda_L1x = 10;
    defmod.lambda_d2XT = 2000;
    defmod.lambda_d2X = 0;
    defmod.lambda_d2T = 2000;
    npq = 2;
    nnq = 0;
    clear kern_types
    kern_types{1} = 'lin';
    for i = 2:(1+npq+nnq)
        kern_types{i} = 'quad';
    end
    init_kerns = 0.01*randn(klen,1+npq+nnq);
    init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
    quad = createGNM(init_kerns,[1 ones(1,npq) -1*ones(1,nnq)],kern_types,defmod,stim_params);
    for i = 2:(1+npq+nnq)
        quad.mods(i).lambda_L1x = 5;
        quad.mods(i).lambda_L2x = 0;
        quad.mods(i).lambda_d2XT = 100;
        quad.mods(i).lambda_d2X = 0;
        quad.mods(i).lambda_d2T = 200;
    end
    quad_fit(cc) = fitGNM_filters(quad,X_tb(tr_inds,:),tr_spkbns,'none',[],1e-4,1e-6,1);
    
    [orig_LL(cc)] = getLL_GNM(quad_fit(cc),X_tb(tr_inds,:),tr_spkbns,'none');
    
    Robs = all_binned_spks(tr_inds,cc);
    avg_rate = mean(Robs);
    null_prate = ones(length(tr_inds),1)*avg_rate;
    null_LL(cc) = -sum(Robs.*log(null_prate) - null_prate)/sum(Robs);
    
    if xv_frac > 0
        xv_Robs = all_binned_spks(xv_inds,use_sus(cc));
        xv_spkbns = convert_to_spikebins(xv_Robs);
        null_prate = ones(length(xv_inds),1)*avg_rate;
        null_xvLL(cc) = -sum(xv_Robs.*log(null_prate) - null_prate)/sum(xv_Robs);
        
        orig_xvLL(cc) = getLL_GNM(quad_fit(cc),X_tb(xv_inds,:),xv_spkbns,'none');
        fprintf('xv Imp: %.3f\n',null_xvLL(cc)-orig_xvLL(cc));
    end
    fprintf('tr Imp: %.3f\n',null_LL(cc)-orig_LL(cc));
end


