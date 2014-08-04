clear all
% close all
cd ~/Data/bruce/G081/
load jbeG081Expts.mat
addpath('~/James_scripts/bruce/G081/');
%%
stim_fs = 100; %in Hz
use_sus = 1:96;
use_lfps = 1:96;

Fs = 3e4;
dsf = 120;Fsd = Fs/dsf;
% dsf = 80;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[1]/niqf,'high');

scales = logspace(log10(15),log10(125),20);
scales = scales*60/dsf;
scales(1:3) = []; %only up to 40 Hz...
% scales = scales*2;
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

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

%%
cd ~/Data/bruce/G081/
all_interp_phasegrams = [];
all_interp_ampgrams = [];
for ee = 1:length(cur_expt_set);
    fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
        
    Vmat = [];
    phasegrams = [];
    ampgrams = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d of %d\n',ll,length(use_lfps));
        filename = sprintf('Expt%d.p%dFullV.mat',cur_expt_set(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
%         V = V + FullV.sumscale*sumv;
        V = V*FullV.intscale(1)/FullV.intscale(2);
        nparts = length(FullV.blklen);
        dVf = [];
        cur_pt = 1;
        for pp = 1:nparts
            cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
            cur_range(cur_range > length(V)) = [];
            dVf = [dVf filtfilt(filt_b,filt_a,decimate(V(cur_range),dsf))]; %do some high-pass filtering
            cur_pt = cur_pt + FullV.blklen(pp);
        end
        Vmat(ll,:) = dVf;
        temp = cwt(dVf,scales,'cmor1-1');
        phasegrams(:,:,ll) = angle(temp)';
        ampgrams(:,:,ll) = abs(temp)';
    end
    t_ax = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        t_ax = [t_ax downsample(cur_t_ax,dsf)];
    end
    t_ax(size(Vmat,2)+1:end) = [];
        
    cur_set = find(all_exptvec == ee);
    
    unwr_phasegram = unwrap(phasegrams);
    interp_phasegrams = interp1(t_ax,unwr_phasegram,all_stim_times(cur_set));
    interp_phasegrams = mod(interp_phasegrams+pi,2*pi)-pi;
    interp_ampgrams = interp1(t_ax,ampgrams,all_stim_times(cur_set));
    
    all_interp_phasegrams = cat(1,all_interp_phasegrams,interp_phasegrams);
    all_interp_ampgrams = cat(1,all_interp_ampgrams,interp_ampgrams);
    clear unwr_phasegram phasegrams interp_phasegrams interp_ampgrams ampgrams Vmat
    
end

% cd ~/James_scripts/bruce/G081/
%%
all_interp_ampgrams = bsxfun(@rdivide,all_interp_ampgrams,nanstd(all_interp_ampgrams));

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
% resh_X = reshape(Xmat(tr_inds,:)',[flen SDIM NT]);
resh_X = reshape(0.5*Xmat(tr_inds,:)',[flen SDIM NT]);

if up_fac > 1
    fprintf('Up-sampling X-matrix\n');
    dresh_X = [];
    for i = 1:SDIM-1
        dresh_X = cat(2,dresh_X,resh_X(:,i,:));
        %for interpolating bar positions
        dresh_X = cat(2,dresh_X,0.5*resh_X(:,i,:) + 0.5*resh_X(:,i+1,:));
        %adding zeros for absent bar positions
%         dresh_X = cat(2,dresh_X,zeros(flen,1,NT));
    end
    dresh_X = cat(2,dresh_X,resh_X(:,end,:));
    
    up_SDIM = up_fac*SDIM-1;
else
    up_SDIM = SDIM;
    dresh_X = resh_X;
end
X_z = reshape(dresh_X,flen*up_SDIM,NT)';

%% LOAD AND INTERPOLATE MODELS
load ./CellList.mat
good_sus = find(all(CellList(:,:,1) > 0));

load ./perfect_fixation_all.mat
% load ./perfect_fixation_all_45.mat
% load ./perfect_fixation_all_90.mat
% load ./perfect_fixation_all_135.mat

klen = flen*SDIM;
lfilt_bank = nan(length(use_sus),klen);
qfilt_bank = nan(length(use_sus),klen);
mod_offset = nan(length(use_sus),1);
for cc = 1:length(use_sus)
    cur_k = get_k_mat(quad_fit(cc));
    lfilt_bank(cc,:) = cur_k(:,1);
    qfilt_bank(cc,:) = cur_k(:,2);
    mod_offset(cc) = quad_fit(cc).spk_theta;
end
resh_lfilt_bank = reshape(lfilt_bank',[flen SDIM 96]);
resh_qfilt_bank = reshape(qfilt_bank',[flen SDIM 96]);
dresh_lf = [];
dresh_qf = [];
for i = 1:SDIM-1
    dresh_lf = cat(2,dresh_lf,resh_lfilt_bank(:,i,:));
    dresh_lf = cat(2,dresh_lf,0.5*resh_lfilt_bank(:,i,:) + 0.5*resh_lfilt_bank(:,i+1,:));
    dresh_qf = cat(2,dresh_qf,resh_qfilt_bank(:,i,:));
    dresh_qf = cat(2,dresh_qf,0.5*resh_qfilt_bank(:,i,:) + 0.5*resh_qfilt_bank(:,i+1,:));
end
dresh_lf = cat(2,dresh_lf,resh_lfilt_bank(:,end,:));
dresh_qf = cat(2,dresh_qf,resh_qfilt_bank(:,end,:));

% for zero Xmatrix
lambdal_L1x = 25;
lambdal_d2X = 1000;
lambdal_d2T = 1000;
lambdal_d2XT = 1000;
lambdaq_L1x = 5;
lambdaq_d2X = 25;
lambdaq_d2T = 25;
lambdaq_d2XT = 25;
% 
% %for interpolated Xmatrix
% lambdal_L1x = 40;
% lambdal_d2X = 5000;
% lambdal_d2T = 5000;
% lambdal_d2XT = 5000;
% lambdaq_L1x = 10;
% lambdaq_d2X = 30;
% lambdaq_d2T = 30;
% lambdaq_d2XT = 30;


% load ./interp_quad_fits_0deg_v4

% REALLY FIT INITIAL MODELS
for cc = 1:96
    interp_quad_fit(cc) = quad_fit(cc);
    interp_quad_fit(cc).stim_params.sdim = up_SDIM;
    interp_quad_fit(cc).stim_params.fsdim = up_SDIM;
    interp_quad_fit(cc).mods(1).k = reshape(dresh_lf(:,:,cc),up_SDIM*flen,1);
    interp_quad_fit(cc).mods(2).k = reshape(dresh_qf(:,:,cc),up_SDIM*flen,1);

    interp_quad_fit(cc).mods(1).lambda_d2T = lambdal_d2T;
    interp_quad_fit(cc).mods(1).lambda_d2X = lambdal_d2X;
    interp_quad_fit(cc).mods(1).lambda_d2XT = lambdal_d2XT;
    interp_quad_fit(cc).mods(1).lambda_L1x = lambdal_L1x;
    interp_quad_fit(cc).mods(2).lambda_d2T = lambdaq_d2T;
    interp_quad_fit(cc).mods(2).lambda_d2X = lambdaq_d2X;
    interp_quad_fit(cc).mods(2).lambda_d2XT = lambdaq_d2XT;
    interp_quad_fit(cc).mods(2).lambda_L1x = lambdaq_L1x;


    Robs = all_binned_spks(tr_inds,cc);
    tr_spkbns = convert_to_spikebins(Robs);

    avg_rate = mean(all_binned_spks(tr_inds,cc));
    null_prate = ones(NT,1)*avg_rate;
    null_LL(cc) = -sum(Robs.*log(null_prate) - null_prate)/sum(Robs);

    interp_quad_fit(cc) = fitGNM_filters(interp_quad_fit(cc),X_z,tr_spkbns,'none',[],1e-4,1e-6,1);

    [old_LL(cc),~,~,~,g] = getLL_GNM(interp_quad_fit(cc),X_z,tr_spkbns,'none');
    interp_mod_theta(cc) = interp_quad_fit(cc).spk_theta;
    fprintf('Unit %d  LL imp: %.4f\n',cc,null_LL(cc)-old_LL(cc));

end

cd ~/Data/bruce/G081
save interp_quad_fits_0deg_interp interp_quad_fit interp_mod_theta old_LL null_LL
%%
% SET UP XV UNITS
NUNITS = 96;

% xv_imp = null_xvLL - orig_xvLL;
% xvLL_thresh = 0.001;
% useable_cells = find(xv_imp >= xvLL_thresh);

% tr_set = useable_cells;
% tr_set(tr_set == 18) = [];
% xv_set = 18;
% tr_set(tr_set == 13) = [];
% tr_set(tr_set == 6) = [];

% xv_set = setdiff(1:96,tr_set);
xv_frac = 0.15;
rperm = randperm(96);
xv_set = rperm(1:round(xv_frac*96));
tr_set = setdiff(1:96,xv_set);

n_tr_chs = length(tr_set);
n_xv_chs = length(xv_set);

fprintf('Using %d TR cells and %d XV cells\n',n_tr_chs,n_xv_chs);

%% DEFINE POINTS TO ALLOW RAPID CHANGES
trial_start_inds = [1; find(diff(all_trialvec(tr_inds)) ~= 0) + 1];

sac_thresh = 10;
% saccade_inds = 1 + find(interp_eye_speed(tr_inds(1:end-1)) < sac_thresh & interp_eye_speed(tr_inds(2:end)) > sac_thresh);
% saccade_inds = unique(saccade_inds);
temp = zeros(length(all_stim_times),1);
temp(interp_saccade_inds) = 1;
temp = temp(tr_inds);
saccade_inds = find(temp==1);
% use_prior(saccade_inds) = 1;

%push the effects of saccades forward in time
sac_shift = round(0.04*stim_fs);

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

shift = round(0.05/dt);
tbmat = [zeros(ntents,tblen-2*shift-1) tbmat];
tent_centers = tent_centers-shift;

rsac_Tmat = zeros(NT,ntents);
for i = 1:ntents
    rsac_Tmat(:,i) = conv(rsac_inds,tbmat(i,:),'same');
end

un_expts = unique(all_exptvec);
expt_Tmat = zeros(NT,length(un_expts));
for i = 1:length(un_expts)
    temp_set = find(all_exptvec(tr_inds)==i);
    expt_Tmat(temp_set,i) = 1;
end

%%
silent = 1;
NT = length(tr_inds);

ampphase_set = [reshape(all_interp_ampgrams(tr_inds,:,:),NT,length(wfreqs)*length(use_lfps)).*reshape(cos(all_interp_phasegrams(tr_inds,:,:)),NT,length(wfreqs)*length(use_lfps)) ...
    reshape(all_interp_ampgrams(tr_inds,:,:),NT,length(wfreqs)*length(use_lfps)).*reshape(sin(all_interp_phasegrams(tr_inds,:,:)),NT,length(wfreqs)*length(use_lfps))];
phase_elec_set = ones(length(wfreqs),1)*(1:length(use_lfps));
phase_elec_set = [phase_elec_set(:); phase_elec_set(:)];

load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
X_pos = ArrayConfig.X;
Y_pos = ArrayConfig.Y;
nearest_lfps = nan(length(use_lfps),1);
for ll = 1:96
    all_dists = sqrt((X_pos-X_pos(ll)).^2 + (Y_pos-Y_pos(ll)).^2);
    %    all_dists(ll) = inf;
    [~,best_loc] = min(all_dists(use_lfps));
    nearest_lfps(ll) = best_loc;
end

%%
NL_type = 0; %exp

reg_params.dl1_ind = 0;
reg_params.dl1_dep = 0;
reg_params.dl2_ind = 1000;
reg_params.dl2_dep = 1000;
reg_params.dl2_freq_ind = 0;
reg_params.dl2_freq_dep = 0;
reg_params.dl2_time_ind = 1000;
reg_params.dl2_time_dep = 0;
reg_params.l1_ind = 0;
reg_params.l1_dep = 0;
reg_params.l1t_ind = 0;
reg_params.l1t_dep = 0;
reg_params.is_phase = 0;

phasemod_out = nan(96,NT);

for ll = 1:length(use_lfps)
    fprintf('Electrode %d of %d\n',ll,length(use_lfps));
    cur_units = find(nearest_lfps==ll);
    use_elecs = ll;
    use_set = find(phase_elec_set==use_elecs);
    
    for cc = 1:length(cur_units)
        fprintf('Cell %d of %d\n',cc,length(cur_units));
        Robs = all_binned_spks(tr_inds,cur_units(cc));

        cur_stim_filts = get_k_mat(interp_quad_fit(cc));
        filt_outs = X_z*cur_stim_filts;
        stim_mod_out = filt_outs(:,1) + filt_outs(:,2).^2;
        
        stim_params = [length(wfreqs),1,ntents];
        X = [ampphase_set(:,use_set) rsac_Tmat expt_Tmat stim_mod_out];
        klen = size(X,2);
        K0 = zeros(klen+1,1);
        [fitp_ampphase] = fit_GLM_phase_model(X, Robs, K0,silent, stim_params, reg_params,1,NL_type);
        ampphase_cfilt(cur_units(cc),:) = fitp_ampphase.k(1:length(wfreqs));
        ampphase_sfilt(cur_units(cc),:) = fitp_ampphase.k((length(wfreqs)+1):length(use_elecs)*length(wfreqs)*2);
        ampphase_sac_tfilt(cur_units(cc),:) = fitp_ampphase.k((length(wfreqs)*2+1):(length(wfreqs)*2+ntents));
        ampphase_expt_filt(cur_units(cc),:) = fitp_ampphase.k((length(wfreqs)*2+ntents+1):(length(wfreqs)*2+ntents+length(un_expts)) );
        ampphase_stimgain(cur_units(cc),:) = fitp_ampphase.k(end-1);
        ampphase_const(cur_units(cc)) = fitp_ampphase.k(end);
        
        phase_kern = [ampphase_cfilt(cur_units(cc),:) ampphase_sfilt(cur_units(cc),:)];
        phasemod_out(cur_units(cc),:) = ampphase_set(:,use_set)*phase_kern';
        
        predrate = X*fitp_ampphase.k(1:end-1) + fitp_ampphase.k(end);
        if NL_type==0
            predrate = log(1+exp(predrate));
        else
            predrate = exp(predrate);
        end
        fullmod_LL(cur_units(cc)) = -sum(Robs.*log(predrate) - predrate)/sum(Robs);
    end
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
% lA = lA + lA_bias; %bias towards returning to 0
lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize

%% ESTIMATE LL for each shift in each stimulus frame
frame_LLs = zeros(NT,n_shifts);
Robs = all_binned_spks(tr_inds,tr_set);

up_klen = up_SDIM*flen;

% lfilt_bank = dresh_lf(:,:,tr_set);
% qfilt_bank = dresh_qf(:,:,tr_set);
lfilt_bank = nan(length(use_sus),up_klen);
qfilt_bank = nan(length(use_sus),up_klen);
for cc = 1:length(use_sus)
    cur_k = get_k_mat(interp_quad_fit(cc));
    lfilt_bank(cc,:) = cur_k(:,1);
    qfilt_bank(cc,:) = cur_k(:,2);
end
lfilt_bank = reshape(lfilt_bank',[flen,up_SDIM,96]);
qfilt_bank = reshape(qfilt_bank',[flen,up_SDIM,96]);
lfilt_bank = lfilt_bank(:,:,tr_set);
qfilt_bank = qfilt_bank(:,:,tr_set);

shifted_lfilt_bank = nan(klen,n_tr_chs);
shifted_qfilt_bank = nan(klen,n_tr_chs);

% rsac_out = rsac_Tmat*init_rsac_kern(tr_set,:)';
% expt_out = expt_Tmat*init_expt_kern(tr_set,:)';
rsac_out = rsac_Tmat*ampphase_sac_tfilt(tr_set,:)';
expt_out = expt_Tmat*ampphase_expt_filt(tr_set,:)';

shift_cnt = 1;
for xx = 1:length(shifts)
    fprintf('Shift %d of %d\n',shift_cnt,n_shifts);
    d2 = dist_shift2d(lfilt_bank,shifts(xx),2,0);
    shifted_lfilt_bank = reshape(d2,up_klen,n_tr_chs);
    d2 = dist_shift2d(qfilt_bank,shifts(xx),2,0);
    shifted_qfilt_bank = reshape(d2,up_klen,n_tr_chs);
    
    lin_out = X_z*shifted_lfilt_bank;
    quad_out = (X_z*shifted_qfilt_bank).^2;
    gfun = lin_out + quad_out;
    
    %weight stimulus term
    gfun = bsxfun(@times,gfun,ampphase_stimgain(tr_set)');
    
    %add contributions from event kernels
    gfun = gfun + rsac_out + expt_out + phasemod_out(tr_set,:)';
    
    gfun = bsxfun(@plus,gfun,ampphase_const(tr_set));
    
    too_large = gfun > 50;
    pred_rate = log(1+exp(gfun));
    pred_rate(too_large) = gfun(too_large);
    
    pred_rate(pred_rate < 1e-20) = 1e-20;
    
    LLs = Robs.*log(pred_rate) - pred_rate;
    frame_LLs(:,shift_cnt) = sum(LLs,2);
    
    shift_cnt = shift_cnt + 1;
end

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
    if mod(t,1000)==0
        fprintf('%d of %d\n',t,NT);
    end
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
    if mod(t,1000)==0
        fprintf('%d of %d\n',t,NT);
    end
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
post_mean = sum(bsxfun(@times,gamma,shifts),2);
post_std = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - post_mean.^2);

%%
resh_X = reshape(Xmat(tr_inds,:)',[flen SDIM NT]);
% resh_X = reshape(0.5*Xmat(tr_inds,:)',[flen SDIM NT]);

if up_fac > 1
    fprintf('Up-sampling X-matrix\n');
    dresh_X = [];
    for i = 1:SDIM-1
        dresh_X = cat(2,dresh_X,resh_X(:,i,:));
        %for interpolating bar positions
%         dresh_X = cat(2,dresh_X,0.5*resh_X(:,i,:) + 0.5*resh_X(:,i+1,:));
        %adding zeros for absent bar positions
        dresh_X = cat(2,dresh_X,zeros(flen,1,NT));
    end
    dresh_X = cat(2,dresh_X,resh_X(:,end,:));
    
    up_SDIM = up_fac*SDIM-1;
else
    up_SDIM = SDIM;
    dresh_X = resh_X;
end
X_z = reshape(dresh_X,flen*up_SDIM,NT)';

%% RECONSTRUCT MAP STIMULUS (OLD METHOD)
[max_post,max_loc] = max(lgamma,[],2);

min_chunk_dur = 0.1;
too_short_chunks = find(chunk_durs < min_chunk_dur);
too_short_inds = find(ismember(chunk_ids,too_short_chunks));
max_loc(too_short_inds) = zero_frame;

shift_cor = shifts(max_loc);
% shift_cor = post_mean;

% shift_cor = zeros(size(post_mean));

% resh_X = reshape(X_zp',[flen pad_SDIM NT]);
resh_X = reshape(X_z',[flen up_SDIM NT]);
resh_X_sh = zeros(size(resh_X));
for ii = 1:NT
    d2 = dist_shift2d(resh_X(:,:,ii), -shift_cor(ii),2,0);
    resh_X_sh(:,:,ii) = d2;
end
% X_sh = reshape(resh_X_sh,flen*pad_SDIM,NT)';
X_sh = reshape(resh_X_sh,flen*up_SDIM,NT)';


%% COMPARE CURRENT INFERRED EYE POSITION TO RECORDED EYE POSITION
% figure
% shadedErrorBar(all_stim_times(tr_inds),shift_cor*.04,post_std*.04,{'b'})
%
% hold on
% plot(all_stim_times(tr_inds),rsac_inds,'r.-')
% % plot(all_stim_times(tr_inds),shift_cor*0.04,'.-')
% plot(all_stim_times(tr_inds),interp_eye_vals(tr_inds,2),'k.-');
% % plot(all_stim_times(tr_inds),interp_eye_vals(tr_inds,4),'g.-');
% % plot(all_stim_times(tr_inds),0.5*interp_eye_vals(tr_inds,4)+0.5*interp_eye_vals(tr_inds,2),'c.-');
%
% figure
% imagesc(1:NT,shifts*0.04,gamma')
% hold on
% plot(inferred_eye_cors*0.04,'r-','linewidth',2)
% plot(interp_eye_vals(tr_inds,2),'w-','linewidth',2)

%% REFIT MODELS
% upklen = flen*pad_SDIM;
% stim_params.spatial_dims = 1;
% stim_params.sdim = pad_SDIM;
% stim_params.flen = flen;
upklen = flen*up_SDIM;
stim_params.spatial_dims = 1;
stim_params.sdim = up_SDIM;
stim_params.flen = flen;
% for cc = 1:length(use_sus)
for cc = 1:96
    % for cc = 13
    if ismember(cc,tr_set)
        temp_str = 'TR';
        fprintf('Fitting training cell %d of %d\n',cc,length(use_sus));
    elseif ismember(cc,xv_set)
        temp_str = 'XV';
        fprintf('Fitting xv cell %d of %d\n',cc,length(use_sus));
    else
        temp_str = 'CR';
        fprintf('Fitting crap cell %d of %d\n',cc,length(use_sus));
    end
    tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,cc));
    
    quad_fit_sh(cc) = fitGNM_filters(interp_quad_fit(cc),X_sh,tr_spkbns,'none',[],1e-4,1e-6,1);
    
    ref_LL(1,cc) = quad_fit_sh(cc).LL;
    
    fprintf('Original: %.3f  New: %.3f\n',null_LL(cc)-old_LL(cc),null_LL(cc)-ref_LL(1,cc));
    
end

% old_imp = null_LL-old_LL;
% new_imp = null_LL - ref_LL(1,:);
% mdl = LinearModel.fit(old_imp(tr_set)',new_imp(tr_set)');
% xx = linspace(0,0.15,100);
% [ypred,pred_errs] = predict(mdl,xx');
% figure;hold on
% plot(xx,ypred,'k','linewidth',2);plot(xx,pred_errs,'r--')
% plot(old_imp(xv_set),new_imp(xv_set),'o')
% line([0 0.2],[0 0.2])
% hold on
% plot(null_LL(tr_set) - old_LL(tr_set),null_LL(tr_set) - ref_LL(1,tr_set),'r*')
% xlim([0 0.2])
%% REFIT EVENT KERNELS
for cc = 1:96
    fprintf('Fitting event kernels for cell %d of %d\n',cc,96);
    
    Robs = all_binned_spks(tr_inds,cc);
    tr_spkbns = convert_to_spikebins(Robs);
    
    [nll, pnll, lpen, prate, g] = getLL_GNM(quad_fit_sh(cc),X_sh,tr_spkbns,'none');
    Xmat = [rsac_Tmat expt_Tmat g];
    klen = size(Xmat,2);
    %     K0 = [init_rsac_kern(cc,:) 1 0];
    K0 = zeros(klen+1,1);
    [fitp_sac_kerns(cc)] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], lamrange2,[], [], [], [klen+1], 0);
    rev_rsac_kern(cc,:) = fitp_sac_kerns(cc).k(1:ntents);
    rev_expt_kern(cc,:) = fitp_sac_kerns(cc).k(ntents+1:end-2);
    rev_gweight(cc) = fitp_sac_kerns(cc).k(end-1);
end


%%
it_shift_cor{1} = shift_cor;
it_quad_fit{1} = quad_fit_sh;
it_rsac_kern{1} = rev_rsac_kern;
it_expt_kern{1} = rev_expt_kern;
it_gweight{1} = rev_gweight;
it_post_mean{1} = post_mean;
it_post_std{1} = post_std;

it_tr_set{1} = tr_set;
it_xv_set{1} = xv_set;
%% NOW ITERATE

n_ITER = 7;
for nn = 2:n_ITER+1
    
%     xv_frac = 0.15;
%     rperm = randperm(96);
%     xv_set = rperm(1:round(xv_frac*96));
%     tr_set = setdiff(1:96,xv_set);
%     
%     n_tr_chs = length(tr_set);
%     n_xv_chs = length(xv_set);
%     
%     it_tr_set{nn} = tr_set;
%     it_xv_set{nn} = xv_set;
    %%
    frame_LLs = zeros(NT,n_shifts);
    
    Robs = all_binned_spks(tr_inds,tr_set);
    fprintf('Iteration %d of %d\n',nn-1,n_ITER);
    lfilt_bank = nan(length(use_sus),upklen);
    qfilt_bank = nan(length(use_sus),upklen);
    new_mod_theta = nan(length(use_sus),1);
    for cc = 1:length(use_sus)
        cur_k = get_k_mat(it_quad_fit{nn-1}(cc));
        lfilt_bank(cc,:) = cur_k(:,1);
        qfilt_bank(cc,:) = cur_k(:,2);
        new_mod_theta(cc) = it_quad_fit{nn-1}(cc).spk_theta;
    end
    lfilt_bank = reshape(lfilt_bank(tr_set,:)',[flen up_SDIM n_tr_chs]);
    qfilt_bank = reshape(qfilt_bank(tr_set,:)',[flen up_SDIM n_tr_chs]);
    shifted_lfilt_bank = nan(upklen,n_tr_chs);
    shifted_qfilt_bank = nan(upklen,n_tr_chs);
    
    rsac_out = rsac_Tmat*it_rsac_kern{nn-1}(tr_set,:)';
    expt_out = expt_Tmat*it_expt_kern{nn-1}(tr_set,:)';
    
    shift_cnt = 1;
    for xx = 1:length(shifts)
        fprintf('Shift %d of %d\n',shift_cnt,n_shifts);
        d2 = dist_shift2d(lfilt_bank,shifts(xx),2,0);
        shifted_lfilt_bank = reshape(d2,upklen,n_tr_chs);
        d2 = dist_shift2d(qfilt_bank,shifts(xx),2,0);
        shifted_qfilt_bank = reshape(d2,upklen,n_tr_chs);
        
        lin_out = X_z*shifted_lfilt_bank;
        quad_out = (X_z*shifted_qfilt_bank).^2;
        gfun = lin_out + quad_out;
        gfun = bsxfun(@plus,gfun,new_mod_theta(tr_set)');
        
        %weight stimulus term
        gfun = bsxfun(@times,gfun,it_gweight{nn-1}(tr_set));
        
        %add contributions from event kernels
        gfun = gfun + rsac_out + expt_out;
        
        too_large = gfun > 50;
        pred_rate = log(1+exp(gfun));
        pred_rate(too_large) = gfun(too_large);
        
        pred_rate(pred_rate < 1e-20) = 1e-20;
        
        LLs = Robs.*log(pred_rate) - pred_rate;
        frame_LLs(:,shift_cnt) = sum(LLs,2);
        
        shift_cnt = shift_cnt + 1;
    end
    
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
        if mod(t,1000)==0
            fprintf('%d of %d\n',t,NT);
        end
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
        if mod(t,1000)==0
            fprintf('%d of %d\n',t,NT);
        end
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
    post_mean = sum(bsxfun(@times,gamma,shifts),2);
    post_std = sqrt(sum(bsxfun(@times,gamma,shifts.^2),2) - post_mean.^2);
    it_post_mean{nn} = post_mean;
    it_post_std{nn} = post_std;
    
    %% RECONSTRUCT MAP STIMULUS (OLD METHOD)
    [max_post,max_loc] = max(lgamma,[],2);
    it_shift_cor{nn} = shifts(max_loc);
    % it_shift_cor{nn} = post_mean;
    % it_shift_std{nn} = post_std;
    
    %     resh_X = reshape(X_zp',[flen pad_SDIM NT]);
    resh_X = reshape(X_z',[flen up_SDIM NT]);
    resh_X_sh = zeros(size(resh_X));
    for ii = 1:NT
        d2 = dist_shift2d(resh_X(:,:,ii), -it_shift_cor{nn}(ii),2,0);
        resh_X_sh(:,:,ii) = d2;
    end
    X_sh = reshape(resh_X_sh,flen*up_SDIM,NT)';
    
    %%
    stim_params.spatial_dims = 1;
    stim_params.sdim = up_SDIM;
    stim_params.flen = flen;
    for cc = 1:length(use_sus)
        %     for cc = 13
        if ismember(cc,tr_set)
            temp_str = 'TR';
            fprintf('Fitting training cell %d of %d\n',cc,length(use_sus));
        elseif ismember(cc,xv_set)
            temp_str = 'XV';
            fprintf('Fitting xv cell %d of %d\n',cc,length(use_sus));
        else
            temp_str = 'CR';
            fprintf('Fitting crap cell %d of %d\n',cc,length(use_sus));
        end
        tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,(cc)));
        
        it_quad_fit{nn}(cc) = fitGNM_filters(it_quad_fit{nn-1}(cc),X_sh,tr_spkbns,'none',[],1e-4,1e-6,1);
        ref_LL(nn,cc) = it_quad_fit{nn}(cc).LL;
        
        fprintf('Original: %.3f  First: %.3f  Previous: %.3f   Current: %.3f\n',...
            null_LL(cc)-old_LL(cc),null_LL(cc)-it_quad_fit{1}(cc).LL,null_LL(cc)-it_quad_fit{nn-1}(cc).LL,null_LL(cc)-it_quad_fit{nn}(cc).LL);
    end
    
    %% REFIT EVENT KERNELS
    % lamrange2 = [300 1 ntents 0; 300 ntents+1 2*ntents 0];
    for cc = 1:96
        fprintf('Fitting event kernels for cell %d of %d\n',cc,96);
        
        Robs = all_binned_spks(tr_inds,cc);
        tr_spkbns = convert_to_spikebins(Robs);
        
        [nll, pnll, lpen, prate, g] = getLL_GNM(it_quad_fit{nn}(cc),X_sh,tr_spkbns,'none');
        Xmat = [rsac_Tmat expt_Tmat g];
        klen = size(Xmat,2);
        %     K0 = [init_rsac_kern(cc,:) init_ssac_kern(cc,:) 1 0];
        K0 = [init_rsac_kern(cc,:) zeros(1,length(un_expts)) 1 0];
        [fitp_sac_kerns(cc)] = GLMsolve_jmm(Xmat, tr_spkbns, K0, 1, [], lamrange2,[], [], [], [klen+1], 0);
        it_rsac_kern{nn}(cc,:) = fitp_sac_kerns(cc).k(1:ntents);
        it_expt_kern{nn}(cc,:) = fitp_sac_kerns(cc).k(ntents+1:end-2);
        it_gweight{nn}(cc) = fitp_sac_kerns(cc).k(end-1);
    end
    
    %%
    save full_eye_correct_0deg_new_interp *_LL it_quad_fit tr_set xv_set it_shift_cor tr_inds interp_quad_fit it_*
    
end

%%
eye_times = all_stim_times(tr_inds);
save full_eye_correct_0deg_new_interp *_LL it_quad_fit tr_set xv_set it_shift_cor tr_inds eye_times interp_quad_fit it_*
%%
pre_improve = (null_LL-old_LL)/log(2);
post_improve = (null_LL - ref_LL(1,:))/log(2);

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
for cc = 1:96
    if ismember(cc,tr_set)
        fprintf('TR Cell %d of %d\n',cc,96);
    else
        fprintf('XV Cell %d of %d\n',cc,96);
    end
    k = get_k_mat(interp_quad_fit(cc));
    subplot(2,2,1)
%     imagesc(up_bar_pos,-(0:flen)*dt,[zeros(flen,zpads) reshape(k(:,1),flen,up_SDIM) zeros(flen,zpads)]);
    imagesc(up_bar_pos,-(0:flen)*dt,reshape(k(:,1),flen,up_SDIM));
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    title('Pre Linear')
    xlabel('Relative bar position (deg)')
    ylabel('Time lag (s)')
    subplot(2,2,2)
%     imagesc(up_bar_pos,-(0:flen)*dt,[zeros(flen,zpads) reshape(k(:,2),flen,up_SDIM) zeros(flen,zpads)])
    imagesc(up_bar_pos,-(0:flen)*dt,reshape(k(:,2),flen,up_SDIM))
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    title('Pre Quad')
    xlabel('Relative bar position (deg)')
    ylabel('Time lag (s)')
    
    k = get_k_mat(it_quad_fit{4}(cc));
    subplot(2,2,3)
    imagesc(up_bar_pos,-(0:flen)*dt,reshape(k(:,1),flen,up_SDIM));
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    title('Post Linear')
    xlabel('Relative bar position (deg)')
    ylabel('Time lag (s)')
    subplot(2,2,1)
    caxis([-ca ca]);colorbar
    
    subplot(2,2,4)
    imagesc(up_bar_pos,-(0:flen)*dt,reshape(k(:,2),flen,up_SDIM))
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    title('Post Quad')
    xlabel('Relative bar position (deg)')
    ylabel('Time lag (s)')
    
    subplot(2,2,2)
    caxis([-ca ca]);colorbar
    orig_imp = (null_LL(cc) - old_LL(cc))/log(2);
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
shadedErrorBar(all_stim_times(tr_inds),it_post_mean{1}*.04,it_post_std{1}*.04,{'r'})
shadedErrorBar(all_stim_times(tr_inds),it_post_mean{2}*.04,it_post_std{2}*.04,{'b'})
% shadedErrorBar(all_stim_times(tr_inds),it_post_mean{end}*.04,it_inferred_std{4}*.04,{'b'},0)
% shadedErrorBar(all_stim_times(tr_inds),inferred_eye_cors1*.04,it_post_std{1}*.04,{'r'})
% shadedErrorBar(all_stim_times(tr_inds),inferred_eye_cors2*.04,it_post_std{end}*.04,{'b'},0)
plot(all_stim_times(tr_inds),it_shift_cor{end}*0.04,'k','linewidth',2)

% plot(all_stim_times(tr_inds),rsac_inds,'r.-')
% plot(all_stim_times(tr_inds),shift_cor*0.04,'.-')
plot(all_stim_times(tr_inds),interp_eye_vals(tr_inds,2),'k.-');
plot(all_stim_times(tr_inds),interp_eye_vals(tr_inds,4),'g.-');
% plot(all_stim_times(tr_inds),0.5*interp_eye_vals(tr_inds,4)+0.5*interp_eye_vals(tr_inds,2),'c.-');

plot(all_stim_times(interp_saccade_inds),0.1*ones(size(interp_saccade_inds)),'k*','markersize',16)
ylim([-0.75 0.75])
yl = ylim();
for i = 1:length(all_trial_start_times)
    line(all_trial_start_times([i i]),yl,'color',[0.2 0.9 0.2],'linewidth',2)
    line(all_trial_end_times([i i]),yl,'color','r','linewidth',2)
end

for i = 500:length(all_trial_start_times)
    xlim([all_trial_start_times(i)-0.1 all_trial_end_times(i)+0.1])
    pause
end
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

%%
% save nocorr_mod_135 nocor_LL test_nocor_mod
