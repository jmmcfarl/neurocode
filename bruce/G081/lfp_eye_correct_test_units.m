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
beg_buffer = round(stim_fs*0.15);
bar_oris = [0];
un_bar_pos = all_un_bar_pos(:,1);

fprintf('Analyzing %d ori expts\n',bar_oris);

%     cur_expt_set = find(is_bar_expt == 1 & expt_image_back == 1 & expt_bar_ori == bar_oris(EE));

%expts with X deg bars and gray back (sim sacs)
cur_expt_set = find(is_bar_expt==1 & expt_image_back == 0 & expt_bar_ori == bar_oris);
cur_expt_set(cur_expt_set<=8) = []; %this expt had continuous bar luminance
cur_expt_set(cur_expt_set == 13) = []; %problem with LFP data?

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
for ee = 1:length(cur_expt_set);
    fprintf('Expt %d of %d\n',ee,length(cur_expt_set));
    cur_expt = cur_expt_set(ee);
    fname = sprintf('Expt%dClusterTimes.mat',cur_expt);
    load(fname);
    
    n_trials = length(Expts{cur_expt}.Trials);
    trial_start_times = [Expts{cur_expt}.Trials(:).TrialStart]/1e4;
    trial_end_times = [Expts{cur_expt}.Trials(:).TrueEnd]/1e4;
    trial_durs = [Expts{cur_expt}.Trials(:).dur]/1e4;
    
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
        
        cur_bar_mat = zeros(length(cur_stim_times),n_bar_pos);
        for bb = 1:n_bar_pos
            cur_set = find(cur_Op==un_bar_pos(bb));
            pset = cur_phase(cur_set) == 0;
            nset = cur_phase(cur_set) == pi;
            cur_bar_mat(cur_set(pset),bb) = 1;
            cur_bar_mat(cur_set(nset),bb) = -1;
            %                 cur_bar_mat(cur_set,bb) = 1;
        end
        
        bar_Xmat = makeStimRows(cur_bar_mat,flen);
        cur_used_inds = ones(length(cur_stim_times),1);
        cur_used_inds(1:flen) = 0;
        cur_used_inds(1:beg_buffer) = 0;
        
        all_stim_times = [all_stim_times; cur_stim_times];
        all_rel_stimes = [all_rel_stimes; cur_stim_times- trial_start_times(tt)];
        all_rel_etimes = [all_rel_etimes; trial_end_times(tt) - cur_stim_times];
        all_Op = [all_Op; cur_Op];
        all_phase = [all_phase; cur_phase'];
        all_used_inds = [all_used_inds; cur_used_inds];
        all_bar_mat = [all_bar_mat; bar_Xmat];
        all_binned_spks = [all_binned_spks; cur_binned_spks];
        all_exptvec = [all_exptvec; ones(size(cur_stim_times))*ee];
        all_trialvec = [all_trialvec; ones(size(cur_stim_times))*tt];
    end
    
end

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
    clear sm_avg_eyepos eye_vel
    sm_avg_eyepos(:,1) = smooth(lEyeXY(:,1),3);
    sm_avg_eyepos(:,2) = smooth(lEyeXY(:,2),3);
    eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/eye_dt;
    eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/eye_dt;
    
    eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
    
    all_eye_vals = [all_eye_vals; lEyeXY];
    all_eye_speed = [all_eye_speed; eye_speed];
    all_eye_ts = [all_eye_ts; eye_ts_interp'];
end
interp_eye_vals = interp1(all_eye_ts,all_eye_vals,all_stim_times);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_stim_times);

orth_eye_pos = interp_eye_vals(:,2);
bad_pts = find(abs(orth_eye_pos) > 1);

%% IF YOU WANT TO DO CROSS-VALIDATION ON INITIAL MODEL FITS, OTHERWISE, 
[c,ia,ic] = unique([all_exptvec all_trialvec],'rows');
n_trials = length(ia);
bad_trials = unique(ic(bad_pts)); %trials with putative blinks

xv_frac = 0;
n_xv_trials = round(n_trials*xv_frac);
xv_set = randperm(n_trials);
xv_set(n_xv_trials+1:end) = [];
tr_set = find(~ismember(1:n_trials,xv_set));

% xv_set(ismember(xv_set,bad_trials)) = [];
% tr_set(ismember(tr_set,bad_trials)) = [];

xv_inds = find(ismember(ic,xv_set));
tr_inds = find(ismember(ic,tr_set));
% tr_inds = find(~ismember(ic,xv_inds))';

tr_inds(all_used_inds(tr_inds) == 0) = [];
xv_inds(all_used_inds(xv_inds) == 0) = [];

Xmat = all_bar_mat;

%% double sample Xmat
bar_dx = 0.08;
up_fac = 2;
bar_dxd = bar_dx/up_fac;
SDIM = n_bar_pos;
NT = length(tr_inds);
resh_X = reshape(Xmat(tr_inds,:)',[flen SDIM NT]);

fprintf('Up-sampling X-matrix\n');
dresh_X = [];
for i = 1:SDIM-1
    dresh_X = cat(2,dresh_X,resh_X(:,i,:));
    dresh_X = cat(2,dresh_X,0.5*resh_X(:,i,:) + 0.5*resh_X(:,i+1,:));
end
dresh_X = cat(2,dresh_X,resh_X(:,end,:));

up_SDIM = up_fac*SDIM-1;

X_z = reshape(dresh_X,flen*up_SDIM,NT)';

if xv_frac > 0
    xvNT = length(xv_inds);
    resh_X = reshape(Xmat(xv_inds,:)',[flen SDIM xvNT]);
    
    dresh_X = [];
    for i = 1:SDIM-1
        dresh_X = cat(2,dresh_X,resh_X(:,i,:));
        dresh_X = cat(2,dresh_X,0.5*resh_X(:,i,:) + 0.5*resh_X(:,i+1,:));
    end
    dresh_X = cat(2,dresh_X,resh_X(:,end,:));
    
    
    xv_X_z = reshape(dresh_X,flen*up_SDIM,xvNT)';        
end
%% FIT models assuming perfect fixation
load ./CellList.mat
good_sus = find(all(CellList(:,:,1) > 0));

klen = flen*up_SDIM;
stim_params.spatial_dims = 1;
stim_params.sdim = up_SDIM;
stim_params.flen = flen;

for cc = 1:length(use_sus)
% for cc = good_sus
    fprintf('Fitting cell %d of %d\n',cc,length(use_sus));
    tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,use_sus(cc)));
    
    clear defmod
    defmod.lambda_L1x = 10;
    defmod.lambda_d2XT =300;
    defmod.lambda_d2X = 600;
    defmod.lambda_d2T = 600;
    npq = 1;
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
        quad.mods(i).lambda_L1x = 1;
        quad.mods(i).lambda_d2XT = 30;
        quad.mods(i).lambda_d2X = 60;
        quad.mods(i).lambda_d2T = 60;
    end
    quad_fit{1}(cc) = fitGNM_filters(quad,X_z,tr_spkbns,'none',[],1e-4,1e-6,1);
    orig_LL(cc) = quad_fit{1}(cc).LL;
    
    Robs = all_binned_spks(tr_inds,use_sus(cc));
    avg_rate = mean(Robs);
    null_prate = ones(length(tr_inds),1)*avg_rate;  
    null_LL(cc) = -sum(Robs.*log(null_prate) - null_prate)/sum(Robs);
    
    if xv_frac > 0
        xv_Robs = all_binned_spks(xv_inds,use_sus(cc));
        xv_spkbns = convert_to_spikebins(xv_Robs);
        null_prate = ones(length(xv_inds),1)*avg_rate;
        null_xvLL(cc) = -sum(xv_Robs.*log(null_prate) - null_prate)/sum(xv_Robs);
        
        orig_xvLL(cc) = getLL_GNM(quad_fit{1}(cc),xv_X_z,xv_spkbns,'none');
        fprintf('xv Imp: %.3f\n',null_xvLL(cc)-orig_xvLL(cc));
    end
end
tr_mean_rate = mean(all_binned_spks(tr_inds,:));
save perfect_fixation_mods3_up2 orig_LL quad_fit null_LL tr_mean_rate defmod stim_params

% load ./perfect_fixation_mods2
%% VIEW MODEL FITS
for cc = 1:96
    fprintf('Cell %d of %d\n',cc,96);
    k = get_k_mat(quad_fit{1}(cc));
    subplot(2,1,1)
    imagesc(reshape(k(:,1),flen,up_SDIM));
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    subplot(2,1,2)
    imagesc(reshape(k(:,2),flen,up_SDIM))
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    fprintf('Imp: %.3f\n',orig_LL(cc)-null_LL(cc));
    pause
    clf
end

%%
% SET UP XV UNITS
NUNITS = 96;
xv_frac = 0;
tr_set = randperm(NUNITS);
tr_set = tr_set(1:round(length(tr_set)*(1-xv_frac)));
xv_set = setdiff(1:NUNITS,tr_set);

ll_imp = null_LL - orig_LL;
% not_useable = find(ll_imp < 0.03);
% tr_set(ismember(tr_set,not_useable)) = [];
% xv_set(ismember(xv_set,not_useable)) = [];
tr_set = find(ll_imp > 0.02);
xv_set = find(ll_imp <= 0.02);

tr_set(tr_set==15) = [];
xv_set = 15;

n_tr_chs = length(tr_set);
n_xv_chs = length(xv_set);
%% PAD XMAT with ZEROS
% zpads = 5;
% resh_X = reshape(X_z',[flen up_SDIM NT]);
% resh_X = cat(2,zeros(flen,zpads,NT),resh_X);
% resh_X = cat(2,resh_X,zeros(flen,zpads,NT));
% pad_SDIM = up_SDIM + 2*zpads;
% X_zp = reshape(resh_X,flen*pad_SDIM,NT)';

%% DEFINE POINTS TO ALLOW RAPID CHANGES
trial_start_inds = find(diff(all_trialvec(tr_inds)) ~= 0) + 1;

sac_thresh = 10;
saccade_inds = 1 + find(interp_eye_speed(tr_inds(1:end-1)) < sac_thresh & interp_eye_speed(tr_inds(2:end)) > sac_thresh);
saccade_inds = unique(saccade_inds);
% use_prior(saccade_inds) = 1;

%push the effects of saccades forward in time
sac_shift = round(0.04*stim_fs);
% saccade_inds = saccade_inds + sac_shift;

use_prior = logical(zeros(size(tr_inds)));
use_prior(trial_start_inds) = 1;
for i = 1:length(saccade_inds)
    next_trial = trial_start_inds(find(trial_start_inds > saccade_inds(i),1,'first'));
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
% 
%%
NT = length(tr_inds);
sp_dx = 0.08/up_fac;
max_shift = 8*up_fac;
dshift = 1;
shifts = -max_shift:dshift:max_shift;
n_shifts = length(shifts);
zero_frame = find(shifts==0);

%overall prior on shifts
eps_prior_sigma = 0.125; %0.2
leps_prior = -(shifts*sp_dx).^2/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
lA_tflip = repmat(leps_prior,n_shifts,1);

eps_prior_sigma2 = 0.25; %0.2
leps_prior2 = -(shifts*sp_dx).^2/(2*eps_prior_sigma2^2);
leps_prior2 = bsxfun(@minus,leps_prior2,logsumexp(leps_prior2)); %normalize
lA_bias = repmat(leps_prior2,n_shifts,1);

cdist = squareform(pdist(shifts'*sp_dx));
deps_sigma = 0.015; %.025
lA = -cdist.^2/(2*deps_sigma^2);
% lA = lA + lA_bias; %bias towards returning to 0
lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize

%% ESTIMATE LL for each shift in each stimulus frame
frame_LLs = zeros(NT,n_shifts);
Robs = all_binned_spks(tr_inds,tr_set);

lfilt_bank = nan(length(use_sus),klen);
qfilt_bank = nan(length(use_sus),klen);
mod_offset = nan(length(use_sus),1);
for cc = 1:length(use_sus)
    cur_k = get_k_mat(quad_fit{1}(cc));
    lfilt_bank(cc,:) = cur_k(:,1);
    qfilt_bank(cc,:) = cur_k(:,2);
    mod_offset(cc) = quad_fit{1}(cc).spk_theta;
end
lfilt_bank = reshape(lfilt_bank(tr_set,:)',[flen up_SDIM n_tr_chs]);
qfilt_bank = reshape(qfilt_bank(tr_set,:)',[flen up_SDIM n_tr_chs]);
shifted_lfilt_bank = nan(klen,n_tr_chs);
shifted_qfilt_bank = nan(klen,n_tr_chs);

shift_cnt = 1;
for xx = 1:length(shifts)
    fprintf('Shift %d of %d\n',shift_cnt,n_shifts);
    d2 = dist_shift2d(lfilt_bank,shifts(xx),2,0);
    shifted_lfilt_bank = reshape(d2,klen,n_tr_chs);
    d2 = dist_shift2d(qfilt_bank,shifts(xx),2,0);
    shifted_qfilt_bank = reshape(d2,klen,n_tr_chs);
    
    lin_out = X_z*shifted_lfilt_bank;
    quad_out = (X_z*shifted_qfilt_bank).^2;
    gfun = lin_out + quad_out;
    gfun = bsxfun(@plus,gfun,mod_offset(tr_set)');
    
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

%% RECONSTRUCT MAP STIMULUS
[max_post,max_loc] = max(lgamma,[],2);

min_chunk_dur = 0.125;
too_short_chunks = find(chunk_durs < min_chunk_dur);
too_short_inds = find(ismember(chunk_ids,too_short_chunks));
max_loc(too_short_inds) = zero_frame;

shift_cor = shifts(max_loc);

% resh_X = reshape(X_zp',[flen pad_SDIM NT]);
resh_X = reshape(X_z',[flen up_SDIM NT]);
resh_X_sh = zeros(size(resh_X));
for ii = 1:NT
    d2 = dist_shift2d(resh_X(:,:,ii), -shift_cor(ii),2,0);
    resh_X_sh(:,:,ii) = d2;
end
% X_sh = reshape(resh_X_sh,flen*pad_SDIM,NT)';
X_sh = reshape(resh_X_sh,flen*up_SDIM,NT)';

%%
figure
hold on
plot(all_stim_times(tr_inds),use_prior,'r.-')
plot(all_stim_times(tr_inds),shift_cor*0.04,'.-')
plot(all_stim_times(tr_inds),interp_eye_vals(tr_inds,2),'k.-')
%%
for cc = 1:96
    cc
    Robs = all_binned_spks(tr_inds,cc);
    tr_spkbns = convert_to_spikebins(Robs);
    avg_rate = mean(all_binned_spks(tr_inds,cc));
    null_prate = ones(NT,1)*avg_rate;
    null_LL(cc) = -sum(Robs.*log(null_prate) - null_prate)/sum(Robs);
    
    [old_LL(cc),~,~,prate] = getLL_GNM(quad_fit{1}(cc),X_z,tr_spkbns,'none');
    new_LL(cc) = getLL_GNM(quad_fit{1}(cc),X_sh,tr_spkbns,'none');
   
end
figure
plot(null_LL - old_LL,null_LL - new_LL,'o')
line([0 0.15],[0 0.15])
hold on
plot(null_LL(tr_set)-old_LL(tr_set),null_LL(tr_set)-new_LL(tr_set),'r*')



%%
to_print = 0;

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
    tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,use_sus(cc)));
    
    clear defmod
    defmod.lambda_L1x = 10;
    defmod.lambda_d2XT = 300;
    defmod.lambda_d2X = 600;
    defmod.lambda_d2T = 600;
    npq = 1;
    nnq = 0;
    clear kern_types
    kern_types{1} = 'lin';
    for i = 2:(1+npq+nnq)
        kern_types{i} = 'quad';
    end
    init_kerns = get_k_mat(quad_fit{1}(cc));
    %     init_kerns = 0.01*randn(upklen,1+npq+nnq);
%     init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
    quad = createGNM(init_kerns,[1 ones(1,npq) -1*ones(1,nnq)],kern_types,defmod,stim_params);
    for i = 2:(1+npq+nnq)
        quad.mods(i).lambda_L1x = 1;
        quad.mods(i).lambda_d2XT = 30;
        quad.mods(i).lambda_d2X = 60;
        quad.mods(i).lambda_d2T = 60;
    end
    quad_fit_sh{1}(use_sus(cc)) = fitGNM_filters(quad,X_sh,tr_spkbns,'none',[],1e-4,1e-6,1);
    
    ref_LL(1,cc) = quad_fit_sh{1}(cc).LL;
    
    fprintf('Original: %.3f  Current: %.3f\n',null_LL(cc)-orig_LL(cc),null_LL(cc)-ref_LL(1,cc));
        
    if to_print==1
    fname = sprintf('~/James_scripts/bruce/G081/cur_mod_fits/Cell%d',cc);
    % VIEW MODEL FITS
    figure('visible','off');
    k = get_k_mat(quad_fit{1}(cc));
    subplot(2,2,1)
    imagesc(reshape(k(:,1),flen,up_SDIM));
    title(sprintf('%s Imp: %.3f',temp_str,null_LL(cc)-orig_LL(cc)));
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    subplot(2,2,2)
    imagesc(reshape(k(:,2),flen,up_SDIM))
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    k = get_k_mat(quad_fit_sh{1}(cc));
    subplot(2,2,3)
    imagesc(reshape(k(:,1),flen,up_SDIM));
    title(sprintf('%s Imp: %.3f',temp_str,null_LL(cc)-ref_LL(1,cc)));
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
%     xlim([zpads pad_SDIM-zpads])
    subplot(2,2,4)
    imagesc(reshape(k(:,2),flen,up_SDIM))
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
%      xlim([zpads pad_SDIM-zpads])
   print(fname,'-dpng');
    close 
    end
end

%% NOW ITERATE
%overall prior on shifts
% eps_prior_sigma = 0.2; %0.2
% leps_prior = -(shifts*sp_dx).^2/(2*eps_prior_sigma^2);
% leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
% lA_tflip = repmat(leps_prior,n_shifts,1);
% 
% cdist = squareform(pdist(shifts'*sp_dx));
% deps_sigma = 0.02;
% lA = -cdist.^2/(2*deps_sigma^2);
% lA = lA + lA_tflip; %bias towards returning to 0
% lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize

% cur_lA = lA;
it_shift_cor{1} = shifts(max_loc);

n_ITER = 1;
for nn = 2:n_ITER+1
    
    fprintf('Iteration %d of %d\n',nn-1,n_ITER);
    lfilt_bank = nan(length(use_sus),upklen);
    qfilt_bank = nan(length(use_sus),upklen);
    mod_offset = nan(length(use_sus),1);
%     for cc = 1:length(use_sus)
    for cc = tr_set
        cur_k = get_k_mat(quad_fit_sh{nn-1}(cc));
        lfilt_bank(cc,:) = cur_k(:,1);
        qfilt_bank(cc,:) = cur_k(:,2);
        mod_offset(cc) = quad_fit_sh{nn-1}(cc).spk_theta;
    end
%     lfilt_bank = reshape(lfilt_bank(tr_set,:)',[flen pad_SDIM n_tr_chs]);
%     qfilt_bank = reshape(qfilt_bank(tr_set,:)',[flen pad_SDIM n_tr_chs]);
    lfilt_bank = reshape(lfilt_bank(tr_set,:)',[flen up_SDIM n_tr_chs]);
    qfilt_bank = reshape(qfilt_bank(tr_set,:)',[flen up_SDIM n_tr_chs]);
    shifted_lfilt_bank = nan(upklen,n_tr_chs);
    shifted_qfilt_bank = nan(upklen,n_tr_chs);
    
    shift_cnt = 1;
    for xx = 1:length(shifts)
        fprintf('Shift %d of %d\n',shift_cnt,n_shifts);
        d2 = dist_shift2d(lfilt_bank,shifts(xx),2,0);
        shifted_lfilt_bank = reshape(d2,upklen,n_tr_chs);
        d2 = dist_shift2d(qfilt_bank,shifts(xx),2,0);
        shifted_qfilt_bank = reshape(d2,upklen,n_tr_chs);
        
%         lin_out = X_zp*shifted_lfilt_bank;
%         quad_out = (X_zp*shifted_qfilt_bank).^2;
        lin_out = X_z*shifted_lfilt_bank;
        quad_out = (X_z*shifted_qfilt_bank).^2;
        gfun = lin_out + quad_out + mod_offset(cc);
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
    
    %% RECONSTRUCT MAP STIMULUS
    [max_post,max_loc] = max(lgamma,[],2);
    
    min_chunk_dur = 0.125;
    too_short_chunks = find(chunk_durs < min_chunk_dur);
    too_short_inds = find(ismember(chunk_ids,too_short_chunks));
    max_loc(too_short_inds) = zero_frame;
    
    it_shift_cor{nn} = shifts(max_loc);
    
%     resh_X = reshape(X_zp',[flen pad_SDIM NT]);
    resh_X = reshape(X_z',[flen up_SDIM NT]);
    resh_X_sh = zeros(size(resh_X));
    for ii = 1:NT
        d2 = dist_shift2d(resh_X(:,:,ii), -it_shift_cor{nn}(ii),2,0);
        resh_X_sh(:,:,ii) = d2;
    end
%     X_sh = reshape(resh_X_sh,flen*pad_SDIM,NT)';
    X_sh = reshape(resh_X_sh,flen*up_SDIM,NT)';
    
    %%
    stim_params.spatial_dims = 1;
    stim_params.sdim = pad_SDIM;
    stim_params.flen = flen;
%     for cc = 1:length(use_sus)
    for cc = tr_set
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
        
        quad_fit_sh{nn}((cc)) = fitGNM_filters(quad_fit_sh{nn-1}(cc),X_sh,tr_spkbns,'none',[],1e-4,1e-6,1);
        ref_LL(nn,cc) = quad_fit_sh{nn}(cc).LL;

        fprintf('Original: %.3f  Previous: %.3f   Current: %.3f\n',quad_fit{1}(cc).LL,quad_fit_sh{nn-1}(cc).LL,quad_fit_sh{nn}(cc).LL);
        
    end
    
end

%%
close all
for cc = xv_set
    fprintf('Cell %d of %d\n',cc,96);
    k = get_k_mat(quad_fit{1}(cc));
    subplot(2,3,1)
    imagesc(reshape(k(:,1),flen,up_SDIM));
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    subplot(2,3,2)
    imagesc(reshape(k(:,2),flen,up_SDIM))
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    subplot(2,3,3)
    imagesc(reshape(k(:,3),flen,up_SDIM))
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar

    k = get_k_mat(quad_fit_sh{end}(cc));
    subplot(2,3,4)
    imagesc(reshape(k(:,1),flen,pad_SDIM));
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    subplot(2,3,5)
    imagesc(reshape(k(:,2),flen,pad_SDIM))
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    
    subplot(2,3,6)
    imagesc(reshape(k(:,3),flen,pad_SDIM))
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    
    quad_fit{1}(cc).LL
    quad_fit_sh{end}(cc).LL
    
    pause
    clf
end

