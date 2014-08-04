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

tr_Xmat = all_bar_mat(tr_inds,:);
xv_Xmat = all_bar_mat(xv_inds,:);

% %%
% SDIM = length(un_bar_pos);
% for cc = 1:96
%     fprintf('Fitting cell %d of %d\n',cc,96);
%     Robs_tr = all_binned_spks(tr_inds,cc);
%     tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,cc));
%     
%     %compute STA
%     sta = mean(tr_Xmat(tr_spkbns,:)) - mean(tr_Xmat);
%     
%     %don't project out STA
%     proj_mat = sta'*inv(sta*sta')*sta;
%     stim_proj = tr_Xmat - tr_Xmat*proj_mat;
%     %     stim_proj = cur_tr_stimemb;
%     stvcv = cov(stim_proj(tr_spkbns,:));
%     utvcv = cov(stim_proj);
%     [evecs,evals] = eig(stvcv-utvcv);
%     evs = diag(evals);
%     STCbvs  = evecs;
%     sta = sta';
%     npos = 4; nneg = 4;
%     stcs_compareset  = evecs(:,[1:nneg,length(evs)-npos+1:end]);
%     stcs_compareset  = stcs_compareset(:,end:-1:1);
%     rstcs = fliplr(stcs_compareset); %reversed STC kernels (suppressive first)
%     used_stcs = [sta stcs_compareset(:,1:npos) rstcs(:,1:nneg)];
%     used_stcs = bsxfun(@rdivide,used_stcs,sqrt(sum(used_stcs.^2)));
%     
%     %FIT NULL MODEL
%     avg_rate = mean(Robs_tr);
%     trpred_rate = ones(1,length(tr_inds))*avg_rate;
%     null_LL(cc) = -sum(Robs_tr.*log(trpred_rate') - trpred_rate')/sum(Robs_tr)
%     
%     subplot(3,4,1)
%     imagesc(reshape(sta,flen,SDIM));
%     cax = caxis(); cm = max(abs(cax));caxis([-cm cm]);
%     subplot(3,4,2)
%     plot(evs(1:30),'.')
%     subplot(3,4,3)
%     plot(evs(end-30:end),'.')
%     for i = 1:8
%         subplot(3,4,4+i)
%         imagesc(reshape(used_stcs(:,i+1),flen,SDIM))
%         cax = caxis(); cm = max(abs(cax));caxis([-cm cm]);
%     end
%     pause
%     clf
% end
% 
%%
SDIM = length(un_bar_pos);
klen = flen*SDIM;
stim_params.spatial_dims = 1;
stim_params.sdim = SDIM;
stim_params.flen = flen;

for cc = 1:96
    fprintf('Fitting cell %d of %d\n',cc,96);
    tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,cc));
    
    clear defmod
    defmod.lambda_L1x = 20;
    defmod.lambda_d2XT = 200;
    defmod.lambda_d2X = 0;
    defmod.lambda_d2T = 200;
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
        quad.mods(i).lambda_L1x = 4;
        quad.mods(i).lambda_d2XT = 15;
        quad.mods(i).lambda_d2X = 0;
        quad.mods(i).lambda_d2T = 30;
    end
    quad_fit(cc) = fitGNM_filters(quad,tr_Xmat,tr_spkbns,'none',[],1e-4,1e-6,1);
    cur_k = get_k_mat(quad_fit(cc));
%     n_zeros = sum(cur_k==0);
%     bad_filts = find(n_zeros > 0.95*klen);
%     cur_k(:,bad_filts) = 0;
%     for ii = 1:length(bad_filts)
%         quad_fit(cc).mods(bad_filts(ii)).k = zeros(klen,1);
%     end
%     quad_fit(cc).LL = getLL_GNM(quad_fit(cc),tr_Xmat,tr_spkbns,'none');
    
    [~,~,~,train_predrate(cc,:),g] = getLL_GNM(quad_fit(cc),tr_Xmat,tr_spkbns,'none');    
%     quad_fit2(cc) = fitGNM_spkNL(quad_fit(cc),g,tr_spkbns,0);

    [orig_LL(cc)] = getLL_GNM(quad_fit(cc),tr_Xmat,tr_spkbns,'none');
%     [orig_LL2(cc)] = getLL_GNM(quad_fit2(cc),tr_Xmat,tr_spkbns,'none');

    Robs = all_binned_spks(tr_inds,cc);
    avg_rate = mean(Robs);
    null_prate = ones(length(tr_inds),1)*avg_rate;
    null_LL(cc) = -sum(Robs.*log(null_prate) - null_prate)/sum(Robs);
    
    if xv_frac > 0
        xv_Robs = all_binned_spks(xv_inds,use_sus(cc));
        xv_spkbns = convert_to_spikebins(xv_Robs);
        null_prate = ones(length(xv_inds),1)*avg_rate;
        null_xvLL(cc) = -sum(xv_Robs.*log(null_prate) - null_prate)/sum(xv_Robs);
        
        orig_xvLL(cc) = getLL_GNM(quad_fit(cc),xv_Xmat,xv_spkbns,'none');
        fprintf('xv Imp: %.3f\n',null_xvLL(cc)-orig_xvLL(cc));
        fprintf('prev xv Imp: %.3f\n',three_imp(cc));
    end
        fprintf('tr Imp: %.3f\n',null_LL(cc)-orig_LL(cc));
end
% save init_mods_new null_xvLL orig_xvLL quad_fit train_predrate null_LL orig_LL

%%
for cc = 1:96
plotfo1d_nopsc(quad_fit(cc),3,'centered');colormap(jet);
set(gcf,'Position',[100 458 550 600])
plotfo1d_nopsc(prev_quad_fit(cc),3,'centered');colormap(jet);
set(gcf,'Position',[650 458 560 600])
fprintf('xv Imp: %.3f\n',null_xvLL(cc)-orig_xvLL(cc));
fprintf('prev xv Imp: %.3f\n',three_imp(cc));
pause
close all
end
%%
[c,ia,ic] = unique([all_exptvec all_trialvec],'rows');
n_trials = length(ia);
bad_trials = unique(ic(bad_pts)); %trials with putative blinks

xv_frac = 0;
xv_set = randperm(n_trials);
xv_set(round(n_trials*xv_frac)+1:end) = [];

tr_set = setdiff(1:n_trials,xv_set);

tr_inds = find(ismember(ic,tr_set));
xv_inds = find(ismember(ic,xv_set));
tr_inds(all_used_inds(tr_inds) == 0) = [];
xv_inds(all_used_inds(xv_inds) == 0) = [];

tr_Xmat = all_bar_mat(tr_inds,:);
xv_Xmat = all_bar_mat(xv_inds,:);

SDIM = length(un_bar_pos);
klen = flen*SDIM;
stim_params.spatial_dims = 1;
stim_params.sdim = SDIM;
stim_params.flen = flen;
load ./init_mods_new
for cc = 1:96
    fprintf('Fitting cell %d of %d\n',cc,96);
    tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,cc));
    
    quad_fit(cc) = fitGNM_filters(quad_fit(cc),tr_Xmat,tr_spkbns,'none',[],1e-4,1e-6,1);
    cur_k = get_k_mat(quad_fit(cc));
    [orig_LL(cc)] = getLL_GNM(quad_fit(cc),tr_Xmat,tr_spkbns,'none');

    Robs = all_binned_spks(tr_inds,cc);
    avg_rate = mean(Robs);
    null_prate = ones(length(tr_inds),1)*avg_rate;
    null_LL(cc) = -sum(Robs.*log(null_prate) - null_prate)/sum(Robs);
    
    if xv_frac > 0
        xv_Robs = all_binned_spks(xv_inds,use_sus(cc));
        xv_spkbns = convert_to_spikebins(xv_Robs);
        null_prate = ones(length(xv_inds),1)*avg_rate;
        null_xvLL(cc) = -sum(xv_Robs.*log(null_prate) - null_prate)/sum(xv_Robs);
        
        orig_xvLL(cc) = getLL_GNM(quad_fit(cc),xv_Xmat,xv_spkbns,'none');
        fprintf('xv Imp: %.3f\n',null_xvLL(cc)-orig_xvLL(cc));
    end
        fprintf('tr Imp: %.3f\n',null_LL(cc)-orig_LL(cc));
end
save init_mods_full_new quad_fit train_predrate null_LL orig_LL












%% double sample Xmat
bar_dx = 0.08;
up_fac = 2;
bar_dxd = bar_dx/up_fac;
SDIM = n_bar_pos;
NT = length(tr_inds);
% resh_X = reshape(Xmat(tr_inds,:)',[flen SDIM NT]);
resh_X = reshape(0.5*tr_Xmat',[flen SDIM NT]);
resh_X_xv = reshape(0.5*xv_Xmat',[flen SDIM length(xv_inds)]);

if up_fac > 1
    fprintf('Up-sampling X-matrix\n');
    dresh_X = [];
    dresh_X_xv = [];
    for i = 1:SDIM-1
        dresh_X = cat(2,dresh_X,resh_X(:,i,:));
        %for interpolating bar positions
        dresh_X = cat(2,dresh_X,0.5*resh_X(:,i,:) + 0.5*resh_X(:,i+1,:));
        %adding zeros for absent bar positions
        %         dresh_X = cat(2,dresh_X,zeros(flen,1,NT));
        
        dresh_X_xv = cat(2,dresh_X_xv,resh_X_xv(:,i,:));
        %for interpolating bar positions
        dresh_X_xv = cat(2,dresh_X_xv,0.5*resh_X_xv(:,i,:) + 0.5*resh_X_xv(:,i+1,:));
        %adding zeros for absent bar positions
        %         dresh_X = cat(2,dresh_X,zeros(flen,1,NT));
    end
    dresh_X = cat(2,dresh_X,resh_X(:,end,:));
    dresh_X_xv = cat(2,dresh_X_xv,resh_X_xv(:,end,:));
    
    up_SDIM = up_fac*SDIM-1;
else
    up_SDIM = SDIM;
    dresh_X = resh_X;
end
X_z = reshape(dresh_X,flen*up_SDIM,NT)';
X_z_xv = reshape(dresh_X_xv,flen*up_SDIM,length(xv_inds))';

%%
% load ./perfect_fixation_all.mat quad_fit
% 
% use_sus = 1:96;
% klen = flen*SDIM;
% lfilt_bank = nan(length(use_sus),klen);
% qfilt_bank = nan(length(use_sus),klen);
% mod_offset = nan(length(use_sus),1);
% for cc = 1:length(use_sus)
%     cur_k = get_k_mat(quad_fit(cc));
%     lfilt_bank(cc,:) = cur_k(:,1);
%     qfilt_bank(cc,:) = cur_k(:,2);
%     mod_offset(cc) = quad_fit(cc).spk_theta;
% end
% resh_lfilt_bank = reshape(0.5*lfilt_bank',[flen SDIM 96]);
% resh_qfilt_bank = reshape(0.5*qfilt_bank',[flen SDIM 96]);
% dresh_lf = [];
% dresh_qf = [];
% for i = 1:SDIM-1
%     dresh_lf = cat(2,dresh_lf,resh_lfilt_bank(:,i,:));
%     dresh_lf = cat(2,dresh_lf,0.5*resh_lfilt_bank(:,i,:) + 0.5*resh_lfilt_bank(:,i+1,:));
%     dresh_qf = cat(2,dresh_qf,resh_qfilt_bank(:,i,:));
%     dresh_qf = cat(2,dresh_qf,0.5*resh_qfilt_bank(:,i,:) + 0.5*resh_qfilt_bank(:,i+1,:));
% end
% dresh_lf = cat(2,dresh_lf,resh_lfilt_bank(:,end,:));
% dresh_qf = cat(2,dresh_qf,resh_qfilt_bank(:,end,:));

stim_params.spatial_dims = 1;
stim_params.sdim = up_SDIM;
stim_params.flen = flen;
spk_nl = 'logexp';
upklen = up_SDIM*flen;
% for zero Xmatrix
lambdal_L1x = 30;
lambdal_d2X = 4000;
lambdal_d2T = 4000;
lambdal_d2XT = 4000;
lambdaq_L1x = 5;
lambdaq_d2X = 25;
lambdaq_d2T = 25;
lambdaq_d2XT = 25;

% REALLY FIT INITIAL MODELS
for cc = 1:96
%     interp_quad_fit(cc) = quad_fit(cc);
init_signs = [1 1];
kern_types{1} = 'lin';
kern_types{2} = 'quad';
init_kerns = 0.01*randn(upklen,2);
interp_quad_fit(cc) = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params,spk_nl);
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

    [interp_LL(cc),~,~,~,g] = getLL_GNM(interp_quad_fit(cc),X_z,tr_spkbns,'none');
    interp_mod_theta(cc) = interp_quad_fit(cc).spk_theta;
    fprintf('Unit %d  LL imp: %.4f\n',cc,null_LL(cc)-interp_LL(cc));

    xv_Robs = all_binned_spks(xv_inds,use_sus(cc));
    xv_spkbns = convert_to_spikebins(xv_Robs);
    interp_xvLL(cc) = getLL_GNM(interp_quad_fit(cc),X_z_xv,xv_spkbns,'none');
     nullxv_prate = ones(length(xv_inds),1)*avg_rate;
    null_xvLL(cc) = -sum(xv_Robs.*log(nullxv_prate) - nullxv_prate)/sum(xv_Robs);
   fprintf('Unit %d  xvLL imp: %.4f\n',cc,null_xvLL(cc)-interp_xvLL(cc));

end
