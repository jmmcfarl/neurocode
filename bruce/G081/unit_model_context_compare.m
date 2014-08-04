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
bar_oris = [135];
un_bar_pos = all_un_bar_pos(:,4); %for 90 deg

fprintf('Analyzing %d ori expts\n',bar_oris);

%     cur_expt_set = find(is_bar_expt == 1 & expt_image_back == 1 & expt_bar_ori == bar_oris(EE));

% %expts with X deg bars and gray back (sim sacs)
% cur_expt_set = find(is_bar_expt==1 & expt_image_back == 0 & expt_bar_ori == bar_oris);

% %expts with X deg bars and image back (real sacs)
% cur_expt_set = find(is_bar_expt==1 & expt_image_back == 1 & expt_bar_ori == bar_oris & expt_sim_sacs==0);

% %expts with X deg bars and image back (real sacs)
% cur_expt_set = find(is_bar_expt==1 & expt_image_back == 1 & expt_bar_ori == bar_oris & expt_sim_sacs>0);

%expts with X deg bars and image back (real sacs)
cur_expt_set = find(is_bar_expt==1 & expt_bar_ori == bar_oris);

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
back_pts = 1 + find(diff(all_eye_ts) <= 0);
double_samples = [];
for i = 1:length(back_pts)
   next_forward = find(all_eye_ts > all_eye_ts(back_pts(i)-1),1,'first'); 
   double_samples = [double_samples back_pts(i):next_forward];
end
all_eye_ts(double_samples) = [];
all_eye_speed(double_samples) = [];
all_eye_vals(double_samples,:) = [];

interp_eye_vals = interp1(all_eye_ts,all_eye_vals,all_stim_times);
interp_eye_speed = interp1(all_eye_ts,all_eye_speed,all_stim_times);

orth_eye_pos = interp_eye_vals(:,2);
bad_pts = find(abs(orth_eye_pos) > 1);

%% IF YOU WANT TO DO CROSS-VALIDATION ON INITIAL MODEL FITS, OTHERWISE,
[c,ia,ic] = unique([all_exptvec all_trialvec],'rows');
n_trials = length(ia);
bad_trials = unique(ic(bad_pts)); %trials with putative blinks

xv_frac = 0.15;
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
up_fac = 1;
bar_dxd = bar_dx/up_fac;
SDIM = n_bar_pos;
NT = length(tr_inds);
resh_X = reshape(Xmat(tr_inds,:)',[flen SDIM NT]);

if up_fac > 1
    fprintf('Up-sampling X-matrix\n');
    dresh_X = [];
    for i = 1:SDIM-1
        dresh_X = cat(2,dresh_X,resh_X(:,i,:));
        dresh_X = cat(2,dresh_X,0.5*resh_X(:,i,:) + 0.5*resh_X(:,i+1,:));
    end
    dresh_X = cat(2,dresh_X,resh_X(:,end,:));
    
    up_SDIM = up_fac*SDIM-1;
else
    up_SDIM = SDIM;
    dresh_X = resh_X;
end
X_z = reshape(dresh_X,flen*up_SDIM,NT)';

if xv_frac > 0
    xvNT = length(xv_inds);
    resh_X = reshape(Xmat(xv_inds,:)',[flen SDIM xvNT]);
    
    if up_fac > 1
        fprintf('Up-sampling X-matrix\n');
        dresh_X = [];
        for i = 1:SDIM-1
            dresh_X = cat(2,dresh_X,resh_X(:,i,:));
            dresh_X = cat(2,dresh_X,0.5*resh_X(:,i,:) + 0.5*resh_X(:,i+1,:));
        end
        dresh_X = cat(2,dresh_X,resh_X(:,end,:));
    else
        dresh_X = resh_X;
    end
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
    defmod.lambda_L1x = 7;
    defmod.lambda_d2XT = 200;
    defmod.lambda_d2X = 200;
    defmod.lambda_d2T = 500;
    %     defmod.lambda_d2XT = 100;
    %     defmod.lambda_d2X = 300;
    %     defmod.lambda_d2T = 300;
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
        quad.mods(i).lambda_L1x = 0.35;
        %         quad.mods(i).lambda_d2XT = 30;
        %         quad.mods(i).lambda_d2X = 70;
        %         quad.mods(i).lambda_d2T = 70;
        quad.mods(i).lambda_d2XT = 15;
        quad.mods(i).lambda_d2X = 15;
        quad.mods(i).lambda_d2T = 30;
    end
    quad_fit(cc) = fitGNM_filters(quad,X_z,tr_spkbns,'none',[],1e-4,1e-6,1);
    cur_k = get_k_mat(quad_fit(cc));
    n_zeros = sum(cur_k==0);
    bad_filts = find(n_zeros > 0.95*klen);
    cur_k(:,bad_filts) = 0;
    for ii = 1:length(bad_filts)
        quad_fit(cc).mods(bad_filts(ii)).k = zeros(klen,1);
    end
    quad_fit(cc).LL = getLL_GNM(quad_fit(cc),X_z,tr_spkbns,'none');
    
    orig_LL(cc) = quad_fit(cc).LL;
    
    Robs = all_binned_spks(tr_inds,use_sus(cc));
    avg_rate = mean(Robs);
    null_prate = ones(length(tr_inds),1)*avg_rate;
    null_LL(cc) = -sum(Robs.*log(null_prate) - null_prate)/sum(Robs);
    
    if xv_frac > 0
        xv_Robs = all_binned_spks(xv_inds,use_sus(cc));
        xv_spkbns = convert_to_spikebins(xv_Robs);
        null_prate = ones(length(xv_inds),1)*avg_rate;
        null_xvLL(cc) = -sum(xv_Robs.*log(null_prate) - null_prate)/sum(xv_Robs);
        
        orig_xvLL(cc) = getLL_GNM(quad_fit(cc),xv_X_z,xv_spkbns,'none');
        fprintf('xv Imp: %.3f\n',null_xvLL(cc)-orig_xvLL(cc));
    end
    
end
% tr_mean_rate = mean(all_binned_spks(tr_inds,:));
% save perfect_fixation_mods_v2 orig_LL quad_fit null_LL tr_mean_rate defmod stim_params
% save perfect_fixation_all2 *xvLL *_LL quad_fit
save perfect_fixation_all_135 *xvLL *_LL quad_fit

%%
load ./perfect_fixation_grayback.mat
gray_xv_imp = null_xvLL-orig_xvLL;
gray_fit = quad_fit;

load ./perfect_fixation_imback.mat
im_xv_imp = null_xvLL-orig_xvLL;
im_fit = quad_fit;

load ./perfect_fixation_imback_simsac.mat
ss_xv_imp = null_xvLL-orig_xvLL;
ss_fit = quad_fit;


%%
for cc = 1:96
    fprintf('Cell %d of %d\n',cc,96);
    k = get_k_mat(gray_fit(cc));
    subplot(2,3,1)
    imagesc(reshape(k(:,1),flen,up_SDIM));
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    title(sprintf('Gray: %.3f',gray_xv_imp(cc)));
    subplot(2,3,4)
    imagesc(reshape(k(:,2),flen,up_SDIM))
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    
    k = get_k_mat(im_fit(cc));
    subplot(2,3,2)
    imagesc(reshape(k(:,1),flen,up_SDIM));
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    title(sprintf('IM: %.3f',im_xv_imp(cc)));
    subplot(2,3,5)
    imagesc(reshape(k(:,2),flen,up_SDIM))
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    
    k = get_k_mat(ss_fit(cc));
    subplot(2,3,3)
    imagesc(reshape(k(:,1),flen,up_SDIM));
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    title(sprintf('SS: %.3f',ss_xv_imp(cc)));
    subplot(2,3,6)
    imagesc(reshape(k(:,2),flen,up_SDIM))
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    pause
    clf
end

%%
xvLL_imp = null_xvLL - orig_xvLL;
for cc = 1:96
    fprintf('Cell %d of %d\n',cc,96);
    k = get_k_mat(quad_fit(cc));
    subplot(2,1,1)
    imagesc(reshape(k(:,1),flen,up_SDIM));
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    title(sprintf('%.3f',xvLL_imp(cc)));
    subplot(2,1,2)
    imagesc(reshape(k(:,2),flen,up_SDIM))
    ca = max(abs(caxis()));
    caxis([-ca ca]);colorbar
    pause
    clf
end