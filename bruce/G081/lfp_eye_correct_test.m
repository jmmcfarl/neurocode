clear all
% close all
cd ~/Data/bruce/G081/
load jbeG081Expts.mat
%%
stim_fs = 100; %in Hz
Fs = 3e4;
dsf = 60;Fsd = Fs/dsf;
niqf = Fsd/2;
[filt_b,filt_a] = butter(2,[2 40]/niqf);
use_lfps = [1:96];
use_sus = 1:96;

%%
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

%%
flen = 20;
beg_buffer = round(stim_fs*0.15);
bar_oris = [0];
un_bar_pos = all_un_bar_pos(:,1);

fprintf('Analyzing %d ori expts\n',bar_oris);

%     cur_expt_set = find(is_bar_expt == 1 & expt_image_back == 1 & expt_bar_ori == bar_oris(EE));

%expts with X deg bars and gray back (sim sacs)
cur_expt_set = find(is_bar_expt==1 & expt_image_back == 0 & expt_bar_ori == bar_oris);
cur_expt_set(cur_expt_set==2) = []; %this expt had continuous bar luminance
cur_expt_set(cur_expt_set == 13) = []; %problem with LFP data?

%%
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

%%
all_Vmat_interp = [];
for ee = 1:length(cur_expt_set);
    cur_expt = cur_expt_set(ee);
    
    Vmat = [];
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d of %d\n',ll,length(use_lfps));
        filename = sprintf('Expt%d.p%dFullV.mat',cur_expt,use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        %     V = V + FullV.sumscale*sumv;
        V = V*FullV.intscale(1)/FullV.intscale(2);
        nparts = length(FullV.blklen);
        dV = [];
        %splice together multiple blocks
        cur_pt = 1;
        for pp = 1:nparts
            cur_range = cur_pt:(cur_pt + FullV.blklen(pp)-1);
            cur_range(cur_range > length(V)) = [];
            curV = decimate(V(cur_range),dsf);
            curV = filtfilt(filt_b,filt_a,curV);
            dV = [dV curV];
            cur_pt = cur_pt + FullV.blklen(pp);
        end
        Vmat(:,ll) = dV;
    end
    
    t_ax = [];
    for pp = 1:nparts
        cur_t_ax = linspace(FullV.blkstart(pp),FullV.blkstart(pp)+FullV.blklen(pp)/Fs,FullV.blklen(pp));
        t_ax = [t_ax downsample(cur_t_ax,dsf)];
    end
    t_ax(size(Vmat,1)+1:end) = [];
    
    fprintf('LFP len: %d\n',range(t_ax));
    
    use_set = find(all_exptvec==ee);
    Vmat_interp = interp1(t_ax,Vmat,all_stim_times(use_set));
    all_Vmat_interp = [all_Vmat_interp; Vmat_interp];
end
%%
[c,ia,ic] = unique([all_exptvec all_trialvec],'rows');
n_trials = length(ia);

xv_frac = 0;
n_xv_trials = round(n_trials*xv_frac);
xv_set = randperm(n_trials);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(ic,xv_set));
tr_inds = find(~ismember(ic,xv_inds))';

tr_inds(all_used_inds(tr_inds) == 0) = [];
xv_inds(all_used_inds(xv_inds) == 0) = [];

Xmat = abs(all_bar_mat);
oXmat = (all_bar_mat);
Yobs = zscore(all_Vmat_interp);

%% double sample Xmat and zero pad
bar_dx = 0.08;
up_fac = 2;
bar_dxd = bar_dx/up_fac;
zpads = 5;
SDIM = n_bar_pos;
NT = length(tr_inds);
resh_X = reshape(Xmat(tr_inds,:)',[flen SDIM NT]);
resh_oX = reshape(oXmat(tr_inds,:)',[flen SDIM NT]);

fprintf('Up-sampling X-matrix\n');
dresh_X = [];
dresh_oX = [];
for i = 1:SDIM-1
    dresh_X = cat(2,dresh_X,resh_X(:,i,:));
    dresh_X = cat(2,dresh_X,0.5*resh_X(:,i,:) + 0.5*resh_X(:,i+1,:));
    
    dresh_oX = cat(2,dresh_oX,resh_oX(:,i,:));
    dresh_oX = cat(2,dresh_oX,0.5*resh_oX(:,i,:) + 0.5*resh_oX(:,i+1,:));
end
dresh_X = cat(2,dresh_X,resh_X(:,end,:));
dresh_X = cat(2,zeros(flen,zpads,NT),dresh_X);
dresh_X = cat(2,dresh_X,zeros(flen,zpads,NT));

dresh_oX = cat(2,dresh_oX,resh_oX(:,end,:));
dresh_oX = cat(2,zeros(flen,zpads,NT),dresh_oX);
dresh_oX = cat(2,dresh_oX,zeros(flen,zpads,NT));

up_SDIM = up_fac*SDIM-1+2*zpads;

X_z = reshape(dresh_X,flen*up_SDIM,NT)';
oX_z = reshape(dresh_oX,flen*up_SDIM,NT)';

%%
stim_params.spatial_dims = 1;
stim_params.sdim = up_SDIM;
stim_params.flen = flen;
klen = up_SDIM*flen;
clear defmod

for cc = 1:length(use_lfps)
    fprintf('Cell %d of %d\n',cc,length(use_lfps));
    
    defmod.lambda_L1x = 0;
    defmod.lambda_d2T = 100;
    defmod.lambda_d2X = 100;
    defmod.lambda_d2XT = 200;
    kern_types{1} = 'lin';
    init_kerns = 0.01*randn(klen,1);
    init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
    glm = createGNM(init_kerns,1,kern_types,defmod,stim_params,'gauss');
    glm_fit(cc) = fitGNM_filters(glm,X_z,Yobs(tr_inds,cc),'none',[],1e-4,1e-6,0);
    
    [glm_sse, ~, ~, prate] = getLL_GNM(glm_fit(cc),X_z,Yobs(tr_inds,cc),'none');
    yobs_var = var(Yobs(tr_inds,cc));
    glm_r2(cc) = 1-glm_sse/yobs_var;
    
    resid(cc,:) = Yobs(tr_inds,cc) - prate;
    
    if xv_frac > 0
        [glm_sse, ~, ~, prate] = getLL_GNM(glm_fit(cc),X_z,Yobs(xv_inds,cc),'none');
        yobs_var = var(Yobs(xv_inds,cc));
        xv_glm_r2(cc) = 1-glm_sse/yobs_var;
    end
end

%%
% for cc = 1:96
%     k = get_k_mat(glm_fit(cc));
%     imagesc(reshape(k,flen,up_SDIM));
%     ca = max(abs(caxis()));
%     caxis([-ca ca]);
%     pause
%     clf
% end
%% PARSE EYE-TRACKING DATA (IGNORING RIGHT EYE BECAUSE OF NOISE)
load ./jbeG081.em.mat
all_eye_vals = [];
all_eye_speed = [];
all_eye_ts = [];
for ee = 1:length(cur_expt_set);
    cur_set = find(all_exptvec==ee);
    [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data(Expt,all_stim_times(cur_set([1 end])));
    
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

sac_thresh = 10;
saccade_inds = 1 + find(interp_eye_speed(tr_inds(1:end-1)) < sac_thresh & interp_eye_speed(tr_inds(2:end)) > sac_thresh);
saccade_inds = unique(saccade_inds);

%%
% SET UP XV CHANNELS
NLFPS = 96;
xv_frac = 0;
tr_set = randperm(NLFPS);
tr_set = tr_set(1:round(length(tr_set)*(1-xv_frac)));
xv_set = setdiff(1:NLFPS,tr_set);
n_tr_chs = length(tr_set);

trial_start_inds = find(diff(all_trialvec(tr_inds)) ~= 0) + 1;

use_prior = logical(zeros(size(tr_inds)));
use_prior(trial_start_inds) = 1;
use_prior(saccade_inds) = 1;
%%
NT = length(tr_inds);
sp_dx = 0.08/up_fac;
max_shift = 8*up_fac;
dshift = 1;
shifts = -max_shift:dshift:max_shift;
n_shifts = length(shifts);
zero_frame = find(shifts==0);

resid_cov = cov(resid');
inv_cov = pinv(resid_cov(tr_set,tr_set));

%overall prior on shifts
eps_prior_sigma = 0.15; %0.2
leps_prior = -(shifts*sp_dx).^2/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
lA_tflip = repmat(leps_prior,n_shifts,1);

cdist = squareform(pdist(shifts'*sp_dx));
deps_sigma = 0.02;
lA = -cdist.^2/(2*deps_sigma^2);
lA = lA + lA_tflip; %bias towards returning to 0
lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize

%% ESTIMATE LL for each shift in each stimulus frame
frame_LLs = zeros(NT,n_shifts);
YY = Yobs(tr_inds,:);
klen = up_SDIM*flen;
filt_bank = nan(length(use_lfps),klen);
mod_offset = nan(length(use_lfps),1);
for cc = 1:length(use_lfps)
    filt_bank(cc,:) = get_k_mat(glm_fit(cc));
    mod_offset(cc) = glm_fit(cc).spk_theta;
end
filt_bank = reshape(filt_bank(tr_set,:)',[flen up_SDIM n_tr_chs]);
shifted_filt_bank = nan(klen,n_tr_chs);

YY = bsxfun(@minus,YY,mod_offset');

shift_cnt = 1;
for xx = 1:length(shifts)
    fprintf('Shift %d of %d\n',shift_cnt,n_shifts);
    d2 = dist_shift2d(filt_bank,shifts(xx),2,0);
    shifted_filt_bank = reshape(d2,klen,n_tr_chs);
    
    y_pred = X_z*shifted_filt_bank;
    cur_resid = YY(:,tr_set) - y_pred;
    
    %     LLs = -2*log(abs(YY(:,tr_set) - y_pred));
    %     frame_LLs(:,shift_cnt) = sum(LLs,2);
    frame_LLs(:,shift_cnt)  = -sum((cur_resid*inv_cov).*cur_resid,2);
    
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
zpads = 5;
[max_post,max_loc] = max(lgamma,[],2);
shift_cor = shifts(max_loc);

resh_X = reshape(X_z',[flen up_SDIM NT]);
resh_X_sh = zeros(size(resh_X));
resh_oX = reshape(oX_z',[flen up_SDIM NT]);
resh_oX_sh = zeros(size(resh_oX));
for ii = 1:NT
    d2 = dist_shift2d(resh_X(:,:,ii), -shift_cor(ii),2,0);
    resh_X_sh(:,:,ii) = d2;
    d2 = dist_shift2d(resh_oX(:,:,ii), -shift_cor(ii),2,0);
    resh_oX_sh(:,:,ii) = d2;
end
X_sh = reshape(resh_X_sh,flen*up_SDIM,NT)';
oX_sh = reshape(resh_oX_sh,flen*up_SDIM,NT)';

%% Try fitting SUs using two Xmat versions
load ./CellList.mat
good_sus = find(all(CellList(:,:,1) > 0));

clear defmod
defmod.lambda_L1x = 10;
defmod.lambda_d2XT = 200;
for cc = 1:length(good_sus)
    
    tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,good_sus(cc)));
    
    npq = 1;
    nnq = 0;
    kern_types{1} = 'lin';
    for i = 2:(1+npq+nnq)
        kern_types{i} = 'quad';
    end
    init_kerns = 0.01*randn(klen,1+npq+nnq);
    init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
    quad = createGNM(init_kerns,[1 ones(1,npq) -1*ones(1,nnq)],kern_types,defmod,stim_params);
    for i = 2:(1+npq+nnq)
        quad.mods(i).lambda_L1x = 2;
        quad.mods(i).lambda_d2XT = 10;
    end
    quad_fit{1}(cc) = fitGNM_filters(quad,oX_z,tr_spkbns,'none',[],1e-4,1e-6,0);
    
    quad_fit_sh{1}(cc) = fitGNM_filters(quad,oX_sh,tr_spkbns,'none',[],1e-4,1e-6,0);
    
end
%%
stim_params.spatial_dims = 1;
stim_params.sdim = up_SDIM;
stim_params.flen = flen;
klen = up_SDIM*flen;
clear defmod
defmod.lambda_L1x = 0;
defmod.lambda_d2T = 200;
defmod.lambda_d2X = 200;
defmod.lambda_d2XT = 400;

for cc = 1:length(use_lfps)
    fprintf('Cell %d of %d\n',cc,length(use_lfps));
    
    init_kerns = 0.01*randn(klen,1);
    init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
    glm = createGNM(init_kerns,1,kern_types,defmod,stim_params,'gauss');
    rev_glm_fit(cc) = fitGNM_filters(glm,X_sh,Yobs(tr_inds,cc),'none',[],1e-4,1e-6,1);
    %     [glm_sse, ~, ~, prate] = getLL_GNM(glm_fit(cc),X_sh,Yobs(tr_inds,cc),'none');
    [rev_glm_sse, ~, ~, prate] = getLL_GNM(rev_glm_fit(cc),X_sh,Yobs(tr_inds,cc),'none');
    yobs_var = var(Yobs(tr_inds,cc));
    %     new_glm_r2(cc) = 1-glm_sse/yobs_var;
    rev_glm_r2(cc) = 1-rev_glm_sse/yobs_var;
    
end
for cc = 1:length(use_lfps)
    [rev_glm_sse, ~, ~, prate] = getLL_GNM(rev_glm_fit(cc),X_sh,Yobs(tr_inds,cc),'none');
    rev_resid(cc,:) = Yobs(tr_inds,cc) - prate;
end

%%
% for cc = 1:96
%     fprintf('Cell %d\n',cc);
%     k = get_k_mat(glm_fit(cc));
%     kr = get_k_mat(rev_glm_fit(cc));
%     subplot(2,1,1)
%     imagesc(reshape(k,flen,up_SDIM));
%     ca = max(abs(caxis()));
%     caxis([-ca ca]);colorbar;
%     subplot(2,1,2)
%     imagesc(reshape(kr,flen,pad_SDIM));
%     ca = max(abs(caxis()));
%     caxis([-ca ca]);colorbar;
%     pause
%     clf
% end


%% NOW ITERATE
%overall prior on shifts
eps_prior_sigma = 0.15; %0.2
leps_prior = -(shifts*sp_dx).^2/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior)); %normalize
lA_tflip = repmat(leps_prior,n_shifts,1);

cdist = squareform(pdist(shifts'*sp_dx));
deps_sigma = 0.02;
lA = -cdist.^2/(2*deps_sigma^2);
lA = lA + lA_tflip; %bias towards returning to 0
lA = bsxfun(@minus,lA,logsumexp(lA,2)); %normalize


n_ITER = 4;
it_glm_fit{1} = rev_glm_fit;
it_glm_r2{1} = rev_glm_r2;
for nn = 2:n_ITER+1
    fprintf('Iteration %d of %d\n',nn-1,n_ITER);
    resid_cov = cov(rev_resid');
    inv_cov = pinv(resid_cov(tr_set,tr_set));
    
    klen = up_SDIM*flen;
    filt_bank = nan(length(use_lfps),klen);
    mod_offset = nan(length(use_lfps),1);
    for cc = 1:length(use_lfps)
        filt_bank(cc,:) = get_k_mat(it_glm_fit{nn-1}(cc));
        mod_offset(cc) = it_glm_fit{nn-1}(cc).spk_theta;
    end
    filt_bank = reshape(filt_bank(tr_set,:)',[flen up_SDIM n_tr_chs]);
    shifted_filt_bank = nan(klen,n_tr_chs);
    
    shift_cnt = 1;
    for xx = 1:length(shifts)
        fprintf('Shift %d of %d\n',shift_cnt,n_shifts);
        d2 = dist_shift2d(filt_bank,shifts(xx),2,0);
        shifted_filt_bank = reshape(d2,klen,n_tr_chs);
        
        y_pred = X_z*shifted_filt_bank;
        cur_resid = YY(:,tr_set) - y_pred;
        
        %     LLs = -2*log(abs(YY(:,tr_set) - y_pred));
        %     frame_LLs(:,shift_cnt) = sum(LLs,2);
        frame_LLs(:,shift_cnt)  = -sum((cur_resid*inv_cov).*cur_resid,2);
        
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
    it_shift_cor{nn} = shifts(max_loc);
    
    resh_X = reshape(X_z',[flen up_SDIM NT]);
    resh_X_sh = zeros(size(resh_X));
    resh_oX = reshape(oX_z',[flen up_SDIM NT]);
    resh_oX_sh = zeros(size(resh_oX));
    for ii = 1:NT
        d2 = dist_shift2d(resh_X(:,:,ii), -it_shift_cor{nn}(ii),2,0);
        resh_X_sh(:,:,ii) = d2;
        d2 = dist_shift2d(resh_oX(:,:,ii), -it_shift_cor{nn}(ii),2,0);
        resh_oX_sh(:,:,ii) = d2;
    end
    X_sh = reshape(resh_X_sh,flen*up_SDIM,NT)';    
    oX_sh = reshape(resh_oX_sh,flen*up_SDIM,NT)';
    %%
    stim_params.spatial_dims = 1;
    stim_params.sdim = up_SDIM;
    stim_params.flen = flen;
    klen = up_SDIM*flen;
    clear defmod
    defmod.lambda_L1x = 0;
    defmod.lambda_d2T = 200;
    defmod.lambda_d2X = 200;
    defmod.lambda_d2XT = 400;
    
    for cc = 1:length(use_lfps)
        fprintf('Cell %d of %d\n',cc,length(use_lfps));
        
        it_glm_fit{nn}(cc) = fitGNM_filters(it_glm_fit{nn-1}(cc),X_sh,Yobs(tr_inds,cc),'none',[],1e-4,1e-6,1);
        [rev_glm_sse, ~, ~, prate] = getLL_GNM(it_glm_fit{nn}(cc),X_sh,Yobs(tr_inds,cc),'none');
        yobs_var = var(Yobs(tr_inds,cc));
        it_glm_r2{nn}(cc) = 1-rev_glm_sse/yobs_var;
        cur_r2 = nan(nn,1);
        for tt = 1:nn
            cur_r2(tt) = it_glm_r2{tt}(cc);
        end
        cur_r2
    end
    for cc = 1:length(use_lfps)
        [rev_glm_sse, ~, ~, prate] = getLL_GNM(it_glm_fit{nn}(cc),X_sh,Yobs(tr_inds,cc),'none');
        rev_resid(cc,:) = Yobs(tr_inds,cc) - prate;
    end
    
    %% Try fitting SUs using two Xmat versions
    load ./CellList.mat
    good_sus = find(all(CellList(:,:,1) > 0));
    
    clear defmod
    defmod.lambda_L1x = 10;
    defmod.lambda_d2XT = 200;
    for cc = 1:length(good_sus)
        
        tr_spkbns = convert_to_spikebins(all_binned_spks(tr_inds,good_sus(cc)));
        
        npq = 1;
        nnq = 0;
        kern_types{1} = 'lin';
        for i = 2:(1+npq+nnq)
            kern_types{i} = 'quad';
        end
        init_kerns = 0.01*randn(klen,1+npq+nnq);
        init_kerns = bsxfun(@rdivide,init_kerns,sqrt(sum(init_kerns.^2)));
        quad = createGNM(init_kerns,[1 ones(1,npq) -1*ones(1,nnq)],kern_types,defmod,stim_params);
        for i = 2:(1+npq+nnq)
            quad.mods(i).lambda_L1x = 2;
            quad.mods(i).lambda_d2XT = 10;
        end
        quad_fit{nn}(cc) = fitGNM_filters(quad,oX_z,tr_spkbns,'none',[],1e-4,1e-6,0);
        
        quad_fit_sh{nn}(cc) = fitGNM_filters(quad,oX_sh,tr_spkbns,'none',[],1e-4,1e-6,0);
        
    end
    
end