clear all
close all
cd
addpath('~/Data/bruce/7_15_12/')

cd ~/Data/bruce/7_15_12/G029/

load ./jbeG029.em.mat
em_data = Expt; clear Expt
load ./CellList.mat
load ./G029Expts.mat
load ./eye_calibration_data

dt = 118/1e4;
frames_per_jump = 44;
ims_per_trial = 4;
dt_per_jump = frames_per_jump*dt;
fst = 1/dt;

Pix2Deg = 0.018837;
dsfrac = 1;
Fsd = 1/Pix2Deg/dsfrac;
Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
% RF_patch = [-0.1 0.9; -0.9 0.1]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [-0.5 1.3; -1.3 0.5]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [-1 2; -2. 1]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));


%%
cd ~/Data/bruce/7_15_12/G029/

Expt_nu = [14 22 27]; %these are the expts
spk_win = 0.1; %in sec
n_allunits = 96;
min_trial_dur = 1.5+spk_win;
all_xo_vec = [];
all_yo_vec = [];
all_image_vec = [];
all_angles = [];
all_jump_sizes = [];
all_expt_vec = [];
all_trial_vec = [];
all_binned_spks = [];
all_stim_num = [];
all_eyepos = [];

single_units = find(CellList(Expt_nu(1),:,1) > 0);
multi_units = setdiff(1:n_allunits,single_units);
for ee = 1:length(Expt_nu)
        
    

    fprintf('Analyzing Expt %d of %d\n',ee,length(Expt_nu));
    n_sus = length(single_units);
    load(sprintf('Expt%dClusterTimes.mat',Expt_nu(ee)));
    Trial_starts = [Expts{Expt_nu(ee)}.Trials(:).Start]/1e4;
    Trial_ends = [Expts{Expt_nu(ee)}.Trials(:).End]/1e4;
    endEvents = [Expts{Expt_nu(ee)}.Trials(:).endevent];
    Trial_durs = (Trial_ends-Trial_starts);
    %     completed_trials = 1:length(Trial_durs); %use all trials
    completed_trials = find(Trial_durs > min_trial_dur);
    Trial_im_nums = [Expts{Expt_nu(ee)}.Trials(:).se];
    Trial_angles = [Expts{Expt_nu(ee)}.Trials(:).Fa];
    
    [eye_vals,eye_ts,eye_insample] = get_eye_tracking_data(em_data,[Trial_starts(1) Trial_ends(end)]);
    [blink_data,in_blink,tot_disp_f] = get_blinks(eye_vals(:,3:4),eye_vals(:,1:2),eye_ts);
    [sac_data,in_sac,eye_speed] = get_saccades(eye_vals(:,3:4),eye_vals(:,1:2),eye_ts);
    
    %correct eye positions
    corrected_left = bsxfun(@plus,eye_vals(:,1:2)*left_gain,left_offset);
    corrected_right = bsxfun(@plus,eye_vals(:,3:4)*right_gain,right_offset);
    eye_vals(:,1:2) = corrected_left;
    eye_vals(:,3:4) = corrected_right;
    
    avg_eyepos = 0.5*eye_vals(:,1:2) + 0.5*eye_vals(:,3:4);
 
%     Trial_angles = mod(Trial_angles + 270,360);
    
    jump_sizes = [Expts{Expt_nu(ee)}.Trials(:).sz];
        
    for i = 1:length(completed_trials)
        cur_images = repmat(Trial_im_nums(completed_trials(i)),ims_per_trial,1);
        
        cur_xo = (0:(ims_per_trial-1))*jump_sizes(completed_trials(i))*cos(degtorad(Trial_angles(completed_trials(i))));
        cur_yo = (0:(ims_per_trial-1))*jump_sizes(completed_trials(i))*sin(degtorad(Trial_angles(completed_trials(i))));
        
        cur_t_starts = Trial_starts(completed_trials(i)) + (0:3)*dt_per_jump;
        cur_t_stops = cur_t_starts + spk_win;
        cur_t_cents = 0.5*cur_t_starts + 0.5*cur_t_stops;
        cur_t_edges = sort([cur_t_starts cur_t_stops]);

        interp_eyepos = interp1(eye_ts,eye_vals,cur_t_cents);
        
        cur_binned_spks = nan(n_allunits,length(cur_t_edges));
        for j = 1:n_allunits
            cur_binned_spks(j,:) = histc(Clusters{(j)}.times,cur_t_edges);
        end
        cur_binned_spks = cur_binned_spks(:,1:2:length(cur_t_edges));
        
        all_image_vec = [all_image_vec; cur_images];
        all_xo_vec = [all_xo_vec; cur_xo'];
        all_yo_vec = [all_yo_vec; cur_yo'];
        all_expt_vec = [all_expt_vec; ones(length(cur_images),1)*Expt_nu(ee)];
        all_trial_vec = [all_trial_vec; ones(length(cur_images),1)*completed_trials(i)];
        all_binned_spks = [all_binned_spks; cur_binned_spks'];
        all_angles = [all_angles; repmat(Trial_angles(completed_trials(i)),ims_per_trial,1)];
        all_jump_sizes = [all_jump_sizes; repmat(jump_sizes(completed_trials(i)),ims_per_trial,1)];
        all_stim_num = [all_stim_num; (1:ims_per_trial)'];
        all_eyepos = [all_eyepos; interp_eyepos];
    end
end

%%

cd ~/James_scripts/data_processing/Images/image_set_A
tot_images = length(unique(all_image_vec));
used_images = unique(all_image_vec);
tot_samps = length(all_image_vec);
all_im_patches = nan(tot_samps,length(ypatch_inds),length(xpatch_inds));
for i = 1:tot_images
    filename = sprintf('%.4d.png',used_images(i));
    IMAGEorg = imread(filename);
    IMAGEorg = double(IMAGEorg); % convert to double format
    IMAGEorg = IMAGEorg - 127.5; %shift range to [-127.5:127.5]
    IMAGE = flipud(IMAGEorg);
    IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
    
    cur_samp_set = find(all_image_vec == used_images(i));
    fprintf('Analyzing image %d, %d samps\n',i,length(cur_samp_set));
    
    for j = 1:length(cur_samp_set)
        ypatch_inds_adj = round(ypatch_inds - all_yo_vec(cur_samp_set(j))*Fsd);
        xpatch_inds_adj = round(xpatch_inds - all_xo_vec(cur_samp_set(j))*Fsd);
        cur_im = IMAGE(ypatch_inds_adj,xpatch_inds_adj);
        all_im_patches(cur_samp_set(j),:,:) = cur_im;
    end
end

fullX = reshape(all_im_patches,size(all_im_patches,1),length(ypatch_inds)*length(xpatch_inds));

%%
xv_frac = 0.;
NT = size(fullX,1);
xv_NT = round(xv_frac*NT);
xv_samp = randperm(NT);
% xv_samp = 1:NT;
xv_samp(xv_NT+1:end) = [];
tr_samp = setdiff(1:NT,xv_samp);

tr_samp(all_stim_num(tr_samp)==1) = [];
xv_samp(all_stim_num(xv_samp)==1) = [];

%%
cd ~/Data/bruce/7_15_12/
cd G029/
load ./grating_mu_data
gr_oris = unique_oris;
load ./oned_fixation_fits_v3.mat

sdim = length(xpatch_inds);
bandwidth = 0.3;
ar = 1.25;
kern_len = 21-1;
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

for t = 1:96
    init_params(1) = mean_x(t); %x0
    init_params(2) = mean_y(t); %y0
    init_params(3) = degtorad(pref_mu_ori(t)); %theta
    init_params(4) = 1/6; %lambda
    init_params(5) = 0.3*init_params(4); %sigma
    init_params(6) = 1; %eccentricity of gaussian ellipse
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,init_params,0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,init_params,pi/2);
    gabor_emp1_filt(t,:) = gabor_emp1(:);
    gabor_emp2_filt(t,:) = gabor_emp2(:);
end

%%
all_gabor_outs1 = gabor_emp1_filt*fullX';
all_gabor_outs2 = gabor_emp2_filt*fullX';
gabor_emp_X = sqrt(all_gabor_outs1.^2 + all_gabor_outs2.^2)';

for t = 1:96
    cur_binned_spks = all_binned_spks(:,t);
    
    [beta(t,:),dev(t),stats(t)] = glmfit(squeeze(gabor_emp_X(tr_samp,t)),cur_binned_spks(tr_samp),'poisson');
    p_vals(t,:) = stats(t).p;
    pred_r = squeeze(gabor_emp_X(tr_samp,t))*beta(t,2) + beta(t,1);
    pred_r = exp(pred_r);
    LL(t) = -sum(cur_binned_spks(tr_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(tr_samp));
    if xv_frac > 0
        pred_r = squeeze(gabor_emp_X(xv_samp,t))*beta(t,2) + beta(t,1);
        pred_r = exp(pred_r);
        xvLL(t) = -sum(cur_binned_spks(xv_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(xv_samp));
    end
    avg_rate(t) = mean(cur_binned_spks(tr_samp));
    null_LL(t) = -sum(bsxfun(@minus,cur_binned_spks(tr_samp)*log(avg_rate(t)),avg_rate(t)))/sum(cur_binned_spks(tr_samp));
    if xv_frac > 0
        null_xvLL(t) = -sum(bsxfun(@minus,cur_binned_spks(xv_samp)*log(avg_rate(t)),avg_rate(t)))/sum(cur_binned_spks(xv_samp));
    end
end
%% COMPUTE BETA FOR EACH MULT FACTOR

gain_mults = [0.25 0.5 0.75 1 1.5 2 3 4 5];
for t = 1:96;
    fprintf('SU %d of %d\n',t,96);
    cur_binned_spks = all_binned_spks(:,t);
    tr_spikebins = convert_to_spikebins(cur_binned_spks(tr_samp));
    
    [fitp,grad] = GLMsolve_jmm(gabor_emp_X(tr_samp,t), tr_spikebins, [0; 0], 1, [], [], [], [], [], [], 0);
    beta(t,:) = fitp.k;
    LL(t) = fitp.LL;
    for i = 1:length(gain_mults)
        [fitp,grad] = GLMsolve_jmm(gabor_emp_X(tr_samp,t), tr_spikebins, [gain_mults(i)*beta(t,1); 0], 1, [], [], [], [], [], [1], 0);
        beta_m(t,i,:) = fitp.k;
        LL_m(t,i) = fitp.LL;
    end
end

%% COMPUTE MASTER LL ARRAY
[stimlen,klen] = size(fullX);

max_shift = 30;
dshift = 2;
x_shifts = -max_shift:dshift:max_shift;
y_shifts = -max_shift:dshift:max_shift;
[Xsh,Ysh] = meshgrid(x_shifts,y_shifts);
SH = [Xsh(:) Ysh(:)];
n_shifts = size(SH,1);

fix_LLs = zeros(stimlen,96,length(gain_mults),n_shifts);

for fix = 1:stimlen
    
    fprintf('Fixation %d of %d\n',fix,stimlen);
    
    cur_im = reshape(fullX(fix,:),sdim,sdim);
    Robs = all_binned_spks(fix,:)';
    
    shifted_ims = nan(length(x_shifts)*length(y_shifts),sdim^2);
    cnt = 1;
    for xx = 1:length(x_shifts)
        for yy = 1:length(y_shifts)
            d2 = dist_shift2d(cur_im, x_shifts(xx), 2,0);
            d2 = dist_shift2d(d2,y_shifts(yy),1,0);
            shifted_ims(cnt,:) = d2(:);
            cnt = cnt + 1;
        end
    end
    
    gabor_outs1 = gabor_emp1_filt*shifted_ims';
    gabor_outs2 = gabor_emp2_filt*shifted_ims';
    gabor_outs = gabor_outs1.^2 + gabor_outs2.^2;
    gabor_outs = sqrt(gabor_outs);
    gabor_outs = reshape(gabor_outs,[96 1 n_shifts]);
    pred_Rs = bsxfun(@times,gabor_outs,beta_m(:,:,1));
    pred_Rs = bsxfun(@plus,pred_Rs,beta_m(:,:,2));
    pred_Rs = log(1+exp(pred_Rs));
    LLs = bsxfun(@times,log(pred_Rs),Robs);
    LLs = LLs - pred_Rs;
    
    fix_LLs(fix,:,:,:) = fix_LLs(fix,:,:,:) + reshape(LLs,[1 size(LLs)]);
end

%%
avg_x_pos = 0.5*all_eyepos(:,1) + 0.5*all_eyepos(:,3);
avg_y_pos = 0.5*all_eyepos(:,2) + 0.5*all_eyepos(:,4);
fix_deltaxs = [0; diff(avg_x_pos)];
fix_deltays = [0; diff(avg_y_pos)];

%%
pvals = [stats(:).p];
sig_fits = find(pvals(2,:) < .01);

NSIG = length(sig_fits);
xv_frac = 0.2;
tr_set = randperm(NSIG);
tr_set = tr_set(1:round(length(tr_set)*(1-xv_frac)));
xv_set = setdiff(1:NSIG,tr_set);
tr_set = sig_fits(tr_set);
xv_set = sig_fits(xv_set);
bad_set = find(~ismember(1:96,sig_fits));

null_shift = find(SH(:,1) == 0 & SH(:,2) == 0);
clear lpeps lptheta
n_iter = 2;
lpeps(1,:,:) = nan(stimlen,n_shifts);

mval = 1e-50;
% lptheta(1,:,:) = mval*ones(96,length(gain_mults));
% lptheta(1,:,6) = 1;
lptheta(1,:,:) = ones(96,length(gain_mults));
lptheta(1,:,:) = log(lptheta(1,:,:));
lptheta(1,:,:) = bsxfun(@minus,lptheta(1,:,:),logsumexp(lptheta(1,:,:),3));

eps_prior_sigma = 0.3;
leps_prior = zeros(n_shifts,1);
leps_prior = -sum((SH/Fsd).^2,2)/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior));

% fix_deltaxs = zeros(stimlen,1);
% fix_deltays = zeros(stimlen,1);

fix_deltax_eps = round(fix_deltaxs*Fsd/dshift);
fix_deltay_eps = round(fix_deltays*Fsd/dshift);
fix_delta_d = sqrt(fix_deltaxs.^2 + fix_deltays.^2);
max_delta_d = 0.5;

deps_sigma = 0.1*Fsd;
NT = stimlen;

for it = 2:n_iter + 1
    
    fprintf('Iteration %d of %d...',it-1,n_iter);
        
    %% COMPUTE POSTERIOR ON EPS GIVEN P(theta)
    temp = reshape(lptheta(it-1,:,:),[1 96 length(gain_mults)]);
    lposterior_eps = bsxfun(@plus,fix_LLs,temp); %multiply L(theta,eps) by p(theta)
    lposterior_eps = squeeze(logsumexp(lposterior_eps,3)); %integrate out theta
    
%     %assume independent errors
%     lposterior_eps = bsxfun(@minus,lposterior_eps,logsumexp(lposterior_eps,3)); %normalize (MAY NOT NEED)
%     lposterior_eps = squeeze(sum(lposterior_eps(:,tr_set,:),2)); %now sum across cells
%     lposterior_eps = bsxfun(@plus,lposterior_eps,leps_prior'); %prior on eps
%     lposterior_eps = bsxfun(@minus,lposterior_eps,logsumexp(lposterior_eps,2));
%     lpeps(it,:,:) = lposterior_eps;
    
    %use HMM
    lB = squeeze(sum(lposterior_eps(:,tr_set,:),2)); %now sum across cells
    if any(isnan(lB))
        error('NANS DETECTED in LL!')
    end
    
    fprintf(' Forward messages ')
    lalpha=zeros(NT,n_shifts);
    lscale=zeros(NT,1); %initialize rescaling parameters
    %compute rescaled forward messages
    lalpha(1,:) = leps_prior' + lB(1,:);
    lscale(1)=logsumexp(lalpha(1,:));
    lalpha(1,:) = lalpha(1,:) - lscale(1);
    for t=2:NT
        if fix_delta_d(t) < max_delta_d
            cdist = pdist2(SH,bsxfun(@plus,SH,[fix_deltax_eps(t) fix_deltay_eps(t)]));
            lA = -cdist.^2/(2*deps_sigma^2);
            lA = bsxfun(@minus,lA,logsumexp(lA,2));
        else
            lA = repmat(leps_prior',n_shifts,1);
        end
        lalpha(t,:) = logmulexp(lalpha(t-1,:),lA) + lB(t,:);
        lscale(t) = logsumexp(lalpha(t,:));
        lalpha(t,:)= lalpha(t,:) - lscale(t);
    end
    
    fprintf(' Backwards messages ')
    %compute rescaled backward messages
    lbeta(NT,:)=log(ones(1,n_shifts)) - lscale(NT);
    for t=NT-1:-1:1
        if fix_delta_d(t) < max_delta_d
            cdist = pdist2(SH,bsxfun(@plus,SH,[fix_deltax_eps(t) fix_deltay_eps(t)]));
            lA = -cdist.^2/(2*deps_sigma^2);
            lA = bsxfun(@minus,lA,logsumexp(lA,2));
        else
            lA = repmat(leps_prior',n_shifts,1);
        end
        %         beta(t,:)=(beta(t+1,:).*B(t+1,:))*(A')/scale(t);
        lf1 = lbeta(t+1,:) + lB(t+1,:);
        lbeta(t,:) = logmulexp(lf1,lA') - lscale(t);
    end
    
    %compute posteriors over hidden states
    lgamma= lalpha + lbeta;
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    lpeps(it,:,:) = lgamma;
    
    %% NOW COMPUTE P(THETA) GIVEN P(EPS)
    fprintf(' Computing p(theta)\n');
    temp = reshape(lpeps(it,:,:),[stimlen 1 1 n_shifts]);
    lposterior_theta = bsxfun(@plus,fix_LLs,temp); %multiply like by current estimates of p(eps)
    lposterior_theta = logsumexp(lposterior_theta,4); %sum over epsilon
    lposterior_theta = squeeze(sum(lposterior_theta)); %sum log posteriors over fixations
    llike_theta(it,:,:) = lposterior_theta;
    lposterior_theta = bsxfun(@minus,lposterior_theta,logsumexp(lposterior_theta,2)); %normalize p(thetas)
    lptheta(it,:,:) = lposterior_theta;
        
end

% FIND MAP SEQUENCE OF ERRORS
clear eps
for i = 1:n_iter
    [max_post,max_loc] = max(lpeps(i+1,:,:),[],3);
    eps(i,:,:) = SH(max_loc,:);
end

%% create new estimate of stimulus
clear beta beta_sh LL_sh

% REFIT GE MODELS
gabor_emp_X1 = fullX*gabor_emp1_filt';
gabor_emp_X2 = fullX*gabor_emp2_filt';
gabor_emp_X = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);
for ii = 1:n_iter
    fullX_sh = zeros(size(fullX));
    for fix = 1:stimlen
        fprintf('Fixation %d of %d\n',fix,stimlen);
        cur_im = reshape(fullX(fix,:),sdim,sdim);
        d2 = dist_shift2d(cur_im, eps(ii,fix,1), 2,0);
        d2 = dist_shift2d(d2,eps(ii,fix,2),1,0);
        fullX_sh(fix,:) = d2(:);
    end
    
    gabor_emp_X1 = fullX_sh*gabor_emp1_filt';
    gabor_emp_X2 = fullX_sh*gabor_emp2_filt';
    gabor_emp_X_sh = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);
    for t = 1:96;
        fprintf('SU %d of %d\n',t,96);
        spikebins = convert_to_spikebins(all_binned_spks(:,t)); 
        avg_rate(t) = mean(all_binned_spks(:,t));
        n_spks = sum(all_binned_spks(:,t));
        null_LL(t) = -sum(bsxfun(@minus,all_binned_spks(:,t)*log(avg_rate(t)),avg_rate(t)))/n_spks;
                
        [fitp,grad] = GLMsolve_jmm(gabor_emp_X_sh(:,t), spikebins, [0; 0], 1, [], [], [], [], [], [], 0);
        beta_sh(ii,t,:) = fitp.k;
        LL_sh(ii,t) = fitp.LL;        
        
        [fitp,grad] = GLMsolve_jmm(gabor_emp_X(:,t), spikebins, [0; 0], 1, [], [], [], [], [], [], 0);
        beta(t,:) = fitp.k;
        LL(t) = fitp.LL;               
    end
end
LL_imp = bsxfun(@minus,LL_sh,LL);

%% FIT RF POSITION MAPs
nsh1 = sqrt(n_shifts);
cur_lpeps = squeeze(lpeps(end,:,:));
lposmap = zeros(n_shifts,96);
for cid = 1:96;
    fprintf('Cell %d of %d\n',cid,96);
    cur_lptheta = squeeze(lptheta(end,cid,:));
    %     cur_lptheta = zeros(1,length(gain_mults));
    %     cur_lptheta(4) = 1;
    %     cur_lptheta = log(cur_lptheta)';
    for t = 1:NT
        %         fprintf('Fixation %d of %d\n',t,NT);
        temp = bsxfun(@plus,squeeze(fix_LLs(t,cid,:,:)),cur_lptheta);
        cur_LL = reshape(logsumexp(temp,1),nsh1,nsh1);
        cur_lprior = reshape(cur_lpeps(t,:),nsh1,nsh1);
        max_L = max(cur_LL(:));
        %         max_p = max(cur_lprior(:));
        cur_L = exp(cur_LL-max_L);
        cur_prior = exp(cur_lprior);
        
        %         cur_prior(cur_prior < max(cur_prior(:))) = 1e-100; cur_prior = bsxfun(@rdivide,cur_prior,sum(cur_prior(:)));
        
        %         cur_L = exp(bsxfun(@minus,cur_LL,max_L));
        %         cur_prior = exp(bsxfun(@minus,cur_lprior,max_p));
        cur_conv = conv2(cur_L,fliplr(flipud(cur_prior)),'same');
        cur_conv_temp = conv2(ones(nsh1,nsh1),fliplr(flipud(cur_prior)),'same');
        cur_conv = cur_conv./cur_conv_temp;
        cur_conv(isnan(cur_conv)) = 0;
        %         lcur_conv = bsxfun(@plus,log(cur_conv),max_L);
        lcur_conv = log(cur_conv);
        lposmap(:,cid) = lposmap(:,cid) + lcur_conv(:);
    end
end
%
prior_sigma = 0.1*Fsd;
poserror_prior = exp(-sum(SH.^2,2)./(2*prior_sigma^2));
poserror_prior = poserror_prior/sum(poserror_prior);

lposmap_wprior = bsxfun(@plus,lposmap,log(poserror_prior));

% lposmap = bsxfun(@minus,lposmap,logsumexp(lposmap,1));
lposmap_wprior = bsxfun(@minus,lposmap_wprior,logsumexp(lposmap_wprior,1));
[a,b] = max(lposmap_wprior);
% [a,b] = max(lposmap);
rf_x = SH(b,1);
rf_y = SH(b,2);
dd = sqrt(rf_x.^2+rf_y.^2)/Fsd;

%%
clear gabor_emp*

for t = 1:96
    init_params(1) = mean_x(t); %x0
    init_params(2) = mean_y(t); %y0
    init_params(3) = degtorad(pref_mu_ori(t)); %theta
    init_params(4) = 1/6; %lambda
    init_params(5) = 0.3*init_params(4); %sigma
    init_params(6) = 1; %eccentricity of gaussian ellipse
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,init_params,0);
    gabor_emp2 = get_pgabor_mask_v2(XX,Y0Y,init_params,pi/2);
    gabor_emp1_filt(t,:) = gabor_emp1(:);
    gabor_emp2_filt(t,:) = gabor_emp2(:);
    
    init_params(1) = init_params(1) - rf_x(t)/Fsd;
    init_params(2) = init_params(2) - rf_y(t)/Fsd;
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,init_params,0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,init_params,pi/2);
    gabor_emp1_filtc(t,:) = gabor_emp1(:);
    gabor_emp2_filtc(t,:) = gabor_emp2(:);
end

% [NT,klen] = size(fullX);

gabor_emp_X1 = fullX*gabor_emp1_filt';
gabor_emp_X2 = fullX*gabor_emp2_filt';
gabor_emp_X = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);

gabor_emp_X1 = fullX*gabor_emp1_filtc';
gabor_emp_X2 = fullX*gabor_emp2_filtc';
gabor_emp_Xc = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);

for t = 1:96;
    fprintf('SU %d of %d\n',t,96);
    cur_binned_spks = all_binned_spks(:,t);
    
    [beta(t,:),dev(t),stats(t)] = glmfit(squeeze(gabor_emp_X(:,t)),cur_binned_spks(:),'poisson');
    [betac(t,:),devc(t),statsc(t)] = glmfit(squeeze(gabor_emp_Xc(:,t)),cur_binned_spks(:),'poisson');
    
end
