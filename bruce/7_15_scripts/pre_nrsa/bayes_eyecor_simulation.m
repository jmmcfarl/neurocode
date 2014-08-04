%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12
cd G029/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

% load ./Expt1_compiled_data_fixedlag.mat
% load ./Expt1_compiled_data_fixedlag.mat
cd ~/Data/bruce/7_15_12/G029/
load ./Expt1_compiled_data_fixedlag
Pix2Deg = 0.018837;

%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

cd ~/Data/bruce/7_15_12/
cd G029/
load ./grating_mu_data
gr_oris = unique_oris;
[~,sf_inds] = max(avg_mu_sf_profile,[],2);
pref_sfs = unique_sfs(sf_inds);
load ./oned_fixation_fits_v2.mat

rf_std = 0.1;
sim_rf_x = mean_x + randn(size(mean_x))*rf_std; 
sim_rf_y = mean_y + randn(size(mean_y))*rf_std; 

clear gabor_emp*
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
    
    init_params(1) = sim_rf_x(t);
    init_params(2) = sim_rf_y(t);
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,init_params,0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,init_params,pi/2);
    gabor_emp1_filt_sim(t,:) = gabor_emp1(:);
    gabor_emp2_filt_sim(t,:) = gabor_emp2(:);    
end

[NT,klen] = size(fullX);

%% PARSE DATA INTO FIXATIONS
min_fix_dur = 0.1;
trial_flips = 1 + find(diff(full_trial_vec) ~= 0);
sac_inds = find(full_insac==1);
all_break_inds = unique([trial_flips; sac_inds]);
fix_start_inds = [1; all_break_inds+1];
fix_stop_inds = [all_break_inds-1; NT];
fix_durs = (fix_stop_inds-fix_start_inds)*dt;
used_fixs = find(fix_durs > min_fix_dur);
fprintf('Using %d of %d fixations\n',length(used_fixs),length(fix_durs));

used_ims = [];
for i = 1:length(used_fixs)
    used_ims = [used_ims fix_start_inds(used_fixs(i)):fix_stop_inds(used_fixs(i))];
end

%use only left eye position
avg_x_pos = full_eyepos(:,1);
avg_y_pos = full_eyepos(:,2);

x_drifts = nan(size(avg_x_pos));
y_drifts = nan(size(avg_y_pos));
medx_pos = nan(size(used_fixs));
medy_pos = nan(size(used_fixs));
for i = 1:length(used_fixs)
    cur_set = fix_start_inds(used_fixs(i)):fix_stop_inds(used_fixs(i));
    medx_pos(i) = median(avg_x_pos(cur_set));
    medy_pos(i) = median(avg_y_pos(cur_set));
    x_drifts(cur_set) = avg_x_pos(cur_set) - medx_pos(i);
    y_drifts(cur_set) = avg_y_pos(cur_set) - medy_pos(i);
end
fix_deltaxs = [0; diff(medx_pos)];
fix_deltays = [0; diff(medy_pos)];

%% CREATE SIMULATED STIMULUS
sigma_deltax = 0.05;
max_deltax = 0.5;
sim_deltaxs = fix_deltaxs + randn(size(fix_deltaxs))*sigma_deltax;
sim_deltays = fix_deltays + randn(size(fix_deltays))*sigma_deltax;
% sim_deltaxs(sim_deltaxs > max_deltax) = max_deltax;
% sim_deltays(sim_deltays > max_deltax) = max_deltax;
% sim_deltaxs(sim_deltaxs < -max_deltax) = -max_deltax;
% sim_deltays(sim_deltays < -max_deltax) = -max_deltax;

[b,a] = butter(2,.025/2,'high');
sim_epsx = cumsum(sim_deltaxs);
sim_epsy = cumsum(sim_deltays);
sim_epsx = filtfilt(b,a,sim_epsx);

max_eps = 0.5;
sim_epsx(sim_epsx > max_eps) = max_eps;
sim_epsx(sim_epsx < -max_eps) = -max_eps;
sim_epsy(sim_epsy > max_eps) = max_eps;
sim_epsy(sim_epsy < -max_eps) = -max_eps;

sim_epsx_pix = round(sim_epsx*Fsd);
sim_epsy_pix = round(sim_epsy*Fsd);

fullX_sh = zeros(size(fullX));
for fix = 1:length(used_fixs)
    fprintf('Fixation %d of %d\n',fix,length(used_fixs));
    cur_ims = fix_start_inds(used_fixs(fix)):fix_stop_inds(used_fixs(fix));
    for jj = 1:length(cur_ims);
        cur_im = reshape(fullX(cur_ims(jj),:),sdim,sdim);
        d2 = dist_shift2d(cur_im, sim_epsx_pix(fix), 2,0);
        d2 = dist_shift2d(d2,sim_epsy_pix(fix),1,0);
        fullX_sh(cur_ims(jj),:) = d2(:);
    end
end

gabor_emp_X1 = fullX_sh*gabor_emp1_filt_sim';
gabor_emp_X2 = fullX_sh*gabor_emp2_filt_sim';
gabor_emp_X_sim = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);

%% CREATE SIMULATED SPIKES
true_gains = rand(96,1)*2e-3 + 6e-4;
poss_bs = linspace(-3,1,250);
true_b = nan(96,1);
sim_spikes_binned = zeros(size(full_binned_spks));
for i = 1:96
    fprintf('Cell %d of %d\n',i,96);
   avg_rate(i) = mean(full_binned_spks(used_ims,i));
   g = gabor_emp_X_sim(:,i)*true_gains(i);
   poss_rates = log(1+exp(bsxfun(@plus,g,poss_bs)));
   implied_mrates = mean(poss_rates);
   [a,b] = min(abs(implied_mrates - avg_rate(i)));
   true_b(i) = poss_bs(b);
   sim_rate = log(1+exp(g + true_b(i)));
   sim_spikes_binned(:,i) = poissrnd(sim_rate);
end

%% COMPUTE BETA FOR EACH MULT FACTOR
gabor_emp_X1 = fullX*gabor_emp1_filt';
gabor_emp_X2 = fullX*gabor_emp2_filt';
gabor_emp_X = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);

trial_gain_mults = [0.5 1 1.5 2 4];
for t = 1:96;
    fprintf('SU %d of %d\n',t,96);
    spikebins = convert_to_spikebins(sim_spikes_binned(used_ims,t));
    
    %for trial models
    [fitp,grad] = GLMsolve_jmm( gabor_emp_X(used_ims,t), spikebins, [0; 0], 1, [], [], [], [], [], [], 0);
    beta(t,:) = fitp.k;
    LL(t) = fitp.LL;
    for i = 1:length(trial_gain_mults)
        [fitp,grad] = GLMsolve_jmm( gabor_emp_X(used_ims,t), spikebins, [trial_gain_mults(i)*beta(t,1); 0], 1, [], [], [], [], [], [1], 0);
        beta_m(t,i,:) = fitp.k;
        LL_m(t,i) = fitp.LL;
    end
end

%% COMPUTE LL ARRAY
max_shift = 20;
dshift = 2;
x_shifts = -max_shift:dshift:max_shift;
y_shifts = -max_shift:dshift:max_shift;
[Xsh,Ysh] = meshgrid(x_shifts,y_shifts);
SH = [Xsh(:) Ysh(:)];
n_shifts = size(SH,1);

fix_LLs = zeros(length(used_fixs),96,length(trial_gain_mults),n_shifts);

for fix = 1:length(used_fixs)
    
    fprintf('Fixation %d of %d\n',fix,length(used_fixs));
    
    tot_LLs = zeros(96,length(trial_gain_mults),n_shifts);
    cur_ims = fix_start_inds(used_fixs(fix)):fix_stop_inds(used_fixs(fix));
    
    for ii = 1:length(cur_ims);
        %         fprintf('Image %d of %d\n',ii,length(cur_ims));
        cur_im = reshape(fullX(cur_ims(ii),:),sdim,sdim);
        Robs = sim_spikes_binned(cur_ims(ii),:)';
      
        tic
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
        toc
        
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
end

%%
NSIG = 96;
xv_frac = 0.3;
tr_set = randperm(NSIG);
tr_set = tr_set(1:round(length(tr_set)*(1-xv_frac)));
xv_set = setdiff(1:NSIG,tr_set);

null_shift = find(SH(:,1) == 0 & SH(:,2) == 0);

clear lpeps lptheta
n_iter = 3;
lpeps(1,:,:) = nan(length(used_fixs),n_shifts);

lptheta(1,:,:) = log(ones(96,length(trial_gain_mults))/length(trial_gain_mults)); %uniform prior
% mval = 1e-50;
% lptheta(1,:,:) = mval*ones(96,length(trial_gain_mults));
% lptheta(1,:,2) = 1;
% lptheta(1,:,:) = log(lptheta(1,:,:));
% lptheta(1,:,:) = bsxfun(@minus,lptheta(1,:,:),logsumexp(lptheta(1,:,:),3));

eps_prior_sigma = 0.2;
leps_prior = zeros(n_shifts,1);
leps_prior = -sum((SH/Fsd).^2,2)/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior));

fix_deltax_eps = round(fix_deltaxs*Fsd/dshift);
fix_deltay_eps = round(fix_deltays*Fsd/dshift);
fix_delta_d = sqrt(fix_deltaxs.^2 + fix_deltays.^2);
max_delta_d = 0.5;

deps_sigma = 1*Fsd;
NT = length(used_fixs);

for it = 2:n_iter + 1
    
    fprintf('Iteration %d of %d: ',it-1,n_iter);
    
    %% COMPUTE POSTERIOR ON EPS GIVEN P(theta)
    %%
    temp = reshape(lptheta(it-1,:,:),[1 96 length(trial_gain_mults)]);
    lposterior_eps = bsxfun(@plus,fix_LLs,temp); %multiply L(theta,eps) by p(theta)
    lposterior_eps_s = squeeze(logsumexp(lposterior_eps,3)); %integrate out theta
    lB = squeeze(sum(lposterior_eps_s(:,tr_set,:),2)); %now sum across cells
    if any(isnan(lB))
        error('NANS DETECTED in LL!')
    end
    
    fprintf(' Forward messages... ')
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
    
    fprintf(' Backwards messages... ')
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
        lf1 = lbeta(t+1,:) + lB(t+1,:);
        lbeta(t,:) = logmulexp(lf1,lA') - lscale(t);
    end
    
    %compute posteriors over hidden states
    lgamma= lalpha + lbeta;
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    lpeps(it,:,:) = lgamma;
    
    %% NOW COMPUTE P(THETA) GIVEN P(EPS)
    fprintf(' Computing p(theta)...\n');
    temp = reshape(lpeps(it,:,:),[length(used_fixs) 1 1 n_shifts]);
    lposterior_theta = bsxfun(@plus,fix_LLs,temp); %multiply like by current estimates of p(eps)
    lposterior_theta = logsumexp(lposterior_theta,4); %sum over epsilon
    %     lposterior_theta = bsxfun(@minus,lposterior_theta,logsumexp(lposterior_theta,3)); %normalize to posteriors on theta
    lposterior_theta = squeeze(sum(lposterior_theta)); %sum log posteriors over fixations
    llike_theta(it,:,:) = lposterior_theta;
    lposterior_theta = bsxfun(@minus,lposterior_theta,logsumexp(lposterior_theta,2)); %normalize p(thetas)
    lptheta(it,:,:) = lposterior_theta;
    
    %%
    
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
    fullX_she = zeros(size(fullX));
    for fix = 1:length(used_fixs)
        fprintf('Fixation %d of %d\n',fix,length(used_fixs));
        cur_ims = fix_start_inds(used_fixs(fix)):fix_stop_inds(used_fixs(fix));
        for jj = 1:length(cur_ims);
            cur_im = reshape(fullX(cur_ims(jj),:),sdim,sdim);
            d2 = dist_shift2d(cur_im, eps(ii,fix,1), 2,0);
            d2 = dist_shift2d(d2,eps(ii,fix,2),1,0);
            fullX_she(cur_ims(jj),:) = d2(:);
        end
    end
    
    gabor_emp_X1 = fullX_she*gabor_emp1_filt';
    gabor_emp_X2 = fullX_she*gabor_emp2_filt';
    gabor_emp_X_sh = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);
    for t = 1:96;
        fprintf('SU %d of %d\n',t,96);
        spikebins = convert_to_spikebins(sim_spikes_binned(used_ims,t)); 
        avg_rate(t) = mean(sim_spikes_binned(used_ims,t));
        n_spks = sum(sim_spikes_binned(used_ims,t));
        null_LL(t) = -sum(bsxfun(@minus,sim_spikes_binned(used_ims,t)*log(avg_rate(t)),avg_rate(t)))/n_spks;
                
        [fitp,grad] = GLMsolve_jmm(gabor_emp_X_sh(used_ims,t), spikebins, [0; 0], 1, [], [], [], [], [], [], 0);
        beta_sh(ii,t,:) = fitp.k;
        LL_sh(ii,t) = fitp.LL;        
        
        [fitp,grad] = GLMsolve_jmm(gabor_emp_X(used_ims,t), spikebins, [0; 0], 1, [], [], [], [], [], [], 0);
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
    for t = 1:NT
        temp = bsxfun(@plus,squeeze(fix_LLs(t,cid,:,:)),cur_lptheta);
        cur_LL = reshape(logsumexp(temp,1),nsh1,nsh1);
        cur_lprior = reshape(cur_lpeps(t,:),nsh1,nsh1);
        max_L = max(cur_LL(:));
        cur_L = exp(cur_LL-max_L);
        cur_prior = exp(cur_lprior);
        
        cur_conv = conv2(cur_L,rot90(cur_prior,2),'same');
        cur_conv_temp = conv2(ones(nsh1,nsh1),rot90(cur_prior,2),'same');
        cur_conv = cur_conv./cur_conv_temp;
        cur_conv(isnan(cur_conv)) = 0;
        lcur_conv = log(cur_conv);
        lposmap(:,cid) = lposmap(:,cid) + lcur_conv(:);
    end
end
%
prior_sigma = 0.25*Fsd;
poserror_prior = exp(-sum(SH.^2,2)./(2*prior_sigma^2));
poserror_prior = poserror_prior/sum(poserror_prior);

lposmap_wprior = bsxfun(@plus,lposmap,log(poserror_prior));

% lposmap = bsxfun(@minus,lposmap,logsumexp(lposmap,1));
lposmap_wprior = bsxfun(@minus,lposmap_wprior,logsumexp(lposmap_wprior,1));
[a,b] = max(lposmap_wprior);
% [a,b] = max(lposmap);
rf_x = -SH(b,1);
rf_y = -SH(b,2);
dd = sqrt(rf_x.^2+rf_y.^2)/Fsd;

%% REFIT GABOR PARAMS to see if RF correction improves fits
clear gabor_emp*

for t = 1:96
    init_params(1) = mean_x(t); %x0
    init_params(2) = mean_y(t); %y0
    init_params(3) = degtorad(pref_mu_ori(t)); %theta
    init_params(4) = 1/6; %lambda
    %     init_params(4) = 1/pref_sfs(t);
    init_params(5) = 0.3*init_params(4); %sigma
    init_params(6) = 1; %eccentricity of gaussian ellipse
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,init_params,0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,init_params,pi/2);
    gabor_emp1_filt(t,:) = gabor_emp1(:);
    gabor_emp2_filt(t,:) = gabor_emp2(:);
    
    init_params(1) = init_params(1) + rf_x(t)/Fsd;
    init_params(2) = init_params(2) + rf_y(t)/Fsd;
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,init_params,0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,init_params,pi/2);
    gabor_emp1_filtc(t,:) = gabor_emp1(:);
    gabor_emp2_filtc(t,:) = gabor_emp2(:);
end

% [NT,klen] = size(fullX);

gabor_emp_X1 = fullX_she*gabor_emp1_filt';
gabor_emp_X2 = fullX_she*gabor_emp2_filt';
gabor_emp_X = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);

gabor_emp_X1 = fullX_she*gabor_emp1_filtc';
gabor_emp_X2 = fullX_she*gabor_emp2_filtc';
gabor_emp_Xc = sqrt(gabor_emp_X1.^2 + gabor_emp_X2.^2);

for t = 1:96;
    fprintf('SU %d of %d\n',t,96);
    cur_binned_spks = sim_spikes_binned(:,t);
    
    [beta(t,:),dev(t),stats(t)] = glmfit(squeeze(gabor_emp_X(used_ims,t)),cur_binned_spks(used_ims),'poisson');
    [betac(t,:),devc(t),statsc(t)] = glmfit(squeeze(gabor_emp_Xc(used_ims,t)),cur_binned_spks(used_ims),'poisson');
    
end

%% CONSTRUCT ADJUSTED LL ARRAY USING NEW RF POSITIONS
fix_LLs_adj = fix_LLs;
adj_set = find(dd ~= 0);
NT = size(fix_LLs,1);
for c = 1:length(adj_set)
    fprintf('Adjusting cell %d of %d\n',c,length(adj_set));
    
    cur_x_adj = round(rf_x(adj_set(c))/dshift);
    cur_y_adj = round(rf_y(adj_set(c))/dshift);
    
    cur_array = reshape(fix_LLs(:,adj_set(c),:,:),[NT length(trial_gain_mults) nsh1 nsh1]);
    nullvals = squeeze(mean(fix_LLs(:,adj_set(c),:,:),4));
    mapper = mod((0:(nsh1-1))-cur_x_adj,nsh1)+1;
    cur_array_sh = cur_array(:,:,:,mapper);
    if cur_x_adj > 0
        cur_array_sh(:,:,:,1:cur_x_adj) = repmat(nullvals,[1 1 nsh1 cur_x_adj]);
    elseif cur_x_adj < 0
        cur_array_sh(:,:,:,(nsh1-cur_x_adj+1:end)) = repmat(nullvals,[1 1 nsh1 cur_x_adj]);
    end
    
    mapper = mod((0:(nsh1-1))-cur_y_adj,nsh1)+1;
    cur_array_sh = cur_array_sh(:,:,mapper,:);
    if cur_y_adj > 0
        cur_array_sh(:,:,1:cur_y_adj,:) = repmat(nullvals,[1 1 cur_y_adj nsh1]);
    elseif cur_y_adj < 0
        cur_array_sh(:,:,(nsh1-cur_y_adj+1:end),:) = repmat(nullvals,[1 1 cur_y_adj nsh1]);
    end
    
    fix_LLs_adj(:,adj_set(c),:,:) = reshape(cur_array_sh,[NT 1 length(trial_gain_mults) n_shifts]);
end

%%
null_shift = find(SH(:,1) == 0 & SH(:,2) == 0);
clear lpeps2 lptheta2
n_iter = 2;
lpeps2(1,:,:) = nan(length(used_fixs),n_shifts);

lptheta2(1,:,:) = lptheta(end,:,:);

eps_prior_sigma = 0.2;
leps_prior = zeros(n_shifts,1);
leps_prior = -sum((SH/Fsd).^2,2)/(2*eps_prior_sigma^2);
leps_prior = bsxfun(@minus,leps_prior,logsumexp(leps_prior));

fix_deltax_pix = round(fix_deltaxs*Fsd);
fix_deltay_pix = round(fix_deltays*Fsd);

deps_sigma = 0.2*Fsd;
NT = length(used_fixs);

for it = 2:n_iter + 1
    
    fprintf('Iteration %d of %d...',it-1,n_iter);
    
    %%
    temp = reshape(lptheta2(it-1,:,:),[1 96 length(trial_gain_mults)]);
    lposterior_eps = bsxfun(@plus,fix_LLs_adj,temp); %multiply L(theta,eps) by p(theta)
    lposterior_eps = squeeze(logsumexp(lposterior_eps,3)); %integrate out theta
    lB = squeeze(sum(lposterior_eps(:,tr_set,:),2)); %now sum across cells
    if any(isnan(lB))
        error('NANS DETECTED in LL!')
    end
    
    
    fprintf(' Forward messages... ')
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
    
    fprintf(' Backwards messages... ')
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
        lf1 = lbeta(t+1,:) + lB(t+1,:);
        lbeta(t,:) = logmulexp(lf1,lA') - lscale(t);
    end
    
        %compute posteriors over hidden states
    lgamma= lalpha + lbeta;
    lgamma = bsxfun(@minus,lgamma,logsumexp(lgamma,2));
    lpeps2(it,:,:) = lgamma;

    %% NOW COMPUTE P(THETA) GIVEN P(EPS)
    fprintf(' Computing p(theta)\n');
    temp = reshape(lpeps2(it,:,:),[length(used_fixs) 1 1 n_shifts]);
    lposterior_theta = bsxfun(@plus,fix_LLs_adj,temp); %multiply like by current estimates of p(eps)
    lposterior_theta = logsumexp(lposterior_theta,4); %sum over epsilon
    lposterior_theta = squeeze(sum(lposterior_theta)); %sum log posteriors over fixations
    llike_theta(it,:,:) = lposterior_theta;
    lposterior_theta = bsxfun(@minus,lposterior_theta,logsumexp(lposterior_theta,2)); %normalize p(thetas)
    lptheta2(it,:,:) = lposterior_theta;
    
    %%
    
end

% FIND MAP SEQUENCE OF ERRORS
[max_post,max_loc] = max(lpeps2(3,:,:),[],3);
eps2 = SH(max_loc,:);
