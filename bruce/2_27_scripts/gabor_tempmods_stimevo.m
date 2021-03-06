clear all
% close all
addpath(genpath('~/James_scripts'));

cd ~/Data/bruce/2_27_12/
load Blocks.mat
accept_window = [-5 5;-5 5];
sac_eyespeed = 10;
thresh_eyespeed = 2.5;
sac_buffer = 0;
min_fix_dur = 0.2;
blink_thresh = 5;
min_blink_dur = 0.05;
% nlags = 4;


Pix2Deg = 0.018837;

% down-sampling fraction for image
dsfrac = 2;
Fsd = 1/Pix2Deg/dsfrac;

Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
% [XX,YY] = meshgrid(xax,yax);

% RF_patch = [3.5 6; -3.5 -1]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [2.5 6.5; -4 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [2.75 6.25; -3.5 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

% patch_inds = find(XX >= RF_patch(1,1) & XX <= RF_patch(1,2) & YY >= RF_patch(2,1) & YY <= RF_patch(2,2));
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_data_v3

cd /Users/James/Data/bruce/2_27_12/stimrecon
load sing_eye_pgabortrack_d2
load sing_eye_pgabortrack_fin_est_d2
cd ~/Data/bruce/2_27_12/ExptA

SDIM = length(xpatch_inds);
cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_image_patches_corrected

NT = length(used_fixs);
X_resh = reshape(X,NT,SDIM^2);
dt = .005;

cellids = [1 2 3 4 5 6 7 8 10];
muaids = [1 2 4 5 6 11 12 13 14];

spk_cnts = [spk_cnts(:,cellids) mua_cnts(:,muaids)];

n_used_cells = size(spk_cnts,2);

%%
fix_stimout = zeros(n_used_cells,length(used_fixs));
all_model_time_axis = [];
all_model_blockids = [];
all_model_fixids = [];
spikes_binned = [];
time_since_fix = [];
for blockid = 1:4;
    fprintf('Block %d of %d\n',blockid,4);
    cur_set = find(blockids(used_fixs) == blockid);
    
    for i = 1:length(cur_set)
        start_T = all_fix_start_times(used_fixs(cur_set(i)));
        end_T = all_fix_stop_times(used_fixs(cur_set(i)));
        cur_tedges = start_T:dt:end_T;
        cur_tcents = 0.5*cur_tedges(1:end-1)+0.5*cur_tedges(2:end);
        all_model_time_axis = [all_model_time_axis cur_tcents];
        all_model_blockids = [all_model_blockids blockid*ones(1,length(cur_tcents))];
        all_model_fixids = [all_model_fixids cur_set(i)*ones(1,length(cur_tcents))];
        temp_binned = zeros(10,length(cur_tcents));
        temp_modouts = zeros(10,length(cur_tcents));
        for c = 1:n_used_cells
            if c <= 9
                temp = histc(Blocks{blockid}.spktimes{cellids(c)},cur_tedges);
            else
                temp = histc(Blocks{blockid}.mutimes{muaids(c-9)},cur_tedges);
            end
            temp_binned(c,:) = temp(1:end-1);
        end
        spikes_binned = [spikes_binned; temp_binned'];
        
        time_since_fix = [time_since_fix cur_tcents-start_T];
        
    end
end

%% Compute TBR time-since fix onset
max_tsf = 0.75; nbins = 30;
used_tsf = time_since_fix;
used_tsf(used_tsf > max_tsf) = [];
tax = prctile(used_tsf,linspace(100/nbins,100-100/nbins,nbins));

% tax = logspace(log10(dt),log10(max_tsf),nbins);
Tmat = tbrep(time_since_fix,tax);

%%
xv_frac = 0.2;
rperm = randperm(length(used_fixs));
xv_fixs = rperm(1:round(xv_frac*length(used_fixs)));
xv_inds = find(ismember(all_model_fixids,xv_fixs));
tr_inds = find(~ismember(all_model_fixids,xv_fixs));

% xv_fixs = find(all_stim_filtered(used_fixs)==0);
% xv_inds = find(ismember(all_model_fixids,xv_fixs));
% tr_inds = find(~ismember(all_model_fixids,xv_fixs));
%% COMPUTE SACCADE TRIG AVG FIRING RATES
sac_trg_rate = zeros(n_used_cells,length(tax));
sac_trg_var = zeros(n_used_cells,length(tax));
taxb = [0 tax];
for i = 1:length(tax)
    curset = find(time_since_fix >= taxb(i) & time_since_fix < taxb(i+1));
    for c = 1:n_used_cells
        sac_trg_rate(c,i) = mean(spikes_binned(curset,c));
        sac_trg_var(c,i) = var(spikes_binned(curset,c));
    end
    n_occ(i) = length(curset);
end
avg_rates = mean(spikes_binned)/dt;
fano_fac = sac_trg_var./sac_trg_rate;
sac_trg_rate = sac_trg_rate/dt;
sac_trg_std = sqrt(sac_trg_var)/dt;
%%
clear mod_*
sac_hist = zeros(n_used_cells,length(tax));
stim_kern = cell(n_used_cells,1);
slope2_lambda = 150;
l1_pen = 20;
stim_mod_predrate = zeros(9,length(time_since_fix));
sac_hist_predrate = zeros(9,length(time_since_fix));
for c = 1:9
    fprintf('Cell %d of %d\n',c,9);
    
    un_spk_cnts = length(unique(spikes_binned(tr_inds,c))) - 1;
    spkbs = [];
    for i = 1:un_spk_cnts
        curset = find(spikes_binned(tr_inds,c) == i);
        spkbs = [spkbs; repmat(curset,i,1)];
    end
    spkbs = sort(spkbs);
    Robs = spikes_binned(tr_inds,c);
    Robs_xv = spikes_binned(xv_inds,c);
    
    %%
    disp('Fitting sac hist model')
    %first fit sachist model
    Xmat = Tmat;
    W0 = zeros(size(Xmat,2),1);
    NT = length(tr_inds);
    
    %null model
        avg_rate = repmat(mean(Robs),length(Robs),1);
    null_LL(c) = sum(-Robs.*log(avg_rate)+avg_rate)/sum(Robs);
    avg_rate_xv = repmat(mean(Robs),length(Robs_xv),1);
    null_xvLL(c) = sum(-Robs_xv.*log(avg_rate_xv)+avg_rate_xv)/sum(Robs_xv);

    %initialize parameters
    silent = 1;
    lamrange = [];
    lamrange2 = [slope2_lambda 1 length(tax)];
    nlxrange = [tax];
    Pcon = [];
    Pmono = [];
    hold_const = [];
    NLtype = 0;
    llist = [];
    
    [fitp,grad] = GLMsolve_jmm( Xmat(tr_inds,:), spkbs, W0, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const, NLtype, 1e-6,nlxrange);
    sac_hist_k(c,:) = fitp.k(1:end-1);
    sac_hist_const(c) = fitp.k(end);
    sac_hist_LP(c) = fitp.LP;
    
    cur_genfun = Tmat(tr_inds,:)*sac_hist_k(c,:)' + sac_hist_const(c);
    cur_predrate = log(1+exp(cur_genfun));
    sac_hist_LL(c) = sum(-Robs.*log(cur_predrate) + cur_predrate)/sum(Robs);
    cur_genfun_xv = Tmat(xv_inds,:)*sac_hist_k(c,:)' + sac_hist_const(c);
    cur_predrate_xv = log(1+exp(cur_genfun_xv));
    sac_hist_xvLL(c) = sum(-Robs_xv.*log(cur_predrate_xv) + cur_predrate_xv)/sum(Robs_xv);
    cur_genfun = Tmat*sac_hist_k(c,:)' + sac_hist_const(c);
    sac_hist_predrate(c,:) = log(1+exp(cur_genfun));
    
    %%
    disp('Fitting total stim model')
    X_resh = reshape(X,size(X,1),SDIM^2);
    cur_mask1 = get_pgabor_mask(gabor_params_fin(c,1:6),0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(gabor_params_fin(c,1:6),pi/2,[SDIM SDIM]);
    mask1_out = X_resh*cur_mask1(:);
    mask2_out = X_resh*cur_mask2(:);
    lin_out = gabor_params_fin(c,8)*mask1_out(all_model_fixids) + gabor_params_fin(c,9)*mask2_out(all_model_fixids);
    energy_out = gabor_params_fin(c,7)*sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
    total_out = lin_out + energy_out;
    total_out = zscore(total_out);
    
    lamrange2 = [slope2_lambda 1 length(tax);
        slope2_lambda length(tax)+1 length(tax)*2];
    nlxrange = [tax tax];
    Xmat = Tmat;
    Xmat = [Xmat bsxfun(@times,Tmat,total_out)];
    
    llist = [l1_pen (length(tax)+1):size(Xmat,2)];
    [fitp,grad] = GLMsolve_jmm( Xmat(tr_inds,:), spkbs, W0, 1, lamrange, lamrange2,Pcon, Pmono, llist, [], NLtype,1e-5,nlxrange);
    stim_mod_const(c) = fitp.k(end);
    stim_mod_sachist(c,:) = fitp.k(1:length(tax));
    stim_mod_tot(c,:) = fitp.k(length(tax)+1:length(tax)*2);
    
    cur_genfun = Xmat(tr_inds,:)*fitp.k(1:end-1) + stim_mod_const(c);
    cur_predrate = log(1+exp(cur_genfun));
    stim_mod_LL(c) = sum(-Robs.*log(cur_predrate) + cur_predrate)/sum(Robs);
    cur_genfun_xv = Xmat(xv_inds,:)*fitp.k(1:end-1) + stim_mod_const(c);
    cur_predrate_xv = log(1+exp(cur_genfun_xv));
    stim_mod_xvLL(c) = sum(-Robs_xv.*log(cur_predrate_xv) + cur_predrate_xv)/sum(Robs_xv);
    cur_genfun = Xmat*fitp.k(1:end-1) + stim_mod_const(c);
    stim_mod_predrate(c,:) = log(1+exp(cur_genfun));

    %%    
    
      
end

%%
close all
% orientations = linspace(0,pi-pi/8,8);
lambdas = [0.25 0.3 0.4 0.5 0.6 0.75 0.9 1.1]*Fsd;
% for c = 1:9;
c = 4;
fprintf('Cell %d of %d\n',c,9);
un_spk_cnts = length(unique(spikes_binned(tr_inds,c))) - 1;
spkbs = [];
for i = 1:un_spk_cnts
    curset = find(spikes_binned(tr_inds,c) == i);
    spkbs = [spkbs; repmat(curset,i,1)];
end
spkbs = sort(spkbs);
Robs = spikes_binned(tr_inds,c);
Robs_xv = spikes_binned(xv_inds,c);

cur_params = gabor_params_fin(c,:);
for ii = 1:length(lambdas)
    fprintf('Orientation %d of %d\n',ii,length(orientations));
    cur_params(4) = lambdas(ii);
    cur_mask1 = get_pgabor_mask(cur_params(1:6),0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(cur_params(1:6),pi/2,[SDIM SDIM]);
    mask1_out = X_resh*cur_mask1(:);
    mask2_out = X_resh*cur_mask2(:);
%     lin_out = cur_params(8)*mask1_out(all_model_fixids) + cur_params(9)*mask2_out(all_model_fixids);
    energy_out = cur_params(7)*sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
%     total_out = lin_out + energy_out;
    total_out = energy_out;
    total_out = zscore(total_out);
    lamrange2 = [slope2_lambda 1 length(tax);
        slope2_lambda length(tax)+1 length(tax)*2];
    nlxrange = [tax tax];
    Xmat = Tmat;
    Xmat = [Xmat bsxfun(@times,Tmat,total_out)];
    llist = [1 (length(tax)+1):size(Xmat,2)];
    [fitp,grad] = GLMsolve_jmm( Xmat(tr_inds,:), spkbs, W0, 1, lamrange, lamrange2,Pcon, Pmono, llist, [], NLtype,1e-5,nlxrange);
    cur_mod_const(ii) = fitp.k(end);
    cur_mod_sachist(ii,:) = fitp.k(1:length(tax));
    cur_mod_tot(ii,:) = fitp.k(length(tax)+1:length(tax)*2);   
end

prate = log(1+exp(repmat(cur_mod_const',1,length(tax)) + cur_mod_sachist + cur_mod_tot));
prate_null = log(1+exp(repmat(cur_mod_const',1,length(tax)) + cur_mod_sachist));
prate_diff = prate-prate_null;
plot(tax,prate_diff)
% end
%%

