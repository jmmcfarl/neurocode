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
load fixation_image_patches_corrected_raw

NT = length(used_fixs);
X_resh = reshape(X,NT,SDIM^2);
dt = .005;

cellids = [1 2 3 4 5 6 7 8 10];
muaids = [1 2 4 5 6 11 12 13 14];

spk_cnts = [spk_cnts(:,cellids) mua_cnts(:,muaids)];

n_used_cells = size(spk_cnts,2);

%%
all_fix_durs = all_fix_stop_times-all_fix_start_times;
min_fix_dur = 0.3;
used_fix_inds = find(all_fix_durs(used_fixs) > min_fix_dur);
X_resh = X_resh(used_fix_inds,:);
used_fixs = used_fixs(used_fix_inds);
%%
spikes_binned_early = [];
spikes_binned_late = [];
for blockid = 1:4
    blockid
    cur_set = find(blockids(used_fixs) == blockid);
    cur_spks_binned_early = nan(n_used_cells,length(cur_set));
    cur_spks_binned_late = nan(n_used_cells,length(cur_set));
    for i = 1:length(cur_set)
        start_T_early = all_fix_start_times(used_fixs(cur_set(i)));
        end_T_early = all_fix_start_times(used_fixs(cur_set(i))) + 0.15;
        start_T_late = all_fix_start_times(used_fixs(cur_set(i))) + 0.15;
        end_T_late = all_fix_start_times(used_fixs(cur_set(i))) + 0.3;
        
        for c = 1:n_used_cells
            if c <= 9
                cur_spk_times = find(Blocks{blockid}.spktimes{cellids(c)} > start_T_early & ...
                    Blocks{blockid}.spktimes{cellids(c)} < end_T_early);
                cur_spks_binned_early(c,i) = length(cur_spk_times);
                cur_spk_times = find(Blocks{blockid}.spktimes{cellids(c)} > start_T_late & ...
                    Blocks{blockid}.spktimes{cellids(c)} < end_T_late);
                cur_spks_binned_late(c,i) = length(cur_spk_times);
            else
                cur_spk_times = find(Blocks{blockid}.mutimes{muaids(c-9)} > start_T_early & ...
                    Blocks{blockid}.mutimes{muaids(c-9)} < end_T_early);
                cur_spks_binned_early(c,i) = length(cur_spk_times);
                cur_spk_times = find(Blocks{blockid}.mutimes{muaids(c-9)} > start_T_late & ...
                    Blocks{blockid}.mutimes{muaids(c-9)} < end_T_late);
                cur_spks_binned_late(c,i) = length(cur_spk_times);
            end
        end
        
    end
    spikes_binned_early = [spikes_binned_early; cur_spks_binned_early'];
    spikes_binned_late = [spikes_binned_late; cur_spks_binned_late'];
    
end

%%
stim_params.spatial_dims = 2;
stim_params.flen = 1;
sparse_lambda = 30;
smooth_lambda = 2500;

clear beta_early beta_late

spatial_scale = 6;
uwin = 15;
xax = -floor(SDIM/2):floor(SDIM/2);
yax = -floor(SDIM/2):floor(SDIM/2);
xax(SDIM+1:end) = [];
yax(SDIM+1:end) = [];
[XX,YY] = meshgrid(xax,yax);
for c = 1:n_used_cells
    c
    cur_params = gabor_params_fin(c,1:6);
    cur_params(5) = spatial_scale;
    uset = find(abs(XX-round(cur_params(1))) <= uwin & abs(YY-round(cur_params(2))) <= uwin);
    
    gauss_mask = get_pgauss_mask(cur_params,[SDIM SDIM]);
    gauss_mask = gauss_mask(uset);
    X_out = bsxfun(@times,X_resh(:,uset),gauss_mask');
    n_dim = length(uset);
    per_dim = sqrt(n_dim);
    X_out = reshape(X_out,[size(X_out,1) per_dim per_dim]);
    fft_out = nan(size(X_out));
    for i = 1:size(X_out,1)
        fft_out(i,:,:) = abs(fftshift(fft2(squeeze(X_out(i,:,:)))));
    end
    pred_mat = reshape(fft_out,size(X_out,1),per_dim^2);
    pred_mat = zscore(pred_mat);
    
    spikebins_early = convert_to_spikebins(spikes_binned_early(:,c));
    spikebins_late = convert_to_spikebins(spikes_binned_late(:,c));
    
    stim_params.sdim = per_dim;
    init_kerns = zeros(n_dim,1);
    init_signs = 1;
    kern_types{1} = 'lin';
    defmod.lambda_L1x = sparse_lambda;
    defmod.lambda_d2X = smooth_lambda;
    gnm0 = createGNM(init_kerns,init_signs,kern_types,defmod,stim_params);
    [gnm] = fitGNM_filters(gnm0,pred_mat,spikebins_early,'none',[],[],[]);
    beta_early(c,:) = gnm.mods(1).k;
    [gnm] = fitGNM_filters(gnm0,pred_mat,spikebins_late,'none',[],[],[]);
    beta_late(c,:) = gnm.mods(1).k;
        
end
%%
 for c = 1:n_used_cells
subplot(2,2,1)
imagesc(reshape(beta_early(c,:),per_dim,per_dim))
colorbar
cax = caxis();
cm1 = max(abs(cax));
caxis([-cm1 cm1]);
subplot(2,2,3)
imagesc(reshape(beta_late(c,:),per_dim,per_dim))
cax = caxis();
cm = max(abs(cax));
caxis([-cm cm]);
colorbar
subplot(2,2,2)
imagesc(reshape(beta_early(c,:),per_dim,per_dim))
colorbar
caxis([-cm1 cm1]);
subplot(2,2,4)
imagesc(reshape(beta_late(c,:),per_dim,per_dim))
caxis([-cm1 cm1]);
colorbar
c
pause
clf
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
% tr_inds = find(~ismember(all_model_fixids,xv_fixs));
tr_inds = 1:length(all_model_fixids);

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
for c = 1:n_used_cells
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
    
    %     %initialize parameters
    silent = 1;
    lamrange = [];
    lamrange2 = [slope2_lambda 1 length(tax)];
    nlxrange = [tax];
    Pcon = [];
    Pmono = [];
    hold_const = [];
    NLtype = 0;
    llist = [];
    
    %%
    disp('Fitting total stim model')
    X_resh = reshape(X,size(X,1),SDIM^2);
    cur_mask1 = get_pgabor_mask(gabor_params_fin(c,1:6),0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(gabor_params_fin(c,1:6),pi/2,[SDIM SDIM]);
    cur_gmask = get_pgauss_mask(gabor_params_fin(c,1:6),[SDIM SDIM]);
    cur_gmask = cur_gmask(:);
    mask1_out = X_resh*cur_mask1(:);
    mask2_out = X_resh*cur_mask2(:);
    gmask_out = bsxfun(@times,X_resh,cur_gmask');
    loc_mean = mean(gmask_out,2);
    loc_cont = std(gmask_out,[],2);
    
    mask1_out = (mask1_out - loc_mean)./loc_cont;
    mask2_out = (mask2_out - loc_mean)./loc_cont;
    
    loc_mean = zscore(loc_mean);
    loc_cont = zscore(loc_cont);
    
    lin_out = gabor_params_fin(c,8)*mask1_out(all_model_fixids) + gabor_params_fin(c,9)*mask2_out(all_model_fixids);
    energy_out = gabor_params_fin(c,7)*sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
    total_out = lin_out + energy_out;
    total_out = zscore(total_out);
    
    lamrange2 = [slope2_lambda 1 length(tax);
        slope2_lambda length(tax)+1 length(tax)*2];
    nlxrange = [tax tax];
    Xmat = Tmat;
    Xmat = [Xmat bsxfun(@times,Tmat,total_out)];
        %     llist = [l1_pen (length(tax)+1):size(Xmat,2)];
    [fitp,grad] = GLMsolve_jmm( Xmat(tr_inds,:), spkbs, W0, 1, lamrange, lamrange2,Pcon, Pmono, llist, [], NLtype,1e-5,nlxrange);
    stim_mod_const(c) = fitp.k(end);
    stim_mod_sachist(c,:) = fitp.k(1:length(tax));
    stim_mod_tot(c,:) = fitp.k(length(tax)+1:length(tax)*2);
end    
