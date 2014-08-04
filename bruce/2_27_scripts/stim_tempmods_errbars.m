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
fix_stimout = zeros(n_used_cells,length(used_fixs));
all_model_time_axis = [];
all_model_blockids = [];
all_model_fixids = [];
spikes_binned = [];
time_since_fix = [];
rel_spk_times = cell(18,1);
spk_fix_inds = cell(18,1);
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
                cur_spk_times = find(Blocks{blockid}.spktimes{cellids(c)} > start_T & ...
                    Blocks{blockid}.spktimes{cellids(c)} < end_T);
                temp = hist(Blocks{blockid}.spktimes{cellids(c)}(cur_spk_times),cur_tcents);
                rel_spk_times{c} = [rel_spk_times{c} Blocks{blockid}.spktimes{cellids(c)}(cur_spk_times)-start_T];
                spk_fix_inds{c} = [spk_fix_inds{c} cur_set(i)*ones(size(cur_spk_times))];
            else
                cur_spk_times = find(Blocks{blockid}.mutimes{muaids(c-9)} > start_T & ...
                    Blocks{blockid}.mutimes{muaids(c-9)} < end_T);
                temp = hist(Blocks{blockid}.mutimes{muaids(c-9)}(cur_spk_times),cur_tcents);
                rel_spk_times{c} = [rel_spk_times{c} Blocks{blockid}.mutimes{muaids(c-9)}(cur_spk_times)-start_T];
                spk_fix_inds{c} = [spk_fix_inds{c} cur_set(i)*ones(size(cur_spk_times))];
            end
            temp_binned(c,:) = temp;
        end
        spikes_binned = [spikes_binned; temp_binned'];
        
        time_since_fix = [time_since_fix cur_tcents-start_T];
        
    end
end

%% Compute TBR time-since fix onset
max_tsf = 0.5; nbins = 40;
used_tsf = time_since_fix;
used_tsf(used_tsf > max_tsf) = [];
% tax = prctile(used_tsf,linspace(100/nbins,100-100/nbins,nbins));
tax = linspace(0,max_tsf,nbins);

% tax = logspace(log10(dt),log10(max_tsf),nbins);
Tmat = tbrep(time_since_fix,tax);

%%
nfold = 10;
rperm = randperm(length(used_fixs));
fixs_per_fold = floor(length(used_fixs)/nfold);
for i = 1:nfold
    cur_fix_set = rperm((i-1)*fixs_per_fold + (1:fixs_per_fold));
    fold_inds{i} = find(ismember(all_model_fixids,cur_fix_set));
end
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
slope2_lambda = 100;
l1_pen = 0;
for c = 1:n_used_cells
    fprintf('Cell %d of %d\n',c,9);
    
        
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

    Xmat = [Tmat bsxfun(@times,Tmat,total_out)];    
    W0 = zeros(size(Xmat,2),1);
    for n = 1:nfold

    Robs = spikes_binned(fold_inds{n},c);
    spkbs = convert_to_spikebins(Robs);
        
    %     %initialize parameters
    
    %%
    silent = 1;
    NLtype = 1;
    Pcon = [];
    Pmono = [];
    llist = [];
    lamrange2 = [slope2_lambda 1 length(tax) 0;
        slope2_lambda length(tax)+1 length(tax)*2 0];
    lamrange = [];
    nlxrange = [tax tax];
    [fitp,grad] = GLMsolve_jmm( Xmat(fold_inds{n},:), spkbs, W0, 1, lamrange, lamrange2,Pcon, Pmono, llist, [], NLtype,1e-5,nlxrange);
    stim_mod_const(c,n) = fitp.k(end);
    stim_mod_sachist(c,n,:) = fitp.k(1:length(tax));
    stim_mod_tot(c,n,:) = fitp.k(length(tax)+1:length(tax)*2);
    end
    
    
end

%%
close all

for c = 1:n_used_cells
    subplot(2,1,1)
    shadedErrorBar(tax,squeeze(mean(stim_mod_sachist(c,:,:),2)),squeeze(std(stim_mod_sachist(c,:,:),[],2))/sqrt(nfold))
    subplot(2,1,2)
    shadedErrorBar(tax,squeeze(mean(stim_mod_tot(c,:,:),2)),squeeze(std(stim_mod_tot(c,:,:),[],2))/sqrt(nfold))
    pause
    clf
end
