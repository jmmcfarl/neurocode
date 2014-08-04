%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/G075/

load ./CellList.mat
% single_units = find(CellList(1,:,1) > 0);

load ./Expt3_fixbased_data.mat
norm_fac = std(resh_all_stims(:));
resh_all_stims = resh_all_stims/norm_fac;
X = resh_all_stims(trial_stimnum,:);

fix_win_dur = 0.15;

Pix2Deg = 0.018837;
NT = length(full_stim_ids);
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

%%
n_trials = length(trial_start_inds);
fix_binned_spks = nan(n_trials,96);
fix_expt_num = nan(n_trials,1);
for i = 1:n_trials
    cur_inds = trial_start_inds(i):trial_stop_winds(i);
    fix_binned_spks(i,:) = sum(full_binned_spks(cur_inds,:));
    
    fix_expt_num(i) = unique(full_expt_vec(cur_inds));
end
% [un_expts,~,full_expt_inds] = unique(fix_expt_num);
% n_un_expts = length(un_expts);
% linX = zeros(n_trials,n_un_expts-1);
% for i = 1:n_un_expts-1
%     linX(fix_expt_num==i,i) = 1;
% end
%%
load ./NS40_gabor_mods.mat gabor*
clear gabor_*filt
for t = 1:96
    
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params_f(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params_f(t,1:6),pi/2);
    
    gabor_emp1_filt(t,:) = gabor_emp1(:);
    gabor_emp2_filt(t,:) = gabor_emp2(:);
end

%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA
all_gabor_out1 = X*gabor_emp1_filt';
all_gabor_out2 = X*gabor_emp2_filt';
spatial_mod_out = sqrt(all_gabor_out1.^2 + all_gabor_out2.^2);

%%
for t = 1:96
    [B,DEV,STATS] = glmfit(spatial_mod_out(:,t),fix_binned_spks(:,t),'poisson');

    weight(t) = B(2);
    weight_se(t) = STATS.se(2);

end

%%
used_trial_set = unique(full_trial_vec);
n_trials = length(used_trial_set);
n_fold = 5;
set_frac =1/n_fold;
trials_per_fold = floor(set_frac*n_trials);

fold_tord = randperm(n_trials);

for i = 1:n_fold
    fold_set = fold_tord((i-1)*trials_per_fold + (1:trials_per_fold));
    fold_inds{i} = find(ismember(full_trial_vec,fold_set));
end

%%

flen = 86;
flen_t = flen*dt;
NT = length(full_expt_vec);

tent_centers = [0:dt:0.47];
tent_centers = round(tent_centers/dt);
if tent_centers(end) >= flen
    tent_centers(end) = [];
end
tbmat = construct_tent_bases(tent_centers,1);
[ntents,tblen] = size(tbmat);

tbmat = [zeros(ntents,tblen-1) tbmat];

trial_inds = zeros(size(full_t));
trial_inds(trial_start_inds) = 1;
trial_Tmat = zeros(length(full_t),ntents);
for i = 1:ntents
    trial_Tmat(:,i) = conv(trial_inds,tbmat(i,:),'same');
end

% [un_expts,~,full_expt_inds] = unique(full_expt_vec);
% n_un_expts = length(un_expts);
% linX = zeros(NT,n_un_expts-1);
% for i = 1:n_un_expts-1
%     linX(full_expt_inds==i,i) = 1;
% end

%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA
NL_type = 1;
all_gabor_out1 = resh_all_stims*gabor_emp1_filt';
all_gabor_out2 = resh_all_stims*gabor_emp2_filt';
spatial_mod_out = sqrt(all_gabor_out1.^2 + all_gabor_out2.^2);
spatial_mean = mean(spatial_mod_out);
spatial_std = std(spatial_mod_out);
spatial_mod_out = bsxfun(@minus,spatial_mod_out,spatial_mean);
spatial_mod_out = bsxfun(@rdivide,spatial_mod_out,spatial_std);
% spatial_mod_out   = zscore(spatial_mod_out);

l2_ind = 400;
l2_dep = 800;
silent = 1;
for t = 1:96
    fprintf('Fitting Cell %d of 96\n',t);
    cur_spatial_mod_out = spatial_mod_out(full_stim_ids,t);
    stim_dep_X = bsxfun(@times,trial_Tmat,cur_spatial_mod_out);
    ov_Xmat1 = [stim_dep_X trial_Tmat];
        ov_Xmat2 = [trial_Tmat];

    for n = 1:n_fold
    Robs = full_binned_spks(fold_inds{n},t);
    spksN = convert_to_spikebins(Robs);
    
    klen = size(ov_Xmat1,2);
    K0 = zeros(klen,1);
    lamrange2 = [l2_dep 1 ntents 0;l2_ind (ntents+1) 2*ntents 0];
    lamrange = [];
    llist = [];
    [fitp,grad] = GLMsolve_jmm(ov_Xmat1(fold_inds{n},:), spksN, K0, silent, lamrange, lamrange2, [], [],llist, [], 1);
    stim_dep_kern(t,n,:) = fitp.k(1:ntents);
    stim_ind_kern(t,n,:) = fitp.k((ntents+1):2*ntents);
    %     block_kern(t,:) = fitp.k((2*ntents+1):end-1);
    ov_const(t,n) = fitp.k(end);
    
%     gfun = ov_Xmat*fitp.k(1:end-1) + fitp.k(end);
%     too_large = find(gfun > 50);
%     predrate = log(1+exp(gfun));
%     predrate(too_large) = gfun(too_large);
%     LL(t) = -sum(Robs.*log(predrate)-predrate)/sum(Robs);
    
    cur_lambda = l2_ind;
    lamrange2 = [cur_lambda 1 ntents 0];
    llist = [];
    [fitp_ns,grad] = GLMsolve_jmm(ov_Xmat2(fold_inds{n},:), spksN, K0, silent, [], lamrange2, [], [],llist, [], 1);
    trig_avg_kern(t,n,:) = fitp_ns.k(1:ntents);
    %     offset = fitp.k(end) + mean(fitp.k(ntents+1:end-1));
    offset = fitp_ns.k(end);
    trig_avg_rate(t,n,:) = log(1+exp(trig_avg_kern(t,n,:) + offset));
    
    
    end
end

%%
cd ~/Data/bruce/G075/
% save expt3_unit_tempmods_norep tent_centers dt stim_* ov_const LL trig_avg_rate spatial_mean spatial_std

%%
close all
for t = 1:96
    subplot(2,1,1)
    shadedErrorBar(tent_centers*dt,squeeze(mean(stim_dep_kern(t,:,:),2)),squeeze(std(stim_dep_kern(t,:,:),[],2))/sqrt(5))
    hold on
    plot(tent_centers*dt,squeeze(stim_dep_kern(t,:,:)),'r')
    xlim([0 0.5])
    subplot(2,1,2)
    shadedErrorBar(tent_centers*dt,squeeze(mean(stim_ind_kern(t,:,:),2)),squeeze(std(stim_ind_kern(t,:,:),[],2))/sqrt(5))
    hold on
    plot(tent_centers*dt,squeeze(stim_ind_kern(t,:,:)),'r')
    xlim([0 0.5])
%     subplot(3,1,3)
%     plot(tent_centers*dt,squeeze(trig_avg_rate(t,:,:))/dt)
%     xlim([0 0.5])
    t
    pause
    clf
end

%% CREATE PREDICTED RESPONSE FOR EACH REPEAT SEQ
cd ~/Data/bruce/Expt_1_8_13_imfolder/
n_images = 8;
repeat_inds = [1e3 2e3 3e3 2.01e5 2.02e5 2.03e5 4.01e5 4.02e5 4.03e5];

samp_fac = 0.5;
rep_t = (1:(320/samp_fac))*dt;
rep_flip_inds = 1:(40/samp_fac):(320/samp_fac);
rep_stim_ids = ceil((1:length(rep_t))/length(rep_t)*8);
for seq = 1:9
    all_im_patches = nan(n_images,length(ypatch_inds),length(xpatch_inds));   
    used_images = repeat_inds(seq) + (1:8);
    for i = 1:n_images
        if used_images(i) < 1e4
            filename = sprintf('IM100%.4d.png',used_images(i));
        elseif used_images(i) < 1e5
            filename = sprintf('IM10%.5d.png',used_images(i));
        else
            filename = sprintf('IM1%.6d.png',used_images(i));
        end
        IMAGEorg = imread(filename);
        IMAGEorg = double(IMAGEorg); % convert to double format
        IMAGEorg = IMAGEorg - 127.5; %shift range to [-127.5:127.5]
        IMAGE = flipud(IMAGEorg); %flip y
        
        IMAGE = imresize(IMAGE,1/dsfrac); %image downsampling
        
        cur_patch = IMAGE(ypatch_inds,xpatch_inds);
        all_im_patches(i,:,:) = cur_patch;
        
    end
    rep_stims{seq} = reshape(all_im_patches,n_images,sdim^2)/norm_fac;      
end

cd ~/Data/bruce/G075/

%create tent-basis stim time rep
rep_trial_inds = zeros(size(rep_t));
rep_trial_inds(rep_flip_inds) = 1;
rep_trial_Tmat = zeros(length(rep_t),ntents);
for i = 1:ntents
    rep_trial_Tmat(:,i) = conv(rep_trial_inds,tbmat(i,:),'same');
end

for seq = 1:9
    all_gabor_out1 = rep_stims{seq}*gabor_emp1_filt';
    all_gabor_out2 = rep_stims{seq}*gabor_emp2_filt';
    spatial_mod_out = sqrt(all_gabor_out1.^2 + all_gabor_out2.^2);
    spatial_mod_out = bsxfun(@minus,spatial_mod_out,spatial_mean);
    spatial_mod_out = bsxfun(@rdivide,spatial_mod_out,spatial_std);

    for t = 1:96
        cur_spatial_mod_out = spatial_mod_out(rep_stim_ids,t);
        stim_dep_X = bsxfun(@times,rep_trial_Tmat,cur_spatial_mod_out);
        ov_Xmat = [stim_dep_X rep_trial_Tmat];
    
        kern = [stim_dep_kern(t,:) stim_ind_kern(t,:)];
        mod_out = ov_Xmat*kern'+ov_const(t);
        rep_pred_r(seq,t,:) = log(1+exp(mod_out));
    end

end

%%
save expt3_repeat_preds rep_t rep_pred_r