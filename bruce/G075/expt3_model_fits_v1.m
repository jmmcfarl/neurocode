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
infix_ids = find(~isnan(full_stim_ids));
n_fixs = length(unique(full_stim_ids(infix_ids)));
xv_frac =0.2;
xv_num = round(xv_frac*n_fixs);
xv_f = randperm(n_fixs);
xv_f(xv_num+1:end) = [];
tr_f = setdiff(1:n_fixs,xv_f);
xv_inds = find(ismember(full_stim_ids,xv_f));
tr_inds = setdiff(find(~isnan(full_stim_ids)),xv_inds);

skip_set = [1000 2000 3000 401000 402000 403000];
rep_inds = find(ismember(full_seof_vec,skip_set));
% tr_inds(rep_inds) = [];
% xv_inds = rep_inds;
tr_inds(ismember(tr_inds,rep_inds)) = [];
xv_inds(ismember(xv_inds,rep_inds)) = [];

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

l2_ind = 1000;
l2_dep = 2000;
silent = 1;
for t = 1:96
    fprintf('Fitting Cell %d of 96\n',t);
    cur_spatial_mod_out = spatial_mod_out(full_stim_ids(tr_inds),t);
    stim_dep_X = bsxfun(@times,trial_Tmat(tr_inds,:),cur_spatial_mod_out);
    %     ov_Xmat = [stim_dep_X trial_Tmat(tr_inds,:) linX(tr_inds,:)];
    ov_Xmat = [stim_dep_X trial_Tmat(tr_inds,:)];
    
    Robs = full_binned_spks(tr_inds,t);
    spksN = convert_to_spikebins(Robs);
    
    klen = size(ov_Xmat,2);
    K0 = zeros(klen,1);
    lamrange2 = [l2_dep 1 ntents 0;l2_ind (ntents+1) 2*ntents 0];
    %     lamrange = [l2_dep/10 1 ntents 0;l2_ind/10 (ntents+1) 2*ntents 0];
    lamrange = [];
    llist = [];
    [fitp,grad] = GLMsolve_jmm(ov_Xmat, spksN, K0, silent, lamrange, lamrange2, [], [],llist, [], 1);
    stim_dep_kern(t,:) = fitp.k(1:ntents);
    stim_ind_kern(t,:) = fitp.k((ntents+1):2*ntents);
    %     block_kern(t,:) = fitp.k((2*ntents+1):end-1);
    ov_const(t) = fitp.k(end);
    
    gfun = ov_Xmat*fitp.k(1:end-1) + fitp.k(end);
    too_large = find(gfun > 50);
    predrate = log(1+exp(gfun));
    predrate(too_large) = gfun(too_large);
    LL(t) = -sum(Robs.*log(predrate)-predrate)/sum(Robs);
    
    cur_lambda = 200;
    ov_Xmat = [trial_Tmat(tr_inds,:)];
    lamrange2 = [cur_lambda 1 ntents 0];
    llist = [];
    [fitp_ns,grad] = GLMsolve_jmm(ov_Xmat, spksN, K0, silent, [], lamrange2, [], [],llist, [], 1);
    trig_avg_kern = fitp_ns.k(1:ntents);
    %     offset = fitp.k(end) + mean(fitp.k(ntents+1:end-1));
    offset = fitp_ns.k(end);
    trig_avg_rate(t,:) = log(1+exp(trig_avg_kern + offset));
    
    
    xv_Robs = full_binned_spks(xv_inds,t);
    cur_spatial_mod_out = spatial_mod_out(full_stim_ids(xv_inds),t);

    Tmat = [trial_Tmat(xv_inds,:)];
    xv_pred_rate = Tmat*fitp_ns.k(1:end-1) + fitp_ns.k(end);
    if NL_type == 0
        xv_pred_rate = log(1+exp(xv_pred_rate));
    else
        xv_pred_rate = exp(xv_pred_rate);
    end
    xv_to_ns_LL(t) = -sum(xv_Robs.*log(xv_pred_rate)-xv_pred_rate)/sum(xv_Robs);
    
    
    Tmat = [trial_Tmat(xv_inds,:) bsxfun(@times,trial_Tmat(xv_inds,:),cur_spatial_mod_out)];
    xv_pred_rate = Tmat*fitp.k(1:end-1) + fitp.k(end);
    if NL_type == 0
            xv_pred_rate = log(1+exp(xv_pred_rate));
        else
            xv_pred_rate = exp(xv_pred_rate);
        end
        xv_to_LL(t) = -sum(xv_Robs.*log(xv_pred_rate)-xv_pred_rate)/sum(xv_Robs);

                avg_rate = mean(Robs);
        null_pred = avg_rate*ones(size(xv_Robs));
        xv_null_LL(t) = -sum(xv_Robs.*log(null_pred)-null_pred)/sum(xv_Robs);
        

end

%%
cd ~/Data/bruce/G075/
save expt3_unit_tempmods_norep tent_centers dt stim_* ov_const LL trig_avg_rate spatial_mean spatial_std

%%
for t = 1:96
    subplot(3,1,1)
    plot(tent_centers*dt,stim_dep_kern(t,:),'r')
    xlim([0 0.5])
    subplot(3,1,2)
    plot(tent_centers*dt,stim_ind_kern(t,:),'r')
    xlim([0 0.5])
    subplot(3,1,3)
    plot(tent_centers*dt,trig_avg_rate(t,:)/dt,'r')
    xlim([0 0.5])
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