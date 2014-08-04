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
    
%     Trial_angles = mod(Trial_angles + 270,360);
    
    jump_sizes = [Expts{Expt_nu(ee)}.Trials(:).sz];
    
    
    for i = 1:length(completed_trials)
        cur_images = repmat(Trial_im_nums(completed_trials(i)),ims_per_trial,1);
        
        cur_xo = (0:(ims_per_trial-1))*jump_sizes(completed_trials(i))*cos(degtorad(Trial_angles(completed_trials(i))));
        cur_yo = (0:(ims_per_trial-1))*jump_sizes(completed_trials(i))*sin(degtorad(Trial_angles(completed_trials(i))));
        
        cur_t_starts = Trial_starts(completed_trials(i)) + (0:3)*dt_per_jump;
        cur_t_stops = cur_t_starts + spk_win;
        cur_t_edges = sort([cur_t_starts cur_t_stops]);
        
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

xv_frac = 0.2;
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
all_gabor_outs = all_gabor_outs1.^2 + all_gabor_outs2.^2;
all_gabor_outs = sqrt(all_gabor_outs)';

for t = 1:96
    cur_binned_spks = all_binned_spks(:,t);
    
    [beta(t,:),dev(t),stats(t)] = glmfit(squeeze(all_gabor_outs(tr_samp,t)),cur_binned_spks(tr_samp),'poisson');
    p_vals(t,:) = stats.p;
    pred_r = squeeze(all_gabor_outs(tr_samp,t))*beta(t,2) + beta(t,1);
    pred_r = exp(pred_r);
    LL(t) = -sum(cur_binned_spks(tr_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(tr_samp));
    pred_r = squeeze(all_gabor_outs(xv_samp,t))*beta(t,2) + beta(t,1);
    pred_r = exp(pred_r);
    xvLL(t) = -sum(cur_binned_spks(xv_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(xv_samp));

    avg_rate(t) = mean(cur_binned_spks(tr_samp));
    null_LL(t) = -sum(bsxfun(@minus,cur_binned_spks(tr_samp)*log(avg_rate(t)),avg_rate(t)))/sum(cur_binned_spks(tr_samp));
    null_xvLL(t) = -sum(bsxfun(@minus,cur_binned_spks(xv_samp)*log(avg_rate(t)),avg_rate(t)))/sum(cur_binned_spks(xv_samp));
end
%%
% rf_shifts = linspace(-0.5,0.5,11);
% for ii = 1:length(rf_shifts)
%     fprintf('Shift X %d of %d\n',ii,length(rf_shifts));
%     for t = 1:96
%         init_params(1) = mean_x(t)+rf_shifts(ii); %x0
%         init_params(2) = mean_y(t); %y0
%         init_params(3) = degtorad(pref_mu_ori(t)); %theta
%         init_params(4) = 1/6; %lambda
%         init_params(5) = 0.3*init_params(4); %sigma
%         init_params(6) = 1; %eccentricity of gaussian ellipse
%         gabor_emp1 = get_pgabor_mask_v2(XX,YY,init_params,0);
%         gabor_emp2 = get_pgabor_mask_v2(XX,YY,init_params,pi/2);
%         gabor_emp1_filt_r(t,:) = gabor_emp1(:);
%         gabor_emp2_filt_r(t,:) = gabor_emp2(:);
%     end
%     
%     all_gabor_outs1 = gabor_emp1_filt_r*fullX';
%     all_gabor_outs2 = gabor_emp2_filt_r*fullX';
%     all_gabor_outs = all_gabor_outs1.^2 + all_gabor_outs2.^2;
%     all_gabor_outs = sqrt(all_gabor_outs)';
%     
%     for t = 1:96
%         cur_binned_spks = all_binned_spks(:,t);
%         
%         [beta_r(t,:),dev(t),stats(t)] = glmfit(squeeze(all_gabor_outs(tr_samp,t)),cur_binned_spks(tr_samp),'poisson');
%         pred_r = squeeze(all_gabor_outs(xv_samp,t))*beta_r(t,2) + beta_r(t,1);
%         pred_r = exp(pred_r);
%         xvLL_r(ii,t) = -sum(cur_binned_spks(xv_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(xv_samp));
%     end
% end
% xvLL_rel_x = bsxfun(@minus,xvLL_r,null_xvLL);
% for ii = 1:length(rf_shifts)
%     fprintf('Shift Y %d of %d\n',ii,length(rf_shifts));
%     for t = 1:96
%         init_params(1) = mean_x(t); %x0
%         init_params(2) = mean_y(t)+rf_shifts(ii); %y0
%         init_params(3) = degtorad(pref_mu_ori(t)); %theta
%         init_params(4) = 1/6; %lambda
%         init_params(5) = 0.3*init_params(4); %sigma
%         init_params(6) = 1; %eccentricity of gaussian ellipse
%         gabor_emp1 = get_pgabor_mask_v2(XX,YY,init_params,0);
%         gabor_emp2 = get_pgabor_mask_v2(XX,YY,init_params,pi/2);
%         gabor_emp1_filt_r(t,:) = gabor_emp1(:);
%         gabor_emp2_filt_r(t,:) = gabor_emp2(:);
%     end
%     
%     all_gabor_outs1 = gabor_emp1_filt_r*fullX';
%     all_gabor_outs2 = gabor_emp2_filt_r*fullX';
%     all_gabor_outs = all_gabor_outs1.^2 + all_gabor_outs2.^2;
%     all_gabor_outs = sqrt(all_gabor_outs)';
%     
%     for t = 1:96
%         cur_binned_spks = all_binned_spks(:,t);
%         
%         [beta_r(t,:),dev(t),stats(t)] = glmfit(squeeze(all_gabor_outs(tr_samp,t)),cur_binned_spks(tr_samp),'poisson');
%         pred_r = squeeze(all_gabor_outs(xv_samp,t))*beta_r(t,2) + beta_r(t,1);
%         pred_r = exp(pred_r);
%         xvLL_r(ii,t) = -sum(cur_binned_spks(xv_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(xv_samp));
%     end
% end
% xvLL_rel_y = bsxfun(@minus,xvLL_r,null_xvLL);
% 
%%
rf_shifts = linspace(-0.5,0.5,11);
for ii = 1:length(rf_shifts)
    for jj = 1:length(rf_shifts)
        fprintf('Shift X %d of %d\n',ii,length(rf_shifts));
        for t = 1:96
            init_params(1) = mean_x(t)+rf_shifts(ii); %x0
            init_params(2) = mean_y(t)+rf_shifts(jj); %y0
            init_params(3) = degtorad(pref_mu_ori(t)); %theta
            init_params(4) = 1/6; %lambda
            init_params(5) = 0.3*init_params(4); %sigma
            init_params(6) = 1; %eccentricity of gaussian ellipse
            gabor_emp1 = get_pgabor_mask_v2(XX,YY,init_params,0);
            gabor_emp2 = get_pgabor_mask_v2(XX,YY,init_params,pi/2);
            gabor_emp1_filt_r(t,:) = gabor_emp1(:);
            gabor_emp2_filt_r(t,:) = gabor_emp2(:);
        end
        
        all_gabor_outs1 = gabor_emp1_filt_r*fullX';
        all_gabor_outs2 = gabor_emp2_filt_r*fullX';
        all_gabor_outs = all_gabor_outs1.^2 + all_gabor_outs2.^2;
        all_gabor_outs = sqrt(all_gabor_outs)';
        
        for t = 1:96
            cur_binned_spks = all_binned_spks(:,t);
            
            [beta_r(t,:),dev(t),stats(t)] = glmfit(squeeze(all_gabor_outs(tr_samp,t)),cur_binned_spks(tr_samp),'poisson');
            pred_r = squeeze(all_gabor_outs(tr_samp,t))*beta_r(t,2) + beta_r(t,1);
            pred_r = exp(pred_r);
            LL_r(ii,jj,t) = -sum(cur_binned_spks(tr_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(tr_samp));
            pred_r = squeeze(all_gabor_outs(xv_samp,t))*beta_r(t,2) + beta_r(t,1);
            pred_r = exp(pred_r);
            xvLL_r(ii,jj,t) = -sum(cur_binned_spks(xv_samp).*log(pred_r) - pred_r)/sum(cur_binned_spks(xv_samp));
        end
    end
end
xvLL_rel = bsxfun(@minus,xvLL_r,shiftdim(null_xvLL,-1));
LL_rel = bsxfun(@minus,LL_r,shiftdim(null_LL,-1));
%%
sdim = sqrt(size(fullX,2));
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

for t = 1:96
    
    fprintf('Cell %d of %d\n',t,96);
    
    init_params(1) = mean_x(t); %x0
    init_params(2) = mean_y(t); %y0
    init_params(3) = degtorad(pref_mu_ori(t)); %theta
    init_params(4) = 1/6; %lambda
    init_params(5) = 0.3*init_params(4); %sigma
    init_params(6) = 1; %eccentricity of gaussian ellipse
    init_params(7) = 0; %weight of quad term
    init_params(8) = 0; %weight of lin term 0 phase
    init_params(9) = 0; %weight of lin term pi/2 phase
    init_params(10) = 0; %const offset
    hold_const = [0 0 1 1 1 1 0 1 1 0];
    LB = [0.1 -0.7 0 0.125 0.02 0.5 -Inf -Inf -Inf -Inf];
    UB = [0.7 -0.1 pi 0.25 0.25 2 Inf Inf Inf Inf];
    
    [gabor_params_emp(t,:),LL_emp(t)] = fit_gabor_params_v3(XX,YY,init_params,fullX(tr_samp,:),all_binned_spks(tr_samp,t),hold_const,LB,UB);
    gabor_bank_emp(t,:,:) = get_pgabor_mask_v2(XX,YY,gabor_params_emp(t,1:6),0);
    xv_LL_emp(t) = get_pgabor_LL_v2(gabor_params_emp(t,:),fullX(xv_samp,:),all_binned_spks(xv_samp,t),XX,YY);    
    
    hold_const = [1 1 1 1 1 1 0 1 1 0];
    [gabor_params_orig(t,:),LL_orig(t)] = fit_gabor_params_v3(XX,YY,init_params,fullX(tr_samp,:),all_binned_spks(tr_samp,t),hold_const,LB,UB);
    xv_LL_orig(t) = get_pgabor_LL_v2(gabor_params_orig(t,:),fullX(xv_samp,:),all_binned_spks(xv_samp,t),XX,YY);
   
    avg_rate(t) = mean(all_binned_spks(tr_samp,t));
%     null_xvLL(t) = -sum(bsxfun(@minus,all_binned_spks(xv_samp,t)*log(avg_rate(t)),avg_rate(t)))/sum(all_binned_spks(xv_samp,t));
    
end


