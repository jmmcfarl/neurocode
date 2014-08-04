%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
cd ~/Data/bruce/7_15_12
cd G029/
load ./eye_calibration_data
% load ./G029Expts.mat
cd ~/Data/bruce/7_15_12/G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt1_newcompiled_data_fixeddelay_d1p25_34_v2.mat
fullX = fullX/std(fullX(:));

Pix2Deg = 0.018837;
[NT,klen] = size(fullX);


%% crop stimulus for the purpose of faster gabor function fitting
new_RF_patch = [-0.11 0.9; -0.9 0.1]; %location of RFs in degrees [x1 x2;y1 y2]
[curXX,curYY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
curXX = curXX(:);
curYY = curYY(:);

xpatch_inds_cropped = find(xax >= new_RF_patch(1,1) & xax <= new_RF_patch(1,2));
ypatch_inds_cropped = find(yax >= new_RF_patch(2,1) & yax <= new_RF_patch(2,2));

new_crop = find(curXX >= new_RF_patch(1,1) & curXX <= new_RF_patch(1,2) & ...
    curYY >= new_RF_patch(2,1) & curYY <= new_RF_patch(2,2));

sdim = length(xpatch_inds);
[XXc,YYc] = meshgrid(xax(xpatch_inds_cropped),yax(ypatch_inds_cropped));
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

cd ~/Data/bruce/7_15_12/G034/
load ./grating_mu_data
gr_oris = unique_oris;
[~,sf_inds] = max(avg_mu_sf_profile,[],2);
pref_sfs = unique_sfs(sf_inds);
cd ~/Data/bruce/7_15_12/G029/
load ./oned_fixation_fits_v3.mat

dt = dt*2;

% PARSE DATA INTO FIXATIONS
diff_used_inds = [1; diff(used_inds)];
rel_fix_start_inds = [1; find(diff_used_inds > 1)];
rel_fix_stop_inds = [(find(diff_used_inds > 1)-1); NT];
n_fixs = length(rel_fix_start_inds);

%%
cd ~/Data/bruce/7_15_12/G034/
load ./gabor_tracking_varmeans

load ./expt1_lfp_alpha_phase_delay_all_v2.mat new* all_interp_phases wfreqs nearest_lfps
% load ./expt1_lfp_alpha_phase_delay_all.mat

%%
fixed_delay = 0.05;
cur_dt = median(diff(full_t));
desired_dt = cur_dt/8;
trial_start_inds = 1+find(diff(new_trial_vec') ~= 0);
trial_stop_inds = trial_start_inds - 1;
trial_start_inds = [1; trial_start_inds];
trial_stop_inds = [trial_stop_inds; length(new_expt_vec)];

% trial_start_inds = trial_start_inds + round(0.1/desired_dt);

n_trials = length(trial_stop_inds);
new_binned_spks = nan(length(new_t_axis),96);
Expt_nu = unique(new_expt_vec);
n_expts = length(Expt_nu);
new_trial_vec = nan(length(new_t_axis),1);
for tt = 1:n_expts
    cur_trial_set = find(new_expt_vec(trial_start_inds) == Expt_nu(tt));
    cur_n_trials = length(cur_trial_set);
    load(sprintf('Expt%dClusterTimes.mat',Expt_nu(tt)));
    for i = 1:cur_n_trials
        cur_inds = trial_start_inds(cur_trial_set(i)):trial_stop_inds(cur_trial_set(i));
        cur_binned = zeros(length(cur_inds),96);
        for c = 1:96
            %             temp = histc(Clusters{c}.times-fixed_delay,new_t_axis(cur_inds));
            temp = histc(Clusters{c}.times,new_t_axis(cur_inds));
            cur_binned(:,c) = temp;
        end
        new_binned_spks(cur_inds,:) = cur_binned;
        new_trial_vec(cur_inds) = cur_trial_set(i);
    end
end
bad_pts = max(isnan(new_binned_spks),[],2);
% uset = find(new_expt_vec ~= 18 & bad_pts' == 0);
uset = find(new_expt_vec ~= 18 & bad_pts' == 0);

%% RECONSTRUCT NEW STIMULUS MATRIX
resh_X = reshape(fullX',[sdim sdim NT]);
resh_X_sh = zeros(size(resh_X));
for ii = 1:NT
    %     if mod(ii,100)==0 fprintf('%d of %d\n',ii,NT); end
    d2 = dist_shift2d(resh_X(:,:,ii), -x_cor{end}(ii), 2,0);
    d2 = dist_shift2d(d2,-y_cor{end}(ii),1,0);
    resh_X_sh(:,:,ii) = d2;
end

%% REFIT GEM PARAMETERS
fullX_sh = reshape(resh_X_sh,sdim^2,NT)';
fullX_cropped = fullX_sh(:,new_crop);

%%
used_trials = unique(new_trial_vec(uset));
n_used_trials = length(used_trials);
xv_frac = 0.2;
n_xv_trials = round(n_used_trials*xv_frac);
xv_set = randperm(n_used_trials);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(new_trial_vec(uset),xv_set));
tr_inds = find(~ismember(uset,xv_inds))';


%%
close all

reg_params.dl1_ind = 1000;
% reg_params.dl1_dep =20000;
reg_params.dl1_dep =5000;
reg_params.dl2_ind = 2000;
% reg_params.dl2_dep = 40000;
reg_params.dl2_dep = 10000;
reg_params.dl2_freq_ind = 500;
reg_params.dl2_freq_dep = 2000;
reg_params.l1_dep = 0;
reg_params.l1_ind = 0;
reg_params.is_phase = 1;

reg_params2 = reg_params;
reg_params2.dl1_ind = reg_params2.dl1_dep;
reg_params2.dl2_ind = reg_params2.dl2_dep;
reg_params2.dl2_freq_ind = reg_params2.dl2_freq_dep;

silent = 1;
nbins = 30;
pax = linspace(-pi,pi,nbins+1);

stim_params = [nbins,length(wfreqs)];
doub_pax = [pax(1:end-1) 2*pi+pax(1:end-1)];

for t = 59:96
    fprintf('Fitting GEM: Cell %d of %d\n',t,96);
    
    Robs = new_binned_spks(uset(tr_inds),t);
    spkbs = convert_to_spikebins(Robs);
    
    gabor_emp1 = get_pgabor_mask_v2(XXc,YYc,gabor_params_f{end}(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XXc,YYc,gabor_params_f{end}(t,1:6),pi/2);
    mask1_out = fullX_cropped*gabor_emp1(:);
    mask2_out = fullX_cropped*gabor_emp2(:);
    tot_mod_out =sqrt(mask1_out.^2 + mask2_out.^2);
    tot_mod_out = tot_mod_out/std(tot_mod_out);
    mod_out = interp1(full_t,tot_mod_out,new_t_axis(uset(tr_inds))-fixed_delay)';
    mod_out(isnan(mod_out)) = nanmean(mod_out);
    
    cur_alpha_phase = squeeze(all_interp_phases(uset(tr_inds),:,nearest_lfps(t)));
    
    klen = 1;
    K0 = 0;
    [fitp_so,grad] = GLMsolve_jmm(mod_out,spkbs,K0,silent,[],[],[],[],[],[],0);
    
    Tmat = [];
    for ww = 1:length(wfreqs)
        cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        Tmat = [Tmat cur_tb];
    end
    
    Xmat =  [Tmat bsxfun(@times,Tmat,mod_out)];
    %     Xmat = [Xmat mod_out];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp] = fit_GLM_phase_model(Xmat,Robs, K0, 1, stim_params, reg_params);
    ov_const(t) = fitp.k(end);
    stim_ind_filt(t,:) = fitp.k(1:nbins*length(wfreqs));
    stim_dep_filt(t,:) = fitp.k((nbins*length(wfreqs)+1):2*nbins*length(wfreqs));
    
    Xmat =  [Tmat];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_po] = fit_GLM_phase_model(Xmat, Robs, K0, 1, stim_params, reg_params);
    ov_const_po(t) = fitp_po.k(end);
    stim_ind_filt_po(t,:) = fitp_po.k(1:nbins*length(wfreqs));
    
    Xmat =  [Tmat mod_out];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_pos] = fit_GLM_phase_model(Xmat, Robs, K0, 1, stim_params, reg_params);
    stim_ind_filt_pos(t,:) = fitp_pos.k(1:nbins*length(wfreqs));
    
%     Xmat =  [bsxfun(@times,Tmat,mod_out) mod_out];
    Xmat =  [bsxfun(@times,Tmat,mod_out)];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_pg] = fit_GLM_phase_model(Xmat, Robs, K0, 1, stim_params, reg_params2);
    stim_dep_filt_pg(t,:) = fitp_pg.k(1:nbins*length(wfreqs));
    
    
    %
    xv_Robs = new_binned_spks(uset(xv_inds),t);
    mod_out = interp1(full_t,tot_mod_out,new_t_axis(uset(xv_inds))-fixed_delay)';
    
    cur_alpha_phase = squeeze(all_interp_phases(uset(xv_inds),:,nearest_lfps(t)));
    
    Tmat = [];
    for ww = 1:length(wfreqs)
        cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        Tmat = [Tmat cur_tb];
    end
    Xmat =  [Tmat bsxfun(@times,Tmat,mod_out)];
    
    xv_pred_rate = Xmat*fitp.k(1:end-1) + fitp.k(end);
    xv_pred_rate = log(1+exp(xv_pred_rate));
    xv_LL(t) = -sum(xv_Robs.*log(xv_pred_rate) - xv_pred_rate)/sum(xv_Robs);
    
    Xmat =  [Tmat];
    xv_pred_po = Xmat*fitp_po.k(1:end-1) + fitp_po.k(end);
    xv_pred_po = log(1+exp(xv_pred_po));
    xv_po_LL(t) = -sum(xv_Robs.*log(xv_pred_po) - xv_pred_po)/sum(xv_Robs);
    
    Xmat =  [Tmat mod_out];
    xv_pred_pos = Xmat*fitp_pos.k(1:end-1) + fitp_pos.k(end);
    xv_pred_pos = log(1+exp(xv_pred_pos));
    xv_pos_LL(t) = -sum(xv_Robs.*log(xv_pred_pos) - xv_pred_pos)/sum(xv_Robs);
    
    Xmat =  [bsxfun(@times,Tmat,mod_out)];
%     Xmat =  [bsxfun(@times,Tmat,mod_out) mod_out];
    xv_pred_pg = Xmat*fitp_pg.k(1:end-1) + fitp_pg.k(end);
    xv_pred_pg = log(1+exp(xv_pred_pg));
    xv_pg_LL(t) = -sum(xv_Robs.*log(xv_pred_pg) - xv_pred_pg)/sum(xv_Robs);
    
    xv_pred_so = mod_out*fitp_so.k(1) + fitp_so.k(end);
    xv_pred_so = log(1+exp(xv_pred_so));
    xv_so_LL(t) = -sum(xv_Robs.*log(xv_pred_so) - xv_pred_so)/sum(xv_Robs);
    
    avg_rate = mean(Robs);
    null_pred = avg_rate*ones(size(xv_Robs));
    xv_null_LL(t) = -sum(xv_Robs.*log(null_pred)-null_pred)/sum(xv_Robs);
    
    % input('')
    % clf
end

%%
save extp1_stim_alpha_models xv_* ind_kern dep_kern wfreqs pax stim*


%%
% norm_dep_kern = zeros(size(dep_kern));
% norm_ind_kern = zeros(size(ind_kern));
clear norm_ind_kern
for t = 1:70
    cur = reshape(stim_ind_filt_pos(t,:),nbins,length(wfreqs))';
    norm_ind_kern(t,:,:) = cur/std(cur(:));
%     norm_ind_kern(t,:,:) = cur;
%     cur = squeeze(ind_kern(i,:,:));
%     norm_ind_kern(i,:,:) = cur/std(cur(:));
end
norm_ind_dom = squeeze(max(norm_ind_kern,[],3)-min(norm_ind_kern,[],3));
% norm_dep_dom = squeeze(max(norm_dep_kern,[],3)-min(norm_dep_kern,[],3));

%%
% for t = 1:96
%     subplot(2,1,1)
%     pcolor(doub_pax,wfreqs,[squeeze(ind_kern(t,:,:)) squeeze(ind_kern(t,:,:))]);shading interp; set(gca,'yscale','log')
%     xlabel('Phase (rad)','fontsize',16)
%     ylabel('Frequency (Hz)','fontsize',16)
%     subplot(2,1,2)
%     pcolor(doub_pax,wfreqs,[squeeze(dep_kern(t,:,:)) squeeze(dep_kern(t,:,:))]);shading interp; set(gca,'yscale','log')
%     xlabel('Phase (rad)','fontsize',16)
%     ylabel('Frequency (Hz)','fontsize',16)
%     %     subplot(2,2,3)
%     %      pcolor(doub_pax,wfreqs,[squeeze(ind_kern_po(t,:,:)) squeeze(ind_kern_po(t,:,:))]);shading interp; set(gca,'yscale','log')
%     set(gca,'fontsize',20)
%     input('')
%     clf
% end

%%
for t =1 :96
    subplot(2,2,1)
    pcolor(pax(1:end-1),wfreqs,reshape(stim_ind_filt_pos(t,:),nbins,length(wfreqs))');shading flat
    subplot(2,2,2)
    pcolor(pax(1:end-1),wfreqs,reshape(stim_dep_filt_pg(t,:),nbins,length(wfreqs))');shading flat
     subplot(2,2,3)
    pcolor(pax(1:end-1),wfreqs,reshape(stim_ind_filt(t,:),nbins,length(wfreqs))');shading flat
    subplot(2,2,4)
    pcolor(pax(1:end-1),wfreqs,reshape(stim_dep_filt(t,:),nbins,length(wfreqs))');shading flat
  t
    input('')
  clf
  
%   temp = [xv_null_LL(t) xv_so_LL(t) xv_po_LL(t) xv_pos_LL(t) xv_pg_LL(t) xv_LL(t)];
%   plot(temp,'o-')
%   input('')
%   clf
end

%%
%ind_kern_po(t,:,:) = reshape(stim_ind_filt_po(t,:),nbins,length(wfreqs))';
%     ind_kern(t,:,:) = reshape(stim_ind_filt(t,:),nbins,length(wfreqs))';
%     dep_kern(t,:,:) = reshape(stim_dep_filt(t,:),nbins,length(wfreqs))';
