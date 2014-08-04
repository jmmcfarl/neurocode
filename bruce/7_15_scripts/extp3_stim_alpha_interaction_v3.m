%% Load Data
cd ~/Data/bruce/7_15_12/G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt3_fixbased_data.mat
fullX = resh_all_stims/std(resh_all_stims(:));

%% crop stimulus for the purpose of faster gabor function fitting
new_RF_patch = [-0.11 0.9; -0.9 0.1]; %location of RFs in degrees [x1 x2;y1 y2]
[curXX,curYY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
curXX = curXX(:);
curYY = curYY(:);

xpatch_inds_cropped = find(xax >= new_RF_patch(1,1) & xax <= new_RF_patch(1,2));
ypatch_inds_cropped = find(yax >= new_RF_patch(2,1) & yax <= new_RF_patch(2,2));

new_crop = find(curXX >= new_RF_patch(1,1) & curXX <= new_RF_patch(1,2) & ...
    curYY >= new_RF_patch(2,1) & curYY <= new_RF_patch(2,2));
fullX_cropped = fullX(:,new_crop);

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

%%
cd ~/Data/bruce/7_15_12/G034/
load ./gabor_tracking_varmeans

load ./expt3_lfp_alpha_phase_delay_all.mat new_* wfreqs nearest_* all_interp_phases

%%
Expt_nu = [7:12]; %expt 3 34

cur_dt = median(diff(full_t));
desired_dt = cur_dt/2;
trial_start_ninds = 1+find(diff(new_trial_vec') ~= 0);
trial_stop_ninds = trial_start_ninds - 1;
trial_start_ninds = [1; trial_start_ninds];
trial_stop_ninds = [trial_stop_ninds; length(new_expt_vec)];

n_trials = length(trial_stop_ninds);
new_binned_spks = nan(length(new_t_axis),96);
Expt_nu = unique(new_expt_vec);
n_expts = length(Expt_nu);
new_trial_vec = nan(length(new_t_axis),1);
for tt = 1:n_expts
    cur_trial_set = find(new_expt_vec(trial_start_ninds) == Expt_nu(tt));
    cur_n_trials = length(cur_trial_set);
    load(sprintf('Expt%dClusterTimes.mat',Expt_nu(tt)));
    for i = 1:cur_n_trials
        cur_inds = trial_start_ninds(cur_trial_set(i)):trial_stop_ninds(cur_trial_set(i));
        cur_binned = zeros(length(cur_inds),96);
        for c = 1:96
             temp = histc(Clusters{c}.times,new_t_axis(cur_inds));
           cur_binned(:,c) = temp;
        end
        new_binned_spks(cur_inds,:) = cur_binned;
        new_trial_vec(cur_inds) = cur_trial_set(i);
    end
end
bad_pts = max(isnan(new_binned_spks),[],2);
% uset = find(new_expt_vec ~= 18 & bad_pts' == 0);
uset = find(bad_pts' == 0);

%%
load ./gabor_tracking_varmeans.mat gabor*
gabor_params = gabor_params_f{end};
clear gabor_*filt
for t = 1:96
    
    gabor_emp1 = get_pgabor_mask_v2(XXc,YYc,gabor_params(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XXc,YYc,gabor_params(t,1:6),pi/2);
    
    gabor_emp1_filt(t,:) = gabor_emp1(:);
    gabor_emp2_filt(t,:) = gabor_emp2(:);
end


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
flen_t = 0.5;
tent_centers = [0:dt:0.15];
cur_sp = dt;
while max(tent_centers) < flen_t
    tent_centers = [tent_centers (tent_centers(end)+cur_sp)];
    cur_sp = cur_sp + dt/6;
end

tent_centers = round(tent_centers/dt);
tbmat = construct_tent_bases(tent_centers,0.5);
[ntents,tblen] = size(tbmat);

tbmat = [zeros(ntents,tblen-1) tbmat];

trial_inds = zeros(size(new_t_axis));
trial_inds(trial_start_ninds) = 1;
trial_Tmat = zeros(length(new_t_axis),ntents);
for i = 1:ntents
    trial_Tmat(:,i) = conv(trial_inds,tbmat(i,:),'same');
end

beg_dur = round(0.1/desired_dt);
late_indicator = zeros(size(new_t_axis));
for i = 1:length(trial_start_ninds)
   cur_inds = (beg_dur+trial_start_ninds(i)):trial_stop_ninds(i);
   late_indicator(cur_inds) = 1;
end
%%
close all
use_freqs = find(wfreqs < 50);

reg_params.dl1_ind = 3000;
reg_params.dl1_dep =2500;
reg_params.dl2_ind = 6000;
reg_params.dl2_dep = 25000;
reg_params.dl2_freq_ind = 500;
reg_params.dl2_freq_dep = 1500;
reg_params.dl2_time_ind = 100;
reg_params.dl2_time_dep = 400;
reg_params.l1_ind = 0;
reg_params.l1_dep = 0;
reg_params.is_phase = 1;

silent = 1;
nbins = 30;
pax = linspace(-pi,pi,nbins+1);
doub_pax = [pax(1:end-1) 2*pi+pax(1:end-1)];

for t = 1:96
    fprintf('Fitting GEM: Cell %d of %d\n',t,96);
    
    Robs = new_binned_spks(uset(tr_inds),t);
    spkbs = convert_to_spikebins(Robs);
    
    mask1_out = fullX_cropped*gabor_emp1_filt(t,:)';
    mask2_out = fullX_cropped*gabor_emp2_filt(t,:)';
    tot_mod_out =sqrt(mask1_out.^2 + mask2_out.^2);
    tot_mod_out = tot_mod_out/std(tot_mod_out);
    tot_mod_out = tot_mod_out(full_stim_ids);
    mod_out = interp1(full_t,tot_mod_out,new_t_axis(uset(tr_inds)))';
    mod_out(isnan(mod_out)) = nanmean(mod_out);
    
    cur_alpha_phase = squeeze(all_interp_phases(uset(tr_inds),use_freqs,nearest_lfps(t)));
    Pmat = [];
    for ww = 1:length(use_freqs)
        cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        Pmat = [Pmat cur_tb];
    end
    
%     aPmat =  [Pmat bsxfun(@times,Pmat,mod_out)];
%     tr_late_indicator = logical(late_indicator(uset(tr_inds)))';
%     aPmat = bsxfun(@times,aPmat,tr_late_indicator);
%     
%     Tmat = [trial_Tmat(uset(tr_inds),:) bsxfun(@times,trial_Tmat(uset(tr_inds),:),mod_out)];
%     Xmat = [aPmat Tmat];
%     
%     stim_params = [nbins,length(wfreqs),ntents,2];
%     klen = size(Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp] = fit_GLM_phase_model(Xmat, Robs, K0,1, stim_params, reg_params,1);
%     stim_ind_pfilt(t,:) = fitp.k(1:nbins*length(wfreqs));
%     stim_dep_pfilt(t,:) = fitp.k((nbins*length(wfreqs)+1):2*nbins*length(wfreqs));
%     stim_ind_tfilt(t,:) = fitp.k((2*nbins*length(wfreqs)+1):(2*nbins*length(wfreqs)+ntents));
%     stim_dep_tfilt(t,:) = fitp.k((2*nbins*length(wfreqs)+ntents+1):end-1);
    
    
    aPmat =  [Pmat];
    tr_late_indicator = logical(late_indicator(uset(tr_inds)))';
    aPmat = bsxfun(@times,aPmat,tr_late_indicator);    
    Tmat = [trial_Tmat(uset(tr_inds),:) bsxfun(@times,trial_Tmat(uset(tr_inds),:),mod_out)];
    Xmat = [aPmat Tmat];   
    stim_params = [nbins,length(use_freqs),ntents,1];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp] = fit_GLM_phase_model(Xmat, Robs, K0,1, stim_params, reg_params,1);
    stim_ind_pfilt(t,:) = fitp.k(1:nbins*length(use_freqs));
    stim_ind_tfilt(t,:) = fitp.k((nbins*length(use_freqs)+1):(nbins*length(use_freqs)+ntents));
    stim_dep_tfilt(t,:) = fitp.k((nbins*length(use_freqs)+ntents+1):(nbins*length(use_freqs)+2*ntents));


    Tmat = [trial_Tmat(uset(tr_inds),:)];
    Xmat = [aPmat Tmat];
    stim_params = [nbins,length(use_freqs),ntents,1];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_ns] = fit_GLM_phase_model(Xmat, Robs, K0,1, stim_params, reg_params,1);
    stim_ind_ns_pfilt(t,:) = fitp_off.k(1:nbins*length(use_freqs));
    stim_ind_ns_tfilt(t,:) = fitp_off.k((nbins*length(use_freqs)+1):(nbins*length(use_freqs)+ntents));
    
    
    stim_params = [0,0,ntents];
    Tmat = [trial_Tmat(uset(tr_inds),:) bsxfun(@times,trial_Tmat(uset(tr_inds),:),mod_out)];
    Xmat = [Tmat];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_to] = fit_GLM_phase_model(Xmat, Robs, K0, 1, stim_params, reg_params,1);
    ov_const(t) = fitp.k(end);
    stim_ind_to_filt(t,:) = fitp_to.k((1):(ntents));
    stim_dep_to_filt(t,:) = fitp_to.k((ntents+1):end-1);

    stim_params = [0,0,ntents];
    Tmat = [trial_Tmat(uset(tr_inds),:)];
    Xmat = [Tmat];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_to_ns] = fit_GLM_phase_model(Xmat, Robs, K0, 1, stim_params, reg_params,1);
    ov_const(t) = fitp.k(end);
    stim_ind_to_ns_filt(t,:) = fitp_to.k((1):(ntents));
    stim_dep_to_ns_filt(t,:) = fitp_to.k((ntents+1):end-1);
    
    
    
    xv_Robs = new_binned_spks(uset(xv_inds),t);
    mod_out = interp1(full_t,tot_mod_out,new_t_axis(uset(xv_inds)))';
    mod_out(isnan(mod_out)) = nanmean(mod_out);
    
     cur_alpha_phase = squeeze(all_interp_phases(uset(xv_inds),use_freqs,nearest_lfps(t)));
    Pmat = [];
    for ww = 1:length(use_freqs)
        cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        Pmat = [Pmat cur_tb];
    end
    cur_late_indicator = logical(late_indicator(uset(xv_inds)))';

    %     aPmat =  [Pmat bsxfun(@times,Pmat,mod_out)];
%     aPmat = bsxfun(@times,aPmat,cur_late_indicator);
%     
%     Tmat = [trial_Tmat(uset(xv_inds),:) bsxfun(@times,trial_Tmat(uset(xv_inds),:),mod_out)];
%     Xmat = [aPmat Tmat];
%    
%     xv_pred_rate = Xmat*fitp.k(1:end-1) + fitp.k(end);
%     Tmat = [trial_Tmat(uset(xv_inds),:) bsxfun(@times,trial_Tmat(uset(xv_inds),:),mod_out)];
%     Xmat = [aPmat Tmat];
%    
%     xv_pred_rate = Xmat*fitp.k(1:end-1) + fitp.k(end);
%     xv_pred_rate = log(1+exp(xv_pred_rate));
%     xv_LL(t) = -sum(xv_Robs.*log(xv_pred_rate) - xv_pred_rate)/sum(xv_Robs);
%     xv_LL_late(t) = -sum(xv_Robs(cur_late_indicator).*log(xv_pred_rate(cur_late_indicator)) - ...
%         xv_pred_rate(cur_late_indicator))/sum(xv_Robs(cur_late_indicator));

aPmat =  [Pmat];
aPmat = bsxfun(@times,aPmat,cur_late_indicator);
Tmat = [trial_Tmat(uset(xv_inds),:) bsxfun(@times,trial_Tmat(uset(xv_inds),:),mod_out)];
Xmat = [aPmat Tmat];
xv_pred_rate = Xmat*fitp.k(1:end-1) + fitp.k(end);
xv_pred_rate = log(1+exp(xv_pred_rate));
xv_LL(t) = -sum(xv_Robs.*log(xv_pred_rate) - xv_pred_rate)/sum(xv_Robs);
xv_LL_late(t) = -sum(xv_Robs(cur_late_indicator).*log(xv_pred_rate(cur_late_indicator)) - ...
    xv_pred_rate(cur_late_indicator))/sum(xv_Robs(cur_late_indicator));
    
    Tmat = [trial_Tmat(uset(xv_inds),:)];
    Xmat = [aPmat Tmat];
    xv_pred_rate = Xmat*fitp_ns.k(1:end-1) + fitp_ns.k(end);
    xv_pred_rate = log(1+exp(xv_pred_rate));
    xv_ns_LL(t) = -sum(xv_Robs.*log(xv_pred_rate) - xv_pred_rate)/sum(xv_Robs);
    xv_ns_LL_late(t) = -sum(xv_Robs(cur_late_indicator).*log(xv_pred_rate(cur_late_indicator)) - ...
        xv_pred_rate(cur_late_indicator))/sum(xv_Robs(cur_late_indicator));

Tmat = [trial_Tmat(uset(xv_inds),:) bsxfun(@times,trial_Tmat(uset(xv_inds),:),mod_out)];
    Xmat = [Tmat];
     xv_pred_rate = Xmat*fitp_to.k(1:end-1) + fitp_to.k(end);
    xv_pred_rate = log(1+exp(xv_pred_rate));
    xv_to_LL(t) = -sum(xv_Robs.*log(xv_pred_rate) - xv_pred_rate)/sum(xv_Robs);
    xv_to_LL_late(t) = -sum(xv_Robs(cur_late_indicator).*log(xv_pred_rate(cur_late_indicator)) - ...
        xv_pred_rate(cur_late_indicator))/sum(xv_Robs(cur_late_indicator));

    Tmat = [trial_Tmat(uset(xv_inds),:)];
        Xmat = [Tmat];
     xv_pred_rate = Xmat*fitp_to_ns.k(1:end-1) + fitp_to_ns.k(end);
    xv_pred_rate = log(1+exp(xv_pred_rate));
    xv_to_ns_LL(t) = -sum(xv_Robs.*log(xv_pred_rate) - xv_pred_rate)/sum(xv_Robs);
    xv_to_ns_LL_late(t) = -sum(xv_Robs(cur_late_indicator).*log(xv_pred_rate(cur_late_indicator)) - ...
        xv_pred_rate(cur_late_indicator))/sum(xv_Robs(cur_late_indicator));

    avg_rate = mean(Robs);
    null_pred = avg_rate*ones(size(xv_Robs));
    xv_null_LL(t) = -sum(xv_Robs.*log(null_pred)-null_pred)/sum(xv_Robs);
     avg_rate = mean(Robs(tr_late_indicator));
    null_pred = avg_rate*ones(size(xv_Robs));
   xv_null_LL_late(t) = -sum(xv_Robs(cur_late_indicator).*log(null_pred(cur_late_indicator)) -...
        null_pred(cur_late_indicator))/sum(xv_Robs(cur_late_indicator));

end

%%
save expt3_stim_alpha_models xv_* stim_* tent_centers wfreqs nbins pax
%%
boxplot([xv_LL_late(:)-xv_null_LL_late(:) xv_to_LL_late(:)-xv_null_LL_late(:)])
figure
boxplot([xv_LL(:)-xv_null_LL(:) xv_to_LL(:)-xv_null_LL(:)])
%%

for t = 1:96
    subplot(2,1,1)
    pcolor(pax(1:end-1),wfreqs,reshape(stim_ind_pfilt(t,:),nbins,length(wfreqs))');shading flat
    subplot(2,1,2)
    pcolor(pax(1:end-1),wfreqs,reshape(stim_dep_pfilt(t,:),nbins,length(wfreqs))');shading flat
    input('')
    clf

end
%%
for t = 1:96
    stim_ind_nmat(t,:,:) = reshape(stim_ind_pfilt(t,:)/std(stim_ind_pfilt(t,:)),nbins,length(wfreqs));
    stim_ind_mat(t,:,:) = reshape(stim_ind_pfilt(t,:),nbins,length(wfreqs));
    stim_dep_nmat(t,:,:) = reshape(stim_dep_pfilt(t,:)/std(stim_dep_pfilt(t,:)),nbins,length(wfreqs));
    stim_dep_mat(t,:,:) = reshape(stim_dep_pfilt(t,:),nbins,length(wfreqs));
end
dom_ind = squeeze(max(stim_ind_nmat,[],2)-min(stim_ind_nmat,[],2));
dom_dep = squeeze(max(stim_dep_nmat,[],2)-min(stim_dep_nmat,[],2));