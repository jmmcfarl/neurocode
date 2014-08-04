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

load ./Expt1_newcompiled_data_fixeddelay_d1p25_34.mat
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

dt = dt*2;

% PARSE DATA INTO FIXATIONS
diff_used_inds = [1; diff(used_inds)];
rel_fix_start_inds = [1; find(diff_used_inds > 1)];
rel_fix_stop_inds = [(find(diff_used_inds > 1)-1); NT];
n_fixs = length(rel_fix_start_inds);

%%
cd ~/Data/bruce/7_15_12/G034/
load ./gabor_tracking_varmeans
% load ./expt1_lfp_alpha_phase_7_15
% load ./expt1_lfp_alpha_phase_7_14_delay_v2.mat
load ./expt1_lfp_alpha_phase_delay_all.mat
% load ./expt1_lfp_alpha_phase_10_20_delay
% load ./expt1_lfp_alpha_phase_5_12_delay
%load ./expt1_lfp_alpha_phase

%%
load ./expt1_eyedata_34_new.mat all_sac_start_times

sac_inds = nan(size(all_sac_start_times));
for i = 1:length(sac_inds)
   [mdist,mloc] = min(abs(all_sac_start_times(i)-new_t_axis)); 
   if mdist < .01
       sac_inds(i) = mloc;
   end
end
sac_inds = sac_inds(~isnan(sac_inds));
new_sac_vec = zeros(size(new_t_axis));
new_sac_vec(sac_inds) = 1;

%%
cur_dt = median(diff(full_t));
desired_dt = cur_dt/4;
trial_start_inds = 1+find(diff(new_trial_vec') ~= 0);
trial_stop_inds = trial_start_inds - 1;
trial_start_inds = [1; trial_start_inds];
trial_stop_inds = [trial_stop_inds; length(new_expt_vec)];

trial_start_inds = trial_start_inds + round(0.1/desired_dt);

n_trials = length(trial_stop_inds);
new_binned_spks = nan(length(new_t_axis),96);
Expt_nu = unique(new_expt_vec);
n_expts = length(Expt_nu);
for tt = 1:n_expts
    cur_trial_set = find(new_expt_vec(trial_start_inds) == Expt_nu(tt));
    cur_n_trials = length(cur_trial_set);
    load(sprintf('Expt%dClusterTimes.mat',Expt_nu(tt)));
    for i = 1:cur_n_trials
        cur_inds = trial_start_inds(cur_trial_set(i)):trial_stop_inds(cur_trial_set(i));
        cur_binned = zeros(length(cur_inds),96);
        for c = 1:96
            temp = histc(Clusters{c}.times,new_t_axis(cur_inds));
            cur_binned(:,c) = temp;
        end
        new_binned_spks(cur_inds,:) = cur_binned;
    end
end
bad_pts = max(isnan(new_binned_spks),[],2);
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

%%
% flen_t = 0.3;
% cNT = length(uset);
% 
% tent_centers = 0:2:70;
% tbmat = construct_tent_bases(tent_centers,1);
% [ntents,tblen] = size(tbmat);
% 
% tbmat = [zeros(ntents,tblen-1) tbmat];
% 
% sac_Tmat = zeros(cNT,ntents);
% for i = 1:ntents
%     sac_Tmat(:,i) = conv(new_sac_vec(uset),tbmat(i,:),'same');
% end

%% REFIT GEM PARAMETERS
fullX_sh = reshape(resh_X_sh,sdim^2,NT)';
fullX_cropped = fullX_sh(:,new_crop);

%%
close all
l2_ind =2000;
l2_dep = 10000;
l1_ind = 2000;
l1_dep = 10000;
silent = 0;

% l2_w = 1000;

clear beta pvals stim_dep_term phase_dep_term emp_stim*
for t = 1:96
    fprintf('Fitting GEM: Cell %d of %d\n',t,96);
    
    spkbs = convert_to_spikebins(new_binned_spks(uset,t));
    
    gabor_emp1 = get_pgabor_mask_v2(XXc,YYc,gabor_params_f{end}(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XXc,YYc,gabor_params_f{end}(t,1:6),pi/2);
    mask1_out = fullX_cropped*gabor_emp1(:);
    mask2_out = fullX_cropped*gabor_emp2(:);
    tot_mod_out =sqrt(mask1_out.^2 + mask2_out.^2);
    tot_mod_out = tot_mod_out/std(tot_mod_out);
    tot_mod_out = interp1(full_t,tot_mod_out,new_t_axis(uset))';
    
     cur_alpha_phase = squeeze(all_interp_phases(uset,:,nearest_lfps(t)));
   nbins = 30;
    pax = linspace(-pi,pi,nbins+1);
    Tmat = [];
    for ww = 1:length(wfreqs)
        cur_tb = tbrep(cur_alpha_phase(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        Tmat = [Tmat cur_tb];
    end
    
    Xmat =  [Tmat bsxfun(@times,Tmat,tot_mod_out)];
    
%     Xmat_sac = [sac_Tmat bsxfun(@times,sac_Tmat,tot_mod_out)];
%     Xmat = [Xmat Xmat_sac];
    
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
%     K0 = 0.01*randn(klen+1,1);
    lamrange2 = [];
    lamrange = [];
    for ww = 1:length(wfreqs)
        lamrange2 = [lamrange2; l2_ind (ww-1)*nbins+1 ww*nbins 1];
        lamrange = [lamrange; l1_ind (ww-1)*nbins+1 ww*nbins 1];
    end
    for ww = 1:length(wfreqs)
        lamrange2 = [lamrange2; l2_dep length(wfreqs)*nbins+(ww-1)*nbins+1 length(wfreqs)*nbins+ww*nbins 1];
        lamrange = [lamrange; l1_dep length(wfreqs)*nbins+(ww-1)*nbins+1 length(wfreqs)*nbins+ww*nbins 1];
    end
    
%     for pp = 1:nbins
%        lamrange2 = [lamrange2; l2_w lamrange2(end,3)+1 lamrange2(end,3)+length(wfreqs) 0];
%     end
    
%     lamrange2 = [lamrange2; l2_ind 2*length(wfreqs)*nbins+1 2*length(wfreqs)*nbins+ntents];
%     lamrange2 = [lamrange2; l2_dep 2*length(wfreqs)*nbins+ntents+1 2*length(wfreqs)*nbins+2*ntents];

    llist = [];
    [fitp,grad] = GLMsolve_jmm(Xmat, spkbs, K0, silent, lamrange, lamrange2, [], [],llist, [], 0);
    ov_const(t) = fitp.k(end);
   stim_ind_filt(t,:) = fitp.k(1:nbins*length(wfreqs));
stim_dep_filt(t,:) = fitp.k((nbins*length(wfreqs)+1):2*nbins*length(wfreqs));
% stim_ind = fitp.k(1801:1800+ntents);
% stim_dep = fitp.k(1800+ntents+1:end-1);
% subplot(2,1,1)
%  pcolor(pax(1:end-1),wfreqs,reshape(temp,nbins,length(wfreqs))');shading flat; set(gca,'yscale','log'); colorbar
%  subplot(2,1,2)
%   pcolor(pax(1:end-1),wfreqs,reshape(temp2,nbins,length(wfreqs))');shading flat;set(gca,'yscale','log'); colorbar
%  subplot(3,1,3)
%  plot(tent_centers*desired_dt,stim_ind)
%  hold on
%  plot(tent_centers*desired_dt,stim_dep,'r')
%  pause
%   clf
end

%%
close all
l2_ind =500;
l2_dep = 1000;
l1_ind = 500;
l1_dep = 1000;
silent = 1;

stim_vals = linspace(0,5,5);
cmap = jet(length(stim_vals));

clear beta pvals stim_dep_term phase_dep_term emp_stim*
for t = 1:96
    fprintf('Fitting GEM: Cell %d of %d\n',t,96);
    
    spkbs = convert_to_spikebins(new_binned_spks(uset,t));
    
    gabor_emp1 = get_pgabor_mask_v2(XXc,YYc,gabor_params_f{end}(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XXc,YYc,gabor_params_f{end}(t,1:6),pi/2);
    mask1_out = fullX_cropped*gabor_emp1(:);
    mask2_out = fullX_cropped*gabor_emp2(:);
    tot_mod_out =sqrt(mask1_out.^2 + mask2_out.^2);
    tot_mod_out = tot_mod_out/std(tot_mod_out);
%     tot_mod_out = zscore(sqrt(mask1_out.^2 + mask2_out.^2));
    tot_mod_out = interp1(full_t,tot_mod_out,new_t_axis(uset))';
    
    stim_val_edges = prctile(tot_mod_out,[0 20 40 60 80 100]);
    stim_val_centers = (stim_val_edges(1:end-1)+stim_val_edges(2:end))/2;
    
    
%     cur_alpha_phase = all_interp_alpha_phase(uset,nearest_lfps(t));
%     cur_alpha_phase = all_interp_alpha(:,nearest_lfps(t));
    %     cur_alpha_phase = 0.5*(cos(all_interp_alpha_phase(:,nearest_lfps(t)))+1);
    
    nbins = 30;
    pax = linspace(-pi,pi,nbins);
%     Tmat = tbrep(cur_alpha_phase,pax);
%     pax = linspace(-3,3,nbins);
    Tmat = tbrep(cur_alpha_phase,pax);
    
%      nebins = 15;
%     epax = linspace(-pi,pi,nebins);
       emp_stim_rate(t,:,:) = nan(5,nbins-1);
        for n = 1:5
           cur_set= find(tot_mod_out >= stim_val_edges(n) & tot_mod_out < stim_val_edges(n+1));
           for b = 1:nbins-1
              cur_cur_set = cur_set(cur_alpha_phase(cur_set) >= pax(b) & cur_alpha_phase(cur_set) < pax(b+1));
              emp_stim_rate(t,n,b) = sum(new_binned_spks(uset(cur_cur_set),t))/length(cur_cur_set);
           end
           emp_stim_rate(t,n,:) = smooth(squeeze(emp_stim_rate(t,n,:)),6,'lowess');
        end
        
    
    %     pred = [cur_alpha_phase cur_alpha_phase.*tot_mod_out tot_mod_out];
    %     pred =  [cur_alpha_phase.*tot_mod_out tot_mod_out];
    %     Xmat =  [Tmat tot_mod_out];
%     Xmat =  [Tmat bsxfun(@times,Tmat,tot_mod_out) tot_mod_out];
    Xmat =  [Tmat bsxfun(@times,Tmat,tot_mod_out)];
    %     [beta(t,:),dev,stats] = glmfit(pred,full_binned_spks(:,t),'poisson');
    %     pvals(t,:) = stats.p;
    
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    %     lamrange2 = [l2_dep 1 nbins];
    lamrange2 = [l2_ind 1 nbins; l2_dep (nbins+1) 2*nbins];
    lamrange = [l1_ind 1 nbins; l1_dep (nbins+1) 2*nbins];
    llist = [];
    [fitp,grad] = GLMsolve_jmm(Xmat, spkbs, K0, silent, lamrange, lamrange2, [], [],llist, [], 0);
    ov_const(t) = fitp.k(end);
    phase_dep_term(t,:) = fitp.k(1:nbins);
    %     stim_dep_term(t) = fitp.k(end-1);
    stim_dep_term(t,:) = fitp.k(nbins+1:2*nbins);
   
    mrate(t) = length(spkbs)/size(Xmat,1);
 
    Xmat =  [Tmat];    
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    lamrange2 = [l2_ind 1 nbins];
    lamrange = [l1_ind 1 nbins];
    llist = [];
    [fitp,grad] = GLMsolve_jmm(Xmat, spkbs, K0, silent, lamrange, lamrange2, [], [],llist, [], 0);
    ov_const_po(t) = fitp.k(end);
    phase_dep_term_po(t,:) = fitp.k(1:nbins);

%     subplot(2,1,1)
% plot(pax,phase_dep_term(t,:))
% xlim([-pi pi])
% subplot(2,1,2)
% plot(pax,stim_dep_term(t,:),'r')
% xlim([-pi pi])
% pause
% clf

%% 
pred_rates = bsxfun(@times,stim_dep_term(t,:)',stim_val_centers);
pred_rates = bsxfun(@plus,pred_rates,phase_dep_term(t,:)');
pred_rates = pred_rates + ov_const(t);
pred_rates = log(1+exp(pred_rates));
null_rate = phase_dep_term(t,:)+ov_const(t);
null_rate = log(1+exp(null_rate));
% norm_rates = bsxfun(@rdivide,pred_rates,mean(pred_rates));
drates = max(pred_rates,[],2)-min(pred_rates,[],2);

null_rate_po = phase_dep_term_po(t,:)+ov_const_po(t);
null_rate_po = log(1+exp(null_rate_po));

% opred_rates = bsxfun(@times,old_stim_dep_term(t,:)',stim_vals);
% opred_rates = bsxfun(@plus,opred_rates,old_phase_dep_term(t,:)');
% opred_rates = opred_rates + old_const(t);
% opred_rates = log(1+exp(opred_rates));
% onull_rate = old_phase_dep_term(t,:)+old_const(t);
% onull_rate = log(1+exp(onull_rate));
% % norm_rates = bsxfun(@rdivide,pred_rates,mean(pred_rates));
% odrates = max(opred_rates,[],2)-min(opred_rates,[],2);
% 
% subplot(2,1,1)
% plot(pax,null_rate/dt)
% hold on
% xlim([-pi pi])
% % plot(pax,onull_rate,'r')
% subplot(2,1,2)
% plot(pax,drates/dt)
% xlim([-pi pi])
% hold on
% % plot(pax,odrates,'r')

subplot(2,2,1)
% plot(pax,null_rate)
plot(pax,phase_dep_term(t,:))
hold on
% plot(pax,onull_rate,'r')
subplot(2,2,2)
% plot(pax,drates)
plot(pax,stim_dep_term(t,:))
hold on
subplot(2,2,3)
plot(pax,pred_rates)
subplot(2,2,4)
plot(pax(1:end-1),squeeze(emp_stim_rate(t,:,:))')
% plot(pax,odrates,'r')
% 
pause
clf

all_null_rates(t,:) = null_rate/mrate(t);
all_null_rates_po(t,:) = null_rate_po/mrate(t);
all_drates(t,:) = drates/mrate(t);
end

%%
norm_drates = bsxfun(@rdivide,all_drates,mean(all_drates,2));

%%
null_mod = (max(all_null_rates,[],2)-min(all_null_rates,[],2))./(max(all_null_rates,[],2)+min(all_null_rates,[],2));
null_mod_po = (max(all_null_rates_po,[],2)-min(all_null_rates_po,[],2))./(max(all_null_rates_po,[],2)+min(all_null_rates_po,[],2));
stim_mod = (max(all_drates,[],2)-min(all_drates,[],2))./(max(all_drates,[],2)+min(all_drates,[],2));
%%

old_stim_dep_term = stim_dep_term;
old_phase_dep_term = phase_dep_term;
old_const = ov_const;
