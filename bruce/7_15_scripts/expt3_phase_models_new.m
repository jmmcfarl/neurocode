clear all
clc

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

load ./expt3_lfp_alpha_phase_new.mat new_* wfreqs nearest_* all_interp_phases all_interp_amps
all_interp_amps = bsxfun(@rdivide,all_interp_amps,std(all_interp_amps));

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
tent_centers = [0:dt:0.5];
cur_sp = dt;
% while max(tent_centers) < flen_t
%     tent_centers = [tent_centers (tent_centers(end)+cur_sp)];
%     cur_sp = cur_sp + dt/6;
% end

tent_centers = round(tent_centers/desired_dt);
tbmat = construct_tent_bases(tent_centers,1);
[ntents,tblen] = size(tbmat);

tbmat = [zeros(ntents,tblen-1) tbmat];

trial_inds = zeros(size(new_t_axis));
trial_inds(trial_start_ninds) = 1;
trial_Tmat = zeros(length(new_t_axis),ntents);
for i = 1:ntents
    trial_Tmat(:,i) = conv(trial_inds,tbmat(i,:),'same');
end

% box_width = 0.03;
% box_edges = 0:box_width:0.5;
% box_cents = 0.5*box_edges(1:end-1)+0.5*box_edges(2:end);
% n_boxs = length(box_cents);
% max_lag = 0.5;
% lags = (round(-max_lag/desired_dt):round(max_lag/desired_dt))*desired_dt;
% boxmat = zeros(n_boxs,length(lags));
% trial_Bmat = zeros(length(new_t_axis),n_boxs);
% for i = 1:n_boxs
%     curset = find(lags >= box_edges(i) & lags < box_edges(i+1));
%     boxmat(i,curset) = 1;
%     trial_Bmat(:,i) = conv(trial_inds,boxmat(i,:),'same');
% end
% trial_Bmat = logical(trial_Bmat);

beg_dur = round(0.15/desired_dt);
late_indicator = zeros(size(new_t_axis));
for i = 1:length(trial_start_ninds)
   cur_inds = (beg_dur+trial_start_ninds(i)):trial_stop_ninds(i);
   late_indicator(cur_inds) = 1;
end
%%
close all
% use_freqs = find(wfreqs < 90);
NL_type = 1; %exp 

reg_params.dl1_ind = 4000;
reg_params.dl1_dep =4000;
reg_params.dl2_ind = 8000;
reg_params.dl2_dep = 8000;
reg_params.dl2_freq_ind = 200;
reg_params.dl2_freq_dep = 200;
reg_params.dl2_time_ind = 200;
reg_params.dl2_time_dep = 600;
reg_params.l1_dep = 0;
reg_params.l1_ind = 0;
reg_params.l1t_ind = 0;
reg_params.l1t_dep = 0;
reg_params.is_phase = 1;


reg_params_sin.dl1_ind = 0;
reg_params_sin.dl1_dep = 0;
reg_params_sin.dl2_ind = 20000;
reg_params_sin.dl2_dep = 20000;
reg_params_sin.dl2_freq_ind = 0;
reg_params_sin.dl2_freq_dep = 0;
reg_params_sin.dl2_time_ind = 200;
reg_params_sin.dl2_time_dep = 600;
reg_params_sin.l1_ind = 0;
reg_params_sin.l1_dep = 0;
reg_params_sin.l1t_ind = 0;
reg_params_sin.l1t_dep = 0;
reg_params_sin.is_phase = 0;


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
%     tot_mod_out = zscore(tot_mod_out);
    tot_mod_out = tot_mod_out(full_stim_ids);
    mod_out = interp1(full_t,tot_mod_out,new_t_axis(uset(tr_inds)))';
    mod_out(isnan(mod_out)) = nanmean(mod_out);
    
    stim_params = [0,0,ntents];
    Tmat = [trial_Tmat(uset(tr_inds),:) bsxfun(@times,trial_Tmat(uset(tr_inds),:),mod_out)];
    Xmat = [Tmat];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_to] = fit_GLM_phase_model(Xmat, Robs, K0, 1, stim_params, reg_params,1,NL_type);
    stim_ind_to_filt(t,:) = fitp_to.k((1):(ntents));
    stim_dep_to_filt(t,:) = fitp_to.k((ntents+1):end-1);

    stim_params = [0,0,ntents];
    Tmat = [trial_Tmat(uset(tr_inds),:)];
    Xmat = [Tmat];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_to_ns] = fit_GLM_phase_model(Xmat, Robs, K0, 1, stim_params, reg_params,1,NL_type);
    stim_ind_to_ns_filt(t,:) = fitp_to_ns.k((1):(ntents));

    cur_phases = squeeze(all_interp_phases(uset(tr_inds),:,nearest_lfps(t)));
    cur_amp = squeeze(all_interp_amps(uset(tr_inds),:,nearest_lfps(t)));    
    Pmat = [];
    for ww = 1:length(wfreqs)
        cur_tb = tbrep(cur_phases(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        Pmat = [Pmat cur_tb];
    end
    tr_late_indicator = logical(late_indicator(uset(tr_inds)))';
    aPmat = bsxfun(@times,Pmat,tr_late_indicator);

    Tmat = [trial_Tmat(uset(tr_inds),:) bsxfun(@times,trial_Tmat(uset(tr_inds),:),mod_out)];
    Xmat = [aPmat Tmat];   
    stim_params = [nbins,length(wfreqs),ntents,1];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp] = fit_GLM_phase_model(Xmat, Robs, K0,1, stim_params, reg_params,1,NL_type);
    stim_ind_pfilt(t,:) = fitp.k(1:nbins*length(wfreqs));
    stim_ind_tfilt(t,:) = fitp.k((nbins*length(wfreqs)+1):(nbins*length(wfreqs)+ntents));
    stim_dep_tfilt(t,:) = fitp.k((nbins*length(wfreqs)+ntents+1):(nbins*length(wfreqs)+2*ntents));


    phase_set = [reshape(cos(cur_phases),length(tr_inds),length(wfreqs)) ...
        reshape(sin(cur_phases),length(tr_inds),length(wfreqs))];
    ampphase_set = [reshape(cur_amp,length(tr_inds),length(wfreqs)).*reshape(cos(cur_phases),length(tr_inds),length(wfreqs)) ...
        reshape(cur_amp,length(tr_inds),length(wfreqs)).*reshape(sin(cur_phases),length(tr_inds),length(wfreqs))];
    
    Xmat = [phase_set Tmat];
    stim_params = [length(wfreqs),1,ntents];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_sinphase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params_sin,1,NL_type);
    res_phase_cfilt(t,:) = fitp_sinphase.k(1:length(wfreqs));
    res_phase_sfilt(t,:) = fitp_sinphase.k((length(wfreqs)+1):length(wfreqs)*2);
    res_phase_ind_tfilt(t,:) = fitp_sinphase.k((length(wfreqs)*2+1):(length(wfreqs)*2+ntents));
    res_phase_dep_tfilt(t,:) = fitp_sinphase.k((length(wfreqs)*2+ntents+1):(length(wfreqs)*2+2*ntents));

        Xmat = [ampphase_set Tmat];
    stim_params = [length(wfreqs),1,ntents];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_sinampphase] = fit_GLM_phase_model(Xmat, Robs, K0,silent, stim_params, reg_params_sin,1,NL_type);
    res_ampphase_cfilt(t,:) = fitp_sinampphase.k(1:length(wfreqs));
    res_ampphase_sfilt(t,:) = fitp_sinampphase.k((length(wfreqs)+1):length(wfreqs)*2);
    res_ampphase_ind_tfilt(t,:) = fitp_sinphase.k((length(wfreqs)*2+1):(length(wfreqs)*2+ntents));
    res_ampphase_dep_tfilt(t,:) = fitp_sinphase.k((length(wfreqs)*2+ntents+1):(length(wfreqs)*2+2*ntents));

    
    Tmat = [trial_Tmat(uset(tr_inds),:)];
    Xmat = [aPmat Tmat];
    stim_params = [nbins,length(wfreqs),ntents,1];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_ns] = fit_GLM_phase_model(Xmat, Robs, K0,1, stim_params, reg_params,1);
    stim_ind_ns_pfilt(t,:) = fitp_ns.k(1:nbins*length(wfreqs));
    stim_ind_ns_tfilt(t,:) = fitp_ns.k((nbins*length(wfreqs)+1):(nbins*length(wfreqs)+ntents));
    
    
    
    xv_Robs = new_binned_spks(uset(xv_inds),t);
    mod_out = interp1(full_t,tot_mod_out,new_t_axis(uset(xv_inds)))';
    mod_out(isnan(mod_out)) = nanmean(mod_out);
    
     cur_phases = squeeze(all_interp_phases(uset(xv_inds),:,nearest_lfps(t)));
    cur_amp = squeeze(all_interp_amps(uset(xv_inds),:,nearest_lfps(t)));    
    Pmat = [];
    for ww = 1:length(wfreqs)
        cur_tb = tbrep(cur_phases(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        Pmat = [Pmat cur_tb];
    end
    cur_late_indicator = logical(late_indicator(uset(xv_inds)))';
    
    aPmat = bsxfun(@times,Pmat,cur_late_indicator);
        
    Tmat = [trial_Tmat(uset(xv_inds),:) bsxfun(@times,trial_Tmat(uset(xv_inds),:),mod_out)];
    Xmat = [aPmat Tmat];
    xv_pred_rate = Xmat*fitp.k(1:end-1) + fitp.k(end);
    if NL_type == 0
        xv_pred_rate = log(1+exp(xv_pred_rate));
    else
        xv_pred_rate = exp(xv_pred_rate);
    end
    xv_LL(t) = -sum(xv_Robs.*log(xv_pred_rate) - xv_pred_rate)/sum(xv_Robs);
    xv_LL_late(t) = -sum(xv_Robs(cur_late_indicator).*log(xv_pred_rate(cur_late_indicator)) - ...
        xv_pred_rate(cur_late_indicator))/sum(xv_Robs(cur_late_indicator));

    
    Tmat = [trial_Tmat(uset(xv_inds),:)];
    Xmat = [aPmat Tmat];
    xv_ns_pred_rate = Xmat*fitp_ns.k(1:end-1) + fitp_ns.k(end);
    if NL_type==0
    xv_ns_pred_rate = log(1+exp(xv_ns_pred_rate));
    else
    xv_ns_pred_rate = exp(xv_ns_pred_rate);        
    end
    xv_ns_LL(t) = -sum(xv_Robs.*log(xv_ns_pred_rate) - xv_ns_pred_rate)/sum(xv_Robs);
    xv_ns_LL_late(t) = -sum(xv_Robs(cur_late_indicator).*log(xv_ns_pred_rate(cur_late_indicator)) - ...
        xv_ns_pred_rate(cur_late_indicator))/sum(xv_Robs(cur_late_indicator));

    Tmat = [trial_Tmat(uset(xv_inds),:) bsxfun(@times,trial_Tmat(uset(xv_inds),:),mod_out)];
    Xmat = [Tmat];
    xv_to_pred_rate = Xmat*fitp_to.k(1:end-1) + fitp_to.k(end);
    if NL_type == 0
    xv_to_pred_rate = log(1+exp(xv_to_pred_rate));
    else
    xv_to_pred_rate = exp(xv_to_pred_rate);        
    end
    xv_to_LL(t) = -sum(xv_Robs.*log(xv_to_pred_rate) - xv_to_pred_rate)/sum(xv_Robs);
    xv_to_LL_late(t) = -sum(xv_Robs(cur_late_indicator).*log(xv_to_pred_rate(cur_late_indicator)) - ...
        xv_to_pred_rate(cur_late_indicator))/sum(xv_Robs(cur_late_indicator));
    
    phase_set = [reshape(cos(cur_phases),length(xv_inds),length(wfreqs)) ...
        reshape(sin(cur_phases),length(xv_inds),length(wfreqs))];
    ampphase_set = [reshape(cur_amp,length(xv_inds),length(wfreqs)).*reshape(cos(cur_phases),length(xv_inds),length(wfreqs)) ...
        reshape(cur_amp,length(xv_inds),length(wfreqs)).*reshape(sin(cur_phases),length(xv_inds),length(wfreqs))];
 
    Xmat = [phase_set Tmat];
    xv_pred_rate = Xmat*fitp_sinphase.k(1:end-1) + fitp_sinphase.k(end);
    if NL_type == 0
    xv_pred_rate = log(1+exp(xv_pred_rate));
    else
    xv_pred_rate = exp(xv_pred_rate);        
    end
    xv_sin_LL(t) = -sum(xv_Robs.*log(xv_pred_rate) - xv_pred_rate)/sum(xv_Robs);
    xv_sin_LL_late(t) = -sum(xv_Robs(cur_late_indicator).*log(xv_pred_rate(cur_late_indicator)) - ...
        xv_pred_rate(cur_late_indicator))/sum(xv_Robs(cur_late_indicator));

    Xmat = [ampphase_set Tmat];
    xv_pred_rate = Xmat*fitp_sinampphase.k(1:end-1) + fitp_sinampphase.k(end);
    if NL_type == 0
    xv_pred_rate = log(1+exp(xv_pred_rate));
    else
    xv_pred_rate = exp(xv_pred_rate);        
    end
    xv_sinamp_LL(t) = -sum(xv_Robs.*log(xv_pred_rate) - xv_pred_rate)/sum(xv_Robs);
    xv_sinamp_LL_late(t) = -sum(xv_Robs(cur_late_indicator).*log(xv_pred_rate(cur_late_indicator)) - ...
        xv_pred_rate(cur_late_indicator))/sum(xv_Robs(cur_late_indicator));
    
    Tmat = [trial_Tmat(uset(xv_inds),:)];
    Xmat = [Tmat];
    xv_to_ns_pred_rate = Xmat*fitp_to_ns.k(1:end-1) + fitp_to_ns.k(end);
    if NL_type == 0
    xv_to_ns_pred_rate = log(1+exp(xv_to_ns_pred_rate));
    else
    xv_to_ns_pred_rate = exp(xv_to_ns_pred_rate);        
    end
    xv_to_ns_LL(t) = -sum(xv_Robs.*log(xv_to_ns_pred_rate) - xv_to_ns_pred_rate)/sum(xv_Robs);
    xv_to_ns_LL_late(t) = -sum(xv_Robs(cur_late_indicator).*log(xv_to_ns_pred_rate(cur_late_indicator)) - ...
        xv_to_ns_pred_rate(cur_late_indicator))/sum(xv_Robs(cur_late_indicator));
             
    
    avg_rate = mean(Robs);
    null_pred = avg_rate*ones(size(xv_Robs));
    xv_null_LL(t) = -sum(xv_Robs.*log(null_pred)-null_pred)/sum(xv_Robs);
    avg_rate = mean(Robs(tr_late_indicator));
    null_pred = avg_rate*ones(size(xv_Robs));
    xv_null_LL_late(t) = -sum(xv_Robs(cur_late_indicator).*log(null_pred(cur_late_indicator)) -...
        null_pred(cur_late_indicator))/sum(xv_Robs(cur_late_indicator));
end

%%
save expt3_stim_alpha_models_new xv_* stim_* tent_centers wfreqs nbins pax res* 
%%
phasemod_kern = sqrt(res_phase_cfilt.^2 + res_phase_sfilt.^2);
ampphasemod_kern = sqrt(res_ampphase_cfilt.^2 + res_ampphase_sfilt.^2);

%%
to_ns_imp = (xv_null_LL - xv_to_ns_LL)/log(2);
to_imp = (xv_null_LL - xv_to_LL)/log(2);
phase_imp = (xv_null_LL - xv_sin_LL)/log(2);
ampphase_imp = (xv_null_LL - xv_sinamp_LL)/log(2);
ucells = single_units;
ucells = setdiff(1:96,single_units);

bad_units = find(to_ns_imp(ucells) < 0);
ucells(bad_units) = [];

figure
boxplot([to_ns_imp(ucells); to_imp(ucells); phase_imp(ucells); ampphase_imp(ucells)]');


%%
ucells = 1:17;
rel_full_LL = xv_full_LL(ucells) - xv_null_LL(ucells);
rel_LL = xv_LL(ucells) - xv_null_LL(ucells);
rel_to_LL = xv_to_LL(ucells)-xv_null_LL(ucells);
rel_to_ns_LL = xv_to_ns_LL(ucells) - xv_null_LL(ucells);
rel_ns_LL = xv_ns_LL(ucells) - xv_null_LL(ucells);
figure
boxplot([-rel_full_LL(:) -rel_LL(:) -rel_ns_LL(:) -rel_to_LL(:) -rel_to_ns_LL(:)])
hold on
plot([-rel_full_LL(:) -rel_LL(:) -rel_ns_LL(:) -rel_to_LL(:) -rel_to_ns_LL(:)]','o-')

rel_full_LL = xv_full_LL_late(ucells) - xv_null_LL_late(ucells);
rel_LL = xv_LL_late(ucells) - xv_null_LL_late(ucells);
rel_to_LL = xv_to_LL_late(ucells)-xv_null_LL_late(ucells);
rel_to_ns_LL = xv_to_ns_LL_late(ucells) - xv_null_LL_late(ucells);
rel_ns_LL = xv_ns_LL_late(ucells) - xv_null_LL_late(ucells);
figure
boxplot([-rel_full_LL(:) -rel_LL(:) -rel_ns_LL(:) -rel_to_LL(:) -rel_to_ns_LL(:)])
hold on
plot([-rel_full_LL(:) -rel_LL(:) -rel_ns_LL(:) -rel_to_LL(:) -rel_to_ns_LL(:)]','o-')

%%
ucells = 1:17;
% ucells = single_units;
% figure
% plot(box_cents,mean(xv_to_box_LL(ucells,:)),'g')
% hold on
% plot(box_cents,mean(xv_ns_box_LL(ucells,:)),'r')
% plot(box_cents,mean(xv_box_LL(ucells,:)),'b')
% plot(box_cents,mean(xv_full_box_LL(ucells,:)),'k')

figure
errorbar(box_cents,mean(xv_to_box_LL(ucells,:)),std(xv_to_box_LL(ucells,:))/sqrt(length(ucells)),'g')
hold on
errorbar(box_cents,mean(xv_ns_box_LL(ucells,:)),std(xv_ns_box_LL(ucells,:))/sqrt(length(ucells)),'r')
errorbar(box_cents,mean(xv_box_LL(ucells,:)),std(xv_box_LL(ucells,:))/sqrt(length(ucells)),'b')
errorbar(box_cents,mean(xv_full_box_LL(ucells,:)),std(xv_full_box_LL(ucells,:))/sqrt(length(ucells)),'k')

%%

% for t = single_units
for t = 1:96
    t
    subplot(3,1,1)
    pcolor(pax(1:end-1),wfreqs,reshape(stim_ind_pfilt(t,:),nbins,length(wfreqs))');shading flat
    subplot(3,1,2)
plot(tent_centers*desired_dt,stim_ind_tfilt(t,:))
hold on
plot(tent_centers*desired_dt,stim_ind_to_filt(t,:),'r')
xlim([0 0.5])
    subplot(3,1,3)
plot(tent_centers*desired_dt,stim_dep_tfilt(t,:))
hold on
plot(tent_centers*desired_dt,stim_dep_to_filt(t,:),'r')
xlim([0 0.5])

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

%%
e3_wfreqs = wfreqs;
e3_stim_ind = stim_ind_pfilt;
e3_nbins = nbins;
% for t = single_units

load ./extp1_stim_alpha_models
for t = 1:96
    subplot(2,1,1)
    pcolor(pax(1:end-1),e3_wfreqs,reshape(e3_stim_ind(t,:),e3_nbins,length(e3_wfreqs))');shading flat

     subplot(2,1,2)
    pcolor(pax(1:end-1),wfreqs,reshape(stim_ind_filt_pos(t,:),nbins,length(wfreqs))');shading flat
   input('')
    clf

end
