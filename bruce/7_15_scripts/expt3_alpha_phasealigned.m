%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
% cd ~/Data/bruce/7_15_12/G029/
cd ~/Data/bruce/7_15_12/G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt3_fixbased_data_buffer.mat
resh_all_stims = resh_all_stims/std(resh_all_stims(:));
used_resh_stims = resh_all_stims(trial_stimnum,:);

fix_win_dur = 0.15;

Pix2Deg = 0.018837;
NT = length(full_stim_ids);
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

Fs = 3e4;
dsf = 200;Fsd = Fs/dsf;
[b_alpha,a_alpha] = butter(2,[5 13]/(Fsd/2));
use_lfps = [1:96];
% alpha_smooth_win = 30;
alpha_smooth_win = 25;

%%
n_trials = length(trial_start_inds);
trial_stop_sacinds = trial_stop_inds;
for i = 1:n_trials
    cur_inds = trial_start_inds(i):trial_stop_inds(i);
    cur_sac_stops = find(full_insac(cur_inds) == 1,1,'first');
    if ~isempty(cur_sac_stops)
        cur_sac_stops = cur_inds(cur_sac_stops);
        trial_stop_sacinds(i) = cur_sac_stops;
    end
end
% trial_stop_inds = trial_stop_sacinds;
trial_start_times = full_t(trial_start_inds);
trial_stop_times = full_t(trial_stop_inds);
trial_durs = trial_stop_times - trial_start_times;

min_trial_dur = 0.4;
used_trials = find(trial_durs >= min_trial_dur);
trial_start_inds = trial_start_inds(used_trials);
trial_stop_inds = trial_stop_inds(used_trials);
trial_imnum = trial_imnum(used_trials);
trial_stimnum = trial_stimnum(used_trials);
trial_start_times = trial_start_times(used_trials);
trial_stop_times = trial_stop_times(used_trials);
trial_stop_winds = trial_start_inds + round(min_trial_dur/dt);
trial_stop_wtimes = full_t(trial_stop_winds);

n_trials = length(trial_start_inds);
fix_expt_num = nan(n_trials,1);
full_intrial = zeros(size(full_t));
for i = 1:n_trials
    cur_inds = trial_start_inds(i):trial_stop_winds(i);
    fix_expt_num(i) = unique(full_expt_vec(cur_inds));
    cur_inds = trial_start_inds(i):trial_stop_inds(i);
    full_intrial(cur_inds) = 1;
end

[un_expts,~,full_expt_inds] = unique(fix_expt_num);
n_un_expts = length(un_expts);


%%
% load ./expt1_eyecor_d1p25_nosac_v2.mat gabor*
load ./gabor_tracking_varmeans.mat gabor*
gabor_params = gabor_params_f{end};
clear gabor_*filt
for t = 1:96
    
    gabor_emp1 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XX,YY,gabor_params(t,1:6),pi/2);
    
    gabor_emp1_filt(t,:) = gabor_emp1(:);
    gabor_emp2_filt(t,:) = gabor_emp2(:);
end

%%
Expt_nu = [7:12]; %expt 3 34
n_expts = length(Expt_nu);
full_alpha_phase = [];
full_alpha_uphase = [];
full_alpha_avg_phase = [];
for ee = 1:n_expts
    fprintf('Expt %d of %d\n',ee,n_expts);
    clear ampgram all_V* V_alpha*
    filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(1));
    load(filename);
    V = double(FullV.V);
    V = decimate(V,dsf);
    V = V(:);
    vlen = length(V);
    V_alpha_all = nan(vlen,length(use_lfps));
    V_alpha_phase = nan(vlen,length(use_lfps));
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d\n',ll);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = decimate(V,dsf);
        V = V(:);
        V_alpha = filtfilt(b_alpha,a_alpha,V);
        V_alpha_all(:,ll) = V_alpha;
        V_alpha_phase(:,ll) = angle(hilbert(V_alpha));
        V_alpha_uphase(:,ll) = unwrap_phase_monotonic(V_alpha_phase(:,ll));
        V_alpha_uphase(:,ll) = smooth(V_alpha_uphase(:,ll),alpha_smooth_win,'lowess');
    end
%     V_alpha_avg = mean(V_alpha_all,2);
%     V_alpha_avg_phase = angle(hilbert(V_alpha));
    %     V_alpha_phase = unwrap_phase_monotonic(V_alpha_avg_phase);
    %     V_alpha_phase = smooth(V_alpha_phase,alpha_smooth_win,'lowess');
    
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
    
    fix_inds = find(full_expt_vec == Expt_nu(ee));
    cur_interp_phase = interp1(t_ax,V_alpha_phase,full_t(fix_inds));
    cur_interp_uphase = interp1(t_ax,V_alpha_uphase,full_t(fix_inds));
%     cur_interp_avg_phase = interp1(t_ax,V_alpha_avg_phase,full_t(fix_inds));
    full_alpha_phase = [full_alpha_phase; cur_interp_phase];
    %     full_alpha_phase = [full_alpha_phase; cur_interp_phase];
    full_alpha_uphase = [full_alpha_uphase; cur_interp_uphase];
%     full_alpha_avg_phase = [full_alpha_avg_phase; cur_interp_avg_phase'];
end

%% phase_since_trial = nan(size(full_alpha_phase));
% for i = 1:n_trials
%     cur_set = trial_start_inds(i):trial_stop_inds(i);
%     %     temp(i,:) = full_alpha_phase(cur_set(1),:);
%     %     phase_since_trial(cur_set) = full_alpha_phase(cur_set)-full_alpha_phase(cur_set(1));
%     phase_since_trial(cur_set,:) = bsxfun(@minus,full_alpha_uphase(cur_set,:),full_alpha_uphase(cur_set(1),:));
% end

for ll = 1:length(use_lfps)
    zero_crossing_inds{ll} = find(full_alpha_phase(1:end-1,ll) < 0 & full_alpha_phase(2:end,ll) > 0);
end
phase_deriv = [zeros(1,length(use_lfps)); diff(full_alpha_phase)];
for ll = 1:length(use_lfps)
    pi_crossing_inds{ll} = find(phase_deriv(1:end-1,ll) > -1 & phase_deriv(2:end,ll) < -1);
end

buffer_t = 0.2;
buffer_inds = round(buffer_t/dt);
phase_start_inds = zeros(length(trial_start_inds),length(use_lfps));
phase_since_trial = nan(length(trial_start_inds),length(use_lfps));
for i = 1:n_trials
    for ll = 1:length(use_lfps)
%         cur_zdiffs = abs(zero_crossing_inds{ll} - trial_start_inds(i));
%         [~,b] = min(cur_zdiffs);
%         phase_start_inds(i,ll) = zero_crossing_inds{ll}(b);
        next_pcross = find(pi_crossing_inds{ll} >= trial_start_inds(i),1,'first');
        phase_start_inds(i,ll) = pi_crossing_inds{ll}(next_pcross);
        cur_set = phase_start_inds(i,ll):trial_stop_inds(i);
        phase_since_trial(cur_set,ll) = full_alpha_uphase(cur_set,ll) - full_alpha_uphase(cur_set(1),ll);
    end
end

full_unique_trial = nan(length(full_t),length(use_lfps));
for i = 1:length(trial_start_inds)
    for ll = 1:length(use_lfps)
        cur_set = phase_start_inds(i,ll):trial_stop_inds(i);
        full_unique_trial(cur_set,ll) = i;
    end
end
%%
% infix_ids = find(full_intrial==1);
% xv_frac =0.;
% xv_num = round(xv_frac*n_trials);
% xv_f = randperm(n_trials);
% xv_f(xv_num+1:end) = [];
% tr_f = setdiff(1:n_trials,xv_f);
% xv_inds = infix_ids(ismember(full_stim_ids(infix_ids),xv_f));
% all_tr_inds = setdiff(infix_ids,xv_inds);

used_trials = 1:length(trial_start_inds);
n_used_trials = length(used_trials);
xv_frac = 0.2;
n_xv_trials = round(n_used_trials*xv_frac);
xv_set = randperm(n_used_trials);
tr_set = xv_set(n_xv_trials+1:end);
xv_set(n_xv_trials+1:end) = [];
for ll = 1:length(use_lfps)
    xv_inds{ll} = find(ismember(full_unique_trial(:,ll),xv_set));
    all_tr_inds{ll} = find(ismember(full_unique_trial(:,ll),tr_set));
end

%%
max_phase = 2*pi*5;
% min_phase = -2*pi;
min_phase = -pi;
% nbins = 62;
nbins = 80;
% tr_inds = all_tr_inds(phase_since_trial(all_tr_inds) <= max_phase & phase_since_trial(all_tr_inds) >= min_phase);
% tax = linspace(0,max_phase,nbins);
tax = linspace(min_phase,max_phase,nbins);
tax_p = tax/(2*pi);
% Tmat = tbrep(phase_since_trial,tax);
% Tmat = Tmat(tr_inds,:);

% NT = length(tr_inds);
[un_expts,~,full_expt_inds] = unique(full_expt_vec);
n_un_expts = length(un_expts);
% linX = zeros(NT,n_un_expts-1);
% for i = 1:n_un_expts-1
%     linX(full_expt_inds(tr_inds)==i,i) = 1;
% end

beg_dur = round(0.15/dt);
late_indicator = zeros(size(full_t),length(use_lfps));
for i = 1:length(trial_start_inds)
    for ll = 1:length(use_lfps)
   cur_inds = (beg_dur+phase_start_inds(i,ll)):trial_stop_inds(i);
   late_indicator(cur_inds,ll) = 1;
    end
end

%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA

all_gabor_out1 = resh_all_stims*gabor_emp1_filt';
all_gabor_out2 = resh_all_stims*gabor_emp2_filt';
% spatial_mod_out = zscore(sqrt(all_gabor_out1.^2 + all_gabor_out2.^2));
spatial_mod_out = sqrt(all_gabor_out1.^2 + all_gabor_out2.^2);
spatial_mod_out = bsxfun(@rdivide,spatial_mod_out,std(spatial_mod_out));

reg_params.dl1_ind = 500;
reg_params.dl1_dep = 3000;
reg_params.dl2_ind = 5000;
reg_params.dl2_dep = 20000;
reg_params.dl2_freq_ind = 0;
reg_params.dl2_freq_dep = 0;
reg_params.dl2_time_ind = 1000;
reg_params.dl2_time_dep = 1000;
reg_params.l1_ind = 2;
reg_params.l1_dep =2;
reg_params.l1t_ind = 1;
reg_params.is_phase = 1;

silent = 1;
npbins = 30;
pax = linspace(-pi,pi,npbins+1);

NL_type = 1;
l2_ind = 1000;
l2_dep = 1000;
l1_ind = 0.4;
l1_dep = 0.4;
silent = 1;
for t = 1:96
    
    tr_inds = all_tr_inds{t}(phase_since_trial(all_tr_inds{t},t) <= max_phase & phase_since_trial(all_tr_inds{t},t) >= min_phase);
    tax = linspace(min_phase,max_phase,nbins);
    Tmat = tbrep(phase_since_trial(:,t),tax);
    
    fprintf('Fitting Cell %d of 96\n',t);
    cur_spatial_mod_out = zeros(length(tr_inds),1);
    uset = find(~isnan(full_stim_ids(tr_inds)));
    cur_spatial_mod_out(uset) = spatial_mod_out(full_stim_ids(tr_inds(uset)),t);
    
    stim_dep_X = bsxfun(@times,Tmat(tr_inds,:),cur_spatial_mod_out);
    ov_Xmat = [stim_dep_X Tmat(tr_inds,:)];
    
    Robs = full_binned_spks(tr_inds,t);
    spksN = convert_to_spikebins(Robs);
    
    klen = size(ov_Xmat,2);
    K0 = zeros(klen+1,1);
    lamrange2 = [l2_dep 1 nbins 0;l2_ind (nbins+1) 2*nbins 0];
    llist = [l1_dep 1:2*nbins];
    [fitp_po,grad] = GLMsolve_jmm(ov_Xmat, spksN, K0, silent, [], lamrange2, [], [],llist, [], NL_type);
    stim_po_dep_kern(t,:) = fitp_po.k(1:nbins);
    stim_po_ind_kern(t,:) = fitp_po.k((nbins+1):2*nbins);
    
    cur_lambda = l2_ind;
    ov_Xmat = [Tmat(tr_inds,:)];
    klen = size(ov_Xmat,2);
    K0 = zeros(klen+1,1);
    lamrange2 = [cur_lambda 1 nbins 0];
    llist = [l1_ind 1:nbins];
    [fitp_avg,grad] = GLMsolve_jmm(ov_Xmat, spksN, K0, silent, [], lamrange2, [], [],llist, [], NL_type);
    trig_avg_kern(t,:) = fitp_avg.k(1:nbins);
    offset_avg(t) = fitp_avg.k(end);
    if NL_type == 0
        trig_avg_rate(t,:) = log(1+exp(trig_avg_kern(t,:) + offset_avg(t)));
    else
        trig_avg_rate(t,:) = exp(trig_avg_kern(t,:) + offset_avg(t));
    end
        
    
    tr_late_indicator = logical(late_indicator(tr_inds,t));
    cur_alpha_phase = squeeze(full_alpha_phase(tr_inds,t));
    
%     late_spks = find(tr_late_indicator(spksN));
%     spk_phases = cur_alpha_phase(spksN(late_spks));
%     mean_phase(t) = circ_mean(spk_phases);
%     phase_term = cos(cur_alpha_phase-mean_phase(t)).*tr_late_indicator;
%     ov_Xmat = [Tmat(tr_inds,:) phase_term];
%     klen = size(ov_Xmat,2);
%     K0 = zeros(klen+1,1);
%     lamrange2 = [cur_lambda 1 nbins 0];
%     llist = [l1_ind 1:nbins];
%     [fitp_np,grad] = GLMsolve_jmm(ov_Xmat, spksN, K0, silent, [], lamrange2, [], [],llist, [], NL_type);
%     trig_np_kern(t,:) = fitp_np.k(1:nbins);
%     phase_dep_np(t) = fitp_np.k(end-1);
%     offset_avg(t) = fitp_np.k(end);
 
%     Pmat = [];
%     cur_tb = tbrep(cur_alpha_phase,pax);
%     cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
%     cur_tb(:,end) = [];
%     Pmat = [Pmat cur_tb];
% %     aPmat = Pmat;
%     aPmat = bsxfun(@times,Pmat,tr_late_indicator);
% 
%     Xmat = [aPmat];
%     stim_params = [npbins,1];
%     klen = size(Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp] = fit_GLM_phase_model(Xmat, Robs, K0,1, stim_params, reg_params,0,NL_type);
%     stim_ponly_pfilt(t,:) = fitp.k(1:npbins);
%     phase_mod_out = aPmat*stim_ponly_pfilt(t,:)';
%     phase_mod_out(tr_late_indicator==1) = zscore(phase_mod_out(tr_late_indicator==1));
%     
%     ov_Xmat = [Tmat(tr_inds,:) phase_mod_out];
%     klen = size(ov_Xmat,2);
%     K0 = zeros(klen+1,1);
%     lamrange2 = [cur_lambda 1 nbins 0];
%     llist = [l1_ind 1:nbins];
%     [fitp_np,grad] = GLMsolve_jmm(ov_Xmat, spksN, K0, silent, [], lamrange2, [], [],llist, [], NL_type);
%     trig_np_kern(t,:) = fitp_np.k(1:nbins);
%     phase_dep_np(t) = fitp_np.k(end-1);
%     offset_np(t) = fitp_np.k(end);
% 
%     ov_Xmat = [stim_dep_X Tmat(tr_inds,:) phase_mod_out];
%     klen = size(ov_Xmat,2);
%     K0 = zeros(klen+1,1);
%     lamrange2 = [l2_dep 1 nbins 0;l2_ind (nbins+1) 2*nbins 0];
%     llist = [l1_dep 1:2*nbins];
%     [fitp_pop,grad] = GLMsolve_jmm(ov_Xmat, spksN, K0, silent, [], lamrange2, [], [],llist, [], NL_type);
%     stim_pop_dep_kern(t,:) = fitp_pop.k(1:nbins);
%     stim_pop_ind_kern(t,:) = fitp_pop.k((nbins+1):2*nbins);
%   
%         ov_Xmat = [stim_dep_X Tmat(tr_inds,:) phase_mod_out phase_mod_out.*cur_spatial_mod_out];
%     klen = size(ov_Xmat,2);
%     K0 = zeros(klen+1,1);
%     lamrange2 = [l2_dep 1 nbins 0;l2_ind (nbins+1) 2*nbins 0];
%     llist = [l1_dep 1:2*nbins];
%     [fitp_pops,grad] = GLMsolve_jmm(ov_Xmat, spksN, K0, silent, [], lamrange2, [], [],llist, [], NL_type);
%     stim_pops_dep_kern(t,:) = fitp_pops.k(1:nbins);
%     stim_pops_ind_kern(t,:) = fitp_pops.k((nbins+1):2*nbins);
    
%     Xmat = [aPmat Tmat(tr_inds,:)];
%     stim_params = [npbins,1,nbins,1];
%     klen = size(Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp] = fit_GLM_phase_model(Xmat, Robs, K0,1, stim_params, reg_params,1,NL_type);
%     stim_ind_pfilt(t,:) = fitp.k(1:npbins);
%     stim_ind_tfilt(t,:) = fitp.k((npbins+1):(npbins+nbins));
%     offset(t) = fitp.k(end);
% %     stim_dep_tfilt(t,:) = fitp.k((npbins+nbins+1):end-1);
     
%     Xmat = [aPmat Tmat(tr_inds,:) stim_dep_X];
%     stim_params = [npbins,1,nbins,1];
%     klen = size(Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp] = fit_GLM_phase_model(Xmat, Robs, K0,1, stim_params, reg_params,1,NL_type);
%     stim_ind_pfilt(t,:) = fitp.k(1:npbins);
%     stim_ind_tfilt(t,:) = fitp.k((npbins+1):(npbins+nbins));
%     stim_dep_tfilt(t,:) = fitp.k((npbins+nbins+1):end-1);

    
%     cur_spatial_mod_out = zeros(length(xv_inds{t}),1);
%     uset = find(~isnan(full_stim_ids(xv_inds{t})));
%     cur_spatial_mod_out(uset) = spatial_mod_out(full_stim_ids(xv_inds{t}(uset)),t);
% 
%     stim_dep_X = bsxfun(@times,Tmat(xv_inds{t},:),cur_spatial_mod_out);
%     xv_Robs = full_binned_spks(xv_inds{t},t);
%     spksN = convert_to_spikebins(xv_Robs);
%     
%     ov_Xmat = [stim_dep_X Tmat(xv_inds{t},:)];
%     xv_po_pred_rate = ov_Xmat*fitp_po.k(1:end-1) + fitp_po.k(end);
%     if NL_type == 0
%         xv_po_pred_rate = log(1+exp(xv_po_pred_rate));
%     else
%         xv_po_pred_rate = exp(xv_po_pred_rate);
%     end
%     xv_po_LL(t) = -sum(xv_Robs.*log(xv_po_pred_rate) - xv_po_pred_rate)/sum(xv_Robs);
%     
%     ov_Xmat = [Tmat(xv_inds{t},:)];
%     xv_avg_pred_rate = ov_Xmat*fitp_avg.k(1:end-1) + fitp_avg.k(end);
%     if NL_type == 0
%         xv_avg_pred_rate = log(1+exp(xv_avg_pred_rate));
%     else
%         xv_avg_pred_rate = exp(xv_avg_pred_rate);
%     end
%     xv_avg_LL(t) = -sum(xv_Robs.*log(xv_avg_pred_rate) - xv_avg_pred_rate)/sum(xv_Robs);
   
%     
%     xv_late_indicator = logical(late_indicator(xv_inds{t},t));
%     cur_alpha_phase = squeeze(full_alpha_phase(xv_inds{t},t));
%     Pmat = [];
%     cur_tb = tbrep(cur_alpha_phase,pax);
%     cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
%     cur_tb(:,end) = [];
%     Pmat = [Pmat cur_tb];
%     %     aPmat = Pmat;
%     aPmat = bsxfun(@times,Pmat,xv_late_indicator);
%     
%     
% %     Xmat = [aPmat Tmat(xv_inds{t},:) stim_dep_X];
%     Xmat = [aPmat Tmat(xv_inds{t},:)];
%     xv_pred_rate = Xmat*fitp.k(1:end-1) + fitp.k(end);
%     if NL_type == 0
%         xv_pred_rate = log(1+exp(xv_pred_rate));
%     else
%         xv_pred_rate = exp(xv_pred_rate);
%     end
%     xv_LL(t) = -sum(xv_Robs.*log(xv_pred_rate) - xv_pred_rate)/sum(xv_Robs);

    
end

%%
ucells = single_units;
rel_po_LL = -(xv_po_LL(ucells) - xv_avg_LL(ucells));
% rel_full_LL = -(xv_full_LL(ucells) - xv_avg_LL(ucells));
% rel_LL = -(xv_LL(ucells) - xv_avg_LL(ucells));
% figure
% boxplot([rel_full_LL(:) rel_LL(:) rel_po_LL(:)])

%%
close all
for t = 1:96
    plot(tax_p,trig_avg_kern(t,:)+offset_avg(t),'b')
    hold on
    plot(tax_p,trig_np_kern(t,:)+offset_np(t),'r')
    plot(tax_p,stim_ind_tfilt(t,:)+offset(t),'k')
    t
    pause
    clf
end

%%
close all
for t = 1:96
    subplot(2,1,1)
    plot(tax_p,stim_po_ind_kern(t,:),'r')
    hold on
    plot(tax_p,trig_avg_kern(t,:),'b')
%     plot(tax_p,stim_pop_ind_kern(t,:),'k')
    subplot(2,1,2)
    plot(tax_p,stim_po_dep_kern(t,:),'r')
    hold on
%     plot(tax_p,stim_pop_dep_kern(t,:),'k')
%      plot(tax_p,stim_pops_dep_kern(t,:),'g')
   t
    pause
    clf
end

%%
ee=1;
full_tsince_tstart = nan(size(full_t));
for i = 1:length(trial_start_inds)
    cur_set = phase_start_inds(i,ee):trial_stop_inds(i);
    cur_ts = full_t(cur_set)-trial_start_times(i);
    full_tsince_tstart(cur_set) = cur_ts;
end
infix_ids = find(~isnan(full_tsince_tstart));
%%

clear phase_dist
% for ee =1:96;
    cur_alpha_phase = squeeze(full_alpha_phase(:,ee));
    ee
    phase_bin_edges = linspace(-pi,pi,25);
    t_bin_edges = linspace(0,0.5,31);
    for t = 1:length(t_bin_edges)-1
        cur_set = infix_ids(find(full_tsince_tstart(infix_ids) >= t_bin_edges(t) & ...
            full_tsince_tstart(infix_ids) < t_bin_edges(t+1)));
        cur_hist = histc(cur_alpha_phase(cur_set),phase_bin_edges);
        cur_hist(end) = [];
        cur_hist = cur_hist/sum(cur_hist);
        phase_dist(ee,t,:) = cur_hist;
        [pval(ee,t),zval(ee,t)] = circ_rtest(cur_alpha_phase(cur_set));
        avgr(ee,t) = circ_r(cur_alpha_phase(cur_set));
    end
% end

%%
load ./abs_temp.mat
close all
for t = 1:96
    subplot(2,2,1)
    plot(tax_p,trig_avg_rate(t,:)/dt);
    hold on
    plot(tax_p,abs_trig_avg_rate(t,:)/dt,'r')
    subplot(2,2,2)
    plot(tax_p,trig_avg_kern(t,:))
hold on
    plot(tax_p,abs_trig_avg_kern(t,:),'r')
    subplot(2,2,3)
    plot(tax_p,stim_po_ind_kern(t,:))
    hold on
    plot(tax_p,abs_stim_po_ind_kern(t,:),'r')
    subplot(2,2,4)
    plot(tax_p,stim_po_dep_kern(t,:))
    hold on
    plot(tax_p,abs_stim_po_dep_kern(t,:),'r')
    t
    pause
    clf
end
