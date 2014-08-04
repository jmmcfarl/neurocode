%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
% cd ~/Data/bruce/7_15_12/G029/
cd ~/Data/bruce/7_15_12/G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt3_fixbased_data.mat
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
% use_lfps = [1:96];
use_lfps = [66];
alpha_smooth_win = 30;
% alpha_smooth_win = 25;

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
trial_stop_inds = trial_stop_sacinds;
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
full_alpha_avg_phase = [];
for ee = 1:n_expts
    fprintf('Expt %d of %d\n',ee,n_expts);
    clear ampgram all_V* alpha_phase
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
%         V_alpha_phase(:,ll) = angle(hilbert(V_alpha));
%         V_alpha_phase(:,ll) = unwrap_phase_monotonic(V_alpha_phase(:,ll));
%         V_alpha_phase(:,ll) = smooth(V_alpha_phase(:,ll),alpha_smooth_win,'lowess');
    end
    V_alpha_avg = mean(V_alpha_all,2);
    V_alpha_avg_phase = angle(hilbert(V_alpha));
    V_alpha_phase = unwrap_phase_monotonic(V_alpha_avg_phase);
    V_alpha_phase = smooth(V_alpha_phase,alpha_smooth_win,'lowess');
    
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
    
    fix_inds = find(full_expt_vec == Expt_nu(ee));
    cur_interp_phase = interp1(t_ax,V_alpha_phase,full_t(fix_inds));
    cur_interp_avg_phase = interp1(t_ax,V_alpha_avg_phase,full_t(fix_inds));
    full_alpha_phase = [full_alpha_phase; cur_interp_phase'];
%     full_alpha_phase = [full_alpha_phase; cur_interp_phase];
    full_alpha_avg_phase = [full_alpha_avg_phase; cur_interp_avg_phase'];
end

%% phase_since_trial = nan(size(full_alpha_phase));
for i = 1:n_trials
    cur_set = trial_start_inds(i):trial_stop_inds(i);
%     temp(i,:) = full_alpha_phase(cur_set(1),:);
    phase_since_trial(cur_set) = full_alpha_phase(cur_set)-full_alpha_phase(cur_set(1));
%     phase_since_trial(cur_set,:) = bsxfun(@minus,full_alpha_phase(cur_set,:),full_alpha_phase(cur_set(1),:));
end

% %relative to initial dip
% for i = 1:n_trials
%     cur_set = trial_start_inds(i):trial_stop_inds(i);
%     temp = find(diff(full_alpha_avg_phase(cur_set)) < -pi,1,'first');
%     if ~isempty(temp)
%         dip_ref(i) = 1 + temp;
%     else
%         dip_ref(i) = 1;
%     end
% %     phase_since_trial(cur_set,:) = bsxfun(@minus,full_alpha_phase(cur_set,:),full_alpha_phase(cur_set(dip_ref(i)),:));
%     phase_since_trial(cur_set) = full_alpha_phase(cur_set)-full_alpha_phase(cur_set(dip_ref(i)));
% end

%%
infix_ids = find(full_intrial==1);
xv_frac =0.;
xv_num = round(xv_frac*n_trials);
xv_f = randperm(n_trials);
xv_f(xv_num+1:end) = [];
tr_f = setdiff(1:n_trials,xv_f);
xv_inds = infix_ids(ismember(full_stim_ids(infix_ids),xv_f));
all_tr_inds = setdiff(infix_ids,xv_inds);

max_phase = 2*pi*5;
% min_phase = -2*pi;
min_phase = 0;
nbins = 62;
% nbins = 80;
tr_inds = all_tr_inds(phase_since_trial(all_tr_inds) <= max_phase & phase_since_trial(all_tr_inds) >= min_phase);
% tax = linspace(0,max_phase,nbins);
tax = linspace(min_phase,max_phase,nbins);
Tmat = tbrep(phase_since_trial,tax);
% Tmat = Tmat(tr_inds,:);

NT = length(tr_inds);
[un_expts,~,full_expt_inds] = unique(full_expt_vec);
n_un_expts = length(un_expts);
linX = zeros(NT,n_un_expts-1);
for i = 1:n_un_expts-1
    linX(full_expt_inds(tr_inds)==i,i) = 1;
end

%% COMPUTE OUTPUTS OF GABOR MODELS BASED ON RF MAPPING DATA

all_gabor_out1 = resh_all_stims*gabor_emp1_filt';
all_gabor_out2 = resh_all_stims*gabor_emp2_filt';
spatial_mod_out = zscore(sqrt(all_gabor_out1.^2 + all_gabor_out2.^2));

l2_ind = 500;
l2_dep = 500;
silent = 1;
% for t = 1:96
for t = 66
%     tr_inds = all_tr_inds(phase_since_trial(all_tr_inds,t) <= max_phase & phase_since_trial(all_tr_inds,t) >= min_phase);
%     tax = linspace(min_phase,max_phase,nbins);
%     Tmat = tbrep(phase_since_trial(:,t),tax);
%     linX = zeros(length(tr_inds),n_un_expts-1);
%     for i = 1:n_un_expts-1
%         linX(full_expt_vec(tr_inds)==i,i) = 1;
%     end
        
    fprintf('Fitting Cell %d of 96\n',t);
    cur_spatial_mod_out = spatial_mod_out(full_stim_ids(tr_inds),t);
    stim_dep_X = bsxfun(@times,Tmat(tr_inds,:),cur_spatial_mod_out);
    ov_Xmat = [stim_dep_X Tmat(tr_inds,:) linX];
    
    Robs = full_binned_spks(tr_inds,t);
    
    poss_spkcnts = unique(Robs(Robs > 0));
    spksN = [];
    for i = 1:length(poss_spkcnts)
        cur_set = find(Robs == poss_spkcnts(i));
        spksN = [spksN; repmat(cur_set,poss_spkcnts(i),1)];
    end
    spksN = sort(spksN);
    
    klen = size(ov_Xmat,2);
    K0 = zeros(klen,1);
    lamrange2 = [l2_dep 1 nbins 0;l2_ind (nbins+1) 2*nbins 0];
    llist = [];
    [fitp,grad] = GLMsolve_jmm(ov_Xmat, spksN, K0, silent, [], lamrange2, [], [],llist, [], 0);
    stim_dep_kern(t,:) = fitp.k(1:nbins);
    stim_ind_kern(t,:) = fitp.k((nbins+1):2*nbins);
    block_kern(t,:) = fitp.k((2*nbins+1):end-1);
    ov_const(t) = fitp.k(end);
    
    gfun = ov_Xmat*fitp.k(1:end-1) + fitp.k(end);
    too_large = find(gfun > 50);
    predrate = log(1+exp(gfun));
    predrate(too_large) = gfun(too_large);
    LL(t) = -sum(Robs.*log(predrate)-predrate)/sum(Robs);
        
    cur_lambda = 400;
    ov_Xmat = [Tmat(tr_inds,:) linX];
    lamrange2 = [cur_lambda 1 nbins 0];
    llist = [];
    [fitp,grad] = GLMsolve_jmm(ov_Xmat, spksN, K0, silent, [], lamrange2, [], [],llist, [], 0);
    trig_avg_kern = fitp.k(1:nbins);
    offset = fitp.k(end) + mean(fitp.k(nbins+1:end-1));
    trig_avg_rate(t,:) = log(1+exp(trig_avg_kern + offset));
    
%%     
n_prc_bins = 5;
dep_avg_rate(t,:,:) = zeros(n_prc_bins,nbins);
bin_edges = prctile(spatial_mod_out(:,t),linspace(0,100,n_prc_bins+1));
for ii = 1:n_prc_bins
    cur_lambda = 50;
    cur_set = find(spatial_mod_out >= bin_edges(ii) & spatial_mod_out < bin_edges(ii+1));
    cur_inds = find(ismember(full_stim_ids(tr_inds),cur_set));
    cur_spk_bns = convert_to_spikebins(Robs(cur_inds));
    
    ov_Xmat = [Tmat(tr_inds(cur_inds),:) linX(cur_inds,:)];
    lamrange2 = [cur_lambda 1 nbins 0];
    llist = [];
    [fitp,grad] = GLMsolve_jmm(ov_Xmat, cur_spk_bns, K0, silent, [], lamrange2, [], [],llist, [], 0);
    trig_avg_kern = fitp.k(1:nbins);
    offset = fitp.k(end) + mean(fitp.k(nbins+1:end-1));
    dep_avg_rate(t,ii,:) = log(1+exp(trig_avg_kern + offset));
end
end

%%
save expt3_unit_tempmods_alpha_v9 tax stim_* ov_const block_* LL trig_avg_rate nbins max_phase min_phase

%%
close all
load ./expt3_unit_tempmods_alpha
e2_block_kern = block_kern;
e2_const = ov_const;
e2_tax = linspace(0,max_phase,62)/(2*pi);
e2_stim_dep_kern = stim_dep_kern;
e2_stim_ind_kern = stim_ind_kern;
e2_trig_avg = trig_avg_rate/dt;
load ./expt3_unit_tempmods_alpha_indalpha
e3_block_kern = block_kern;
e3_const = ov_const;
e3_tax = tax/(2*pi);
e3_stim_dep_kern = stim_dep_kern;
e3_stim_ind_kern = stim_ind_kern;
e3_trig_avg = trig_avg_rate/dt;
load ./expt3_unit_tempmods_alpha_indalpha_relpeak_sm1.mat
for t = 1:96
    subplot(3,3,1)
    plot(tax/(2*pi),stim_dep_kern(t,:));
    xlim([-1 5])
    subplot(3,3,4)
    plot(tax/(2*pi),stim_ind_kern(t,:));
    xlim([-1 5])
    subplot(3,3,7)
    plot(tax/(2*pi),trig_avg_rate(t,:)/dt);
    xlim([-1 5])
     subplot(3,3,2)
    plot(e3_tax,e3_stim_dep_kern(t,:));
    xlim([0 5])
    subplot(3,3,5)
    plot(e3_tax,e3_stim_ind_kern(t,:));
    xlim([0 5])
    subplot(3,3,8)
    plot(e3_tax,e3_trig_avg(t,:)/dt);
    xlim([0 5])
     subplot(3,3,3)
    plot(e2_tax,e2_stim_dep_kern(t,:));
    xlim([0 5])
    subplot(3,3,6)
    plot(e2_tax,e2_stim_ind_kern(t,:));
    xlim([0 5])
    subplot(3,3,9)
    plot(e2_tax,e2_trig_avg(t,:)/dt);
    xlim([0 5])
   t
    pause
    clf
end


%%
close all
load ./expt3_unit_tempmods_v2
e3_block_kern = block_kern;
e3_const = ov_const;
e3_tents = tent_centers*dt;
e3_stim_dep_kern = stim_dep_kern;
e3_stim_ind_kern = stim_ind_kern;
e3_trig_avg = trig_avg_rate/dt;
d3_dt = dt;
% load ./expt3_unit_tempmods_alpha
load ./expt3_unit_tempmods_alpha_indalpha
for t = 1:96
    subplot(3,2,1)
    plot(tax/(2*pi),stim_dep_kern(t,:));
    xlim([0 5])
    subplot(3,2,3)
    plot(tax/(2*pi),stim_ind_kern(t,:));
    xlim([0 5])
    subplot(3,2,5)
    plot(tax/(2*pi),trig_avg_rate(t,:)/dt);
    xlim([0 5])
     subplot(3,2,2)
    plot(e3_tents,e3_stim_dep_kern(t,:));
    xlim([0 .5])
    subplot(3,2,4)
    plot(e3_tents,e3_stim_ind_kern(t,:));
    xlim([0 .5])
    subplot(3,2,6)
    plot(e3_tents,e3_trig_avg(t,:)/d3_dt);
    xlim([0 .5])
   t
    pause
    clf
end

%%
close all
load ./expt3_unit_tempmods_alpha_v2
old_stim_dep = stim_dep_kern;
old_stim_ind = stim_ind_kern;
load ./expt3_unit_tempmods_alpha_v7.mat

for t = 1:96
    subplot(2,2,1)
    plot(tax/(2*pi),stim_dep_kern(t,:));
    xlim([0 5])
    subplot(2,2,2)
    plot(tax/(2*pi),stim_ind_kern(t,:));
    xlim([0 5])
     subplot(2,2,3)
    plot(tax/(2*pi),old_stim_dep(t,:));
    xlim([0 5])
    subplot(2,2,4)
    plot(tax/(2*pi),old_stim_ind(t,:));
    xlim([0 5])
   t
    pause
    clf
end

%%
cd /home/james/Data/bruce/7_15_12/G034
close all
load ./expt3_unit_tempmods_v2
e3_block_kern = block_kern;
e3_const = ov_const;
e3_tents = tent_centers*dt;
e3_stim_dep_kern = stim_dep_kern;
e3_stim_ind_kern = stim_ind_kern;
e3_trig_avg = trig_avg_rate/dt;
d3_dt = dt;

load ./expt3_unit_tempmods_alpha_v9.mat
stim_outs = [-2 0 2];
cmap = jet(3);
% ucells = [9 18 21 27 29 35 37 50 52 54 55 61 62 64 71 77 83];
% ucells = [54];
ucells = 1:96;
% ucells = [3 27 35 67 70 71 76 83 93];
ucells = [9 17 30 34 71 77 66 78];
% ucells = [9 18 21 27 30 31 33 35 37 52 61 62 71 77 81 83 88 93];
% ucells = 71
cd ~/Desktop/   
for t = 1:length(ucells)
    cur_cell = ucells(t);
    
    phase_offset = mean(block_kern(cur_cell,:),2) + ov_const(cur_cell);
    phase_preds = repmat(stim_dep_kern(cur_cell,:),3,1);
    phase_preds = bsxfun(@times,phase_preds,stim_outs');
    phase_preds = bsxfun(@plus,phase_preds,stim_ind_kern(cur_cell,:));
    phase_preds = log(1+exp(phase_preds + ov_const(cur_cell)))/dt;

        e3_offset = mean(e3_block_kern(cur_cell,:),2) + ov_const(cur_cell);
    e3_preds = repmat(e3_stim_dep_kern(cur_cell,:),3,1);
    e3_preds = bsxfun(@times,e3_preds,stim_outs');
    e3_preds = bsxfun(@plus,e3_preds,e3_stim_ind_kern(cur_cell,:));
    e3_preds = log(1+exp(e3_preds + e3_const(cur_cell)))/dt;

    %     figure
    subplot(4,2,1)
    plot(tax/(2*pi),stim_dep_kern(cur_cell,:));
     grid on
   xlim([0 4.5])
    xlabel('Alpha cycles','fontsize',16)
    title('Stim-dependent','fontsize',16)
    subplot(4,2,3)
    plot(tax/(2*pi),stim_ind_kern(cur_cell,:));
     grid on
   xlim([0 4.5])
     title('Stim-Independent','fontsize',16)
   xlabel('Alpha cycles','fontsize',16)
%     subplot(3,2,5)
%     plot(tax/(2*pi),trig_avg_rate(cur_cell,:)/dt);
%     xlim([0 4.5])
%      title('Average Rate','fontsize',16)
%    xlabel('Alpha cycles','fontsize',16)
    subplot(4,2,7)
    plot(tax/(2*pi),phase_preds');
    grid on
    xlim([0 4.5])
    xlabel('Alpha cycles','fontsize',16)
    title('Predicted Responses','fontsize',16)
    
       subplot(4,2,2)
    plot(e3_tents,e3_stim_dep_kern(cur_cell,:));
    grid on
    xlim([0 0.45])
    xlabel('Time (s)','fontsize',16)
    title('Stim-dependent','fontsize',16)
    subplot(4,2,4)
    plot(e3_tents,e3_stim_ind_kern(cur_cell,:));
    xlim([0 0.45])
    xlabel('Time (s)','fontsize',16)
     title('Stim-Independent','fontsize',16)
     subplot(4,2,8)
    plot(e3_tents,e3_preds');
    grid on
    xlim([0 0.45])
    xlabel('Alpha cycles','fontsize',16)
    title('Predicted Responses','fontsize',16)

    subplot(4,2,6)
    plot(e3_tents,e3_trig_avg(cur_cell,:))
    xlim([0 0.45])
    subplot(4,2,5)
    plot(tax/(2*pi),trig_avg_rate(cur_cell,:)/dt)
    xlim([0 4.5])
%     fillPage(gcf,'PaperSize',[12 14])
%     figname = sprintf('New_Warped_Cell_%d',cur_cell);
%     print('-dpdf',figname);
%     close
   ucells(t)
    pause
    clf
end


%%
cd /home/james/Data/bruce/7_15_12/G034
close all
% load ./expt3_unit_tempmods_alpha_relpeak.mat
load ./expt3_unit_tempmods_alpha_v9.mat
stim_outs = [-2 0 2];
cmap = jet(3);
% ucells = [9 18 21 27 35 37 52 61 62 71 77];
ucells = [1:96];
for t = 1:length(ucells)
    cur_cell = ucells(t);
    
    e3_offset = mean(block_kern(cur_cell,:),2) + ov_const(t);
    e3_preds = repmat(stim_dep_kern(cur_cell,:),3,1);
    e3_preds = bsxfun(@times,e3_preds,stim_outs');
    e3_preds = bsxfun(@plus,e3_preds,stim_ind_kern(cur_cell,:));
    e3_preds = log(1+exp(e3_preds + ov_const(cur_cell)));
%     figure
    subplot(2,2,1)
    plot(tax/(2*pi),stim_dep_kern(cur_cell,:));
    title('Stim-dependent','fontsize',16)
    xlim([0 4.5])
    xlabel('Alpha cycles','fontsize',16)
    subplot(2,2,2)
    plot(tax/(2*pi),stim_ind_kern(cur_cell,:));
    title('Stim-Independent','fontsize',16)
    xlim([0 4.5])
    xlabel('Alpha cycles','fontsize',16)
    subplot(2,2,3)
    plot(tax/(2*pi),trig_avg_rate(cur_cell,:)/dt);
    xlim([0 4.5])
    title('Average Rate','fontsize',16)
    xlabel('Alpha cycles','fontsize',16)
    subplot(2,2,4)
    plot(tax/(2*pi),e3_preds');
    grid on
    xlim([0 4.5])
%     title('Predicted Responses','fontsize',16)
%     xlabel('Alpha cycles','fontsize',16)
%     figname = sprintf('Warped_Cell_%d',cur_cell);
%     print('-dpdf',figname);
%     close
   t
    pause
    clf
end

%%

%%
cd /home/james/Data/bruce/7_15_12/G034
close all
load ./expt3_unit_tempmods_v2
e3_block_kern = block_kern;
e3_const = ov_const;
e3_tents = tent_centers*dt;
e3_stim_dep_kern = stim_dep_kern;
e3_stim_ind_kern = stim_ind_kern;
e3_trig_avg = trig_avg_rate/dt;
d3_dt = dt;

load ./expt3_unit_tempmods_alpha_v9.mat
stim_outs = [-2 0 2];
cmap = jet(3);
% ucells = [9 18 21 27 29 35 37 50 52 54 55 61 62 64 71 77 83];
% ucells = [54];
% ucells = 1:96;
% ucells = [3 27 35 67 70 71 76 83 93];
% ucells = [9 17 30 34 71 77 66 78];
ucells = [9 18 21 27 30 31 33 35 37 52 61 62 71 77 81 83 88 93];
% ucells = 71
cd ~/Desktop/presentation_fig_files/   
for t = 1:length(ucells)
    cur_cell = ucells(t);
    
    phase_offset = mean(block_kern(cur_cell,:),2) + ov_const(cur_cell);
    phase_preds = repmat(stim_dep_kern(cur_cell,:),3,1);
    phase_preds = bsxfun(@times,phase_preds,stim_outs');
    phase_preds = bsxfun(@plus,phase_preds,stim_ind_kern(cur_cell,:));
    phase_preds = log(1+exp(phase_preds + ov_const(cur_cell)))/dt;

        e3_offset = mean(e3_block_kern(cur_cell,:),2) + ov_const(cur_cell);
    e3_preds = repmat(e3_stim_dep_kern(cur_cell,:),3,1);
    e3_preds = bsxfun(@times,e3_preds,stim_outs');
    e3_preds = bsxfun(@plus,e3_preds,e3_stim_ind_kern(cur_cell,:));
    e3_preds = log(1+exp(e3_preds + e3_const(cur_cell)))/dt;

    %     figure
    subplot(3,2,1)
     plot(tax/(2*pi),trig_avg_rate(cur_cell,:)/dt)
    xlim([0 4.5])
    box off
    subplot(3,2,3)
    plot(tax/(2*pi),stim_dep_kern(cur_cell,:));
   xlim([0 4.5])
    box off
    xlabel('Alpha cycles','fontsize',16)
    title('Stim-dependent','fontsize',16)
    subplot(3,2,5)
     plot(tax/(2*pi),phase_preds');
    xlim([0 4.5])
    xlabel('Alpha cycles','fontsize',16)
    title('Predicted Responses','fontsize',16)
    box off
   
    subplot(3,2,2)
        plot(e3_tents,e3_trig_avg(cur_cell,:))
    xlim([0 0.45])
    box off
        subplot(3,2,4)
    plot(e3_tents,e3_stim_dep_kern(cur_cell,:));
    xlim([0 0.45])
    xlabel('Time (s)','fontsize',16)
    title('Stim-dependent','fontsize',16)
     box off
    subplot(3,2,6)
    plot(e3_tents,e3_preds');
    xlim([0 0.45])
    xlabel('Alpha cycles','fontsize',16)
    title('Predicted Responses','fontsize',16)
    box off

    fillPage(gcf,'PaperSize',[12 14])
    figname = sprintf('New_Warped_Cell_%d',cur_cell);
    print('-dpdf',figname);
    close
   ucells(t)
%     pause
%     clf
end

