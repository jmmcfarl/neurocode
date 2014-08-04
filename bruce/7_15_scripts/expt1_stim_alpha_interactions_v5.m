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

load ./Expt1_newcompiled_data_fixeddelay_d1p25_34.mat full_t
old_full_t = full_t;

load ./Expt1_newcompiled_data_fixeddelay_d1p25_34_ht.mat
fullX = fullX/std(fullX(:));

Pix2Deg = 0.018837;
[Nims,klen] = size(fullX);
NT = length(full_im_ids);

%% crop stimulus for the purpose of faster gabor function fitting
% new_RF_patch = [-0.11 0.9; -0.9 0.1]; %location of RFs in degrees [x1 x2;y1 y2]
new_RF_patch = [-0.4 1.2; -1.2 0.42]; %location of RFs in degrees [x1 x2;y1 y2]
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

%%
cd ~/Data/bruce/7_15_12/G034/
load ./gabor_tracking_varmeans

load ./expt1_lfp_alpha_phase_delay_all_sameelec.mat all_interp_phases nearest_lfps new_* desired_dt wfreqs frac
% load ./expt1_lfp_alpha_phase_delay_all_nearestelec
% load ./expt1_lfp_alpha_phase_delay_all.mat

%%
% frac = 4;
cur_dt = median(diff(full_t));
desired_dt = cur_dt/frac;
trial_start_inds = 1+find(diff(new_trial_vec) ~= 0);
trial_stop_inds = trial_start_inds - 1;
trial_start_inds = [1; trial_start_inds(:)];
trial_stop_inds = [trial_stop_inds(:); length(new_trial_vec)];

n_trials = length(trial_stop_inds);
new_binned_spks = nan(length(new_t_axis),96);
Expt_nu = unique(new_expt_vec);
n_expts = length(Expt_nu);
% new_trial_vec = nan(length(new_t_axis),1);
new_early_inds = nan(length(new_t_axis),1);
for tt = 1:n_expts
    fprintf('Expt %d of %d\n',tt,n_expts);
    cur_trial_set = find(new_expt_vec(trial_start_inds) == Expt_nu(tt));
    cur_n_trials = length(cur_trial_set);
    load(sprintf('Expt%dClusterTimes.mat',Expt_nu(tt)));
    for i = 1:cur_n_trials
        cur_inds = trial_start_inds(cur_trial_set(i)):trial_stop_inds(cur_trial_set(i));
        cur_binned = zeros(length(cur_inds),96);
        for c = 1:96
            cur_binned(:,c) = histc(Clusters{c}.times,new_t_axis(cur_inds));
        end
        new_binned_spks(cur_inds,:) = cur_binned;
        %         new_trial_vec(cur_inds) = cur_trial_set(i);
        cur_inds((round(0.125/desired_dt)+1):end) = [];
        new_early_inds(cur_inds) = 1;
    end
end

im_inds_up = round(interp1(full_t,full_im_ids,new_t_axis));

nan_spks = max(isnan(new_binned_spks),[],2);
uset = find(nan_spks'==0 | ~isnan(im_inds_up));

im_inds_up(isnan(im_inds_up)) = 1;
%% RECONSTRUCT NEW STIMULUS MATRIX
im_flip_inds = [1; 1+find(diff(full_im_ids)~=0)];
im_flip_cents = 1 + im_flip_inds;
im_flip_times = full_t(im_flip_cents);
x_cor_interp = round(interp1(old_full_t,-x_cor{end},im_flip_times));
y_cor_interp = round(interp1(old_full_t,-y_cor{end},im_flip_times));
x_cor_interp(isnan(x_cor_interp)) = 0;
y_cor_interp(isnan(y_cor_interp)) = 0;

resh_X = reshape(fullX',[sdim sdim Nims]);
resh_X_sh = zeros(size(resh_X));
for ii = 1:Nims
    d2 = dist_shift2d(resh_X(:,:,ii), x_cor_interp(ii), 2,0);
    d2 = dist_shift2d(d2,y_cor_interp(ii),1,0);
    resh_X_sh(:,:,ii) = d2;
end

%% REFIT GEM PARAMETERS
fullX_sh = reshape(resh_X_sh,sdim^2,Nims)';
fullX_sh_cropped = fullX_sh(:,new_crop);
fullX_cropped = fullX(:,new_crop);
clear fullX fullX_sh resh_X_sh resh_X

%%
used_trials = unique(new_trial_vec(uset));
n_used_trials = length(used_trials);
xv_frac = 0.2;
n_xv_trials = round(n_used_trials*xv_frac);
xv_set = randperm(n_used_trials);
xv_set(n_xv_trials+1:end) = [];
xv_inds = find(ismember(new_trial_vec(uset),xv_set));
tr_inds = find(~ismember(uset,xv_inds))';

early_inds = find(new_early_inds==1);
xv_inds(ismember(xv_inds,early_inds)) = [];
tr_inds(ismember(tr_inds,early_inds));
%%
% tent_centers = [0:desired_dt:0.125];
tent_centers = [0:.006:0.125];
tent_centers = round(tent_centers/desired_dt);
tbmat = construct_tent_bases(tent_centers,1);
[ntents,tblen] = size(tbmat);

% tbmat = [tbmat zeros(ntents,tblen-1)];
tbmat = [zeros(ntents,tblen-1) tbmat ];
tbmat = bsxfun(@rdivide,tbmat,sum(tbmat,2));

%%
close all
NL_type = 1; %exp
% NL_type = 0; %logexp

reg_params.dl1_ind = 4000;
% reg_params.dl1_dep =20000;
reg_params.dl1_dep =4000;
reg_params.dl2_ind = 8000;
% reg_params.dl2_dep = 40000;
reg_params.dl2_dep = 8000;
reg_params.dl2_freq_ind = 200;
reg_params.dl2_freq_dep = 200;
reg_params.l1_dep = 0;
reg_params.l1_ind = 0;
reg_params.is_phase = 1;
% 
% reg_params.dl1_ind = 4000;
% % reg_params.dl1_dep =20000;
% reg_params.dl1_dep =4000;
% reg_params.dl2_ind = 8000;
% % reg_params.dl2_dep = 40000;
% reg_params.dl2_dep = 8000;
% reg_params.dl2_freq_ind = 100;
% reg_params.dl2_freq_dep = 100;
% reg_params.l1_dep = 0;
% reg_params.l1_ind = 0;
% reg_params.is_phase = 1;

reg_params2 = reg_params;
reg_params2.dl1_ind = reg_params2.dl1_dep;
reg_params2.dl2_ind = reg_params2.dl2_dep;
reg_params2.dl2_freq_ind = reg_params2.dl2_freq_dep;

silent =1;
nbins = 30;
pax = linspace(-pi,pi,nbins+1);

stim_params = [nbins,length(wfreqs)];
doub_pax = [pax(1:end-1) 2*pi+pax(1:end-1)];

prctile_bins = [1 5 10 25 50 75 90 95 99];
mid_prctile = find(prctile_bins==50);

NT_new = length(new_t_axis);
for t = 1:96
    fprintf('Fitting GEM: Cell %d of %d\n',t,96);
    
    cur_tr_inds = tr_inds;
    
    cur_trial_vec = new_trial_vec(uset(cur_tr_inds));
    cur_trial_set = unique(cur_trial_vec);
    rperm = randperm(length(cur_trial_set));
    rand_trial_order = cur_trial_set(rperm);
    rand_tr_inds = nan(size(cur_trial_vec));
    for i = 1:length(rperm)
       cur_from_set = find(cur_trial_vec == cur_trial_set(i));
       cur_to_set = find(cur_trial_vec == rand_trial_order(i));
       mlen = min(length(cur_from_set),length(cur_to_set));
       rand_tr_inds(cur_from_set(1:mlen)) = cur_tr_inds(cur_to_set(1:mlen));
    end
     cur_tr_inds(isnan(rand_tr_inds)) = [];
   rand_tr_inds(isnan(rand_tr_inds)) = [];
    
    Robs = new_binned_spks(uset(cur_tr_inds),t);    
%     Robs = new_binned_spks(uset(rand_tr_inds),t);

    spkbs = convert_to_spikebins(Robs);
    
    gabor_emp1 = get_pgabor_mask_v2(XXc,YYc,gabor_params_f{end}(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XXc,YYc,gabor_params_f{end}(t,1:6),pi/2);
    %     mask1_out = fullX_cropped*gabor_emp1(:);
    %     mask2_out = fullX_cropped*gabor_emp2(:);
    mask1_out = fullX_sh_cropped*gabor_emp1(:);
    mask2_out = fullX_sh_cropped*gabor_emp2(:);
    
    tot_mod_out =sqrt(mask1_out.^2 + mask2_out.^2);
    tot_mod_out = tot_mod_out/std(tot_mod_out);
    
    tot_mod_out = tot_mod_out(im_inds_up);
    
    trial_Tmat = zeros(length(new_t_axis),ntents);
    for i = 1:ntents
        trial_Tmat(:,i) = conv(tot_mod_out,tbmat(i,:),'same');
    end
    %     trial_Tmat = fliplr(trial_Tmat);
    cur_trial_Tmat = trial_Tmat(uset(cur_tr_inds),:);
    
    cur_phases = squeeze(all_interp_phases(uset(cur_tr_inds),:,nearest_lfps(t)));
    
    klen = ntents;
    K0 = zeros(klen+1,1);
%     lamrange2 = [2000 1 ntents 0];
    lamrange2 = [5000 1 ntents 0];
    [fitp_so,grad] = GLMsolve_jmm(cur_trial_Tmat,spkbs,K0,silent,[],lamrange2,[],[],[],[],NL_type);
    stim_kern(t,:) = fitp_so.k(1:end-1);
    
    cur_stim_mod_out = cur_trial_Tmat*stim_kern(t,:)';
    
    Pmat = [];
    for ww = 1:length(wfreqs)
        cur_tb = tbrep(cur_phases(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        Pmat = [Pmat cur_tb];
    end
    
    %     Xmat =  [Pmat bsxfun(@times,Pmat,mod_out)];
    %     %     Xmat = [Xmat mod_out];
    %     klen = size(Xmat,2);
    %     K0 = zeros(klen+1,1);
    %     [fitp] = fit_GLM_phase_model(Xmat,Robs, K0, 1, stim_params, reg_params);
    %     ov_const(t) = fitp.k(end);
    %     stim_ind_filt(t,:) = fitp.k(1:nbins*length(wfreqs));
    %     stim_dep_filt(t,:) = fitp.k((nbins*length(wfreqs)+1):2*nbins*length(wfreqs));
    
    Xmat =  [Pmat];
%     Xmat =  [Pmat(:,(23-1)*nbins + (1:nbins))];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_po] = fit_GLM_phase_model(Xmat, Robs, K0, 1, stim_params, reg_params,0,NL_type);
    ov_const_po(t) = fitp_po.k(end);
    stim_ind_filt_po(t,:) = fitp_po.k(1:nbins*length(wfreqs));
    
    Xmat =  [Pmat cur_stim_mod_out];
    klen = size(Xmat,2);
    K0 = zeros(klen+1,1);
    [fitp_pos] = fit_GLM_phase_model(Xmat, Robs, K0, 1, stim_params, reg_params,0,NL_type);
    stim_ind_filt_pos(t,:) = fitp_pos.k(1:nbins*length(wfreqs));
    
    temp = Pmat*stim_ind_filt_pos(t,:)';
     temp2 = fitp_pos.k(end-1)*cur_stim_mod_out;
     phase_prctiles(t,:) = prctile(exp(temp),prctile_bins);
     stim_prctiles(t,:) = prctile(exp(temp2),prctile_bins);
     phase_prctiles(t,:) = phase_prctiles(t,:)/phase_prctiles(t,mid_prctile);
     stim_prctiles(t,:) = stim_prctiles(t,:)/stim_prctiles(t,mid_prctile);
     
%     Xmat =  [bsxfun(@times,Pmat,cur_stim_mod_out)];
%     klen = size(Xmat,2);
%     K0 = zeros(klen+1,1);
%     [fitp_pg] = fit_GLM_phase_model(Xmat, Robs, K0, 1, stim_params, reg_params2,0,NL_type);
%     stim_dep_filt_pg(t,:) = fitp_pg.k(1:nbins*length(wfreqs));
    
    %
    xv_Robs = new_binned_spks(uset(xv_inds),t);
    cur_trial_Tmat = trial_Tmat(uset(xv_inds),:);
    cur_stim_mod_out = cur_trial_Tmat*stim_kern(t,:)';
    
    cur_phases = squeeze(all_interp_phases(uset(xv_inds),:,nearest_lfps(t)));
    
    Pmat = [];
    for ww = 1:length(wfreqs)
        cur_tb = tbrep(cur_phases(:,ww),pax);
        cur_tb(:,1) = cur_tb(:,1) + cur_tb(:,end);
        cur_tb(:,end) = [];
        Pmat = [Pmat cur_tb];
    end
    %     Xmat =  [Tmat bsxfun(@times,Tmat,mod_out)];
    
    %     xv_pred_rate = Xmat*fitp.k(1:end-1) + fitp.k(end);
    %     xv_pred_rate = log(1+exp(xv_pred_rate));
    %     xv_LL(t) = -sum(xv_Robs.*log(xv_pred_rate) - xv_pred_rate)/sum(xv_Robs);
    %
    Xmat =  [Pmat];
    xv_pred_po = Xmat*fitp_po.k(1:end-1) + fitp_po.k(end);
    if NL_type == 0
        xv_pred_po = log(1+exp(xv_pred_po));
    else
        xv_pred_po = exp(xv_pred_po);
    end
    xv_po_LL(t) = -sum(xv_Robs.*log(xv_pred_po) - xv_pred_po)/sum(xv_Robs);
    
    Xmat =  [Pmat cur_stim_mod_out];
    xv_pred_pos = Xmat*fitp_pos.k(1:end-1) + fitp_pos.k(end);
    if NL_type == 0
        xv_pred_pos = log(1+exp(xv_pred_pos));
    else
        xv_pred_pos = exp(xv_pred_pos);
    end
    xv_pos_LL(t) = -sum(xv_Robs.*log(xv_pred_pos) - xv_pred_pos)/sum(xv_Robs);
    
%     Xmat =  [bsxfun(@times,Pmat,cur_stim_mod_out)];
%     xv_pred_pg = Xmat*fitp_pg.k(1:end-1) + fitp_pg.k(end);
%     if NL_type == 0
%         xv_pred_pg = log(1+exp(xv_pred_pg));
%     else
%         xv_pred_pg = exp(xv_pred_pg);
%     end
%     xv_pg_LL(t) = -sum(xv_Robs.*log(xv_pred_pg) - xv_pred_pg)/sum(xv_Robs);
    
    xv_pred_so = cur_trial_Tmat*fitp_so.k(1:end-1) + fitp_so.k(end);
    if NL_type == 0
        xv_pred_so = log(1+exp(xv_pred_so));
    else
        xv_pred_so = exp(xv_pred_so);
    end
    xv_so_LL(t) = -sum(xv_Robs.*log(xv_pred_so) - xv_pred_so)/sum(xv_Robs);
    
    avg_rate = mean(Robs);
    null_pred = avg_rate*ones(size(xv_Robs));
    xv_null_LL(t) = -sum(xv_Robs.*log(null_pred)-null_pred)/sum(xv_Robs);
    
    % input('')
    % clf
end

%%
save extp1_stim_alpha_models_same_elec xv_* wfreqs pax stim_* *_prctiles prctile*


%%
% ucells = 1:96;
ucells = single_units;
% ucells = find(xv_so_LL < xv_null_LL);
% ucells(ismember(ucells,single_units)) = [];
po_rel_LL = xv_po_LL(ucells) - xv_null_LL(ucells);
pos_rel_LL = xv_pos_LL(ucells) - xv_null_LL(ucells);
% pg_rel_LL = xv_pg_LL(ucells)-xv_null_LL(ucells);
so_rel_LL = xv_so_LL(ucells) - xv_null_LL(ucells);

figure
boxplot([-pos_rel_LL(:) -so_rel_LL(:) -po_rel_LL(:)],{'Phase + Stim','Stim','Phase'})
hold on
plot([-pos_rel_LL(:) -so_rel_LL(:) -po_rel_LL(:)]','o-')
%%
% ucells = 1:96;
ucells = find(xv_so_LL < xv_null_LL);
ucells(ismember(ucells,single_units)) = [];
% ucells = single_units;

avg_phase_pos = reshape(mean(stim_ind_filt_pos(ucells,:)),nbins,length(wfreqs))';
figure
pcolor(doub_pax*180/pi,wfreqs,[avg_phase_pos avg_phase_pos]);shading flat
% set(gca,'yscale','log')

avg_phase_po = reshape(mean(stim_ind_filt_po(ucells,:)),nbins,length(wfreqs));
figure
pcolor(pax(1:end-1),wfreqs,avg_phase_po');shading flat

figure
shadedErrorBar(tent_centers*desired_dt,mean(stim_kern(ucells,:)),std(stim_kern(ucells,:))/sqrt(length(ucells)))
shadedErrorBar(tent_centers*desired_dt,mean(stim_kern(ucells,:)),std(stim_kern(ucells,:))/sqrt(length(ucells)),{'color','r'})


%%
ucells = single_units;

figure
% shadedErrorBar(prctile_bins,mean(stim_prctiles(ucells,:)),std(stim_prctiles(ucells,:))/sqrt(length(ucells)))
hold on
shadedErrorBar(prctile_bins,mean(phase_prctiles(ucells,:)),std(phase_prctiles(ucells,:))/sqrt(length(ucells)),{'color','r'})
% ucells = 1:96;
% shadedErrorBar(prctile_bins,mean(stim_prctiles(ucells,:)),std(stim_prctiles(ucells,:))/sqrt(length(ucells)),{'color','b'})
% shadedErrorBar(prctile_bins,mean(phase_prctiles(ucells,:)),std(phase_prctiles(ucells,:))/sqrt(length(ucells)),{'color','g'})
% plot(prctile_bins,stim_prctiles(ucells,:),'k')
% plot(prctile_bins,phase_prctiles(ucells,:),'r')


%%
% norm_dep_kern = zeros(size(dep_kern));
% norm_ind_kern = zeros(size(ind_kern));
clear norm_ind_kern
for t = 1:96
    cur = reshape(stim_ind_filt_pos(t,:),nbins,length(wfreqs))';
%     norm_ind_kern(t,:,:) = cur/std(cur(:));
        norm_ind_kern(t,:,:) = cur;
    %     cur = squeeze(ind_kern(i,:,:));
    %     norm_ind_kern(i,:,:) = cur/std(cur(:));
end
norm_ind_dom = squeeze(max(norm_ind_kern,[],3)-min(norm_ind_kern,[],3));
% norm_dep_dom = squeeze(max(norm_dep_kern,[],3)-min(norm_dep_kern,[],3));

% ucells = 1:96;
ucells = single_units;
figure
shadedErrorBar(wfreqs,mean(norm_ind_dom(ucells,:)),std(norm_ind_dom(ucells,:))/sqrt(length(ucells)))
figure
plot(wfreqs,norm_ind_dom(ucells,:))
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
close all
for t =22
    gabor_emp1 = get_pgabor_mask_v2(XXc,YYc,gabor_params_f{end}(t,1:6),0);
    gabor_emp2 = get_pgabor_mask_v2(XXc,YYc,gabor_params_f{end}(t,1:6),pi/2);
    phase_kern = reshape(stim_ind_filt_pos(t,:),nbins,length(wfreqs))';
    subplot(2,2,[1 2])
    pcolor(doub_pax*180/pi,wfreqs,[phase_kern phase_kern]);shading flat
    subplot(2,2,3)
    plot(tent_centers*desired_dt,stim_kern(t,:))
    axis tight
    subplot(2,2,4)
    imagesc(xax(xpatch_inds_cropped),yax(ypatch_inds_cropped),gabor_emp1);set(gca,'ydir','normal');%colormap(gray)
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
