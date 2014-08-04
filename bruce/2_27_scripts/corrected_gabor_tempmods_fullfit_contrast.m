clear all
% close all
addpath(genpath('~/James_scripts'));

cd ~/Data/bruce/2_27_12/
load Blocks.mat
accept_window = [-5 5;-5 5];
sac_eyespeed = 10;
thresh_eyespeed = 2.5;
sac_buffer = 0;
min_fix_dur = 0.2;
blink_thresh = 5;
min_blink_dur = 0.05;
% nlags = 4;


Pix2Deg = 0.018837;

% down-sampling fraction for image
dsfrac = 2;
Fsd = 1/Pix2Deg/dsfrac;

Nyp = 1024;
Nxp = 1280;
Ny = round(Nyp/dsfrac);
Nx = round(Nxp/dsfrac);
xax = linspace(-Nx/2,Nx/2,Nx)/Fsd; yax = linspace(-Ny/2,Ny/2,Ny)/Fsd;
% [XX,YY] = meshgrid(xax,yax);

% RF_patch = [3.5 6; -3.5 -1]; %location of RFs in degrees [x1 x2;y1 y2]
% RF_patch = [2.5 6.5; -4 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch = [2.75 6.25; -3.5 0]; %location of RFs in degrees [x1 x2;y1 y2]
RF_patch_pix = Fsd*RF_patch;
RF_patch_cent = mean(RF_patch,2);
RF_patch_width = diff(RF_patch,[],2);

% patch_inds = find(XX >= RF_patch(1,1) & XX <= RF_patch(1,2) & YY >= RF_patch(2,1) & YY <= RF_patch(2,2));
xpatch_inds = find(xax >= RF_patch(1,1) & xax <= RF_patch(1,2));
ypatch_inds = find(yax >= RF_patch(2,1) & yax <= RF_patch(2,2));

cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_data_v3

cd /Users/James/Data/bruce/2_27_12/stimrecon
load sing_eye_pgabortrack_d2
load sing_eye_pgabortrack_fin_est_d2
cd ~/Data/bruce/2_27_12/ExptA

SDIM = length(xpatch_inds);
cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_image_patches_corrected

load ./avg_lfp_hfpow

NT = length(used_fixs);
X_resh = reshape(X,NT,SDIM^2);
dt = .005;

cellids = [1 2 3 4 5 6 7 8 10];
muaids = [1 2 4 5 6 11 12 13 14];

spk_cnts = [spk_cnts(:,cellids) mua_cnts(:,muaids)];

n_used_cells = size(spk_cnts,2);

%%
fix_stimout = zeros(n_used_cells,length(used_fixs));
all_model_time_axis = [];
all_model_blockids = [];
all_model_fixids = [];
spikes_binned = [];
time_since_fix = [];
rel_spk_times = cell(18,1);
spk_fix_inds = cell(18,1);
for blockid = 1:4;
    fprintf('Block %d of %d\n',blockid,4);
    cur_set = find(blockids(used_fixs) == blockid);
    
    for i = 1:length(cur_set)
        start_T = all_fix_start_times(used_fixs(cur_set(i)));
        end_T = all_fix_stop_times(used_fixs(cur_set(i)));
        cur_tedges = start_T:dt:end_T;
        cur_tcents = 0.5*cur_tedges(1:end-1)+0.5*cur_tedges(2:end);
        all_model_time_axis = [all_model_time_axis cur_tcents];
        all_model_blockids = [all_model_blockids blockid*ones(1,length(cur_tcents))];
        all_model_fixids = [all_model_fixids cur_set(i)*ones(1,length(cur_tcents))];
        temp_binned = zeros(10,length(cur_tcents));
        temp_modouts = zeros(10,length(cur_tcents));
        for c = 1:n_used_cells
            if c <= 9
                cur_spk_times = find(Blocks{blockid}.spktimes{cellids(c)} > start_T & ...
                    Blocks{blockid}.spktimes{cellids(c)} < end_T);
                temp = hist(Blocks{blockid}.spktimes{cellids(c)}(cur_spk_times),cur_tcents);
                rel_spk_times{c} = [rel_spk_times{c} Blocks{blockid}.spktimes{cellids(c)}(cur_spk_times)-start_T];
                spk_fix_inds{c} = [spk_fix_inds{c} cur_set(i)*ones(size(cur_spk_times))];
            else
                cur_spk_times = find(Blocks{blockid}.mutimes{muaids(c-9)} > start_T & ...
                    Blocks{blockid}.mutimes{muaids(c-9)} < end_T);
                temp = hist(Blocks{blockid}.mutimes{muaids(c-9)}(cur_spk_times),cur_tcents);
                rel_spk_times{c} = [rel_spk_times{c} Blocks{blockid}.mutimes{muaids(c-9)}(cur_spk_times)-start_T];
                spk_fix_inds{c} = [spk_fix_inds{c} cur_set(i)*ones(size(cur_spk_times))];
            end
            temp_binned(c,:) = temp;
        end
        spikes_binned = [spikes_binned; temp_binned'];
        
        time_since_fix = [time_since_fix cur_tcents-start_T];
        
    end
end

%% Compute TBR time-since fix onset
max_tsf = 0.75; nbins = 30;
used_tsf = time_since_fix;
used_tsf(used_tsf > max_tsf) = [];
tax = prctile(used_tsf,linspace(100/nbins,100-100/nbins,nbins));

% tax = logspace(log10(dt),log10(max_tsf),nbins);
Tmat = tbrep(time_since_fix,tax);

%%
xv_frac = 0;
rperm = randperm(length(used_fixs));
xv_fixs = rperm(1:round(xv_frac*length(used_fixs)));
xv_inds = find(ismember(all_model_fixids,xv_fixs));
% tr_inds = find(~ismember(all_model_fixids,xv_fixs));
tr_inds = 1:length(all_model_fixids);

% xv_fixs = find(all_stim_filtered(used_fixs)==0);
% xv_inds = find(ismember(all_model_fixids,xv_fixs));
% tr_inds = find(~ismember(all_model_fixids,xv_fixs));
%% COMPUTE SACCADE TRIG AVG FIRING RATES
sac_trg_rate = zeros(n_used_cells,length(tax));
sac_trg_var = zeros(n_used_cells,length(tax));
taxb = [0 tax];
for i = 1:length(tax)
    curset = find(time_since_fix >= taxb(i) & time_since_fix < taxb(i+1));
    for c = 1:n_used_cells
        sac_trg_rate(c,i) = mean(spikes_binned(curset,c));
        sac_trg_var(c,i) = var(spikes_binned(curset,c));
    end
    n_occ(i) = length(curset);
end
avg_rates = mean(spikes_binned)/dt;
fano_fac = sac_trg_var./sac_trg_rate;
sac_trg_rate = sac_trg_rate/dt;
sac_trg_std = sqrt(sac_trg_var)/dt;
%%
clear mod_*
sac_hist = zeros(n_used_cells,length(tax));
stim_kern = cell(n_used_cells,1);
slope2_lambda = 150;
l1_pen = 20;
stim_mod_predrate = zeros(9,length(time_since_fix));
sac_hist_predrate = zeros(9,length(time_since_fix));
for c = 1:n_used_cells
    fprintf('Cell %d of %d\n',c,9);
    
    un_spk_cnts = length(unique(spikes_binned(tr_inds,c))) - 1;
    spkbs = [];
    for i = 1:un_spk_cnts
        curset = find(spikes_binned(tr_inds,c) == i);
        spkbs = [spkbs; repmat(curset,i,1)];
    end
    spkbs = sort(spkbs);
    Robs = spikes_binned(tr_inds,c);
    Robs_xv = spikes_binned(xv_inds,c);
    
    %%
    disp('Fitting sac hist model')
    %first fit sachist model
    Xmat = Tmat;
    W0 = zeros(size(Xmat,2),1);
    NT = length(tr_inds);
    
    %null model
    avg_rate = repmat(mean(Robs),length(Robs),1);
    null_LL(c) = sum(-Robs.*log(avg_rate)+avg_rate)/sum(Robs);
    avg_rate_xv = repmat(mean(Robs),length(Robs_xv),1);
%     null_xvLL(c) = sum(-Robs_xv.*log(avg_rate_xv)+avg_rate_xv)/sum(Robs_xv);
    
    %     %initialize parameters
    silent = 1;
    lamrange = [];
    lamrange2 = [slope2_lambda 1 length(tax)];
    nlxrange = [tax];
    Pcon = [];
    Pmono = [];
    hold_const = [];
    NLtype = 0;
    llist = [];
    %
    %     [fitp,grad] = GLMsolve_jmm( Xmat(tr_inds,:), spkbs, W0, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const, NLtype, 1e-6,nlxrange);
    %     sac_hist_k(c,:) = fitp.k(1:end-1);
    %     sac_hist_const(c) = fitp.k(end);
    %     sac_hist_LP(c) = fitp.LP;
    %
    %     cur_genfun = Tmat(tr_inds,:)*sac_hist_k(c,:)' + sac_hist_const(c);
    %     cur_predrate = log(1+exp(cur_genfun));
    %     sac_hist_LL(c) = sum(-Robs.*log(cur_predrate) + cur_predrate)/sum(Robs);
    %     cur_genfun_xv = Tmat(xv_inds,:)*sac_hist_k(c,:)' + sac_hist_const(c);
    %     cur_predrate_xv = log(1+exp(cur_genfun_xv));
    %     sac_hist_xvLL(c) = sum(-Robs_xv.*log(cur_predrate_xv) + cur_predrate_xv)/sum(Robs_xv);
    %     cur_genfun = Tmat*sac_hist_k(c,:)' + sac_hist_const(c);
    %     sac_hist_predrate(c,:) = log(1+exp(cur_genfun));
    
    %%
    disp('Fitting total stim model')
    X_resh = reshape(X,size(X,1),SDIM^2);
    cur_mask1 = get_pgabor_mask(gabor_params_fin(c,1:6),0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(gabor_params_fin(c,1:6),pi/2,[SDIM SDIM]);
    cur_gmask = get_pgauss_mask(gabor_params_fin(c,1:6),[SDIM SDIM]);
    cur_gmask = cur_gmask(:);
    mask1_out = X_resh*cur_mask1(:);
    mask2_out = X_resh*cur_mask2(:);
    gmask_out = bsxfun(@times,X_resh,cur_gmask');
    loc_mean = mean(gmask_out,2);
    loc_cont = std(gmask_out,[],2);
    
    mask1_out = (mask1_out - loc_mean)./loc_cont;
    mask2_out = (mask2_out - loc_mean)./loc_cont;
    
    loc_mean = zscore(loc_mean);
    loc_cont = zscore(loc_cont);
    
    lin_out = gabor_params_fin(c,8)*mask1_out(all_model_fixids) + gabor_params_fin(c,9)*mask2_out(all_model_fixids);
    energy_out = gabor_params_fin(c,7)*sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
    total_out = lin_out + energy_out;
    total_out = zscore(total_out);
    
    lamrange2 = [slope2_lambda 1 length(tax);
        slope2_lambda length(tax)+1 length(tax)*2];
    nlxrange = [tax tax];
    Xmat = Tmat;
    Xmat = [Xmat bsxfun(@times,Tmat,total_out)];
           
    [fitp,grad] = GLMsolve_jmm( Xmat(tr_inds,:), spkbs, W0, 1, lamrange, lamrange2,Pcon, Pmono, llist, [], NLtype,1e-5,nlxrange);
    stim_mod_const(c) = fitp.k(end);
    stim_mod_sachist(c,:) = fitp.k(1:length(tax));
    stim_mod_tot(c,:) = fitp.k(length(tax)+1:length(tax)*2);
    
    %%
    disp('Fitting separate stim model')
    ori_vals = linspace(0,pi,13);
    ori_vals(end) = [];
    
    ori_outs = zeros(length(total_out),length(ori_vals));
    cur_params = gabor_params_fin(c,1:6);
    for oo = 1:length(ori_vals)
        cur_params(3) = ori_vals(oo);
        cur_mask1 = get_pgabor_mask(cur_params,0,[SDIM SDIM]);
        cur_mask2 = get_pgabor_mask(cur_params,pi/2,[SDIM SDIM]);
        mask1_out = X_resh*cur_mask1(:);
        mask2_out = X_resh*cur_mask2(:);
%         lin_out = gabor_params_fin(c,8)*mask1_out(all_model_fixids) + gabor_params_fin(c,9)*mask2_out(all_model_fixids);
%         energy_out = gabor_params_fin(c,7)*sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
        energy_out = sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
%         ori_outs(:,oo) = lin_out + energy_out;
        ori_outs(:,oo) = energy_out;
    end
    cur_mask1 = get_pgabor_mask(gabor_params_fin(c,1:6),0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(gabor_params_fin(c,1:6),pi/2,[SDIM SDIM]);
    mask1_out = X_resh*cur_mask1(:);
    mask2_out = X_resh*cur_mask2(:);
%     lin_out = gabor_params_fin(c,8)*mask1_out(all_model_fixids) + gabor_params_fin(c,9)*mask2_out(all_model_fixids);
    energy_out = sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
    total_out = energy_out;
    
    cur_probe = Blocks{1}.suprobes(cellids(c));
    adj_probes = [(cur_probe-1) (cur_probe+1)];
    adj_probes(adj_probes > 24 | adj_probes < 1) = [];
    net_out = mean(lfp_avg_pow(adj_probes,:),1);
    net_out = zscore(net_out(all_model_fixids))';
    
    contrast_out = sum(ori_outs,2);
    total_out = total_out./contrast_out;
    
    contrast_out = zscore(contrast_out);
    total_out = zscore(total_out);
    
    lamrange2 = [slope2_lambda 1 length(tax);
        slope2_lambda length(tax)+1 length(tax)*2;
        slope2_lambda length(tax)*2+1 length(tax)*3];
    nlxrange = [tax tax tax];
    
    Xmat = Tmat;
    Xmat = [Xmat bsxfun(@times,Tmat,total_out)];
    Xmat = [Xmat bsxfun(@times,Tmat,contrast_out)];
    llist = [];
    [fitp,grad] = GLMsolve_jmm( Xmat(tr_inds,:), spkbs, W0, 1, lamrange, lamrange2,Pcon, Pmono, llist, [], NLtype,1e-5,nlxrange);
    sep_mod_const(c) = fitp.k(end);
    sep_mod_sachist(c,:) = fitp.k(1:length(tax));
    sep_mod_stim(c,:) = fitp.k(length(tax)+1:length(tax)*2);
    sep_mod_cont(c,:) = fitp.k(length(tax)*2+1:length(tax)*3);
    cur_genfun = Xmat(tr_inds,:)*fitp.k(1:end-1) + sep_mod_const(c);
    cur_predrate = log(1+exp(cur_genfun));
    sep_mod_LL(c) = sum(-Robs.*log(cur_predrate) + cur_predrate)/sum(Robs);

    lamrange2 = [slope2_lambda 1 length(tax);
        slope2_lambda length(tax)+1 length(tax)*2;
        slope2_lambda length(tax)*2+1 length(tax)*3;
        slope2_lambda length(tax)*3+1 length(tax)*4];
    nlxrange = [tax tax tax tax];
    
    Xmat = Tmat;
    Xmat = [Xmat bsxfun(@times,Tmat,total_out)];
    Xmat = [Xmat bsxfun(@times,Tmat,contrast_out)];
    Xmat = [Xmat bsxfun(@times,Tmat,net_out)];
    llist = [];
    [fitp,grad] = GLMsolve_jmm( Xmat(tr_inds,:), spkbs, W0, 1, lamrange, lamrange2,Pcon, Pmono, llist, [], NLtype,1e-5,nlxrange);
    asep_mod_const(c) = fitp.k(end);
    asep_mod_sachist(c,:) = fitp.k(1:length(tax));
    asep_mod_stim(c,:) = fitp.k(length(tax)+1:length(tax)*2);
    asep_mod_cont(c,:) = fitp.k(length(tax)*2+1:length(tax)*3);
    asep_mod_net(c,:) = fitp.k(length(tax)*3+1:length(tax)*4);
    cur_genfun = Xmat(tr_inds,:)*fitp.k(1:end-1) + sep_mod_const(c);
    cur_predrate = log(1+exp(cur_genfun));
    asep_mod_LL(c) = sum(-Robs.*log(cur_predrate) + cur_predrate)/sum(Robs);
   
    
end

%%
cd /Users/James/James_scripts/bruce/modelfits
for c = 1:9
% clf
[~,peakloc] = max(stim_mod_sachist(c,:));
subplot(4,1,1)
plot(tax,stim_mod_sachist(c,:))
axis tight
xlim([0 0.4])
yl = ylim();
line([tax(peakloc) tax(peakloc)],yl,'color','k','linestyle','--')
title('Stimulus-Independent')
subplot(4,1,2)
plot(tax,stim_mod_tot(c,:),'r')
axis tight
xlim([0 0.4])
yl = ylim();
line([tax(peakloc) tax(peakloc)],yl,'color','k','linestyle','--')
title('Stimulus Only')
subplot(4,1,3)
plot(tax,sep_mod_stim(c,:),'r')
hold on
plot(tax,sep_mod_cont(c,:),'k')
axis tight
xlim([0 0.4])
yl = ylim();
line([tax(peakloc) tax(peakloc)],yl,'color','k','linestyle','--')
title('Stimulus + Contrast')
subplot(4,1,4)
plot(tax,asep_mod_stim(c,:),'r')
hold on
plot(tax,asep_mod_cont(c,:),'k')
plot(tax,asep_mod_net(c,:),'g')
axis tight
xlim([0 0.4])
yl = ylim();
line([tax(peakloc) tax(peakloc)],yl,'color','k','linestyle','--')
title('Stimulus + Contrast + HF Power')
fillPage(gcf,'margins',[0 0 0 0],'papersize',[6 14]);
fname = sprintf('Network_contrast_compare_cell%d',cellids(c));
print('-dpdf','-painters',fname)
close 
end


%%
% for c = 1:9;
%     fprintf('Cell %d of %d\n',c,9);
% un_spk_cnts = length(unique(spikes_binned(tr_inds,c))) - 1;
% spkbs = [];
% for i = 1:un_spk_cnts
%     curset = find(spikes_binned(tr_inds,c) == i);
%     spkbs = [spkbs; repmat(curset,i,1)];
% end
% spkbs = sort(spkbs);
% Robs = spikes_binned(tr_inds,c);
% Robs_xv = spikes_binned(xv_inds,c);
%
% init_params = gabor_params_fin(c,:);
% cur_gain = Tmat*stim_mod_tot(c,:)';
% cur_offset = Tmat*stim_mod_sachist(c,:)';
% LB = [-5 -5 0 6 2 0.5 -Inf -Inf -Inf -Inf];
% UB = [10 10 pi SDIM(1)/2 SDIM(1)/6 2 Inf Inf Inf Inf];
% hold_const = [1 1 1 1 1 1 1 1 1 0];
% [gabor_params,LL] = fit_temporal_gabor_params(init_params,X_resh,cur_gain(tr_inds),cur_offset(tr_inds),...
%     Robs,all_model_fixids(tr_inds),[SDIM SDIM],hold_const,LB,UB);
%
% % hold_const = [1 1 1 1 1 1 0 0 0 0];
% % [gabor_params,LL] = fit_temporal_gabor_params(gabor_params,X_resh,cur_gain(tr_inds),cur_offset(tr_inds),...
% %     Robs,all_model_fixids(tr_inds),[SDIM SDIM],hold_const,LB,UB);
%
% hold_const = [1 1 0 0 0 0 0 0 0 0];
% [new_gabor_params(c,:),LL] = fit_temporal_gabor_params(gabor_params,X_resh,cur_gain(tr_inds),cur_offset(tr_inds),...
%     Robs,all_model_fixids(tr_inds),[SDIM SDIM],hold_const,LB,UB);
%
% cur_mask1 = get_pgabor_mask(new_gabor_params(c,1:6),0,[SDIM SDIM]);
% cur_mask2 = get_pgabor_mask(new_gabor_params(c,1:6),pi/2,[SDIM SDIM]);
% mask1_out = X_resh*cur_mask1(:);
% mask2_out = X_resh*cur_mask2(:);
% lin_out = new_gabor_params(c,8)*mask1_out(all_model_fixids) + new_gabor_params(c,9)*mask2_out(all_model_fixids);
% energy_out = new_gabor_params(c,7)*sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
% total_out = lin_out + energy_out;
% total_out = zscore(total_out);
% lamrange2 = [slope2_lambda 1 length(tax);
%     slope2_lambda length(tax)+1 length(tax)*2];
% nlxrange = [tax tax];
% Xmat = Tmat;
% Xmat = [Xmat bsxfun(@times,Tmat,total_out)];
% llist = [l1_pen (length(tax)+1):size(Xmat,2)];
% [fitp,grad] = GLMsolve_jmm( Xmat(tr_inds,:), spkbs, W0, 1, lamrange, lamrange2,Pcon, Pmono, llist, [], NLtype,1e-5,nlxrange);
% rstim_mod_const(c) = fitp.k(end);
% rstim_mod_sachist(c,:) = fitp.k(1:length(tax));
% rstim_mod_tot(c,:) = fitp.k(length(tax)+1:length(tax)*2);
%
% cur_genfun = Xmat(tr_inds,:)*fitp.k(1:end-1) + rstim_mod_const(c);
% cur_predrate = log(1+exp(cur_genfun));
% rstim_mod_LL(c) = sum(-Robs.*log(cur_predrate) + cur_predrate)/sum(Robs);
% cur_genfun_xv = Xmat(xv_inds,:)*fitp.k(1:end-1) + stim_mod_const(c);
% cur_predrate_xv = log(1+exp(cur_genfun_xv));
% rstim_mod_xvLL(c) = sum(-Robs_xv.*log(cur_predrate_xv) + cur_predrate_xv)/sum(Robs_xv);
% cur_genfun = Xmat*fitp.k(1:end-1) + rstim_mod_const(c);
% rstim_mod_predrate(c,:) = log(1+exp(cur_genfun));
% end
%%
stim_pred_diff = zeros(n_used_cells,length(tax));
obs_diff = zeros(n_used_cells,length(tax));
taxb = [0 tax];
for i = 1:length(tax)
    curset = find(time_since_fix >= taxb(i) & time_since_fix < taxb(i+1));
    for c = 1:n_used_cells
        %         stim_pred_diff(c,i) =mean(abs(stim_mod_predrate(c,curset)-sac_hist_predrate(c,curset)));
        stim_pred_diff(c,i) =sqrt(mean((stim_mod_predrate(c,curset)-sac_hist_predrate(c,curset)).^2));
        %         obs_diff(c,i) = mean(abs(spikes_binned(curset,c)'-sac_hist_predrate(c,curset)));
        obs_diff(c,i) = sqrt(mean((spikes_binned(curset,c)'-sac_hist_predrate(c,curset)).^2));
    end
    n_occ(i) = length(curset);
end
stim_pred_diff = stim_pred_diff/dt;
obs_diff = obs_diff/dt;

%%
cd /Users/James/James_scripts/bruce/modelfits/
% save gabor_tempmodfits_en_lin stim_mod* sac_* null* tax stim_mod_predrate all_model_* fano* time_since_fix sac_trg*
save gabor_tempmodfits_sepmod stim_mod* sep_mod* sac_* null* tax stim_mod_predrate all_model_* fano* time_since_fix sac_trg*

%%
c = 14;
close all
plot(tax,sep_mod_stim(c,:),'r')
hold on
% plot(tax,sep_mod_mean(c,:),'k')
plot(tax,sep_mod_cont(c,:),'b')
% plot(tax,sep_mod_locmean(c,:),'k--')
plot(tax,sep_mod_loccont(c,:),'b--')
plot(tax,stim_mod_tot(c,:),'k')

% figure
% plot(tax,sep_mod_sachist(c,:),'b')
% hold on
% plot(tax,stim_mod_sachist(c,:),'k')
%%
cd /Users/James/James_scripts/bruce/modelfits
close all
% for cur_cell = 1:9;
cur_cell = 6;
X_resh = reshape(X,size(X,1),SDIM^2);
cur_mask1 = get_pgabor_mask(gabor_params_fin(cur_cell,1:6),0,[SDIM SDIM]);
cur_mask2 = get_pgabor_mask(gabor_params_fin(cur_cell,1:6),pi/2,[SDIM SDIM]);
mask1_out = X_resh*cur_mask1(:);
mask2_out = X_resh*cur_mask2(:);
lin_out = gabor_params_fin(cur_cell,8)*mask1_out + gabor_params_fin(cur_cell,9)*mask2_out;
energy_out = gabor_params_fin(cur_cell,7)*sqrt(mask1_out.^2 + mask2_out.^2);
total_out = lin_out + energy_out;
total_out = zscore(total_out);

n_fixs = length(total_out);
maxlag = round(0.5/dt);
fix_trg_mat = nan(n_fixs,maxlag);
for i = 1:n_fixs
    cur_set = find(all_model_fixids==i);
    cur_set(maxlag+1:end) = [];
    fix_trg_mat(i,1:length(cur_set)) = spikes_binned(cur_set,cur_cell);
end

n_stims = 10;
cmap = jet(n_stims);
% stim_vals = linspace(-2,2,n_stims);
stim_val_edges = prctile(total_out,linspace(0,100,n_stims+1));
stim_vals = 0.5*stim_val_edges(1:end-1) + 0.5*stim_val_edges(2:end);
stim_vals(end) = stim_vals(end-1) + median(diff(stim_vals));
clear cur_rate
cur_tax = linspace(0,0.5,200);
test_tax = tbrep(cur_tax,tax);
cur_rate = zeros(n_stims,length(cur_tax));
cur_trg_avg = zeros(n_stims,maxlag);
for i = 1:n_stims
    cur_gen = stim_mod_sachist(cur_cell,:)*test_tax';
    cur_gen = cur_gen + stim_mod_tot(cur_cell,:)*test_tax'*stim_vals(i);
    cur_gen = cur_gen + stim_mod_const(cur_cell);
    cur_rate(i,:) = log(1+exp(cur_gen));
    
    cur_set = find(total_out >= stim_val_edges(i) & total_out < stim_val_edges(i+1));
    cur_trg_avg(i,:) = nanmean(fix_trg_mat(cur_set,:));
    cur_trg_avg(i,:) = smooth(cur_trg_avg(i,:),8,'lowess');
end
norm_crate = bsxfun(@rdivide,cur_rate,max(cur_rate,[],2));
norm_trg_rate = bsxfun(@rdivide,cur_trg_avg,max(cur_trg_avg,[],2));
cur_rate = cur_rate/dt;
cur_trg_avg = cur_trg_avg/dt;

[~,stim_ord] = sort(total_out);
ordered_spikes = [];
for i = 1:length(stim_ord)
    cur_set = find(spk_fix_inds{cur_cell} == stim_ord(i));
    cur_set = cur_set(rel_spk_times{cur_cell}(cur_set) <= 0.5);
    ordered_spikes = [ordered_spikes rel_spk_times{cur_cell}(cur_set)];
end
% for i = 1:n_stims
%     cur_set = find(total_out >= stim_val_edges(i) & total_out < stim_val_edges(i+1));
%     cc(i) = length(cur_set);
%     cur_spk_set = find(ismember(spk_fix_inds{cur_cell},cur_set));
%     cur_spk_set = cur_spk_set(rel_spk_times{cur_cell}(cur_spk_set) <= 0.5);
%     ordered_spikes = [ordered_spikes rel_spk_times{cur_cell}(cur_spk_set)];
% end

for i = 1:n_stims
    subplot(2,2,3)
        plot(cur_tax,norm_crate(i,:),'color',cmap(i,:),'linewidth',1)
        hold on
        subplot(2,2,4)
        plot((1:maxlag)*dt,norm_trg_rate(i,:),'color',cmap(i,:),'linewidth',1)
        hold on
    
    subplot(2,2,1)
    plot(cur_tax,cur_rate(i,:),'color',cmap(i,:),'linewidth',1)
    hold on
    subplot(2,2,2)
    plot((1:maxlag)*dt,cur_trg_avg(i,:),'color',cmap(i,:),'linewidth',1)
    hold on
end
subplot(2,2,1)
title('Predicted','fontsize',14)
xlim([0 0.4])
xlabel('Time since fixation onset (s)','fontsize',14)
ylabel('Average Rate (Hz)','fontsize',14)
subplot(2,2,2)
title('Observed','fontsize',14)
xlim([0 0.4])
xlabel('Time since fixation onset (s)','fontsize',14)
ylabel('Average Rate (Hz)','fontsize',14)
subplot(2,2,3)
title('Predicted','fontsize',14)
xlim([0 0.2])
xlabel('Time since fixation onset (s)','fontsize',14)
ylabel('Average Rate (relative)','fontsize',14)
subplot(2,2,4)
title('Observed','fontsize',14)
xlim([0 0.2])
xlabel('Time since fixation onset (s)','fontsize',14)
ylabel('Average Rate (relative)','fontsize',14)

% fillPage(gcf,'margins',[0 0 0 0],'papersize',[10 8]);
% fname = sprintf('Cell%d_modelcompare',cellids(cur_cell));
% print('-dpdf','-painters',fname);
% close 

figure
raster(ordered_spikes)
xlim([0 0.4])
xlabel('Time since fixation onset (s)','fontsize',14)
ylabel('Fixation','fontsize',14)
fillPage(gcf,'margins',[0 0 0 0],'papersize',[10 8]);
% fname = sprintf('Cell%d_raster',cellids(cur_cell));
% print('-dpdf','-painters',fname);
% close 
% end


%% CREATE SU STIM MODEL FIGURES
close all
for cur_cell = 1:9;
    % cur_cell = 6;
    cur_tsf = linspace(0,0.5,100);
    cur_Tmat = tbrep(cur_tsf,tax);
    x_out = ones(length(cur_tsf),1)*-1;
    cur_X = bsxfun(@times,cur_Tmat,x_out);
    x_out2 = ones(length(cur_tsf),1)*1;
    cur_X2 = bsxfun(@times,cur_Tmat,x_out2);
    x_out3 = zeros(length(cur_tsf),1);
    cur_X3 = bsxfun(@times,cur_Tmat,x_out3);
    
    cur_genfun = cur_Tmat*stim_mod_sachist(cur_cell,:)' + cur_X*stim_mod_tot(cur_cell,:)' + stim_mod_const(cur_cell);
    cur_predrate = log(1+exp(cur_genfun))/dt;
    cur_genfun2 = cur_Tmat*stim_mod_sachist(cur_cell,:)' + cur_X2*stim_mod_tot(cur_cell,:)' + stim_mod_const(cur_cell);
    cur_predrate2 = log(1+exp(cur_genfun2))/dt;
    cur_genfun3 = cur_Tmat*stim_mod_sachist(cur_cell,:)' + cur_X3*stim_mod_tot(cur_cell,:)' + stim_mod_const(cur_cell);
    cur_predrate3 = log(1+exp(cur_genfun3))/dt;
    
    cur_genfun4 = cur_Tmat*sac_hist(cur_cell,:)'+ sac_hist_const(cur_cell);
    cur_predrate4 = log(1+exp(cur_genfun4))/dt;
    
    stim_dep = cur_predrate2-cur_predrate;
    [~,peakloc] = max(cur_predrate3);
    
    % sachist_var = var(sac_hist_predrate(cur_cell,:));
    % stim_var = zeros(size(tax));
    % for i = 1:length(tax)
    %     cur = stim_mod_predrate(cur_cell,:)'.*Tmat(:,i);
    %     stim_var(i) = var(cur);
    % end
    
    figure
    subplot(3,3,1)
    set(gca,'fontsize',14)
    plot(cur_tsf,cur_predrate3,'k','linewidth',1)
    hold on
    plot(tax,sac_trg_rate(cur_cell,:),'b','linewidth',1)
    plot(cur_tsf,cur_predrate2-cur_predrate,'r','linewidth',1)
    xlim([0 0.4])
    yl = ylim();
    line(cur_tsf([peakloc peakloc]),yl,'color','k','linestyle','--')
    xlabel('Time since fixation onset (s)','fontsize',18)
    ylabel('Firing Rate (Hz)','fontsize',18)
    
    subplot(3,3,4)
    set(gca,'fontsize',14)
    plot(cur_tsf,stim_dep./cur_predrate3,'b','linewidth',1)
    hold on
    % plot(tax,stim_var/sachist_var,'r')
    xlim([0 0.4])
    yl = ylim();
    line(cur_tsf([peakloc peakloc]),yl,'color','k','linestyle','--')
    xlabel('Time since fixation onset (s)','fontsize',18)
    ylabel('Stimulus Dependence Ratio','fontsize',18)
    
    
    xl = [0 0.4];
    n_stims = 10;
    cmap = jet(n_stims);
    stim_vals = linspace(-2,2,n_stims);
    clear cur_rate
    cur_tax = linspace(0,0.5,200);
    test_tax = tbrep(cur_tax,tax);
    for i = 1:n_stims
        cur_gen = stim_mod_sachist(cur_cell,:)*test_tax';
        cur_gen = cur_gen + stim_mod_tot(cur_cell,:)*test_tax'*stim_vals(i);
        cur_gen = cur_gen + stim_mod_const(cur_cell);
        cur_rate(i,:) = log(1+exp(cur_gen));
    end
    norm_crate = bsxfun(@rdivide,cur_rate,max(cur_rate,[],2));
    cur_rate = cur_rate/dt;
    
    for i = 1:n_stims
        subplot(3,3,2)
        plot(cur_tax,norm_crate(i,:),'color',cmap(i,:),'linewidth',1)
        hold on
        subplot(3,3,5)
        plot(cur_tax,cur_rate(i,:),'color',cmap(i,:),'linewidth',1)
        hold on
    end
    % subplot(2,1,1)
    % plot(tax,sac_trg_rate(cur_cell,:)/max(sac_trg_rate(cur_cell,:)),'color','k','linewidth',2)
    
    subplot(3,3,2)
    set(gca,'fontsize',14)
    xlim(xl);
    xlabel('Time since fixation onset (s)','fontsize',18)
    ylabel('Relative predicted rate','fontsize',18)
    yl = ylim();
    line(cur_tsf([peakloc peakloc]),yl,'color','k','linestyle','--')
    subplot(3,3,5)
    set(gca,'fontsize',14)
    xlim(xl);
    xlabel('Time since fixation onset (s)','fontsize',18)
    ylabel('Predicted rate (Hz)','fontsize',18)
    yl = ylim();
    line(cur_tsf([peakloc peakloc]),yl,'color','k','linestyle','--')
    
    subplot(3,3,7)
    set(gca,'fontsize',14)
    plot(tax,stim_mod_sachist(cur_cell,:),'linewidth',1)
    xlabel('Time since fixation onset (s)','fontsize',18)
    ylabel('Constant Term','fontsize',18)
    xlim([0 0.4])
    yl = ylim();
    line(cur_tsf([peakloc peakloc]),yl,'color','k','linestyle','--')
    subplot(3,3,8)
    set(gca,'fontsize',14)
    plot(tax,stim_mod_tot(cur_cell,:),'linewidth',1)
    xlabel('Time since fixation onset (s)','fontsize',18)
    ylabel('Stimulus Gain','fontsize',18)
    xlim([0 0.4])
    yl = ylim();
    line(cur_tsf([peakloc peakloc]),yl,'color','k','linestyle','--')
    
    cur_mask1 = get_pgabor_mask(gabor_params_fin(cur_cell,1:6),0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(gabor_params_fin(cur_cell,1:6),pi/2,[SDIM SDIM]);
    cur_ovmask = gabor_params_fin(cur_cell,8)*cur_mask1 + gabor_params_fin(cur_cell,9)*cur_mask2;
    subplot(3,3,3)
    imagesc(xax(xpatch_inds),yax(ypatch_inds),cur_ovmask);
    set(gca,'ydir','normal')
    colormap(gray)
    %     caxis([-20 20]*1e-4);
    X_resh = reshape(X,size(X,1),SDIM^2);
    mask1_out = X_resh*cur_mask1(:);
    mask2_out = X_resh*cur_mask2(:);
    lin_out = gabor_params_fin(cur_cell,8)*mask1_out + gabor_params_fin(cur_cell,9)*mask2_out;
    energy_out = gabor_params_fin(cur_cell,7)*sqrt(mask1_out.^2 + mask2_out.^2);
    total_out = lin_out + energy_out;
    total_out = zscore(total_out);
    
    cur_spk_cnts = spk_cnts(used_fixs,cur_cell);
    T = tabulate(cur_spk_cnts);
    un_cnts = T(:,1);
    un_cnts(T(:,2) < 5) = [];
    avg_stim = zeros(size(un_cnts));
    for i = 1:length(un_cnts)
        cur_set = find(cur_spk_cnts==un_cnts(i));
        avg_stim(i) = mean(total_out(cur_set));
    end
    subplot(3,3,6)
    plot(spk_cnts(used_fixs,cur_cell),total_out,'.','markersize',8)
    hold on
    plot(un_cnts,avg_stim,'r','linewidth',2)
    axis tight
    xlabel('Spike counts','fontsize',18)
    ylabel('Stimulus Strength','fontsize',18)
    
    % end
    subplot(3,3,9)
    hist(total_out,100)
    axis tight
    xlabel('Stimulus Strength','fontsize',18)
    
    
    fillPage(gcf,'margins',[0 0 0 0],'papersize',[30 25]);
    fname = sprintf('Cell%d_model_char',cellids(cur_cell));
    print('-dpdf','-painters',fname)
    close
end

%%
close all
    cur_cell = 7;
    [~,peakloc] = max(sac_trg_rate(cur_cell,:));
        
    figure
    subplot(3,1,1)
    set(gca,'fontsize',14)
    hold on
    plot(tax,sac_trg_rate(cur_cell,:),'o-','linewidth',1)
    axis tight
    xlim([0 0.4])
    yl = ylim();
    ylim([0 yl(2)]);
    yl = ylim();
    line(tax([peakloc peakloc]),yl,'color','k','linestyle','--')
    xlabel('Time since fixation onset (s)','fontsize',18)
    ylabel('Firing Rate (Hz)','fontsize',18)
    box('off')
    
    subplot(3,1,2)
    set(gca,'fontsize',14)
    plot(tax,stim_mod_sachist(cur_cell,:),'o-','linewidth',1)
    xlabel('Time since fixation onset (s)','fontsize',18)
    ylabel('Constant Term','fontsize',18)
    axis tight
    xlim([0 0.4])
    yl = ylim();
    line(tax([peakloc peakloc]),yl,'color','k','linestyle','--')
    box('off')
    
    subplot(3,1,3)
    set(gca,'fontsize',14)
    plot(tax,stim_mod_tot(cur_cell,:),'o-','linewidth',1)
    xlabel('Time since fixation onset (s)','fontsize',18)
    ylabel('Stimulus Gain','fontsize',18)
    axis tight
    xlim([0 0.4])
    yl = ylim();
    line(tax([peakloc peakloc]),yl,'color','k','linestyle','--')
    box('off')
    fillPage(gcf,'margins',[0 0 0 0],'papersize',[7 12]);
    fname = sprintf('Cell%d_new_modelfit',cellids(cur_cell));
    print('-dpdf','-painters',fname)
    close
    
%% PLOT ALL SU PROPS
close all
simp_tuning = sqrt(gabor_params_fin(1:9,8).^2 + gabor_params_fin(1:9,9).^2);
comp_tuning = gabor_params_fin(1:9,7);
comp_ind = log2(comp_tuning./simp_tuning);
ov_tuning = sqrt(simp_tuning.^2 + comp_tuning.^2);
for i = 1:9
    subplot(3,3,i)
    cur_mask1 = get_pgabor_mask(gabor_params_fin(i,1:6),0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(gabor_params_fin(i,1:6),pi/2,[SDIM SDIM]);
    cur_ovmask = gabor_params_fin(i,8)*cur_mask1 + gabor_params_fin(i,9)*cur_mask2;
    imagesc(xax(xpatch_inds),yax(ypatch_inds),cur_ovmask);
    colormap(gray)
    set(gca,'ydir','normal')
    title(sprintf('Cell%d',cellids(i)),'fontsize',18);
end
fillPage(gcf,'margins',[0 0 0 0],'papersize',[15 15])

figure
imagesc(reshape(ov_tuning(:),3,3)');%set(gca,'ydir','normal')
colorbar
figure
imagesc(reshape(comp_ind(:),3,3)');%set(gca,'ydir','normal')
colorbar

%% PLOT ALL MU PROPS
% close all

figure
simp_tuning = sqrt(gabor_params_fin(10:end,8).^2 + gabor_params_fin(10:end,9).^2);
comp_tuning = gabor_params_fin(10:end,7);
comp_ind = log2(comp_tuning./simp_tuning);
ov_tuning = sqrt(simp_tuning.^2 + comp_tuning.^2);
for i = 1:9
    subplot(3,3,i)
    cur_mask1 = get_pgabor_mask(gabor_params_fin(i+9,1:6),0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(gabor_params_fin(i+9,1:6),pi/2,[SDIM SDIM]);
    cur_ovmask = gabor_params_fin(i+9,8)*cur_mask1 + gabor_params_fin(i+9,9)*cur_mask2;
    imagesc(xax(xpatch_inds),yax(ypatch_inds),cur_ovmask);
    colormap(gray)
    set(gca,'ydir','normal')
    title(sprintf('MUA%d',muaids(i)),'fontsize',18);
end
fillPage(gcf,'margins',[0 0 0 0],'papersize',[15 15])

figure
imagesc(reshape(ov_tuning(:),3,3)');%set(gca,'ydir','normal')
colorbar
figure
imagesc(reshape(comp_ind(:),3,3)');%set(gca,'ydir','normal')
colorbar

%%
clf
cur_cell = 6;
subplot(3,1,1)
plot(tax,sac_trg_rate(cur_cell,:))
xlim([0 0.5])
subplot(3,1,2)
plot(tax,sac_hist_k(cur_cell,:),'r')
hold on
plot(tax,stim_mod_sachist(cur_cell,:),'k')
title('Constant')
xlim([0 0.5])
subplot(3,1,3)
plot(tax,stim_mod_tot(cur_cell,:),'k')
hold on
% plot(tax,sep_mod_lin(cur_cell,:),'b')
% plot(tax,sep_mod_en(cur_cell,:),'r')
hold on
xlim([0 0.5])

%%
used = 1:n_used_cells;
sac_LL_imp = (null_LL(used) - sac_hist_LL(used))/log(2);
sac_xvLL_imp = (null_xvLL(used) - sac_hist_xvLL(used))/log(2);

stim_LL_imp = (null_LL(used) - stim_mod_LL(used))/log(2);
stim_xvLL_imp = (null_xvLL(used) - stim_mod_xvLL(used))/log(2);

% sep_LL_imp = null_LL(used) - sep_mod_LL(used);
% sep_xvLL_imp = null_xvLL(used) - sep_mod_xvLL(used);

%%
% taxold = tax;
% max_tsf = 0.75; nbins = 25;
% used_tsf = time_since_fix;
% used_tsf(used_tsf > max_tsf) = [];
% tax = prctile(used_tsf,linspace(100/nbins,100-100/nbins,nbins));
stim_mod_predrate_hz = stim_mod_predrate/dt;
taxb = [0 tax];
bin_widths = taxb(2:end)-taxb(1:end-1);
rel_spike_binned = nan(n_used_cells,length(used_fixs),length(tax));
combined_spike_binned = nan(n_used_cells,length(used_fixs),length(tax));
sac_spike_binned = nan(n_used_cells,length(used_fixs),length(tax));
for blockid = 1:4;
    fprintf('Block %d of %d\n',blockid,4);
    cur_set = find(blockids(used_fixs) == blockid);
    cur_model_set = find(all_model_blockids==blockid);
    for i = 1:length(cur_set)
        start_T = all_fix_start_times(used_fixs(cur_set(i)));
        end_T = all_fix_stop_times(used_fixs(cur_set(i)));
        cur_tedges = start_T + taxb;
        model_t_interp = round(interp1(all_model_time_axis(cur_model_set),1:length(cur_model_set),start_T+tax));
        cur_tcents = 0.5*cur_tedges(1:end-1)+0.5*cur_tedges(2:end);
        cur_bad = find(cur_tedges > end_T);
        cur_bad2 = find(cur_tcents > end_T | isnan(model_t_interp));
        for c = 1:n_used_cells
            if c <= 9
                temp = histc(Blocks{blockid}.spktimes{cellids(c)},cur_tedges);
            else
                temp = histc(Blocks{blockid}.mutimes{muaids(c-9)},cur_tedges);
            end
            temp(cur_bad) = nan;
            %             rel_spike_binned(c,cur_set(i),:) = temp(1:end-1);
            rel_spike_binned(c,cur_set(i),:) = temp(1:end-1);
            
            model_t_interp(cur_bad2) = 1;
            combined_pred_cnts = stim_mod_predrate(c,cur_model_set(model_t_interp));
            combined_pred_cnts(cur_bad2) = nan;
            combined_spike_binned(c,cur_set(i),:) = combined_pred_cnts.*bin_widths;
            %             combined_spike_binned(c,cur_set(i),:) = combined_pred_cnts;
            
            %             sac_pred_cnts = sac_pred_r(cur_model_set(model_t_interp),c);
            %             sac_pred_cnts(cur_bad2) = nan;
            %             sac_spike_binned(c,cur_set(i),:) = sac_pred_cnts.*bin_widths';
        end
    end
end
% combined_residual = rel_spike_binned - combined_spike_binned;
% combined_resvar = squeeze(nanvar(combined_residual,[],2));
% orivar = squeeze(nanvar(rel_spike_binned,[],2));
% combined_mean = squeeze(nanmean(combined_spike_binned,2));
% ori_mean = squeeze(nanmean(rel_spike_binned,2));
% combined_fano = combined_resvar./ori_mean;
% ori_fano = orivar./ori_mean;

cur_corr = zeros(n_used_cells,length(tax));
for i = 1:length(tax)
    uset = find(~isnan(combined_spike_binned(1,:,i)));
    for c = 1:n_used_cells
        [a,b] = corrcoef(rel_spike_binned(c,uset,i),combined_spike_binned(c,uset,i));
        cur_corr(c,i) = a(2,1);
        cur_p(c,i) = b(2,1);
    end
end

%%
for cur_cell = 1:n_used_cells;
    % cur_cell = 8
    % clf
    [~,peakloc] = max(sac_trg_rate(cur_cell,:));
    
    subplot(3,2,1)
    plot(tax,sac_trg_rate(cur_cell,:))
    axis tight
    xlim([0 0.5])
    yl = ylim();
    line(tax([peakloc peakloc]),yl,'color','k','linestyle','--')
    ylabel('Avreage rate (Hz)','fontsize',14)
    title('Fixation triggered average rate','fontsize',14)
    
    subplot(3,2,3)
    plot(tax,sac_hist_k(cur_cell,:)+sac_hist_const(cur_cell),'r')
    hold on
    plot(tax,stim_mod_sachist(cur_cell,:)+stim_mod_const(cur_cell),'k')
    title('Offset','fontsize',14)
    axis tight
    xlim([0 0.5])
    ylabel('Offset term','fontsize',14)
    legend('Without Stim','With Stim')
    yl = ylim();
    line(tax([peakloc peakloc]),yl,'color','k','linestyle','--')
    
    subplot(3,2,5)
    plot(tax,stim_mod_tot(cur_cell,:),'k')
    title('Stimulus tuning gain','fontsize',14)
    axis tight
    xlim([0 0.5])
    ylabel('Stimulus filter gain (1/z)','fontsiz',14)
    xlabel('Time since fixation onset (s)','fontsize',14)
    yl = ylim();
    line(tax([peakloc peakloc]),yl,'color','k','linestyle','--')
    
    subplot(3,2,2)
    % plot(tax,obs_diff(cur_cell,:))
    hold on
    plot(tax,stim_pred_diff(cur_cell,:),'b')
    % legend('Observed','Model')
    axis tight
    xlim([0 0.5])
    ylabel('RMS rate difference (Hz)','fontsize',14)
    title('Stim-Sac Model Difference','fontsize',14)
    yl = ylim();
    line(tax([peakloc peakloc]),yl,'color','k','linestyle','--')
    
    subplot(3,2,4)
    plot(tax,fano_fac(cur_cell,:),'ko-')
    hold on
    plot(tax,smooth(fano_fac(cur_cell,:),3),'r')
    % plot(tax,combined_fano(cur_cell,:),'r')
    % plot(tax,ori_fano(cur_cell,:),'k')
    axis tight
    xlim([0 0.5])
    ylabel('Fano Factor','fontsize',14)
    title('Fano Factor','fontsize',14)
    yl = ylim();
    line(tax([peakloc peakloc]),yl,'color','k','linestyle','--')
    
    subplot(3,2,6)
    plot(tax,cur_corr(cur_cell,:).^2,'ko-')
    hold on
    plot(tax,smooth(cur_corr(cur_cell,:).^2,3),'r')
    title('Stimulus prediction R2','fontsize',14)
    hold on
    % plot(tax,smooth(cur_corr(cur_cell,:).^2,3,'moving'),'r')
    axis tight
    xlim([0 0.5])
    xlabel('Time since fixation onset (s)','fontsize',14)
    ylabel('Spike count R^2','fontsize',14)
    yl = ylim();
    line(tax([peakloc peakloc]),yl,'color','k','linestyle','--')
    
    
    fillPage(gcf,'margins',[0 0 0 0],'papersize',[12 14]);
    if cur_cell < 10
        print(sprintf('Cell%d_modfit',cellids(cur_cell)),'-dpdf','-painters');
    else
        print(sprintf('MUA%d_modfit',muaids(cur_cell-9)),'-dpdf','-painters');
    end
    close
end

%%
