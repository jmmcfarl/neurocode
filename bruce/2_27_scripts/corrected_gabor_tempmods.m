clear all
% close all

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

%%
% for n = 1:9
%
%     cur_gabor_params = gabor_params_fin(n,:);
%     comp_ratio(n) = log(cur_gabor_params(7)/norm(cur_gabor_params(8:9)));
%
% %     cur_gabor_params(1:5) = ip(n,:);
%     cur_mask1 = get_pgabor_mask(cur_gabor_params(1:6),0,[SDIM SDIM]);
%     cur_mask2 = get_pgabor_mask(cur_gabor_params(1:6),pi/2,[SDIM SDIM]);
%
%     cur_mask = cur_gabor_params(8)*cur_mask1 + cur_gabor_params(9)*cur_mask2;
%     subplot(3,3,n)
%     imagesc(xax(xpatch_inds),yax(ypatch_inds),cur_mask); colormap(gray);
%     axis square
%     set(gca,'ydir','normal');
%     set(gca,'fontsize',12)
%     xlabel('Horizontal position (degrees)','fontsize',14)
%     ylabel('Vertical position (degrees)','fontsize',14)
%     title(sprintf('Cell%d, CI: %.2f PO: %.1f',cellids(n),comp_ratio(n),180/pi*cur_gabor_params(3)),'fontsize',14);
% %     fname = sprintf('Fit_Gabor_filter_cell%d',cellids(n));
% %     print('-dpdf',fname);
% %     close
% end
% fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [11 11]);
% %%
% muaids = [1 2 4 5 6 11 12 13 14];
%
% for n = 10:18
%
%     cur_gabor_params = gabor_params_fin(n,:);
%     comp_ratio(n) = log(cur_gabor_params(7)/norm(cur_gabor_params(8:9)));
%
% %     cur_gabor_params(1:5) = ip(n,:);
%     cur_mask1 = get_pgabor_mask(cur_gabor_params(1:6),0,[SDIM SDIM]);
%     cur_mask2 = get_pgabor_mask(cur_gabor_params(1:6),pi/2,[SDIM SDIM]);
%
%     cur_mask = cur_gabor_params(8)*cur_mask1 + cur_gabor_params(9)*cur_mask2;
%     subplot(3,3,n-9)
%     imagesc(xax(xpatch_inds),yax(ypatch_inds),cur_mask); colormap(gray);
%     axis square
%     set(gca,'ydir','normal');
%     set(gca,'fontsize',12)
%     xlabel('Horizontal position (degrees)','fontsize',14)
%     ylabel('Vertical position (degrees)','fontsize',14)
%     title(sprintf('MUA%d, CI: %.2f PO: %.1f',muaids(n-9),comp_ratio(n),180/pi*cur_gabor_params(3)),'fontsize',14);
% %     fname = sprintf('Fit_Gabor_filter_cell%d',cellids(n));
% %     print('-dpdf',fname);
% %     close
% end
% fp = fillPage(gcf, 'margins', [0 0 0 0], 'papersize', [11 11]);

%%
NT = length(used_fixs);
X_resh = reshape(X,NT,SDIM^2);
dt = .005;
spikes_binned = [];
time_since_fix = [];
stim_outs = [];

%%
cellids = [1 2 3 4 5 6 7 8 10];
% n_used_cells = length(cellids);
% for c = 1:n_used_cells
%     fprintf('Fitting cell %d of %d\n',c,n_used_cells);
%     hold_const = [0 0 1 1 1 1 1 1 1 1];
%     LB = [-10 -10 6 2 0.5 -Inf -Inf -Inf -Inf];
%     UB = [10 10 pi SDIM(1)/2 SDIM(1)/6 2 Inf Inf Inf Inf];
%     [gabor_params_cor(c,:),LL_cor(c)] = fit_gabor_params(gabor_params_fin(c,:),X_resh,spk_cnts(used_fixs,cellids(c)),[SDIM SDIM],hold_const,LB,UB);
%     R_avg(c) = mean(spk_cnts(used_fixs,cellids(c)));
%     LL_nul(c) = sum(-bsxfun(@times,spk_cnts(used_fixs,cellids(c)),log(R_avg(c))) + R_avg(c))/sum(spk_cnts(used_fixs,cellids(c)));
% end

muaids = [1 2 4 5 6 11 12 13 14];

spk_cnts = [spk_cnts(:,cellids) mua_cnts(:,muaids)];

n_used_cells = size(spk_cnts,2);

%%
gabor_bank = zeros(n_used_cells,2,SDIM,SDIM);
for n = 1:n_used_cells
    gabor_bank(n,1,:,:) = get_pgabor_mask(gabor_params_fin(n,1:6),0,[SDIM SDIM]);
    gabor_bank(n,2,:,:) = get_pgabor_mask(gabor_params_fin(n,1:6),pi/2,[SDIM SDIM]);
end
gabor_bank_resh = reshape(gabor_bank,[size(gabor_bank,1),2,SDIM^2]);

w0 = gabor_params_fin(1:n_used_cells,8)';
w90 = gabor_params_fin(1:n_used_cells,9)';
wen = gabor_params_fin(1:n_used_cells,7)';
%%
fix_stimout = zeros(n_used_cells,length(used_fixs));
for blockid = 1:4;
    fprintf('Block %d of %d\n',blockid,4);
    cur_set = find(blockids(used_fixs) == blockid);
    
    for i = 1:length(cur_set)
        start_T = all_fix_start_times(used_fixs(cur_set(i)));
        end_T = all_fix_stop_times(used_fixs(cur_set(i)));
        cur_tedges = start_T:dt:end_T;
        cur_tcents = 0.5*cur_tedges(1:end-1)+0.5*cur_tedges(2:end);
        temp_binned = zeros(10,length(cur_tcents));
        temp_modouts = zeros(10,length(cur_tcents));
        for c = 1:n_used_cells
            if c <= 9
                temp = histc(Blocks{blockid}.spktimes{cellids(c)},cur_tedges);
            else
                temp = histc(Blocks{blockid}.mutimes{muaids(c-9)},cur_tedges);
            end
            temp_binned(c,:) = temp(1:end-1);
        end
        spikes_binned = [spikes_binned; temp_binned'];
        
        time_since_fix = [time_since_fix cur_tcents-start_T];
        
        cur_gout1 = X_resh(cur_set(i),:)*squeeze(gabor_bank_resh(:,1,:))';
        cur_gout2 = X_resh(cur_set(i),:)*squeeze(gabor_bank_resh(:,2,:))';
        cur_energy = sqrt(cur_gout1.^2+cur_gout2.^2);
        cur_stim_out = w0.*cur_gout1 + w90.*cur_gout2 + wen.*cur_energy;
        
        fix_stimout(:,cur_set(i)) = cur_stim_out;
        stim_outs = [stim_outs; repmat(cur_stim_out,length(cur_tcents),1)];
    end
end

%%
all_model_time_axis = [];
all_model_blockids = [];
for blockid = 1:4
    fprintf('Block %d of %d\n',blockid,4);
    cur_set = find(blockids(used_fixs) == blockid);
    
    for i = 1:length(cur_set)
        start_T = all_fix_start_times(used_fixs(cur_set(i)));
        end_T = all_fix_stop_times(used_fixs(cur_set(i)));
        cur_tedges = start_T:dt:end_T;
        cur_tcents = 0.5*cur_tedges(1:end-1)+0.5*cur_tedges(2:end);
        all_model_time_axis = [all_model_time_axis cur_tcents];
        all_model_blockids = [all_model_blockids blockid*ones(1,length(cur_tcents))];
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
sac_hist = zeros(n_used_cells,length(tax));
sac_hist_only = zeros(n_used_cells,length(tax));
stim_kern = zeros(n_used_cells,length(tax));
stim_kern_only = zeros(n_used_cells,length(tax));
combined_LL = zeros(n_used_cells,1);
stimonly_LL = zeros(n_used_cells,1);
saconly_LL = zeros(n_used_cells,1);
LL_nul = zeros(n_used_cells,1);
for c = 1:9;
    fprintf('Cell %d of %d\n',c,n_used_cells);
    un_spk_cnts = length(unique(spikes_binned(:,c))) - 1;
    spkbs = [];
    for i = 1:un_spk_cnts
        curset = find(spikes_binned(:,c) == i);
        spkbs = [spkbs; repmat(curset,i,1)];
    end
    spkbs = sort(spkbs);
    
    Xmat = [Tmat bsxfun(@times,Tmat,stim_outs(:,c))];
    
    %%
    W0 = zeros(size(Xmat,2),1);
    
    %initialize parameters
    silent = 1;
    slope2_lambda = 50;
    lamrange = [];
    lamrange2 = [slope2_lambda 1 length(tax); slope2_lambda (length(tax)+1) (length(W0)-1)];
    Pcon = [];
    Pmono = [];
    hold_const = [];
    NLtype = 0;
    llist = [];
    
    [fitp,grad] = GLMsolve_jmm( Xmat, spkbs, W0, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const, NLtype, 1e-6, [tax'; tax'] );
    sac_hist(c,:) = fitp.k(1:length(tax));
    stim_kern(c,:) = fitp.k((length(tax)+1):end-1);
    combined_LP(c) = fitp.LL;
    combined_const(c) = fitp.k(end);
    gout = Xmat*fitp.k(1:end-1)+fitp.k(end);
    combined_rpred(c,:) = log(1+exp(gout));
    combined_LL(c) = -sum(spikes_binned(:,c)'.*log(combined_rpred(c,:))-combined_rpred(c,:))/sum(spikes_binned(:,c));

    [fitp,grad] = GLMsolve_jmm(Xmat(:,1:length(tax)), spkbs, W0(1:length(tax)), silent, lamrange, lamrange2(1,:), Pcon, Pmono, llist, hold_const, NLtype, 1e-6, tax');
    sac_hist_only(c,:) = fitp.k(1:length(tax));
    saconly_LP(c) = fitp.LL;
    saconly_const(c) = fitp.k(end);
    gout = Tmat*fitp.k(1:end-1)+fitp.k(end);
    saconly_rpred(c,:) = log(1+exp(gout));
    saconly_LL(c) = -sum(spikes_binned(:,c)'.*log(saconly_rpred(c,:))-saconly_rpred(c,:))/sum(spikes_binned(:,c));
    
    [fitp,grad] = GLMsolve_jmm(Xmat(:,(length(tax)+1):end), spkbs, W0(1:length(tax)), silent, lamrange, lamrange2(1,:), Pcon, Pmono, llist, hold_const, NLtype, 1e-6, tax');
    stim_kern_only(c,:) = fitp.k(1:length(tax));
    stimonly_LP(c) = fitp.LL;
    stimonly_const(c) = fitp.k(end);
    gout = Xmat(:,(length(tax)+1):end)*fitp.k(1:end-1)+fitp.k(end);
    stimonly_rpred(c,:) = log(1+exp(gout));
    stimonly_LL(c) = -sum(spikes_binned(:,c)'.*log(stimonly_rpred(c,:))-stimonly_rpred(c,:))/sum(spikes_binned(:,c));
    
    R_avg = mean(spikes_binned(:,c));
    LL_nul(c) = sum(-bsxfun(@times,spikes_binned(:,c),log(R_avg)) + R_avg)/sum(spikes_binned(:,c));
end
%%

sac_trg_rate = zeros(n_used_cells,length(tax));
taxb = [0 tax];
for i = 1:length(tax)
    curset = find(time_since_fix >= taxb(i) & time_since_fix < taxb(i+1));
    for c = 1:n_used_cells
        sac_trg_rate(c,i) = mean(spikes_binned(curset,c))/dt;
    end
    n_occ(i) = length(curset);
end
avg_rates = mean(spikes_binned)/dt;
combined_imp = combined_LL - LL_nul;
saconly_imp = saconly_LL - LL_nul;
stimonly_imp = stimonly_LL - LL_nul;
stimspecific_imp = combined_LL - saconly_LL;

%%
sac_pred_r = Tmat*sac_hist_only';
sac_pred_r = bsxfun(@plus,sac_pred_r,saconly_const);
sac_pred_r = log(1+exp(sac_pred_r))/dt;

for c = 1:n_used_cells
    cur_X = bsxfun(@times,Tmat,stim_outs(:,c));
    combined_pred_r(:,c) = Tmat*sac_hist(c,:)' + cur_X*stim_kern(c,:)';
    combined_pred_r(:,c) = bsxfun(@plus,combined_pred_r(:,c),combined_const(c));
end
combined_pred_r = log(1+exp(combined_pred_r))/dt;

squared_diff = (combined_pred_r - sac_pred_r).^2;

clear X
% cd /Users/James/James_scripts/bruce/modelfits/
% save fixation_gabor_models

%%
rms_tax = zeros(length(tax),n_used_cells);
std_tax = zeros(length(tax),n_used_cells);
info_sachist_tax = zeros(length(tax),n_used_cells);
info_combined_tax = zeros(length(tax),n_used_cells);
taxb = [0 tax];
for i = 1:length(tax)
    curset = find(time_since_fix >= taxb(i) & time_since_fix < taxb(i+1));
    n_occ(i) = length(curset);
    rms_tax(i,:) = mean(squared_diff(curset,:));
    std_tax(i,:) = std(spikes_binned(curset,:)/dt);
    for c = 1:n_used_cells
        %         uset = curset(spikes_binned(curset,c) > 0);
        %         info_sachist_tax(i,c) = mean(spikes_binned(uset,c)./sac_pred_r(uset,c).*log2(spikes_binned(uset,c)./sac_pred_r(uset,c)));
        %         info_combined_tax(i,c) = mean(spikes_binned(uset,c)./combined_pred_r(uset,c).*log2(spikes_binned(uset,c)./combined_pred_r(uset,c)));
        info_sachist_tax(i,c) = mean(sac_pred_r(curset,c).*log2(bsxfun(@rdivide,sac_pred_r(curset,c),mean(spikes_binned(curset,c)))));
        info_combined_tax(i,c) = mean(combined_pred_r(curset,c).*log2(bsxfun(@rdivide,combined_pred_r(curset,c),mean(spikes_binned(curset,c)))));
    end
end
rms_tax = sqrt(rms_tax);
% std_tax = std_tax/dt;

%%
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
            rel_spike_binned(c,cur_set(i),:) = temp(1:end-1);
            
            model_t_interp(cur_bad2) = 1;
            combined_pred_cnts = combined_pred_r(cur_model_set(model_t_interp),c);
            combined_pred_cnts(cur_bad2) = nan;
            combined_spike_binned(c,cur_set(i),:) = combined_pred_cnts.*bin_widths';
            sac_pred_cnts = sac_pred_r(cur_model_set(model_t_interp),c);
            sac_pred_cnts(cur_bad2) = nan;
            sac_spike_binned(c,cur_set(i),:) = sac_pred_cnts.*bin_widths';
        end
    end
end

%%
% cur_var = squeeze(nanvar(rel_spike_binned(c,:,:)));
% resvar_sac =
% squeeze(nanvar(rel_spike_binned(c,:,:)-sac_spike_binned(c,:,:))); 
combined_residual = rel_spike_binned - combined_spike_binned;
resvar_stim = squeeze(nanvar(rel_spike_binned-combined_spike_binned,[],2))./squeeze(nanvar(rel_spike_binned,[],2));

for i = 1:length(tax)
    uset = find(~isnan(combined_spike_binned(1,:,i)));
    for c = 1:n_used_cells
        [a,b] = corrcoef(rel_spike_binned(c,uset,i),combined_spike_binned(c,uset,i));
        cur_corr(c,i) = a(2,1);
        cur_p(c,i) = b(2,1);
    end
end
%%
xl = [0 0.4];
c = 4;
clf
% figure
subplot(2,2,1)
plot(tax,sac_trg_rate(c,:))
axis tight
xlim(xl)
title(sprintf('Avg Rate: %.2f',avg_rates(c)));
subplot(2,2,2)
plot(tax(1:end-1),sac_hist(c,1:end-1))
hold on
plot(tax(1:end-1),sac_hist_only(c,1:end-1),'k')
% set(gca,'xscale','log')
xlim(tax([1 end]))
axis tight
xlim(xl)
subplot(2,2,3)
plot(tax(1:end-1),stim_kern(c,1:end-1))
hold on
% plot(tax,stim_kern_only(c,:),'k')
% set(gca,'xscale','log')
xlim(tax([1 end]))
axis tight
xlim(xl)
subplot(2,2,4)
% plot(tax,rms_tax(:,c))
hold on
% plot(tax,std_tax(:,c),'r')
% plot(tax,rms_tax(:,c)./std_tax(:,c),'k')
plot(tax,cur_corr(c,:))
xlim(xl)
% subplot(3,2,5)
% plot(tax,info_sachist_tax(:,c))
% hold on
% plot(tax,info_combined_tax(:,c),'r')
% xlim(xl)
% subplot(3,2,6)
% plot(tax,info_combined_tax(:,c)-info_sachist_tax(:,c))
% xlim(xl)

%%
