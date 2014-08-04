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

NT = length(used_fixs);
X_resh = reshape(X,NT,SDIM^2);
dt = .005;

cellids = [1 2 3 4 5 6 7 8 10];
muaids = [1 2 4 5 6 11 12 13 14];

spk_cnts = [spk_cnts(:,cellids) mua_cnts(:,muaids)];

n_used_cells = size(spk_cnts,2);

%%
% orientations = linspace(0,pi-pi/12,12);
% phases = [0 pi/2];
% lambdas = [0.25 0.3 0.35 0.4 0.5 0.65 0.85 1.1 1.4]*Fsd;
% sigmas = lambdas*0.5;
% dxs = [-0.05:0.05:0.35]*Fsd;
% dys = [-0.05 -0.1 -0.15 -0.2]*Fsd;
% n_gabors = length(orientations)*length(phases)*length(lambdas)*length(dxs)*length(dys);
% gabor_bank = zeros(n_gabors,SDIM^2);
% % gabor_bank_params.orientations = zeros(n_gabors,1);
% % gabor_bank_params.phases = zeros(n_gabors,1);
% % gabor_bank_params.lambdas = zeros(n_gabors,1);
% % gabor_bank_params.sigmas = zeros(n_gabors,1);
% % gabor_bank_params.dxs = zeros(n_gabors,1);
% % gabor_bank_params.dys = zeros(n_gabors,1);
% cnt = 1;
% cur_params = zeros(1,6);
% cur_params(6) = 1;
% for jj = 1:length(phases)
%     for ii = 1:length(orientations)
%         cur_params(3) = orientations(ii);
%         for kk = 1:length(lambdas)
%             cur_params(4) = lambdas(kk);
%             cur_params(5) = sigmas(kk);
%             for ll = 1:length(dxs)
%                 cur_params(1) = dxs(ll);
%                 for mm = 1:length(dys)
%                     cur_params(2) = dys(mm);
%                     temp = get_pgabor_mask(cur_params,phases(jj),[SDIM SDIM]);
%                     gabor_bank(cnt,:) = temp(:);
%                     gabor_bank_params(cnt).orientations = orientations(ii);
%                     gabor_bank_params(cnt).phases = phases(jj);
%                     gabor_bank_params(cnt).lambdas = lambdas(kk);
%                     gabor_bank_params(cnt).sigmas = sigmas(kk);
%                     gabor_bank_params(cnt).dxs = dxs(ll);
%                     gabor_bank_params(cnt).dys = dys(mm);
%                     cnt = cnt + 1;
%                 end
%             end
%         end
%     end
% end
% gabor_bank_output = X_resh*gabor_bank';
% n_unique_gabors = n_gabors/2;
% gabor_bank_energy = sqrt(gabor_bank_output(:,1:n_unique_gabors).^2 + gabor_bank_output(:,n_unique_gabors+1:end).^2);
% gabor_bank_output = [gabor_bank_output -gabor_bank_output gabor_bank_energy];
% gabor_bank_output = gabor_bank_energy;
% gabor_bank = gabor_bank(1:n_unique_gabors,:);
%
% %feed through threshold linear internal NL
% gabor_bank_output(gabor_bank_output < 0) = 0;
%
% %normalize gabor bank outputs
% gabor_bank_means = mean(gabor_bank_output);
% gabor_bank_stds = std(gabor_bank_output);
% gabor_bank_output = bsxfun(@minus,gabor_bank_output,gabor_bank_means);
% gabor_bank_output = bsxfun(@rdivide,gabor_bank_output,gabor_bank_stds);

%%
clear gabor_bank* gabor_bank_params 
for c = 1:9
    c
    orientations = linspace(0,pi-pi/15,15);
    phases = [0 pi/2 pi 3*pi/2];
    lambdas = [0.25 0.3 0.35 0.4 0.5 0.6 0.7 0.85 1.1 1.2 1.4]*Fsd;
    sigmas = [0.5];
    n_gabors = length(orientations)*length(phases)*length(lambdas)*length(sigmas);
    gabor_bank{c} = zeros(n_gabors,SDIM^2);
    gabor_bank_orientations = zeros(n_gabors,1);
    gabor_bank_phases = zeros(n_gabors,1);
    gabor_bank_lambdas = zeros(n_gabors,1);
    gabor_bank_sigmas = zeros(n_gabors,1);
    gabor_bank_dxs(:,c) = ones(n_gabors,1)*gabor_params_fin(c,1);
    gabor_bank_dys(:,c) = ones(n_gabors,1)*gabor_params_fin(c,2);
    cnt = 1;
    cur_params = zeros(1,6);
    cur_params(6) = 1;
    cur_params(1) = gabor_params_fin(c,1);
    cur_params(2) = gabor_params_fin(c,2);
    for jj = 1:length(phases)
        for ii = 1:length(orientations)
            cur_params(3) = orientations(ii);
            for kk = 1:length(lambdas)
                for ll = 1:length(sigmas)
                    cur_params(4) = lambdas(kk);
                    cur_params(5) = sigmas(ll)*lambdas(kk);
                    temp = get_pgabor_mask(cur_params,phases(jj),[SDIM SDIM]);
                    gabor_bank{c}(cnt,:) = temp(:);
                    gabor_bank_orientations(cnt) = orientations(ii);
                    gabor_bank_phases(cnt) = phases(jj);
                    gabor_bank_lambdas(cnt) = lambdas(kk);
                    gabor_bank_sigmas(cnt) = sigmas(ll)*lambdas(kk);
                    cnt = cnt + 1;
                end
            end
        end
    end
    %%
    gabor_bank_output{c} = X_resh*gabor_bank{c}';
    n_unique_gabors = n_gabors/4;
%     gabor_bank_energy = gabor_bank_output{c}(:,1:n_unique_gabors).^2;
%     gabor_bank_energy2 = gabor_bank_output{c}(:,(n_unique_gabors+1):(2*n_unique_gabors)).^2;
%     gabor_bank_energy = sqrt(gabor_bank_energy + gabor_bank_energy2);
%     
%     gabor_bank_output{c} = [gabor_bank_output{c} gabor_bank_energy];
%     gabor_bank{c} = [gabor_bank{c}; gabor_bank{c}(1:n_unique_gabors,:)];
    
    
    %normalize gabor bank outputs
    gabor_bank_means = mean(gabor_bank_output{c});
    gabor_bank_stds = std(gabor_bank_output{c});
    gabor_bank_output{c} = bsxfun(@minus,gabor_bank_output{c},gabor_bank_means);
    gabor_bank_output{c} = bsxfun(@rdivide,gabor_bank_output{c},gabor_bank_stds);
    
    %feed through threshold linear internal NL
    %     gabor_bank_output{c}(gabor_bank_output{c} < 0) = 0;
    
    %feed through threshold quadratic internal NL
    gabor_bank_output{c} = gabor_bank_output{c}.^2;
    gabor_bank_output{c}(gabor_bank_output{c} < 0) = 0;

end

% %account for added energy filters
% is_energy_filt = [zeros(n_gabors,1); ones(n_unique_gabors,1)];
% gabor_bank_orientations = [gabor_bank_orientations; gabor_bank_orientations(1:n_unique_gabors)];
% gabor_bank_phases = [gabor_bank_phases; gabor_bank_phases(1:n_unique_gabors)];
% gabor_bank_lambdas = [gabor_bank_lambdas; gabor_bank_lambdas(1:n_unique_gabors)];
% gabor_bank_sigmas = [gabor_bank_sigmas; gabor_bank_sigmas(1:n_unique_gabors)];
% gabor_bank_dxs = [gabor_bank_dxs; gabor_bank_dxs(1:n_unique_gabors,:)];
% gabor_bank_dys = [gabor_bank_dys; gabor_bank_dys(1:n_unique_gabors,:)];

%%
gabor_bank_old = zeros(n_used_cells,2,SDIM,SDIM);
for n = 1:n_used_cells
    gabor_bank_old(n,1,:,:) = get_pgabor_mask(gabor_params_fin(n,1:6),0,[SDIM SDIM]);
    gabor_bank_old(n,2,:,:) = get_pgabor_mask(gabor_params_fin(n,1:6),pi/2,[SDIM SDIM]);
end
gabor_bank_old_resh = reshape(gabor_bank_old,[size(gabor_bank_old,1),2,SDIM^2]);
w0 = gabor_params_fin(1:n_used_cells,8)';
w90 = gabor_params_fin(1:n_used_cells,9)';
wen = gabor_params_fin(1:n_used_cells,7)';

%%
fix_stimout = zeros(n_used_cells,length(used_fixs));
all_model_time_axis = [];
all_model_blockids = [];
all_model_fixids = [];
spikes_binned = [];
time_since_fix = [];
stim_outs = zeros(n_used_cells,length(used_fixs));
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
        for c = 1:9
            if c <= 9
                temp = histc(Blocks{blockid}.spktimes{cellids(c)},cur_tedges);
            else
                temp = histc(Blocks{blockid}.mutimes{muaids(c-9)},cur_tedges);
            end
            temp_binned(c,:) = temp(1:end-1);
        end
        spikes_binned = [spikes_binned; temp_binned'];
        
        time_since_fix = [time_since_fix cur_tcents-start_T];
        
        cur_gout1 = X_resh(cur_set(i),:)*squeeze(gabor_bank_old_resh(:,1,:))';
        cur_gout2 = X_resh(cur_set(i),:)*squeeze(gabor_bank_old_resh(:,2,:))';
        cur_energy = sqrt(cur_gout1.^2+cur_gout2.^2);
        stim_outs(:,cur_set(i)) = w0.*cur_gout1 + w90.*cur_gout2 + wen.*cur_energy;
        
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
stim_kern = cell(n_used_cells,1);
sac_hist_only = zeros(n_used_cells,length(tax));
lambdaW = 1000;
slope2_lambda = 70;
n_gabors = size(gabor_bank_output{2},2);
for c = 1:9
    fprintf('Cell %d of %d\n',c,9);
    
    Robs = spk_cnts(used_fixs,c)';
    w_init = zeros(n_gabors,1);
    const = 0;
    W0 = [w_init; const];
    lambda = ones(1,n_gabors)*lambdaW;
    lambda = [lambda 0]/sum(Robs);
    [params LL] = L1General2_PSSas(@(K) gabor_weights_LLinternal(K, Robs, gabor_bank_output{c}),W0,lambda',[],1,500);
    gabor_weights = params(1:end-1);
    
    un_spk_cnts = length(unique(spikes_binned(:,c))) - 1;
    spkbs = [];
    for i = 1:un_spk_cnts
        curset = find(spikes_binned(:,c) == i);
        spkbs = [spkbs; repmat(curset,i,1)];
    end
    spkbs = sort(spkbs);
    
    
        %% old combined model
        Xmat = [Tmat bsxfun(@times,Tmat,stim_outs(c,all_model_fixids)')];
        W0 = zeros(size(Xmat,2),1);
    
        %initialize parameters
        silent = 1;
        lamrange = [];
        lamrange2 = [slope2_lambda 1 length(tax)];
        nlxrange = [tax];
        lamrange2 = [lamrange2; slope2_lambda length(tax)+1 length(tax)*2];
        nlxrange = [nlxrange tax];
        Pcon = [];
        Pmono = [];
        hold_const = [];
        NLtype = 0;
        llist = [];
    
        [fitp,grad] = GLMsolve_jmm(Xmat, spkbs, W0, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const, NLtype, 1e-6,nlxrange);
        old_sac_hist(c,:) = fitp.k(1:length(tax));
        old_stim_kern(c,:) = fitp.k((length(tax)+1):end-1);
        old_combined_LP(c) = fitp.LL;
        old_combined_const(c) = fitp.k(end);
        gout = Xmat*fitp.k(1:end-1)+fitp.k(end);
        old_combined_rpred(c,:) = log(1+exp(gout));
        old_combined_LL(c) = -sum(spikes_binned(:,c)'.*log(old_combined_rpred(c,:))-old_combined_rpred(c,:))/sum(spikes_binned(:,c));
    
        %% sac only model
        [fitp,grad] = GLMsolve_jmm( Xmat(:,1:length(tax)), spkbs, W0(1:length(tax)), silent, lamrange, lamrange2(1,:), Pcon, Pmono, llist, hold_const, NLtype, 1e-6,nlxrange(1:length(tax)));
    
        sac_hist_only(c,:) = fitp.k(1:length(tax));
        sachist_LP(c) = fitp.LL;
        sachist_const(c) = fitp.k(end);
        gout = Tmat*fitp.k(1:end-1)+fitp.k(end);
        sachist_rpred(c,:) = log(1+exp(gout));
        sachist_LL(c) = -sum(spikes_binned(:,c)'.*log(sachist_rpred(c,:))-sachist_rpred(c,:))/sum(spikes_binned(:,c));
    
    %% new combined model
    Xmat = Tmat;
    used_gabors{c} = find(gabor_weights ~= 0);
    if ~isempty(used_gabors{c})
        used_gabor_outs = gabor_bank_output{c}(:,used_gabors{c});
        cur_n_gabors = length(used_gabors{c});
        for i = 1:cur_n_gabors
            Xmat = [Xmat bsxfun(@times,Tmat,used_gabor_outs(all_model_fixids,i))];
        end
        
        W0 = zeros(size(Xmat,2),1);
        
        %initialize parameters
        silent = 1;
        lamrange = [];
        lamrange2 = [slope2_lambda 1 length(tax)];
        nlxrange = [tax];
        for i = 1:cur_n_gabors
            lamrange2 = [lamrange2; slope2_lambda length(tax)*i+1 length(tax)*(i+1)];
            nlxrange = [nlxrange tax];
        end
        Pcon = [];
        Pmono = [];
        hold_const = [];
        NLtype = 0;
        llist = [];
        
        [fitp,grad] = GLMsolve_jmm( Xmat, spkbs, W0, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const, NLtype, 1e-6,nlxrange);
        
        sac_hist(c,:) = fitp.k(1:length(tax));
        stim_kern{c} = zeros(cur_n_gabors,length(tax));
        for i = 1:cur_n_gabors
            stim_kern{c}(i,:) = fitp.k((length(tax)*i)+(1:length(tax)));
        end
        combined_LP(c) = fitp.LL;
        combined_const(c) = fitp.k(end);
        gout = Xmat*fitp.k(1:end-1)+fitp.k(end);
        combined_rpred(c,:) = log(1+exp(gout));
        combined_LL(c) = -sum(spikes_binned(:,c)'.*log(combined_rpred(c,:))-combined_rpred(c,:))/sum(spikes_binned(:,c));
        
    end
    
    R_avg = mean(spikes_binned(:,c));
    LL_nul(c) = sum(-bsxfun(@times,spikes_binned(:,c),log(R_avg)) + R_avg)/sum(spikes_binned(:,c));
    
end

%%
cd /Users/James/James_scripts/bruce/modelfits
save tempmodfits_v2 *LL LL_nul combined_* old_* sac* tax time_since_fix
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

%%
close all
cur_c = 9;
tempmat = gabor_bank{cur_c}(used_gabors{cur_c},:);
cur_ng = length(used_gabors{cur_c});
for i = 1:cur_ng
    
    subplot(cur_ng,2,(i-1)*2+1)
    imagesc(reshape(tempmat(i,:),SDIM,SDIM))
    set(gca,'ydir','normal')
    subplot(cur_ng,2,(i-1)*2+2)
    plot(tax,stim_kern{cur_c}(i,:))
    axis tight
    xlim([0 0.5])
end

figure
subplot(3,1,1)
plot(tax,old_sac_hist(cur_c,:))
xlim([0 0.5])
subplot(3,1,2)
plot(tax,old_stim_kern(cur_c,:))
xlim([0 0.5])
subplot(3,1,3)
imagesc(squeeze(gabor_bank_old(cur_c,1,:,:)))
set(gca,'ydir','normal')
