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
dt = .01;

cellids = [1 2 3 4 5 6 7 8 10];
muaids = [1 2 4 5 6 11 12 13 14];

spk_cnts = [spk_cnts(:,cellids) mua_cnts(:,muaids)];

n_used_cells = size(spk_cnts,2);

%%
clear gabor_bank* gabor_bank_params
for c = 1:9
    c
    phases = [0 pi/2 pi 3*pi/2];
    lambdas = [0.25 0.3 0.4 0.5 0.6 0.8 1.1]*Fsd;
    sigmas = [0.5];
    dxs = [0]*Fsd;
    dys = [0]*Fsd;
    n_gabors = length(phases)*length(lambdas)*length(sigmas)*length(dxs)*length(dys);
    gabor_bank{c} = zeros(n_gabors,SDIM^2);
    gabor_bank_phases = zeros(n_gabors,1);
    gabor_bank_lambdas = zeros(n_gabors,1);
    gabor_bank_sigmas = zeros(n_gabors,1);
    gabor_bank_dxs(:,c) = ones(n_gabors,1)*gabor_params_fin(c,1);
    gabor_bank_dys(:,c) = ones(n_gabors,1)*gabor_params_fin(c,2);
    gabor_bank_orientations(:,c) = ones(n_gabors,1)*gabor_params_fin(c,3);
    cnt = 1;
    cur_params = zeros(1,6);
    cur_params(6) = 1;
    cur_params(1) = gabor_params_fin(c,1);
    cur_params(2) = gabor_params_fin(c,2);
    for jj = 1:length(phases)
        cur_params(3) = gabor_params_fin(c,3);
        for kk = 1:length(lambdas)
            for ll = 1:length(sigmas)
                for mm = 1:length(dxs)
                    for nn = 1:length(dys)
                        cur_params(4) = lambdas(kk);
                        cur_params(5) = sigmas(ll)*lambdas(kk);
                        temp = get_pgabor_mask(cur_params,phases(jj),[SDIM SDIM]);
                        gabor_bank{c}(cnt,:) = temp(:);
                        gabor_bank_orientations(cnt,c) = cur_params(3);
                        gabor_bank_phases(cnt) = phases(jj);
                        gabor_bank_lambdas(cnt) = lambdas(kk);
                        gabor_bank_sigmas(cnt) = sigmas(ll)*lambdas(kk);
                        gabor_bank_dxs(cnt,c) = gabor_bank_dxs(cnt,c) + dxs(mm);
                        gabor_bank_dys(cnt,c) = gabor_bank_dys(cnt,c) + dys(nn);
                        cnt = cnt + 1;
                    end
                end
            end
        end
%         cur_params(3) = gabor_params_fin(c,3)+pi/2;
%         for kk = 1:length(lambdas)
%             for ll = 1:length(sigmas)
%                 for mm = 1:length(dxs)
%                     for nn = 1:length(dys)
%                         cur_params(4) = lambdas(kk);
%                         cur_params(5) = sigmas(ll)*lambdas(kk);
%                         temp = get_pgabor_mask(cur_params,phases(jj),[SDIM SDIM]);
%                         gabor_bank{c}(cnt,:) = temp(:);
%                         gabor_bank_orientations(cnt,c) = cur_params(3);
%                         gabor_bank_phases(cnt) = phases(jj);
%                         gabor_bank_lambdas(cnt) = lambdas(kk);
%                         gabor_bank_sigmas(cnt) = sigmas(ll)*lambdas(kk);
%                         gabor_bank_dxs(cnt,c) = gabor_bank_dxs(cnt,c) + dxs(mm);
%                         gabor_bank_dys(cnt,c) = gabor_bank_dys(cnt,c) + dys(nn);
%                         cnt = cnt + 1;
%                     end
%                 end
%             end
%         end
    end
    %%
    gabor_bank_output{c} = X_resh*gabor_bank{c}';
    n_unique_gabors = n_gabors/4;
    
    %normalize gabor bank outputs
    gabor_bank_means = mean(gabor_bank_output{c});
    gabor_bank_stds = std(gabor_bank_output{c});
    gabor_bank_output{c} = bsxfun(@minus,gabor_bank_output{c},gabor_bank_means);
    gabor_bank_output{c} = bsxfun(@rdivide,gabor_bank_output{c},gabor_bank_stds);
    
    gabor_bank_energy = gabor_bank_output{c}(:,1:n_unique_gabors).^2;
    gabor_bank_energy2 = gabor_bank_output{c}(:,(n_unique_gabors+1):(2*n_unique_gabors)).^2;
    gabor_bank_energy = sqrt(gabor_bank_energy + gabor_bank_energy2);
    
    %     %feed through threshold quadratic internal NL
    %     gabor_bank_output{c} = gabor_bank_output{c}.^2;
    %     gabor_bank_output{c}(gabor_bank_output{c} < 0) = 0;
    
    %feed through threshold linear internal NL
    gabor_bank_output{c}(gabor_bank_output{c} < 0) = 0;
    
    gabor_bank_output{c} = [gabor_bank_output{c} gabor_bank_energy];
    gabor_bank{c} = [gabor_bank{c}; gabor_bank{c}(1:n_unique_gabors,:)];
    
    
end

%account for added energy filters
is_energy_filt = [zeros(n_gabors,1); ones(n_unique_gabors,1)];
gabor_bank_phases = [gabor_bank_phases; gabor_bank_phases(1:n_unique_gabors)];
gabor_bank_lambdas = [gabor_bank_lambdas; gabor_bank_lambdas(1:n_unique_gabors)];
gabor_bank_sigmas = [gabor_bank_sigmas; gabor_bank_sigmas(1:n_unique_gabors)];
gabor_bank_dxs = [gabor_bank_dxs; gabor_bank_dxs(1:n_unique_gabors,:)];
gabor_bank_dys = [gabor_bank_dys; gabor_bank_dys(1:n_unique_gabors,:)];
gabor_bank_orientations = [gabor_bank_orientations; gabor_bank_orientations(1:n_unique_gabors,:)];


%%
fix_stimout = zeros(n_used_cells,length(used_fixs));
all_model_time_axis = [];
all_model_blockids = [];
all_model_fixids = [];
spikes_binned = [];
time_since_fix = [];
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
        
    end
end

%% Compute TBR time-since fix onset
max_tsf = 0.75; nbins = 25;
used_tsf = time_since_fix;
used_tsf(used_tsf > max_tsf) = [];
tax = prctile(used_tsf,linspace(100/nbins,100-100/nbins,nbins));

% tax = logspace(log10(dt),log10(max_tsf),nbins);
Tmat = tbrep(time_since_fix,tax);

%%
xv_frac = 0.2;
rperm = randperm(length(used_fixs));
xv_fixs = rperm(1:round(xv_frac*length(used_fixs)));
xv_inds = find(ismember(all_model_fixids,xv_fixs));
tr_inds = find(~ismember(all_model_fixids,xv_fixs));

%%
clear mod_*
n_iter = 6;
sac_hist = zeros(n_used_cells,length(tax));
stim_kern = cell(n_used_cells,1);
slope2_lambda = 250;
n_gabors = size(gabor_bank_output{2},2);
for c = 1:9
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
    %first fit sachist model
    Xmat = Tmat;
    W0 = zeros(size(Xmat,2),1);
    NT = length(tr_inds);
    
    %initialize parameters
    silent = 1;
    lamrange = [];
    lamrange2 = [slope2_lambda 1 length(tax)];
    nlxrange = [tax];
    Pcon = [];
    Pmono = [];
    hold_const = [];
    NLtype = 0;
    llist = [];
    
    [fitp,grad] = GLMsolve_jmm( Xmat(tr_inds,:), spkbs, W0, silent, lamrange, lamrange2, Pcon, Pmono, llist, hold_const, NLtype, 1e-6,nlxrange);
    sac_hist_k(c,:) = fitp.k(1:end-1);
    sac_hist_const(c) = fitp.k(end);
    sac_hist_LP(c) = fitp.LP;
    
    cur_stim_model_out = zeros(NT,1);
    cur_genfun = Tmat(tr_inds,:)*sac_hist_k(c,:)' + sac_hist_const(c);
    cur_predrate = log(1+exp(cur_genfun));
    sac_hist_LL(c) = sum(-Robs.*log(cur_predrate) + cur_predrate)/sum(Robs);
    cur_genfun_xv = Tmat(xv_inds,:)*sac_hist_k(c,:)' + sac_hist_const(c);
    cur_predrate_xv = log(1+exp(cur_genfun_xv));
    sac_hist_xvLL(c) = sum(-Robs_xv.*log(cur_predrate_xv) + cur_predrate_xv)/sum(Robs_xv);
    
    %%
%     cnt = 1;
%     cur_used_mods = [];
%     clear cur_mod_* ss_err best_loc cur_real_LL cur_xv_LL
%     while cnt <= n_iter
%         fprintf('Iteration %d of %d\n',cnt,n_iter);
%         cur_residual = (Robs./cur_predrate - 1).*(exp(cur_genfun)./(1+exp(cur_genfun)));
%         ss_err(cnt,:) = zeros(1,n_gabors);
%         for i = 1:n_gabors;
%             fprintf('Gabor %d of %d\n',i,n_gabors);
%             Xmat = [bsxfun(@times,Tmat(tr_inds,:),gabor_bank_output{c}(all_model_fixids(tr_inds),i))];
%             cur_beta = Xmat\cur_residual;
%             ss_err(cnt,i) = sum((Xmat*cur_beta-cur_residual).^2);
%             pause(0);
%         end
%         unused = setdiff(1:n_gabors,cur_used_mods);
%         [best_sserr(cnt),best_loc(cnt)] = min(ss_err(cnt,unused));
%         best_loc(cnt) = unused(best_loc(cnt));
%         cur_used_mods = [cur_used_mods best_loc(cnt)];
%         
%         Xmat = Tmat;
%         lamrange2 = [slope2_lambda 1 length(tax)];
%         nlxrange = [tax];
%         for i = 1:cnt
%             Xmat = [Xmat bsxfun(@times,Tmat,gabor_bank_output{c}(all_model_fixids,best_loc(i)))];
%             lamrange2 = [lamrange2; slope2_lambda length(tax)*i+1 length(tax)*(i+1)];
%             nlxrange = [nlxrange tax];
%         end
%         [fitp,grad] = GLMsolve_jmm( Xmat(tr_inds,:), spkbs, W0, 1, lamrange, lamrange2,Pcon, Pmono, llist, hold_const, NLtype,1e-5,nlxrange);
%         cur_mod_k{cnt} = fitp.k(1:end-1);
%         cur_mod_const(cnt) = fitp.k(end);
%         cur_mod_LP(cnt) = fitp.LP;
%         cur_mod_LL(cnt) = fitp.LL;
%         
%         cur_genfun = Xmat(tr_inds,:)*cur_mod_k{cnt} + cur_mod_const(cnt);
%         cur_predrate = log(1+exp(cur_genfun));
%         cur_real_LL(cnt) = sum(-Robs.*log(cur_predrate) + cur_predrate)/sum(Robs);
%         cur_genfun_xv = Xmat(xv_inds,:)*cur_mod_k{cnt} + cur_mod_const(cnt);
%         cur_predrate_xv = log(1+exp(cur_genfun_xv));
%         cur_xv_LL(cnt) = sum(-Robs_xv.*log(cur_predrate_xv) + cur_predrate_xv)/sum(Robs_xv);
%         
%         cnt = cnt + 1;
%     end

Xmat = Tmat;
lamrange2 = [slope2_lambda 1 length(tax)];
nlxrange = [tax];
for i = 1:n_gabors
    i
    Xmat = [Xmat bsxfun(@times,Tmat,gabor_bank_output{c}(all_model_fixids,i))];
    lamrange2 = [lamrange2; slope2_lambda length(tax)*i+1 length(tax)*(i+1)];
    nlxrange = [nlxrange tax];
    pause(0);
end
[fitp,grad] = GLMsolve_jmm( Xmat(tr_inds,:), spkbs, W0, 0, lamrange, lamrange2,Pcon, Pmono, llist, hold_const, NLtype,1e-5,nlxrange);

%%
stim_mod_const = fitp.k(end);
stim_mod_sac = fitp.k(1:length(tax));
for i = 1:n_gabors
   stim_mod_k(i,:) = fitp.k(length(tax)*i+1:length(tax)*(i+1)); 
end

for i = 1:length(lambdas)
   cur_set = find(gabor_bank_lambdas==lambdas(i)); 
   cur_en = cur_set(is_energy_filt(cur_set)==1);
   cur_phase = cur_set(is_energy_filt(cur_set)==0);
   lambda_en_k(i,:) = stim_mod_k(cur_en,:);
   lambda_phase1_k(i,:) = stim_mod_k(cur_phase(1),:) - stim_mod_k(cur_phase(3),:); 
   lambda_phase2_k(i,:) = stim_mod_k(cur_phase(2),:) - stim_mod_k(cur_phase(4),:); 
   lambda_phase_k(i,:) = sqrt(lambda_phase1_k(i,:).^2+lambda_phase2_k(i,:).^2); 
end

%%
    avg_rate = repmat(mean(Robs),length(Robs),1);
    null_LL(c) = sum(-Robs.*log(avg_rate)+avg_rate)/sum(Robs);
    avg_rate_xv = repmat(mean(Robs),length(Robs_xv),1);
    null_xvLL(c) = sum(-Robs_xv.*log(avg_rate_xv)+avg_rate_xv)/sum(Robs_xv);
    
    mod_sequence{c} = cur_mod_k;
    mod_LL(c,:) = cur_real_LL;
    mod_xvLL(c,:) = cur_xv_LL;
    mod_const(c,:) = cur_mod_const;
    mod_best_gabors(c,:) = best_loc;
    
end

%%
cd /Users/James/James_scripts/bruce/modelfits/
save gabor_tempmodfits_new mod_* gabor_params_fin tax gabor_bank SDIM is_energy_filt
%%
complexity_inds = log(abs(gabor_params_fin(:,7)) ./ sqrt(sum(gabor_params_fin(:,8:9).^2,2)));
close all

cur_cell = 7;
subplot(2,1,1)
orig_gabor = get_pgabor_mask(gabor_params_fin(cur_cell,1:6),0,[SDIM SDIM]);
imagesc(orig_gabor);set(gca,'ydir','normal');
title(sprintf('CI: %.3f',complexity_inds(cur_cell)));
subplot(2,1,2)
plot(tax,mod_sequence{cur_cell}{1}(1:length(tax)));
axis tight

figure
plot(mod_LL(cur_cell,:))
hold on
plot(mod_xvLL(cur_cell,:),'r')

% figure
% for i = 1:3
%     subplot(3,4,4*(i-1)+1)
%     cur_gabor = gabor_bank{cur_cell}(mod_best_gabors(cur_cell,i),:);
%     imagesc(reshape(cur_gabor,SDIM,SDIM));set(gca,'ydir','normal');
%     subplot(3,4,4*(i-1)+2)
%     plot(tax,mod_sequence{cur_cell}{i}(length(tax)*i+1:length(tax)*(i+1)));
%     hold on
%     plot(tax,mod_sequence{cur_cell}{end}(length(tax)*i+1:length(tax)*(i+1)),'r');
%     axis tight
%     if is_energy_filt(mod_best_gabors(cur_cell,i)) == 1
%         title('Energy Filt');
%     end
% end
% for i = 4:6
%     subplot(3,4,4*(i-4)+3)
%     cur_gabor = gabor_bank{cur_cell}(mod_best_gabors(cur_cell,i),:);
%     imagesc(reshape(cur_gabor,SDIM,SDIM));set(gca,'ydir','normal');
%     subplot(3,4,4*(i-4)+4)
%     plot(tax,mod_sequence{cur_cell}{i}(length(tax)*i+1:length(tax)*(i+1)));
%     hold on
%     plot(tax,mod_sequence{cur_cell}{end}(length(tax)*i+1:length(tax)*(i+1)),'r');
%     axis tight
%     if is_energy_filt(mod_best_gabors(cur_cell,i)) == 1
%         title('Energy Filt');
%     end
% end
% % 
figure
for i = 1:5
    subplot(5,4,4*(i-1)+1)
    cur_gabor = gabor_bank{cur_cell}(mod_best_gabors(cur_cell,i),:);
    imagesc(reshape(cur_gabor,SDIM,SDIM));set(gca,'ydir','normal');
    subplot(5,4,4*(i-1)+2)
    plot(tax,mod_sequence{cur_cell}{i}(length(tax)*i+1:length(tax)*(i+1)));
    hold on
    plot(tax,mod_sequence{cur_cell}{end}(length(tax)*i+1:length(tax)*(i+1)),'r');
    axis tight
%     if is_energy_filt(mod_best_gabors(cur_cell,i)) == 1
%         title('Energy Filt');
%     end
end
for i = 6:10
    subplot(5,4,4*(i-6)+3)
    cur_gabor = gabor_bank{cur_cell}(mod_best_gabors(cur_cell,i),:);
    imagesc(reshape(cur_gabor,SDIM,SDIM));set(gca,'ydir','normal');
    subplot(5,4,4*(i-6)+4)
    plot(tax,mod_sequence{cur_cell}{i}(length(tax)*i+1:length(tax)*(i+1)));
    hold on
    plot(tax,mod_sequence{cur_cell}{end}(length(tax)*i+1:length(tax)*(i+1)),'r');
    axis tight
%     if is_energy_filt(mod_best_gabors(cur_cell,i)) == 1
%         title('Energy Filt');
%     end
end


