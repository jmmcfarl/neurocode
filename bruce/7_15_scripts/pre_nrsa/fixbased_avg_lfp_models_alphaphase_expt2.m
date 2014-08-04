%% Load Data
clear all;
addpath(genpath('~/Code/James_scripts'));
% cd ~/Data/bruce/7_15_12/G029/
cd ~/Data/bruce/7_15_12/G034/

load ./CellList.mat
single_units = find(CellList(1,:,1) > 0);

load ./Expt2_compiled_windata_d1p5.mat
resh_all_stims = resh_all_stims/std(resh_all_stims(:));

fix_win_dur = 0.15;

Pix2Deg = 0.018837;
NT = length(full_fix_ids);
sdim = length(xpatch_inds);
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

Fs = 3e4;
dsf = 30;Fsd = Fs/dsf;
[b_alpha,a_alpha] = butter(2,[6 12]/(Fsd/2));
scales = logspace(log10(3.2),log10(60),20);
scales = [scales 70 80 100 120 150];
% scales = [5:10 12 14 16 20 25 30 35 40 45 50 55 60];
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);
use_lfps = [1:96];

alpha_smooth_win = 50;

dir_fit = 0;
load ./expt2_lfp_amp_models_full_alphaphase lfp_gabor_params fft_beta
% lfp_gabor_params = lfp_gabor_params(1:2:end,:,:);
% fft_beta = fft_beta(1:2:end,:);
% load ./expt2_lfp_amp_models_full_v2

%% eliminate image flip fixations
%used_fixs = 1:length(fix_is_sac);
used_fixs = find(fix_is_sac==1);
flip_fixs = find(fix_is_sac==0);
full_fix_ids(ismember(full_fix_ids,flip_fixs)) = nan;
full_fix_starts = full_fix_starts(used_fixs);
full_fix_wends = full_fix_wends(used_fixs);
full_fix_ends = full_fix_ends(used_fixs);
fix_start_inds = fix_start_inds(used_fixs);
fix_start_times = fix_start_times(used_fixs);
fix_end_inds = fix_end_inds(used_fixs);
fix_end_times = fix_end_times(used_fixs);
%%
n_fixs = length(full_fix_wends);
fix_expt_num = nan(n_fixs,1);
for i = 1:n_fixs
    cur_inds = full_fix_starts(i):full_fix_wends(i);
    fix_expt_num(i) = unique(full_expt_vec(cur_inds));
end
%%
Expt_nu = [13 14 15 16 25 28 29];
n_expts = length(Expt_nu);
fixavg_lfp_amps = nan(length(fix_start_inds),length(wfreqs),length(use_lfps));

fix_start_log = zeros(size(full_t));
fix_end_log = fix_start_log;
fix_start_log(full_fix_starts) = 1;
fix_end_log(full_fix_ends) = 1;

full_ampgrams = [];
full_alpha_phase = [];
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
    %     alpha_phase = nan(vlen,length(use_lfps));
    ampgram = nan(vlen,length(wfreqs),length(use_lfps));
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d\n',ll);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = decimate(V,dsf);
        V = V(:);
        V_alpha = filtfilt(b_alpha,a_alpha,V);
        V_alpha_all(:,ll) = V_alpha;
                V_alpha_phase = angle(hilbert(V_alpha));
                V_alpha_phase = unwrap_phase_monotonic(V_alpha_phase);
                alpha_phase(:,ll) = V_alpha_phase;

                ampgram(:,:,ll) = abs(cwt(V,scales,'cmor1-1'))';
    end
%     avg_alpha = mean(V_alpha_all,2);
%     alpha_phase = angle(hilbert(avg_alpha));
%     alpha_phase = unwrap_phase_monotonic(alpha_phase);
%     alpha_phase = smooth(alpha_phase,alpha_smooth_win,'lowess');
    
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
    
    if dir_fit == 1
        fix_set = find(fix_expt_num == ee);
        for n = 1:length(fix_set)
            cur_inds = find(t_ax >= fix_start_times(fix_set(n)) & t_ax < fix_end_times(fix_set(n)));
            fixavg_lfp_amps(fix_set(n),:,:) = mean(ampgram(cur_inds,:,:));
        end
    end
    fix_inds = find(full_expt_vec == ee);
    ufix_set = fix_inds(~isnan(full_fix_ids(fix_inds)));
    cur_interp_amps = interp1(t_ax,ampgram,full_t(ufix_set));
    cur_interp_phase = interp1(t_ax,alpha_phase,full_t(ufix_set));
    full_ampgrams = [full_ampgrams; cur_interp_amps];
    full_alpha_phase = [full_alpha_phase; cur_interp_phase];
end
if dir_fit==1
    fixavg_lfp_amps = nanzscore(fixavg_lfp_amps);
end
full_ampgrams = nanzscore(full_ampgrams);

%%
use_inds = find(~isnan(full_fix_ids));
rel_fix_start_inds = find(fix_start_log(use_inds) == 1);
rel_fix_end_inds = find(fix_end_log(use_inds) == 1);
phase_since_trial = nan(size(full_alpha_phase));
for i = 1:n_fixs
    cur_set = rel_fix_start_inds(i):rel_fix_end_inds(i);
       phase_since_trial(cur_set,:) = bsxfun(@minus,full_alpha_phase(cur_set,:),full_alpha_phase(cur_set(1),:));
%     phase_since_trial(cur_set) = full_alpha_phase(cur_set)-full_alpha_phase(cur_set(1));
end

%% RECONSTRUCT MAP STIMULUS
new_RF_patch = [-0.3 1.2; -1.2 0.3]; %location of RFs in degrees [x1 x2;y1 y2]
[curXX,curYY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));
curXX = curXX(:);
curYY = curYY(:);

xpatch_inds_cropped = find(xax >= new_RF_patch(1,1) & xax <= new_RF_patch(1,2));
ypatch_inds_cropped = find(yax >= new_RF_patch(2,1) & yax <= new_RF_patch(2,2));

new_crop = find(curXX >= new_RF_patch(1,1) & curXX <= new_RF_patch(1,2) & ...
    curYY >= new_RF_patch(2,1) & curYY <= new_RF_patch(2,2));
% fullX_cropped = fullX(:,new_crop);

sdim = length(xpatch_inds);
[XXc,YYc] = meshgrid(xax(xpatch_inds_cropped),yax(ypatch_inds_cropped));
[XX,YY] = meshgrid(xax(xpatch_inds),yax(ypatch_inds));

load ./fixation_based_corrections_34
resh_X = reshape(resh_all_stims',[sdim sdim length(fix_is_sac)]);
resh_X_sh = zeros(size(resh_X));
for ii = 1:n_fixs
    d2 = dist_shift2d(resh_X(:,:,ii), -x_cor{end}(ii), 2,0);
    d2 = dist_shift2d(d2,-y_cor{end}(ii),1,0);
    resh_X_sh(:,:,ii) = d2;
end
fullX_sh = reshape(resh_X_sh,sdim^2,length(fix_is_sac))';
fullX_sh_cropped = fullX_sh(:,new_crop);
clear fullX_sh

if dir_fit==1
    %% INTERPOLATE RF CENTERS OF LFP TUNING BASED ON UNIT GABOR FITS
    % load ./gabor_initfits_d1p25_varrate
    load ./gabor_tracking_varmeans
    gabor_params_used = gabor_params_f{2};
    load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
    X_pos = ArrayConfig.X;
    Y_pos = ArrayConfig.Y;
    
    %fit smoothed retinotopic surface
    orientations = linspace(0,pi-pi/12,12);
    for i = 1:96
        tempx(Y_pos(i),X_pos(i)) = gabor_params_used(i,1);
        tempy(Y_pos(i),X_pos(i)) = gabor_params_used(i,2);
        tempinds(Y_pos(i),X_pos(i)) = i;
    end
    weights = ones(10,10);
    weights(tempx==0) = 0;
    xpos_interp = smoothn(tempx,weights,'robust');
    ypos_interp = smoothn(tempy,weights,'robust');
    used_inds = find(weights == 1);
    tempinds = tempinds(used_inds);
    interp_x(tempinds) = xpos_interp(used_inds);
    interp_y(tempinds) = ypos_interp(used_inds);
    gabor_params_used(:,1) = interp_x;
    gabor_params_used(:,2) = interp_y;
    
    %% estimate preferred orientations
    rsq = zeros(length(use_lfps),length(wfreqs),12);
    en_out = zeros(12,n_fixs);
    for cc = 1:length(use_lfps)
        fprintf('Ch %d of %d\n',cc,length(use_lfps));
        init_params(1:5) = gabor_params_used(use_lfps(cc),1:5);
        init_params(6) = 1;
        for i = 1:12
            init_params(3) = orientations(i);
            cur_mask1 = get_pgabor_mask_v2(XXc,YYc,init_params,0);
            cur_mask2 = get_pgabor_mask_v2(XXc,YYc,init_params,pi/2);
            
            mask1_out = fullX_sh_cropped(used_fixs,:)*cur_mask1(:);
            mask2_out = fullX_sh_cropped(used_fixs,:)*cur_mask2(:);
            
            en_out(i,:) = sqrt(mask1_out.^2+mask2_out.^2);
        end
        en_out = bsxfun(@minus,en_out,mean(en_out,2));
        en_out = bsxfun(@rdivide,en_out,std(en_out,[],2));
        
        for ww = 1:length(wfreqs)
            %     for ww = 9
            for i = 1:12
                temp = [en_out(i,:)' ones(n_fixs,1)];
                [B,BINT,R,RINT,STATS] = regress(squeeze(fixavg_lfp_amps(:,ww,cc)),temp);
                rsq(cc,ww,i) = STATS(1);
            end
        end
    end
    
    %% FIT LFP MODEL PARAMS
    [max_rsq,maxloc] = max(rsq,[],3);
    
    hold_const = [1 1 0 0 0 1 0 0];
    LB = [-0.1 -0.8 0 0.125 0.025 0.2 0 -Inf];
    UB = [0.8 0.1 pi 0.5 0.4 6 Inf Inf];
    lfp_gabor_params = zeros(length(use_lfps),length(wfreqs),9);
    LL = zeros(length(use_lfps),length(wfreqs));
    clear init_params
    for cc = 1:length(use_lfps)
        fprintf('CC %d of %d\n',cc,length(use_lfps));
        for ww = 1:length(wfreqs)
            fprintf('freq %d of %d\n',ww,length(wfreqs));
            init_params(1:5) = gabor_params_used(use_lfps(cc),1:5);
            init_params(6) = 1;
            init_params(3) = orientations(maxloc(cc,ww));
            cur_init_params = [init_params zeros(1,3)];
            [lfp_gabor_params(cc,ww,:),LL(cc,ww)] = fit_gabor_params_lfp_v2(XXc,YYc,cur_init_params,fullX_sh_cropped(used_fixs,:)...
                ,squeeze(fixavg_lfp_amps(:,ww,cc)),hold_const,LB,UB);            
        end
    end
end
%% FIT FFT MODEL PARAMS
spatial_scale = 0.2;
sparse_lambda = 500;
smooth_lambda = 10000;
for cc = 1:length(use_lfps)
    fprintf('Ch %d of %d\n',cc,length(use_lfps));
    for ww = 1:length(wfreqs)
        cur_params = squeeze(lfp_gabor_params(cc,ww,1:6));
        cur_params(5) = spatial_scale;
        uset = (abs(XXc-cur_params(1)) <= 2*cur_params(5) & ...
            abs(YYc-cur_params(2)) <= 2*cur_params(5));
        xdims = max(uset);
        ydims = max(uset,[],2);
        if sum(xdims) > sum(ydims)
            last_col = find(xdims == 1,1,'first');
            uset(:,last_col) = 0;
        end
        if sum(ydims) > sum(xdims)
            last_col = find(ydims == 1,1,'first');
            uset(last_col,:) = 0;
        end
        uset = find(uset==1);
        
        gauss_mask = get_pgauss_mask_v2(XXc,YYc,cur_params);
        gauss_mask = gauss_mask(uset);
        X_out = bsxfun(@times,fullX_sh_cropped(:,uset),gauss_mask');
        n_dim = length(uset);
        per_dim = sqrt(n_dim);
        X_out = reshape(X_out,[size(X_out,1) per_dim per_dim]);
        fft_out = nan(size(X_out));
        for i = 1:size(X_out,1)
            fft_out(i,:,:) = abs(fftshift(fft2(squeeze(X_out(i,:,:)))));
        end
        pred_mat = reshape(fft_out,size(X_out,1),per_dim^2);
        pred_mat = zscore(pred_mat);
        
        if dir_fit==1
            init_beta = zeros(n_dim,1);
            [fft_beta{cc,ww},best_LL(cc,ww)] = smooth_regress_2d(squeeze(fixavg_lfp_amps(:,ww,cc)),pred_mat(used_fixs,:),init_beta,smooth_lambda,sparse_lambda);
        end
        fft_pred_out(cc,ww,:) = pred_mat*squeeze(fft_beta{cc,ww});
    end
end

%%
infix_ids = find(~isnan(full_fix_ids));
n_fixs = length(unique(full_fix_ids(infix_ids)));
NT = length(full_expt_vec);

max_phase = 2*pi*5;
nbins = 62;
% max_phase = 2*pi*4;
% nbins = 50;
uset = find(phase_since_trial(:,1) <= max_phase);
tax = prctile(phase_since_trial(uset,1),linspace(100/nbins,100-100/nbins,nbins));
% tax = linspace(0,max_phase,nbins);
% Tmat = tbrep(phase_since_trial,tax);
% Tmat = Tmat(uset,:);
%%
% nat_set = 1:685;
% white_set = 686:913;
% sim_set = 914:1141;
% obj_set = 1142:1369;
% noise_set = [white_set sim_set];
% uinds = find(~isnan(full_fix_ids));
% uinds = uinds(uset);

% nat_inds = find(ismember(full_image_vec(uinds),nat_set));
% noise_inds = find(ismember(full_image_vec(uinds),noise_set));
% obj_inds = find(ismember(full_image_vec(uinds),obj_set));

% clear B_stimmod*
for cc =1:length(use_lfps)
    fprintf('Fitting Cell %d of %d\n',cc,length(use_lfps));
    
    uset = find(phase_since_trial(:,cc) <= max_phase);
Tmat = tbrep(phase_since_trial(:,cc),tax);
Tmat = Tmat(uset,:);
uinds = infix_ids(uset);

    for ww = 1:length(wfreqs)
        fprintf('Frequency %d of %d\n',ww,length(wfreqs));
        
%          gabor_emp1 = get_pgabor_mask_v2(XXc,YYc,lfp_gabor_params(cc,ww,1:6),0);
%         gabor_emp2 = get_pgabor_mask_v2(XXc,YYc,lfp_gabor_params(cc,ww,1:6),pi/2);
%         all_gabor_out1 = fullX_sh_cropped*gabor_emp1(:);
%         all_gabor_out2 = fullX_sh_cropped*gabor_emp2(:);
%         spatial_mod_out = sqrt(all_gabor_out1.^2 + all_gabor_out2.^2);
%         cur_spatial_mod_out = zscore(spatial_mod_out(full_fix_ids(uinds)));
%  
%         stim_dep_X = bsxfun(@times,Tmat,cur_spatial_mod_out);
%         ov_Xmat = [stim_dep_X Tmat];
%         yobs = full_ampgrams(uset,ww,cc);
%         init_params = randn(2*nbins,1);
%         [B_stimmod_gabor(cc,ww,:)] = smoothed_regression(yobs,ov_Xmat,init_params,[10e3 2e3],[1 nbins; nbins+1 2*nbins],0);

        total_out = zscore(squeeze(fft_pred_out(cc,ww,(full_fix_ids(uinds))))');
        
%         stim_dep_X = bsxfun(@times,Tmat,total_out');
        stim_dep_X = bsxfun(@times,Tmat,total_out');
        ov_Xmat = [stim_dep_X Tmat];
        yobs = full_ampgrams(uset,ww,cc);
        init_params = randn(2*nbins,1);
%         [B_stimmod_fft(cc,ww,:)] = smoothed_regression(yobs,ov_Xmat,init_params,[10e3 2e3],[1 nbins; nbins+1 2*nbins],0);
        [B_stimmod_fft(cc,ww,:)] = smoothed_regression(yobs,ov_Xmat,init_params,[50e3 2e3],[1 nbins; nbins+1 2*nbins],0);

%             total_out = zscore(squeeze(fft_pred_out(cc,ww,(full_fix_ids(uinds(nat_inds)))))');
%         stim_dep_X = bsxfun(@times,Tmat(nat_inds,:),total_out');
%         ov_Xmat = [stim_dep_X Tmat(nat_inds,:)];
%         yobs = full_ampgrams(uset(nat_inds),ww,cc);
%         init_params = randn(2*nbins,1);
%         [B_stimmod_fft_nat(cc,ww,:)] = smoothed_regression(yobs,ov_Xmat,init_params,[10e3 2e3],[1 nbins; nbins+1 2*nbins],0);
%          
%         total_out = zscore(squeeze(fft_pred_out(cc,ww,(full_fix_ids(uinds(noise_inds)))))');
%         stim_dep_X = bsxfun(@times,Tmat(noise_inds,:),total_out');
%         ov_Xmat = [stim_dep_X Tmat(noise_inds,:)];
%         yobs = full_ampgrams(uset(noise_inds),ww,cc);
%         init_params = randn(2*nbins,1);
%         [B_stimmod_fft_noise(cc,ww,:)] = smoothed_regression(yobs,ov_Xmat,init_params,[10e3 2e3],[1 nbins; nbins+1 2*nbins],0);
% 
%          total_out = zscore(squeeze(fft_pred_out(cc,ww,(full_fix_ids(uinds(obj_inds)))))');
%         stim_dep_X = bsxfun(@times,Tmat(obj_inds,:),total_out');
%         ov_Xmat = [stim_dep_X Tmat(obj_inds,:)];
%         yobs = full_ampgrams(uset(obj_inds),ww,cc);
%         init_params = randn(2*nbins,1);
%         [B_stimmod_fft_obj(cc,ww,:)] = smoothed_regression(yobs,ov_Xmat,init_params,[10e3 2e3],[1 nbins; nbins+1 2*nbins],0);

   end
end

%%
cd ~/Data/bruce/7_15_12/G034/
save expt2_lfp_amp_models_full_alphaphase_indalpha B_stimmod* wfreqs tax use_lfps scales fft* lfp_gabor_params LL 

%%
% cur_freq = 20;
% figure
% cnt = 1;
% for i = 1:10
%     for j = 1:10
%         cur = find(X_pos == i & Y_pos == j);
%         if ~isempty(cur)
%             subplot(10,10,cnt)
%             hold on
%             plot(tent_centers,squeeze(B_stimmod(cur,cur_freq,1:ntents)),'b');
% %             hold on
% %             plot(tent_centers,squeeze(B_stimmod(cur,cur_freq,ntents+1:end)),'r');
%
%             axis tight;
%             xl = xlim();
%             line(xl,[0 0],'color','k')
%         end
%         cnt = cnt + 1;
%     end
% end
%
%%
load ./expt2_lfp_amp_models_full_alphaphase_v2
close all
for i = 1:length(use_lfps)
    fprintf('Ch %d of %d\n',i,length(use_lfps));
    subplot(2,1,1)
    pcolor(tax/(2*pi),wfreqs,squeeze(B_stimmod_fft(i,:,1:nbins)));shading flat
    set(gca,'yscale','log')
    caxis([-0.02 0.15])
    colorbar
    subplot(2,1,2)
    pcolor(tax/(2*pi),wfreqs,squeeze(B_stimmod_fft(i,:,nbins+1:end)));shading flat
    set(gca,'yscale','log')
    caxis([-0.3 0.3])
    colorbar
    pause
    clf
end

%%
close all
for i = 1:length(use_lfps)
    fprintf('Ch %d of %d\n',i,length(use_lfps));
    subplot(2,2,1)
    pcolor(tax/(2*pi),wfreqs,squeeze(B_stimmod_fft(i,:,1:nbins)));shading flat
    set(gca,'yscale','log')
    caxis([-0.02 0.15])
    colorbar
    subplot(2,2,2)
    pcolor(tax/(2*pi),wfreqs,squeeze(B_stimmod_gabor(i,:,1:nbins)));shading flat
    set(gca,'yscale','log')
    caxis([-0.02 0.15])
    colorbar
    subplot(2,2,3)
    pcolor(tax/(2*pi),wfreqs,squeeze(B_stimmod_fft(i,:,nbins+1:end)));shading flat
    set(gca,'yscale','log')
    caxis([-0.3 0.3])
    colorbar
    subplot(2,2,4)
    pcolor(tax/(2*pi),wfreqs,squeeze(B_stimmod_gabor(i,:,nbins+1:end)));shading flat
    set(gca,'yscale','log')
    caxis([-0.3 0.3])
    colorbar
    pause
    clf
end

%%
cd ~/Data/bruce/7_15_12/G034/
close all
% for cc = 1:length(use_lfps)
for cc = 27
    fprintf('Ch %d of %d\n',cc,length(use_lfps));
    for ww = 1:length(wfreqs);
        figure(1)
        subplot(5,5,ww)
        cur_pdim = sqrt(length(fft_beta{cc,ww}));
        imagesc(reshape(fft_beta{cc,ww},cur_pdim,cur_pdim))
        cax = caxis();
        cm = max(abs(cax));
        caxis([-cm cm])
        title(sprintf('F: %.1f Hz',wfreqs(ww)),'fontsize',20);
            set(gca,'xticklabel',[],'ytickLabel',[])
    end
%     figure(2)
%     plot(wfreqs,best_LL(cc,:))
%     hold on
%     plot(wfreqs,LL(cc,:),'r')
%     pause
% %     figure(2)
%     clf
end
cd ~/Desktop/
fillPage(gcf,'PaperSize',[20 20])
print('-dpdf','gamma_fft_models_ch27')

%%
cd ~/Data/bruce/7_15_12/G034/
cur_freq = 13;
    load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
    X_pos = ArrayConfig.X;
    Y_pos = ArrayConfig.Y;
    
cnt = 1;
for i = 1:10
    for j = 1:10
        cur = find(X_pos == i & Y_pos == j);
        if ~isempty(cur)
            subplot(10,10,cnt)
            cur_pdim = sqrt(length(fft_beta{cur,cur_freq}));
            imagesc(reshape(fft_beta{cur,cur_freq},cur_pdim,cur_pdim))
            cax = caxis();
            cm = max(abs(cax));
            caxis([-cm cm])
            set(gca,'xticklabel',[],'ytickLabel',[])
        end
        cnt = cnt + 1;
    end
end

cd ~/Desktop/
fillPage(gcf,'PaperSize',[20 20])
print('-dpdf','gamma_fft_models_array')
%%
    load ./gabor_tracking_varmeans
    gabor_params_used = gabor_params_f{2};
    load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
    X_pos = ArrayConfig.X;
    Y_pos = ArrayConfig.Y;
    
    %fit smoothed retinotopic surface
    orientations = linspace(0,pi-pi/12,12);
    for i = 1:96
        tempx(Y_pos(i),X_pos(i)) = gabor_params_used(i,1);
        tempy(Y_pos(i),X_pos(i)) = gabor_params_used(i,2);
        tempinds(Y_pos(i),X_pos(i)) = i;
    end
    weights = ones(10,10);
    weights(tempx==0) = 0;
    xpos_interp = smoothn(tempx,weights,'robust');
    ypos_interp = smoothn(tempy,weights,'robust');
    used_inds = find(weights == 1);
    tempinds = tempinds(used_inds);
    interp_x(tempinds) = xpos_interp(used_inds);
    interp_y(tempinds) = ypos_interp(used_inds);

%%
figure
    load ~/Data/bruce/7_15_12/G029/ArrayConfig.mat
    X_pos = ArrayConfig.X;
    Y_pos = ArrayConfig.Y;
cnt = 1;
for i = 1:10
    for j = 1:10
        cur = find(X_pos == i & Y_pos == j);
        if ~isempty(cur)
            subplot(10,10,cnt)
            pcolor(tax/(2*pi),wfreqs,squeeze(B_stimmod_fft(cur,:,1:nbins)));shading flat
            set(gca,'yscale','log')
            caxis([-0.02 0.15])
        end
        cnt = cnt + 1;
    end
end

%%
load ./expt2_lfp_amp_models_full_alphaphase_v2
B_stimmod_fft_2 = B_stimmod_fft;
load ./expt2_lfp_amp_models_full_v2.mat
ntents = length(tent_centers);

close all
% for i = 25:length(use_lfps)
for i = [1 11 13 14 33 34 48]
    
% i = 14
i
    2*(i-1)+1;
%     fprintf('Ch %d of %d\n',i,length(use_lfps));
%     subplot(2,2,1)
figure
pcolor(tax/(2*pi),wfreqs,squeeze(B_stimmod_fft_2(2*(i-1)+1,:,1:nbins)));shading interp
    set(gca,'yscale','log')
    caxis([-0.0 0.16])
    xlabel('Alpha Cycles','fontsize',24)
    ylabel('Frequency (Hz)','fontsize',24)
    title('Stimulus-Dependent','fontsize',24)
    colorbar
%     subplot(2,2,3)
figure
    pcolor(tax/(2*pi),wfreqs,squeeze(B_stimmod_fft_2(2*(i-1)+1,:,nbins+1:end)));shading interp
    set(gca,'yscale','log')
 caxis([-0.3 0.75])
    xlabel('Alpha Cycles','fontsize',24)
    ylabel('Frequency (Hz)','fontsize',24)
    title('Stimulus-Independent','fontsize',24)
     colorbar
%     subplot(2,2,2)
figure
    pcolor(tent_centers,wfreqs,squeeze(B_stimmod_fft(i,:,1:ntents)));shading interp
    set(gca,'yscale','log')
%     caxis([-0.0 0.16])
    caxis([-0.02 0.075])
    xlabel('Time (s)','fontsize',24)
    ylabel('Frequency (Hz)','fontsize',24)
    title('Stimulus-Dependent','fontsize',24)
    colorbar
%     subplot(2,2,4)
figure
    pcolor(tent_centers,wfreqs,squeeze(B_stimmod_fft(i,:,ntents+1:end)));shading interp
    set(gca,'yscale','log')
     xlabel('Time (s)','fontsize',24)
    ylabel('Frequency (Hz)','fontsize',24)
    title('Stimulus-Independent','fontsize',24)
    colorbar
caxis([-0.3 0.75])
 pause
%     clf
close all
end

%%
for i = 1:48
    i
    2*(i-1)+1;

stim_vals = [-2 0 2];
wrange = find(wfreqs >= 30 & wfreqs <= 60);
avg_stim_dep_p = squeeze(mean(B_stimmod_fft_2(2*(i-1)+1,wrange,1:nbins),2));
avg_stim_ind_p = squeeze(mean(B_stimmod_fft_2(2*(i-1)+1,wrange,nbins+1:end),2));
avg_stim_dep = squeeze(mean(B_stimmod_fft(i,wrange,1:ntents),2));
avg_stim_ind = squeeze(mean(B_stimmod_fft(i,wrange,ntents+1:end),2));

pred_rates_p = bsxfun(@times,avg_stim_dep_p,stim_vals);
pred_rates_p = bsxfun(@plus,pred_rates_p,avg_stim_ind_p);
pred_rates = bsxfun(@times,avg_stim_dep,stim_vals);
pred_rates = bsxfun(@plus,pred_rates,avg_stim_ind);

subplot(2,2,1)
plot(tax/(2*pi),squeeze(mean(B_stimmod_fft_2(2*(i-1)+1,wrange,1:nbins),2)))
xlabel('Alpha cycles')
subplot(2,2,2)
plot(tent_centers,squeeze(mean(B_stimmod_fft(i,wrange,1:ntents),2)),'b')
xlabel('Time (s)')
xlim([0 0.5])
subplot(2,2,3)
plot(tax/(2*pi),pred_rates_p)
xlabel('Alpha cycles')
subplot(2,2,4)
plot(tent_centers,pred_rates)
xlim([0 0.5])
xlabel('Time (s)')

pause
clf
end
%%
close all
for i = 1:length(use_lfps)
    fprintf('Ch %d of %d\n',i,length(use_lfps));
    subplot(2,2,1)
    pcolor(tax/(2*pi),wfreqs,squeeze(B_stimmod_fft(i,:,1:nbins)));shading flat
    set(gca,'yscale','log')
    caxis([-0.02 0.15])
%     colorbar
    subplot(2,2,2)
    pcolor(tax/(2*pi),wfreqs,squeeze(B_stimmod_fft_nat(i,:,1:nbins)));shading flat
    set(gca,'yscale','log')
    caxis([-0.02 0.15])
    subplot(2,2,3)
    pcolor(tax/(2*pi),wfreqs,squeeze(B_stimmod_fft_noise(i,:,1:nbins)));shading flat
    set(gca,'yscale','log')
    caxis([-0.02 0.15])
%     colorbar
    subplot(2,2,4)
    pcolor(tax/(2*pi),wfreqs,squeeze(B_stimmod_fft_obj(i,:,1:nbins)));shading flat
    set(gca,'yscale','log')
    caxis([-0.02 0.15])
    pause
    clf
end
