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
scales = logspace(log10(2.5),log10(40),25);
% scales = [5:10 12 14 16 20 25 30 35 40 45 50 55 60];
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);
use_lfps = [1:2:96];

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
earlyavg_lfp_amps = nan(length(fix_start_inds),length(wfreqs),length(use_lfps));
lateavg_lfp_amps = nan(length(fix_start_inds),length(wfreqs),length(use_lfps));

full_ampgrams = [];
for ee = 1:n_expts
    fprintf('Expt %d of %d\n',ee,n_expts);
    clear ampgram all_V*
    for ll = 1:length(use_lfps)
        fprintf('Electrode %d\n',ll);
        filename = sprintf('Expt%d.p%dFullV.mat',Expt_nu(ee),use_lfps(ll));
        load(filename);
        V = double(FullV.V);
        V = decimate(V,dsf);
        V = V(:);
        ampgram(:,:,ll) = abs(cwt(V,scales,'cmor1-1'))';
    end
    t_ax = linspace(FullV.start,FullV.start+FullV.blklen/Fs,FullV.blklen);
    t_ax = downsample(t_ax,dsf);
    
    fix_set = find(fix_expt_num == ee);
    for n = 1:length(fix_set)
        cur_inds = find(t_ax >= fix_start_times(fix_set(n)) & t_ax < fix_end_times(fix_set(n)));
        fixavg_lfp_amps(fix_set(n),:,:) = mean(ampgram(cur_inds,:,:));
        cur_early_inds = find(t_ax >= fix_start_times(fix_set(n)) & t_ax < (fix_start_times(fix_set(n)) + 0.15));
        earlyavg_lfp_amps(fix_set(n),:,:) = mean(ampgram(cur_early_inds,:,:));
        cur_late_inds = find(t_ax >= (0.15 + fix_start_times(fix_set(n))) & t_ax < (fix_end_times(fix_set(n))));
        if ~isempty(cur_late_inds)
            lateavg_lfp_amps(fix_set(n),:,:) = mean(ampgram(cur_late_inds,:,:),1);
        end
    end
    fix_inds = find(full_expt_vec == ee);
    ufix_set = fix_inds(~isnan(full_fix_ids(fix_inds)));
    cur_interp_amps = interp1(t_ax,ampgram,full_t(ufix_set));
    full_ampgrams = [full_ampgrams; cur_interp_amps];
end

fixavg_lfp_amps = nanzscore(fixavg_lfp_amps);
earlyavg_lfp_amps = nanzscore(earlyavg_lfp_amps);
lateavg_lfp_amps = nanzscore(lateavg_lfp_amps);
full_ampgrams = nanzscore(full_ampgrams);
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
resh_X = reshape(resh_all_stims',[sdim sdim n_fixs]);
resh_X_sh = zeros(size(resh_X));
for ii = 1:n_fixs
    d2 = dist_shift2d(resh_X(:,:,ii), -x_cor{end}(ii), 2,0);
    d2 = dist_shift2d(d2,-y_cor{end}(ii),1,0);
    resh_X_sh(:,:,ii) = d2;
end
fullX_sh = reshape(resh_X_sh,sdim^2,n_fixs)';
fullX_sh_cropped = fullX_sh(:,new_crop);
% fullX = reshape(resh_X,sdim^2,n_fixs)';
% fullX_cropped = fullX(:,new_crop);
clear fullX_sh fullX
%%
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


%% FIT LFP MODEL PARAMS
cd ~/Data/bruce/7_15_12/G034/
load ./expt2_lfp_amp_models_v3 lfp_gabor_params LL rsq
lfp_gabor_params = lfp_gabor_params(use_lfps,:,:);
LL = LL(use_lfps,:);

%%
spatial_scale = 0.2;
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
        
        sparse_lambda = 500;
        smooth_lambda = 5000;
        init_beta = zeros(n_dim,1);
        [beta{cc,ww},best_LL(cc,ww)] = smooth_regress_2d(squeeze(fixavg_lfp_amps(:,ww,cc)),pred_mat,init_beta,smooth_lambda,sparse_lambda);
        fft_pred_out(cc,ww,:) = pred_mat*squeeze(beta{cc,ww});
    end
end

%%
infix_ids = find(~isnan(full_fix_ids));
n_fixs = length(unique(full_fix_ids(infix_ids)));
flen = 80;
flen_t = flen*dtd;
NT = length(full_expt_vec);

tent_centers = [0:dtd:0.1];
cur_sp = dtd;
while max(tent_centers) < flen_t
    tent_centers = [tent_centers (tent_centers(end)+cur_sp)];
    cur_sp = cur_sp + dtd/4;
end
tent_centers = round(tent_centers/dtd)*dtd;
if tent_centers(end) >= flen_t
    tent_centers(end) = [];
end
tbmat = construct_tent_bases(tent_centers,dtd);
[ntents,tblen] = size(tbmat);

tbmat = [zeros(ntents,tblen-1) tbmat];

sac_inds = zeros(size(full_t));
sac_inds(full_sac_inds(:,2)) = 1;

sac_Tmat = zeros(length(full_t),ntents);
for i = 1:ntents
    sac_Tmat(:,i) = conv(sac_inds,tbmat(i,:),'same');
end

%%
nat_set = 1:685;
white_set = 686:913;
sim_set = 914:1141;
obj_set = 1142:1369;
noise_set = [white_set sim_set];
uinds = find(~isnan(full_fix_ids));
nat_inds = find(ismember(full_image_vec(uinds),nat_set));
noise_inds = find(ismember(full_image_vec(uinds),noise_set));

for cc = 1:length(use_lfps)
    fprintf('Fitting Cell %d of %d\n',cc,length(use_lfps));
    for ww = 1:length(wfreqs)
        fprintf('Frequency %d of %d\n',ww,length(wfreqs));
        
        
%         gabor_emp1 = get_pgabor_mask_v2(XXc,YYc,lfp_gabor_params(cc,ww,1:6),0);
%         gabor_emp2 = get_pgabor_mask_v2(XXc,YYc,lfp_gabor_params(cc,ww,1:6),pi/2);
%         all_gabor_out1 = fullX_sh_cropped*gabor_emp1(:);
%         all_gabor_out2 = fullX_sh_cropped*gabor_emp2(:);
%         spatial_mod_out = sqrt(all_gabor_out1.^2 + all_gabor_out2.^2);
%         cur_spatial_mod_out = zscore(spatial_mod_out(full_fix_ids(uinds)));
        total_out = zscore(squeeze(fft_pred_out(cc,ww,(full_fix_ids(uinds))))');
        
%         stim_dep_X = bsxfun(@times,sac_Tmat(uinds,:),cur_spatial_mod_out);
        stim_dep_X = bsxfun(@times,sac_Tmat(uinds,:),total_out');
        ov_Xmat = [stim_dep_X sac_Tmat(uinds,:)];
        yobs = full_ampgrams(:,ww,cc);
        init_params = randn(2*ntents,1);
        [B_stimmod(cc,ww,:)] = smoothed_regression(yobs,ov_Xmat,init_params,[20e3 2e3],[1 ntents; ntents+1 2*ntents],0);
        
        %         cur_spatial_mod_out = zscore(spatial_mod_out(full_fix_ids(uinds(nat_inds))));
        %         stim_dep_X = bsxfun(@times,sac_Tmat(uinds(nat_inds),:),cur_spatial_mod_out);
        %         ov_Xmat = [stim_dep_X sac_Tmat(uinds(nat_inds),:)];
        %         yobs = full_ampgrams(nat_inds,ww,cc);
        %         init_params = randn(2*ntents,1);
        %         [B_stimmod_nat(cc,ww,:)] = smoothed_regression(yobs,ov_Xmat,init_params,[20e3 2e3],[1 ntents; ntents+1 2*ntents],0);
        %
        %         cur_spatial_mod_out = zscore(spatial_mod_out(full_fix_ids(uinds(noise_inds))));
        %         stim_dep_X = bsxfun(@times,sac_Tmat(uinds(noise_inds),:),cur_spatial_mod_out);
        %         ov_Xmat = [stim_dep_X sac_Tmat(uinds(noise_inds),:)];
        %         yobs = full_ampgrams(noise_inds,ww,cc);
        %         init_params = randn(2*ntents,1);
        %         [B_stimmod_noise(cc,ww,:)] = smoothed_regression(yobs,ov_Xmat,init_params,[20e3 2e3],[1 ntents; ntents+1 2*ntents],0);
    end
end

%%
cd ~/Data/bruce/7_15_12/G034/
save expt2_lfp_amp_models_fft B_stimmod wfreqs tent_centers lfp_gabor_params best_LL beta rsq* orientations

%%

for cc = 1:length(use_lfps)
        fprintf('Ch %d of %d\n',cc,length(use_lfps));
    for ww = 1:length(wfreqs);
        subplot(5,5,ww)
        cur_pdim = sqrt(length(beta{cc,ww}));
        imagesc(reshape(beta{cc,ww},cur_pdim,cur_pdim))
        cax = caxis();
        cm = max(abs(cax));
        caxis([-cm cm])
        title(sprintf('%.3f',wfreqs(ww)));
    end
    
    pause
    clf
end

%%
load ./expt2_lfp_amp_models_fft
B_stimmod_fft = B_stimmod;
wfreqs_fft = wfreqs;
use_lfps_fft = use_lfps;

load ./expt2_lfp_amp_models_v3

for i = 1:length(use_lfps)
    subplot(2,2,1)
pcolor(tent_centers,wfreqs,squeeze(B_stimmod(i,:,1:ntents)));shading flat;colorbar
set(gca,'yscale','log')
    subplot(2,2,3)
pcolor(tent_centers,wfreqs,squeeze(B_stimmod(i,:,ntents+1:end)));shading flat;colorbar
set(gca,'yscale','log')
    subplot(2,2,2)
pcolor(tent_centers,wfreqs_fft,squeeze(B_stimmod_fft(i,:,1:ntents)));shading flat;colorbar
set(gca,'yscale','log')
    subplot(2,2,4)
pcolor(tent_centers,wfreqs_fft,squeeze(B_stimmod_fft(i,:,ntents+1:end)));shading flat;colorbar
set(gca,'yscale','log')
pause
clf
end
%%
cur_freq = 20;
figure
cnt = 1;
for i = 1:10
    for j = 1:10
        cur = find(X_pos == i & Y_pos == j);
        if ~isempty(cur)
            subplot(10,10,cnt)
            hold on
            plot(tent_centers,squeeze(B_stimmod(cur,cur_freq,1:ntents)),'b');
%             hold on
%             plot(tent_centers,squeeze(B_stimmod(cur,cur_freq,ntents+1:end)),'r');
                     
            axis tight;
            xl = xlim();
            line(xl,[0 0],'color','k')
        end
        cnt = cnt + 1;
    end
end

