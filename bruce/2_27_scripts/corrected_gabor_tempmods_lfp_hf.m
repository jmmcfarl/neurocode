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

cd /Users/James/James_scripts/bruce/modelfits
load pref_oris
cd /Users/James/James_scripts/bruce/modelfits/
% load fixation_gabor_models
load gabor_tempmodfits_en_lin

%%
Fs = 1000;
dsf = 4;
Fsd = Fs/dsf;
niqf = Fs/2;
% [b,a] = butter(4,[45 60]/niqf);
[b,a] = butter(2,100/niqf,'high');
scales = logspace(log10(2.5),log10(50),30);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

%%
all_lfp_pow = [];
within_fix_avgs = [];
fix_nums = [];
all_used_inds = [];
all_amps = [];
for blockid = 1:4;
    fprintf('Block %d of %d\n',blockid,3);
    
    %%
    cd /Users/James/Data/bruce/2_27_12
%     load(sprintf('lemM232A.5%d.lfp.mat',blockid));
%     %get start times of each LFP trial
%     n_lfp_trials = length(LFP.Trials);
%     lfp_trial_start = nan(n_lfp_trials,1);
%     lfp_trial_stop = nan(n_lfp_trials,1);
%     for i = 1:n_lfp_trials
%         lfp_trial_start(i) = LFP.Trials(i).Start/1e4; %originally in tenths of ms
%         lfp_trial_stop(i) = LFP.Trials(i).End/1e4;
%         lfp_dur(i) = size(LFP.Trials(i).LFP,1)/1000;
%     end
%     
%     lfp_time = [];
%     lfp_samps = [];
%     for i = 1:n_lfp_trials
%         lfp_time = [lfp_time linspace(lfp_trial_start(i),lfp_trial_start(i)+lfp_dur(i),size(LFP.Trials(i).LFP,1))];
%         lfp_samps = [lfp_samps; LFP.Trials(i).LFP];
%     end
%     
%     lfp_samps = filtfilt(b,a,lfp_samps);
%     lfp_samps = abs(hilbert(lfp_samps));
%     lfp_sampsd = downsample(lfp_samps,dsf);
%     
%     lfp_timed = downsample(lfp_time,dsf);

    cd ~/Data/bruce/2_27_12/M232/
    load(sprintf('Expt5%dhfPow.mat',blockid));
    hf_pow = zscore(hf_pow);
    
    cur_all_model = find(all_model_blockids==blockid);
%     interp_amps = interp1(lfp_timed,lfp_sampsd,all_model_time_axis(cur_all_model));
    interp_amps = interp1(t,hf_pow,all_model_time_axis(cur_all_model));
    
    %%
    used_tpoints = find(~isnan(interp_amps(:,1)));
    all_amps = cat(1,all_amps,interp_amps(used_tpoints,:));
    all_used_inds = [all_used_inds; cur_all_model(used_tpoints)'];
    
end
%%
cur_used_fixs = unique(all_model_fixids(all_used_inds));
n_fixs = length(cur_used_fixs);
fix_avg_amps = zeros(24,n_fixs);
for i = 1:n_fixs
    cur_set = find(all_model_fixids(all_used_inds) == cur_used_fixs(i));
    fix_avg_amps(:,i) = squeeze(mean(all_amps(cur_set,:)));
end
fix_avg_amps = bsxfun(@minus,fix_avg_amps,mean(fix_avg_amps,2));
fix_avg_amps = bsxfun(@rdivide,fix_avg_amps,std(fix_avg_amps,[],2));

%%
X_resh = reshape(X,NT,SDIM^2);
orientations = linspace(0,pi-pi/12,12);
init_params = [5 -5 0 15 6 1];
for i = 1:12
    init_params(3) = orientations(i);
    cur_mask1 = get_pgabor_mask(init_params,0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(init_params,pi/2,[SDIM SDIM]);
    
    mask1_out = X_resh(cur_used_fixs,:)*cur_mask1(:);
    mask2_out = X_resh(cur_used_fixs,:)*cur_mask2(:);
    
    en_out(i,:) = sqrt(mask1_out.^2+mask2_out.^2);
end
en_out = bsxfun(@minus,en_out,mean(en_out,2));
en_out = bsxfun(@rdivide,en_out,std(en_out,[],2));

rsq = zeros(24,12);
% for cc = 1:24
cc = 17
for i = 1:12
    temp = [en_out(i,:)' ones(length(cur_used_fixs),1)];
    [B,BINT,R,RINT,STATS] = regress(squeeze(fix_avg_amps(cc,:))',temp);
    rsq(cc,i) = STATS(1);
end
% end

[max_rsq,maxloc] = max(rsq,[],2);

% used_freq = 8;
hold_const = [0 0 0 0 0 0 0 1 1 0];
LB = [-5 -5 0 6 2 0.5 -Inf -Inf -Inf -Inf];
UB = [10 10 pi SDIM(1)/2 SDIM(1)/6 2 Inf Inf Inf Inf];
lfp_gabor_params = zeros(24,10);
LL_fin = zeros(24,length(wfreqs));
for cc = 1:24
fprintf('CC %d of %d\n',cc,24);
init_params(3) = orientations(maxloc(cc));
cur_init_params = [init_params zeros(1,4)];
[lfp_gabor_params(cc,:),LL_fin(cc)] = fit_gabor_params_lfp(cur_init_params,X_resh(cur_used_fixs,:),...
    squeeze(fix_avg_amps(cc,:)),[SDIM SDIM],hold_const,LB,UB);
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
rperm = randperm(length(cur_used_fixs));
xv_fixs = rperm(1:round(xv_frac*length(cur_used_fixs)));
xv_inds = find(ismember(all_model_fixids,xv_fixs));
% tr_inds = find(~ismember(all_model_fixids,xv_fixs));
tr_inds = 1:length(all_model_fixids);
tr_inds(all_model_blockids(tr_inds) > 3) = [];

% xv_fixs = find(all_stim_filtered(cur_used_fixs)==0);
% xv_inds = find(ismember(all_model_fixids,xv_fixs));
% tr_inds = find(~ismember(all_model_fixids,xv_fixs));
%% COMPUTE SACCADE TRIG AVG LFP AMPS
sac_trg_amp = zeros(24,length(tax));
taxb = [0 tax];
for i = 1:length(tax)
    curset = find(time_since_fix(all_used_inds) >= taxb(i) & time_since_fix(all_used_inds) < taxb(i+1));
    sac_trg_amp(:,i) = mean(all_amps(curset,:));
    n_occ(i) = length(curset);
end

%%
% used_freq = 8;
% sac_thresh = 1;
% min_sac_amp = 0;
% micro_sacs = find(all_sac_amps(used_fixs(cur_used_fixs)) < sac_thresh & all_sac_amps(used_fixs(cur_used_fixs)) > min_sac_amp);
% macro_sacs = find(all_sac_amps(used_fixs(cur_used_fixs)) > sac_thresh);
% 
% all_sacs = randperm(length(cur_used_fixs));
% micro_sacs = all_sacs(1:length(micro_sacs));
% macro_sacs = all_sacs(end-length(micro_sacs)+1:end);
% 
% % macro_sacs = macro_sacs(randperm(length(macro_sacs)));
% % macro_sacs(length(micro_sacs)+1:end) = [];
% 
% micro_set = find(ismember(all_model_fixids(all_used_inds),micro_sacs));
% macro_set = find(ismember(all_model_fixids(all_used_inds),macro_sacs));

%         micro_sacs = find(all_stim_filtered(used_fixs(cur_used_fixs)) == 1);
%         macro_sacs = find(all_stim_filtered(used_fixs(cur_used_fixs)) == 0);
%         micro_set = find(ismember(all_model_fixids(all_used_inds),micro_sacs));
%         macro_set = find(ismember(all_model_fixids(all_used_inds),macro_sacs));

% used_freq = 8;
for cc = 1:24
% cc = 14;

    %         disp('Fitting sac hist model')
    %         %first fit sachist model
    %         Xmat = Tmat;
    %
    %         Y = all_ampgrams(cc,:,used_freq)';
    %         XX = [Tmat(all_used_inds,:)];
    %         %     [B_sachist(cc,:),~] = REGRESS(Y,XX);
    %
    %         init_params = randn(length(tax),1);
    %         [B_sachist(cc,:)] = smoothed_regression(Y,XX,init_params,2e3,[1 length(tax)]);
    
    disp('Fitting total stim model')
    X_resh = reshape(X,size(X,1),SDIM^2);
    cur_mask1 = get_pgabor_mask(lfp_gabor_params(cc,1:6),0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(lfp_gabor_params(cc,1:6),pi/2,[SDIM SDIM]);
    mask1_out = X_resh*cur_mask1(:);
    mask2_out = X_resh*cur_mask2(:);
    energy_out = lfp_gabor_params(cc,7)*sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
    total_out = zscore(energy_out);
    
    
    
    %         %         sac_thresh = 1;
    % %         min_sac_amp = 0;
    % %         micro_sacs = cur_used_fixs(all_sac_amps(cur_used_fixs) < sac_thresh & all_sac_amps(cur_used_fixs) > min_sac_amp);
    % %         macro_sacs = cur_used_fixs(all_sac_amps(cur_used_fixs) > sac_thresh);
    % %         micro_set = find(ismember(all_model_fixids(all_used_inds),micro_sacs));
    % %         macro_set = find(ismember(all_model_fixids(all_used_inds),macro_sacs));
    %
    %         micro_sacs = cur_used_fixs(all_stim_filtered(cur_used_fixs) == 1);
    %         macro_sacs = cur_used_fixs(all_stim_filtered(cur_used_fixs) == 0);
    %         micro_set = find(ismember(all_model_fixids(all_used_inds),micro_sacs));
    %         macro_set = find(ismember(all_model_fixids(all_used_inds),macro_sacs));
    %
    %     Y = all_amps(micro_set,cc);
    %     XX = [Tmat(all_used_inds(micro_set),:) bsxfun(@times,Tmat(all_used_inds(micro_set),:),total_out(all_used_inds(micro_set)))];
    %     init_params = randn(2*length(tax),1);
    %     [B_stimmod_micro(cc,:)] = smoothed_regression(Y,XX,init_params,[2e3 5e3],[1 length(tax); length(tax)+1 2*length(tax)]);
    %
    %     Y = all_amps(macro_set,cc);
    %     XX = [Tmat(all_used_inds(macro_set),:) bsxfun(@times,Tmat(all_used_inds(macro_set),:),total_out(all_used_inds(macro_set)))];
    %     init_params = randn(2*length(tax),1);
    %     [B_stimmod_macro(cc,:)] = smoothed_regression(Y,XX,init_params,[2e3 5e3],[1 length(tax); length(tax)+1 2*length(tax)]);
    
        Y = all_amps(:,cc);
        XX = [Tmat(all_used_inds,:) bsxfun(@times,Tmat(all_used_inds,:),total_out(all_used_inds))];
        %     [B_stimmod(cc,:),~] = REGRESS(Y,XX);
        init_params = randn(2*length(tax),1);
        [B_stimmod(cc,:)] = smoothed_regression(Y,XX,init_params,[2e3 5e3],[1 length(tax); length(tax)+1 2*length(tax)]);
    %     [B_stimmod(cc,:)] = smoothed_regression(Y,XX,init_params,[2e3 5e3],[1 length(tax); length(tax)+1 2*length(tax)]);
    
%     
%     for bb = 1:3
%         set = find(blockids(used_fixs(cur_used_fixs)) == bb);
%         set = find(ismember(all_model_fixids(all_used_inds),set));
%         
%         Y = all_amps(set,cc);
%         XX = [Tmat(all_used_inds(set),:) bsxfun(@times,Tmat(all_used_inds(set),:),total_out(all_used_inds(set)))];
%         init_params = randn(2*length(tax),1);
%         [B_stimmod_set(cc,bb,:)] = smoothed_regression(Y,XX,init_params,[2e3 5e3],[1 length(tax); length(tax)+1 2*length(tax)]);
%     end
    
                disp('Fitting separate stim model')
    ori_vals = linspace(0,pi,13);
    ori_vals(end) = [];
    
    ori_outs = zeros(length(total_out),length(ori_vals));
    cur_params = lfp_gabor_params(cc,1:6);
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
    cur_mask1 = get_pgabor_mask(lfp_gabor_params(cc,1:6),0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(lfp_gabor_params(cc,1:6),pi/2,[SDIM SDIM]);
    mask1_out = X_resh*cur_mask1(:);
    mask2_out = X_resh*cur_mask2(:);
%     lin_out = gabor_params_fin(c,8)*mask1_out(all_model_fixids) + gabor_params_fin(c,9)*mask2_out(all_model_fixids);
    energy_out = sqrt(mask1_out(all_model_fixids).^2 + mask2_out(all_model_fixids).^2);
    total_out = energy_out;
    
    contrast_out = sum(ori_outs,2);
    total_out = total_out./contrast_out;
    
    contrast_out = zscore(contrast_out);
    total_out = zscore(total_out);
        
    Xmat = Tmat(all_used_inds,:);
    Xmat = [Xmat bsxfun(@times,Tmat(all_used_inds,:),total_out(all_used_inds))];
    Xmat = [Xmat bsxfun(@times,Tmat(all_used_inds,:),contrast_out(all_used_inds))];
    init_params = randn(3*length(tax),1);
    [sep_mod] = smoothed_regression(Y,Xmat,init_params,[2e3 5e3 5e3],[1 length(tax); length(tax)+1 2*length(tax); 2*length(tax)+1 3*length(tax)]);
    sep_mod_sachist(cc,:) = sep_mod(1:length(tax));
    sep_mod_stim(cc,:) = sep_mod(length(tax)+1:length(tax)*2);
    sep_mod_cont(cc,:) = sep_mod(length(tax)*2+1:length(tax)*3);

end

%%
figure
pcolor(tax,1:24,B_stimmod(:,31:end));shading flat
xlabel('Time (s)','fontsize',14)
ylabel('Channel Number','fontsize',14)

figure
pcolor(tax,1:24,B_stimmod(:,1:30));shading flat
xlabel('Time (s)','fontsize',14)
ylabel('Channel Number','fontsize',14)

%%
close all
ch = 16;
figure
pcolor(tax,wfreqs,squeeze(B_stimmod(ch,:,31:end)));shading flat
xlabel('Time (s)','fontsize',14)
ylabel('Frequency (Hz)','fontsize',14)
caxis([0 0.35])

figure
pcolor(tax,wfreqs,squeeze(B_stimmod(ch,:,1:30)));shading flat
xlabel('Time (s)','fontsize',14)
ylabel('Frequency (Hz)','fontsize',14)

%%
close all
cfreq = 8;
figure
pcolor(tax,1:24,squeeze(B_stimmod(:,cfreq,31:end)));shading flat
xlabel('Time (s)','fontsize',14)
ylabel('Frequency (Hz)','fontsize',14)
caxis([0 0.35])
colorbar

figure
pcolor(tax,1:24,squeeze(B_stimmod(:,cfreq,1:30)));shading flat
xlabel('Time (s)','fontsize',14)
ylabel('Frequency (Hz)','fontsize',14)
colorbar
