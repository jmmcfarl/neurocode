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
cd /Users/James/James_scripts/bruce/modelfits/
% load fixation_gabor_models
load gabor_tempmodfits_en_lin
%%
Fs = 1000;
dsf = 2;
Fsd = Fs/dsf;
niqf = Fs/2;
[b,a] = butter(2,[1 200]/niqf);

scales = logspace(log10(3),log10(150),50);
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

%%
cd ~/Data/bruce/2_27_12


all_ampgrams = [];
all_used_inds = [];

for blockid = 1:3;
    
    fprintf('Processing block %d of %d...\n',blockid,4);
    
    block_times = Blocks{blockid}.blocktimes;
    stim_times = Blocks{blockid}.stimtime;
    stimID = Blocks{blockid}.stimids;
    
    %%
    cd /Users/James/Data/bruce/2_27_12
    load(sprintf('lemM232A.5%d.lfp.mat',blockid));
    %get start times of each LFP trial
    n_lfp_trials = length(LFP.Trials);
    lfp_trial_start = nan(n_lfp_trials,1);
    lfp_trial_stop = nan(n_lfp_trials,1);
    for i = 1:n_lfp_trials
        lfp_trial_start(i) = LFP.Trials(i).Start/1e4; %originally in tenths of ms
        lfp_trial_stop(i) = LFP.Trials(i).End/1e4;
        lfp_dur(i) = size(LFP.Trials(i).LFP,1)/1000;
    end
    
    lfp_time = [];
    lfp_samps = [];
    for i = 1:n_lfp_trials
        lfp_time = [lfp_time linspace(lfp_trial_start(i),lfp_trial_start(i)+lfp_dur(i),size(LFP.Trials(i).LFP,1))];
        lfp_samps = [lfp_samps; LFP.Trials(i).LFP];
    end
    % lfp_samps = zscore(lfp_samps);
    
    lfp_samps = filtfilt(b,a,lfp_samps);
    
    lfp_sampsd = downsample(lfp_samps,dsf);
    lfp_sampsd = zscore(lfp_sampsd);
    lfp_timed = downsample(lfp_time,dsf);
    
    %%
    cur_all_model = find(all_model_blockids==blockid);
    interp_ampgram = zeros(24,length(cur_all_model),length(wfreqs));
    for cc = 1:24
        %     cc = 17;
        fprintf('Channel %d of %d\n',cc,24);
        temp = cwt(lfp_sampsd(:,cc),scales,'cmor1-1');
        tampgram = abs(temp);
        tampgram = bsxfun(@minus,tampgram,mean(tampgram,2));
        tampgram = bsxfun(@rdivide,tampgram,std(tampgram,[],2));
        tampgram = tampgram';       
        interp_ampgram(cc,:,:) = interp1(lfp_timed,tampgram,all_model_time_axis(cur_all_model));
   end
    clear temp tampgram
    %%
    used_tpoints = find(~isnan(interp_ampgram(1,:,1)));    
    all_ampgrams = cat(2,all_ampgrams,interp_ampgram(:,used_tpoints,:));
    all_used_inds = [all_used_inds; cur_all_model(used_tpoints)'];
    
end


%%
sac_trg_amp = zeros(24,length(wfreqs),length(tax));
taxb = [0 tax];
for i = 1:length(tax)
    curset = find(time_since_fix(all_used_inds) >= taxb(i) & time_since_fix(all_used_inds) < taxb(i+1));
    sac_trg_amp(:,:,i) = mean(all_ampgrams(:,curset,:),2);
    n_occ(i) = length(curset);
end

%%
sac_pred_amps = zeros(24,length(all_used_inds),length(wfreqs));
for cc = 1:24
    sac_pred_amps(cc,:,:) = Tmat(all_used_inds,:)*squeeze(sac_trg_amp(cc,:,:))';
end
sac_resid_amps = all_ampgrams - sac_pred_amps;

%%
cellid = 7;
for cc = 1:24
    [lfp_stimout_corr(cc,:),~] = corr(squeeze(sac_resid_amps(cc,:,:)),stim_mod_stimout(cellid,all_used_inds)');
end

%% COMPUTE FIXATION AVERAGE SPECTRA
cd /Users/James/Data/bruce/2_27_12/stimrecon
load fixation_image_patches_corrected
[NT,SDIM,~] = size(X);
Xmatr = reshape(X,NT,SDIM^2);

used_fixs = unique(all_model_fixids(all_used_inds));
n_fixs = length(used_fixs);
fix_avg_amps = zeros(24,n_fixs,length(wfreqs));
for i = 1:n_fixs
    cur_set = find(all_model_fixids(all_used_inds) == used_fixs(i));
    fix_avg_amps(:,i,:) = squeeze(mean(all_ampgrams(:,cur_set,:),2));
end
fix_avg_amps = bsxfun(@minus,fix_avg_amps,mean(fix_avg_amps,2));
fix_avg_amps = bsxfun(@rdivide,fix_avg_amps,std(fix_avg_amps,[],2));

%%

c = 8;
cur_lfp = Blocks{1}.suprobes(c);
hold_const = [1 1 0 0 0 0 0 1 1 0];
LB = [-5 -5 0 6 2 0.5 -Inf -Inf -Inf -Inf];
UB = [10 10 pi SDIM(1)/2 SDIM(1)/6 2 Inf Inf Inf Inf];
init_params = gabor_params{end}(c,:);
init_params(7:10) = 0;
for i = 1:length(wfreqs)
    fprintf('Freq %d of %d\n',i,length(wfreqs));
    [lfp_gabor_params(i,:),LL_fin(i)] = fit_gabor_params_lfp(init_params,Xmatr(used_fixs,:),...
        squeeze(fix_avg_amps(cur_lfp,:,i)),[SDIM SDIM],hold_const,LB,UB);
end

orientations = linspace(0,pi-pi/12,12);
init_params = [5 -5 0 15 6 1];
for i = 1:12
    i
    init_params(3) = orientations(i);
    cur_mask1 = get_pgabor_mask(init_params,0,[SDIM SDIM]);
    cur_mask2 = get_pgabor_mask(init_params,pi/2,[SDIM SDIM]);
    
    mask1_out = Xmatr(used_fixs,:)*cur_mask1(:);
    mask2_out = Xmatr(used_fixs,:)*cur_mask2(:);
    
    en_out(i,:) = sqrt(mask1_out.^2+mask2_out.^2);
end
en_out = bsxfun(@minus,en_out,mean(en_out,2));
en_out = bsxfun(@rdivide,en_out,std(en_out,[],2));
for ww = 1:length(wfreqs)
    for cc = 1:24
        for i = 1:12
            temp = [en_out(i,:)' ones(n_fixs,1)];
            [B,BINT,R,RINT,STATS] = regress(squeeze(fix_avg_amps(cc,:,ww))',temp);
            rsq(cc,ww,i) = STATS(1);
        end
    end
end

for i = 1:24
    fprintf('CC %d of %d\n',i,24);
    init_params(3) = orientations(maxloc(i,7));
    cur_init_params = [init_params zeros(1,4)];
    [lfp_gabor_params2(i,:),LL_fin2(i)] = fit_gabor_params_lfp(init_params,Xmatr(used_fixs,:),...
        squeeze(fix_avg_amps(i,:,17)),[SDIM SDIM],hold_const,LB,UB);
end
