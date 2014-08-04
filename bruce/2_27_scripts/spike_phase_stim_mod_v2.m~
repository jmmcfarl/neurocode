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
% dt = .01;
dt = .0025;

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
for blockid = 1:3;
    fprintf('Block %d of %d\n',blockid,3);
    cur_set = find(blockids(used_fixs) == blockid);
    
    for i = 1:length(cur_set)
        start_T = all_fix_start_times(used_fixs(cur_set(i)));
        end_T = all_fix_stop_times(used_fixs(cur_set(i)));
        cur_tedges = start_T:dt:end_T;
        cur_tcents = 0.5*cur_tedges(1:end-1)+0.5*cur_tedges(2:end);
        all_model_time_axis = [all_model_time_axis cur_tcents];
        all_model_blockids = [all_model_blockids blockid*ones(1,length(cur_tcents))];
        all_model_fixids = [all_model_fixids cur_set(i)*ones(1,length(cur_tcents))];
%         temp_binned = zeros(10,length(cur_tcents));
%         temp_modouts = zeros(10,length(cur_tcents));
%         for c = 1:24
%             if c <= 10
%                 temp = histc(Blocks{blockid}.spktimes{c},cur_tedges);
%             else
%                 temp = histc(Blocks{blockid}.mutimes{c-10},cur_tedges);
%             end
%             temp_binned(c,:) = temp(1:end-1);
%         end
%         spikes_binned = [spikes_binned; temp_binned'];
        
        time_since_fix = [time_since_fix cur_tcents-start_T];
        
    end
end

%%
Fs = 1000;
dsf = 2;
Fsd = Fs/dsf;
niqf = Fs/2;
[b,a] = butter(2,[30 50]/niqf);


all_lfp_amp = [];
all_lfp_phase = [];
all_used_inds = [];
spk_tsf = cell(n_used_cells,1);
spk_stim = cell(n_used_cells,1);
spk_phases = cell(n_used_cells,1);
spk_fixids = cell(n_used_cells,1);

for blockid = 1:3;
    fprintf('Block %d of %d\n',blockid,3);
    
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
    
    lfp_samps = filtfilt(b,a,lfp_samps);
    lfp_sampsd = downsample(lfp_samps,dsf);
    lfp_h = hilbert(lfp_sampsd);
    lfp_phase = angle(lfp_h);
    lfp_amp = abs(lfp_h);
    clear lfp_samps lfp_h
    
    lfp_timed = downsample(lfp_time,dsf)';
    
    cur_all_model = find(all_model_blockids==blockid);
    model_lfp_inds = round(interp1(all_model_time_axis(cur_all_model),all_model_fixids(cur_all_model),lfp_timed));
    model_tsf = interp1(all_model_time_axis(cur_all_model),time_since_fix(cur_all_model),lfp_timed);
    bad = find(isnan(model_lfp_inds) | isnan(model_tsf));
    model_lfp_inds(bad) = 1;
    model_tsf(bad) = 1;
    X_resh = reshape(X,size(X,1),SDIM^2);
    for c = 1:n_used_cells
        if c <= 9
            cur_spk_times = Blocks{blockid}.spktimes{cellids(c)};
            cur_probe = Blocks{1}.suprobes(cellids(c));
        else
            cur_spk_times = Blocks{blockid}.mutimes{muaids(c-9)};
            cur_probe = Blocks{1}.muprobes(muaids(c-9));
        end
        cur_spk_inds = round(interp1(lfp_timed,1:length(lfp_timed),cur_spk_times));
        cur_spk_inds(isnan(cur_spk_inds)) = [];
        cur_spk_inds(ismember(cur_spk_inds,bad)) = [];

        cur_spk_phases = lfp_phase(cur_spk_inds,cur_probe);
        cur_spk_fixids = model_lfp_inds(cur_spk_inds);
        
        cur_mask1 = get_pgabor_mask(gabor_params_fin(c,1:6),0,[SDIM SDIM]);
        cur_mask2 = get_pgabor_mask(gabor_params_fin(c,1:6),pi/2,[SDIM SDIM]);
        mask1_out = X_resh*cur_mask1(:);
        mask2_out = X_resh*cur_mask2(:);
        lin_out = gabor_params_fin(c,8)*mask1_out(model_lfp_inds) + gabor_params_fin(c,9)*mask2_out(model_lfp_inds);
        energy_out = gabor_params_fin(c,7)*sqrt(mask1_out(model_lfp_inds).^2 + mask2_out(model_lfp_inds).^2);
        total_out = lin_out + energy_out;
        total_out = zscore(total_out);
        cur_spk_stim = total_out(cur_spk_inds);
        
        cur_spk_tsf = model_tsf(cur_spk_inds);
        
        spk_fixids{c} = [spk_fixids{c}; cur_spk_fixids];
        spk_tsf{c} = [spk_tsf{c}; cur_spk_tsf];
        spk_stim{c} = [spk_stim{c}; cur_spk_stim];
        spk_phases{c} = [spk_phases{c}; cur_spk_phases];
    end
  
    
end

%%
% n_fixs = length(unique(all_model_fixids(all_model_blockids < 4)));
% for c = 1:10
%     spk_isfirst{c} = zeros(size(spk_fixids{c}));
%     for i = 1:n_fixs
%         temp = find(spk_fixids{c}==i,1,'first');
%         spk_isfirst{c}(temp) = 1;
%     end
% end

%%
all_sac_dirs = atan2(all_sac_dy,all_sac_dx);

%%
% Now divide spikes into quantiles based on time since fixation onset, and look for circular-lin correlation between spike phase and stimulus
% Also, divide spikes into quantiles based on stimulus (maybe separately for early vs late spikes) and look for relationship between phase-locking strength and stimulus
n_sbins = 10;
clear stim_mean* stim_kappa*
for c = 1:n_used_cells
    cur_set = find(spk_tsf{c} > 0.15);
    mean_spkphase(c) = circ_mean(spk_phases{c}(cur_set));
    kappa_spkphase(c) = circ_kappa(spk_phases{c}(cur_set));
    p_val(c) = circ_rtest(spk_phases{c}(cur_set));
    
    [stim_spkphase(c),stim_spkphase_p(c)] = circ_corrcl(spk_phases{c}(cur_set),spk_stim{c}(cur_set));
    
    sbin_edges = prctile(spk_stim{c}(cur_set),linspace(0,100,n_sbins+1));
    sbin_cents = 0.5*sbin_edges(1:end-1)+0.5*sbin_edges(2:end);
    
    for s = 1:length(sbin_cents)
        cc_set = cur_set(find(spk_stim{c}(cur_set) >= sbin_edges(s) & spk_stim{c}(cur_set) < sbin_edges(s+1)));
        stim_mean_spkphase(c,s) = circ_mean(spk_phases{c}(cc_set));
        stim_kappa_spkphase(c,s) = circ_kappa(spk_phases{c}(cc_set));
        
    end
end





%%
c = 10;
n_bins = 15;
bin_edges = [0.001 0.025 0.05 0.075 0.1 0.125 0.15 0.25 1 10];

bin_cents = 0.5*bin_edges(1:end-1)+0.5*bin_edges(2:end);
n_bins = length(bin_cents);

spk_sac_amps = all_sac_amps(used_fixs(spk_fixids{c}));    
spk_sac_dirs = all_sac_dirs(used_fixs(spk_fixids{c}));    

ov_avg_sac_dir = circ_mean(all_sac_dirs);

clear cur_set mean_spkphase kappa_spkphase stim_spk* amp_spk* dir_spk*
for i = 1:n_bins
    cur_set{i} = find(spk_tsf{c} >= bin_edges(i) & spk_tsf{c} < bin_edges(i+1));
    n_spks(i) = length(cur_set{i});
%     avg_sac_dir(i) = circ_mean(spk_sac_dirs(cur_set{i}));
%     avg_sac_dist(i) = circ_dist(avg_sac_dir(i),ov_avg_sac_dir);
%     sac_cont_p(i) = circ_mtest(spk_sac_dirs(cur_set{i}),ov_avg_sac_dir);
         [stim_spktsf(i),stim_spktsf_p(i)] = corr(spk_tsf{c}(cur_set{i}),spk_stim{c}(cur_set{i}));
   for ww = 1:length(wfreqs)
        mean_spkphase(i,ww) = circ_mean(spk_phases{c}(cur_set{i},ww));
        kappa_spkphase(i,ww) = circ_kappa(spk_phases{c}(cur_set{i},ww));
        [stim_spkphase(i,ww),stim_spkphase_p(i,ww)] = circ_corrcl(spk_phases{c}(cur_set{i},ww),spk_stim{c}(cur_set{i}));
        [amp_spkphase(i,ww),amp_spkphase_p(i,ww)] = circ_corrcl(spk_phases{c}(cur_set{i},ww),spk_sac_amps(cur_set{i}));
        [dir_spkphase(i,ww),dir_spkphase_p(i,ww)] = circ_corrcc(spk_phases{c}(cur_set{i},ww),spk_sac_dirs(cur_set{i}));
    end
end
clf
pcolor(bin_cents,wfreqs,stim_spkphase_p');shading flat
xlim([0 0.4])
shg
caxis([0 0.05])

% clear cur_set mean_spkphase kappa_spkphase stim_spk* amp_spk* dir_spk* n_spks
% n_bins = 15;
% uspks = find(spk_tsf{c} < 0.1);
% bin_edges = prctile(spk_stim{c}(uspks),linspace(0,100,n_bins+1));
% bin_cents = 0.5*bin_edges(1:end-1)+0.5*bin_edges(2:end);
% for i = 1:n_bins
%     cur_set{i} = uspks(spk_stim{c}(uspks) >= bin_edges(i) & spk_stim{c}(uspks) < bin_edges(i+1));
%     n_spks(i) = length(cur_set{i});
%     for ww = 1:length(wfreqs)
%         mean_spkphase(i,ww) = circ_mean(spk_phases{c}(cur_set{i},ww));
%         kappa_spkphase(i,ww) = circ_kappa(spk_phases{c}(cur_set{i},ww));
%     end
% end

%%
c = 6;
spk_sac_amps = all_sac_amps(used_fixs(spk_fixids{c}));    
spk_sac_dirs = all_sac_dirs(used_fixs(spk_fixids{c}));    

cur_used_fixs = unique(all_model_fixids(all_model_blockids < 4));
used_sac_dirs = all_sac_dirs(used_fixs(cur_used_fixs));
used_sac_amps = all_sac_amps(used_fixs(cur_used_fixs));

% n_tbins = 10;
% bin_edges = prctile(spk_tsf{c},linspace(0,100,n_tbins+1));
% bin_cents = 0.5*bin_edges(1:end-1)+0.5*bin_edges(2:end);
bin_edges = [0.001 0.025 0.05 0.075 0.1 0.125 0.15 0.25 1 10];
bin_cents = 0.5*bin_edges(1:end-1)+0.5*bin_edges(2:end);
n_tbins = length(bin_cents);

n_bins = 8;
dir_bin_edges = linspace(-pi,pi,n_bins+1);
dir_bin_cents = 0.5*dir_bin_edges(1:end-1)+0.5*dir_bin_edges(2:end);

n_ampbins = 8;
amp_bin_edges = prctile(used_sac_amps,linspace(0,100,n_ampbins+1));
amp_bin_cents = 0.5*amp_bin_edges(1:end-1)+0.5*amp_bin_edges(2:end);
% 
% sac_dir_spk_rates = zeros(length(n_tbins),length(n_bins));
% for j = 1:n_bins
%     temp = find(used_sac_dirs >= dir_bin_edges(j) & used_sac_dirs < dir_bin_edges(j+1));
%     for i = 1:n_tbins
%         spk_set = find(spk_tsf{c} >= bin_edges(i) & spk_tsf{c} < bin_edges(i+1));
%         sac_dir_spk_rates(i,j) = sum(ismember(spk_fixids{c}(spk_set),temp))/length(temp);
%     end
% end
% 
% sac_dir_spk_rates = bsxfun(@rdivide,sac_dir_spk_rates,sum(sac_dir_spk_rates,2));

sac_amp_spk_rates = zeros(length(n_tbins),length(n_ampbins));
for j = 1:n_ampbins
    temp = find(used_sac_amps >= amp_bin_edges(j) & used_sac_amps < amp_bin_edges(j+1));
    for i = 1:n_tbins
        spk_set = find(spk_tsf{c} >= bin_edges(i) & spk_tsf{c} < bin_edges(i+1));
        ss = spk_set(ismember(spk_fixids{c}(spk_set),temp));
        sac_amp_spk_rates(i,j) = length(ss)/length(temp);
        for ww = 1:length(wfreqs)
            sac_amp_spk_mphase(i,j,ww) = circ_mean(spk_phases{c}(ss,ww));
            sac_amp_spk_kphase(i,j,ww) = circ_kappa(spk_phases{c}(ss,ww));
        end
    end
end
sac_amp_spk_kphase = bsxfun(@rdivide,sac_amp_spk_kphase,mean(sac_amp_spk_kphase,2));
for i = 1:n_tbins
    for ww = 1:length(wfreqs)
        sac_mphase(i,ww) = circ_mean(squeeze(sac_amp_spk_mphase(i,:,ww))');
    end
end
 for j = 1:n_ampbins
    for i = 1:n_tbins
        for ww = 1:length(wfreqs)
            sac_amp_spk_mphase(i,j,ww) = circ_dist(sac_amp_spk_mphase(i,j,ww),sac_mphase(i,ww));
        end
    end
end
   
sac_amp_spk_rates = bsxfun(@rdivide,sac_amp_spk_rates,sum(sac_amp_spk_rates,2));
    
    
 
%%
clc
c = 2;
bin_edges = [0.001 0.05 0.1 0.15 0.25 0.75];
bin_cents = 0.5*bin_edges(1:end-1)+0.5*bin_edges(2:end);
n_tbins = length(bin_cents);

n_sbins = 20;
sbin_edges = prctile(spk_stim{c},linspace(0,100,n_sbins+1));
sbin_cents = 0.5*sbin_edges(1:end-1)+0.5*sbin_edges(2:end);

clear stim_n_spikes stim_t_*
for t = 1:n_tbins
    spk_set = find(spk_tsf{c} >= bin_edges(t) & spk_tsf{c} < bin_edges(t+1));
    for s = 1:n_sbins
        cur_set = spk_set(spk_stim{c}(spk_set) >= sbin_edges(s) & spk_stim{c}(spk_set) < sbin_edges(s+1));
        for ww = 1:length(wfreqs)
            stim_t_kappa(t,s,ww) = circ_kappa(spk_phases{c}(cur_set,ww));
            stim_t_mean(t,s,ww) = circ_mean(spk_phases{c}(cur_set,ww));
            stim_n_spikes(t,s,ww) = length(cur_set);
        end
    end
end

% clf
% cmap = jet(n_tbins);
% subplot(2,1,1)
% for i = 1:n_tbins
% plot(sbin_cents,squeeze(stim_t_kappa(i,:,4)),'o-','color',cmap(i,:),'linewidth',2)
% hold on
% end
% subplot(2,1,2)
% for i = 1:n_tbins
% plot(sbin_cents,squeeze(stim_n_spikes(i,:,4)),'o-','color',cmap(i,:),'linewidth',2)
% hold on
% end
% shg

close all
figure
imagesc(squeeze(stim_t_kappa(:,:,4)));
figure
imagesc(squeeze(stim_t_mean(:,:,4)));
% figure
% imagesc(squeeze(stim_n_spikes(:,:,4)));
% colorbar
% 
for i = 1:n_tbins
    [a,b] = corrcoef(1:n_sbins,squeeze(stim_t_kappa(i,:,4)));
    [cca,ccb] = circ_corrcl(squeeze(stim_t_mean(i,:,4)),1:n_sbins);
    fprintf('Bin %d C: %.3f  p: %.3f  CCp: %.3f\n',i,a(2,1),b(2,1),ccb);
end

