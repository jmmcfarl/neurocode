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

cellids = [1 2 3 4 5 6 7 8 9 10];
muaids = [1 2 3 4 5 6 7 8 9 10 11 12 13 14];

spk_cnts = [spk_cnts(:,cellids) mua_cnts(:,muaids)];

n_used_cells = size(spk_cnts,2);
su_probes = [Blocks{1}.suprobes(cellids)];
mu_probes = [Blocks{1}.muprobes(muaids)];
all_probes = [su_probes mu_probes];


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
dsf = 3;
Fsd = Fs/dsf;
niqf = Fs/2;
[bb,aa] = butter(2,120/niqf,'low');

scales = logspace(log10(4),log10(150),50);
% scales = [scales 70 80 90];
wfreqs = scal2frq(scales,'cmor1-1',1/Fsd);
nwfreqs = length(wfreqs);

all_lfp_amp = [];
all_lfp_phase = [];
all_used_inds = [];
spk_tsf = cell(n_used_cells,1);
spk_stim = cell(n_used_cells,1);
spk_phases = cell(n_used_cells,1);
spk_fixids = cell(n_used_cells,1);


%%
for blockid = 1:4;
    fprintf('Block %d of %d\n',blockid,4);
    
    %%
    cd /Users/James/Data/bruce/2_27_12
    load(sprintf('lemM232A.5%d.lfp.mat',blockid));
    if blockid == 4
        LFP.Trials = LFP.Trials(1:5);
    end
    
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
    
    lfp_samps = filtfilt(bb,aa,lfp_samps);
    lfp_sampsd = zscore(downsample(lfp_samps,dsf));
    clear lfp_samps 
    
    lfp_timed = downsample(lfp_time,dsf)';
    
    cur_all_model = find(all_model_blockids==blockid);
    
    all_phasegrams = zeros(24,length(lfp_timed),length(wfreqs));
    for cc = 1:24
        fprintf('Channel %d of %d\n',cc,24);
        temp = cwt(lfp_sampsd(:,cc),scales,'cmor1-1');
        temp = temp';
        tampgram = abs(temp);
        tphasegram = angle(temp);
        all_phasegrams(cc,:,:) = tphasegram;
    end
    
    model_lfp_inds = round(interp1(all_model_time_axis(cur_all_model),all_model_fixids(cur_all_model),lfp_timed));
    model_tsf = interp1(all_model_time_axis(cur_all_model),time_since_fix(cur_all_model),lfp_timed);
    for c = 1:n_used_cells
        if c <= 10
            cur_spk_times = Blocks{blockid}.spktimes{cellids(c)};
            cur_probe = Blocks{1}.suprobes(cellids(c));
        else
            cur_spk_times = Blocks{blockid}.mutimes{muaids(c-10)};
            cur_probe = Blocks{1}.muprobes(muaids(c-10));
        end
        cur_spk_inds = round(interp1(lfp_timed,1:length(lfp_timed),cur_spk_times));
        cur_spk_inds(isnan(cur_spk_inds)) = [];

        cur_spk_phases = all_phasegrams(:,cur_spk_inds,:);
        cur_spk_tsf = model_tsf(cur_spk_inds);
        
        spk_phases{c} = cat(2,spk_phases{c},cur_spk_phases);
            spk_tsf{c} = [spk_tsf{c}; cur_spk_tsf];
end
  
    
end

%%
for c = 1:n_used_cells
   use_spks = find(spk_tsf{c} >= 0.15); 
   
   late_phase_locking(c,:,:) = squeeze(nansum(exp(1i*spk_phases{c}(:,use_spks,:)),2));
   late_phase_locking(c,:,:) = abs(late_phase_locking(c,:,:))/sum(~isnan(spk_phases{c}(1,use_spks,1)),2);

   all_phase_locking(c,:,:) = squeeze(nansum(exp(1i*spk_phases{c}),2));
   all_phase_locking(c,:,:) = abs(all_phase_locking(c,:,:))/sum(~isnan(spk_phases{c}(1,:,1)),2);
  
   use_spks = find(spk_tsf{c} < 0.15); 

   early_phase_locking(c,:,:) = squeeze(nansum(exp(1i*spk_phases{c}(:,use_spks,:)),2));
   early_phase_locking(c,:,:) = abs(early_phase_locking(c,:,:))/sum(~isnan(spk_phases{c}(1,use_spks,1)),2);
end

%%
close all
[~,pord] = sort(all_probes);
for c = 1:n_used_cells
   all_probes(c)
   subplot(3,1,1)
   pcolor(wfreqs,1:24,squeeze(late_phase_locking(pord(c),:,:)));shading flat;colorbar
   cm = caxis();
   subplot(3,1,2)
   pcolor(wfreqs,1:24,squeeze(early_phase_locking(pord(c),:,:)));shading flat;colorbar
   caxis(cm);
   subplot(3,1,3)
   pcolor(wfreqs,1:24,squeeze(all_phase_locking(pord(c),:,:)));shading flat;colorbar
   caxis(cm);
   pause
   clf
end

%%
avg_late_phaselocking = squeeze(mean(late_phase_locking,2));
avg_early_phaselocking = squeeze(mean(early_phase_locking,2));
avg_all_phaselocking = squeeze(mean(all_phase_locking,2));

% use_probes = all_probes+1;
% use_probes(use_probes > 24) = 23;
use_probes = all_probes;
best_late_phaselocking = nan(size(avg_late_phaselocking));
best_early_phaselocking = nan(size(avg_late_phaselocking));
best_all_phaselocking = nan(size(avg_late_phaselocking));
for c = 1:n_used_cells
    c
    cur_use_channels = use_probes(c) + [-1 1];
    cur_use_channels(cur_use_channels < 1 | cur_use_channels > 24) = [];
    best_late_phaselocking(c,:) = squeeze(mean(late_phase_locking(c,cur_use_channels,:),2));
    best_early_phaselocking(c,:) = squeeze(mean(early_phase_locking(c,cur_use_channels,:),2));
    best_all_phaselocking(c,:) = squeeze(mean(all_phase_locking(c,cur_use_channels,:),2));
end


%%
save allunit_fv_phaselocking wfreqs all_probes *phase_locking 

%%
sus = 1:10;
mus = 11:24;
% figure
% shadedErrorBar(wfreqs,mean(best_late_phaselocking(sus,:)),std(best_late_phaselocking(sus,:))/sqrt(length(sus)));
% hold on
% shadedErrorBar(wfreqs,mean(best_late_phaselocking(mus,:)),std(best_late_phaselocking(mus,:))/sqrt(length(mus)),{'color','r'});
% 
% figure
% shadedErrorBar(wfreqs,mean(best_all_phaselocking(sus,:)),std(best_all_phaselocking(sus,:))/sqrt(length(sus)));
% hold on
% shadedErrorBar(wfreqs,mean(best_all_phaselocking(mus,:)),std(best_all_phaselocking(mus,:))/sqrt(length(mus)),{'color','r'});

figure
shadedErrorBar(wfreqs,mean(best_late_phaselocking),std(best_late_phaselocking)/sqrt(24));
hold on
shadedErrorBar(wfreqs,mean(best_all_phaselocking),std(best_all_phaselocking)/sqrt(24),{'color','r'});
xlim(wfreqs([end 1]));
set(gca,'xscale','log');

%%
% for c = 1:10
c = 10;
figure
    pcolor(wfreqs,1:24,squeeze(late_phase_locking(c,:,:)));shading interp;
    set(gca,'xscale','log');
    set(gca,'fontsize',16,'fontname','arial')
%     pause
%     clf
% end