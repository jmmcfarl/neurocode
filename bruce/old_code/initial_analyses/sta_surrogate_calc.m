clear all
% close all

cd ~/Data/bruce/2_27_12/
load Blocks.mat


block_id = 1;
% cell_id = 6;
% spk_times = Blocks{block_id}.spktimes{cell_id};
block_times = Blocks{block_id}.blocktimes;
stim_times = Blocks{block_id}.stimtime;
n_stims = length(stim_times);

cd ~/Data/bruce/2_27_12/saccades/
load(sprintf('lemM232.5%d.em.sac.mat',block_id))
sac_times = Expt.Trials.EyeMsac.sacT;

cd ~/Data/bruce/2_27_12/stimrecon/
load ./surrogate_spikes2.mat

% load stimrecon_t
load(['BLOCK' num2str(block_id) 'IMAGE' num2str(1) 'r2.mat'], 'STIMrec');

stimres = 0.025; %in s

stim_NY = size(STIMrec,2);
stim_NX = size(STIMrec,3);
% [X,Y] = meshgrid(1:2:stim_NX,1:2:stim_NY);
% [Xi,Yi] = meshgrid(1:stim_NX,1:stim_NY);

%%
n_units = length(Blocks{block_id}.spktimes);
delta_bin_vals = [-8:0];

ov_avg_stim = zeros(stim_NY,stim_NX);
ov_spktrg_avg = zeros(stim_NY,stim_NX);
ov_spktrg_avg2 = zeros(stim_NY,stim_NX);
n_tot_spks = 0;
n_tot_spks2 = 0;
n_used_samps = zeros(n_stims,1);
n_tot_samps = zeros(n_stims,1);

% sac_win = [0.05 0.15];
sac_win = [0 0];

target_window = [-5 5; -5 5];
% target_window = [-8 8; -8 8];

eyepos = [Expt.Trials.Eyevals.rh Expt.Trials.Eyevals.rv];
Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
eye_tx = Expt.Trials.Start/1e4:Eyedt:Expt.Trials.End/1e4;

delta_bin = 0; %number of time bins to shift stim relative to spike data
for i = 1:n_stims
    i
    %     load(['BLOCK' num2str(block_id) 'IMAGE' num2str(i) 'hres.mat'], 'STIMrec');
    load(['BLOCK' num2str(block_id) 'IMAGE' num2str(i) 'r2.mat']);
    n_recon_samps(i) = size(STIMrec,1);
    
    cur_stim_set = find(stim_num==i);
    cur_stim_taxis = recon_t(cur_stim_set);
    cur_bin_edges = [(cur_stim_taxis-stimres/2) (cur_stim_taxis(end)+stimres/2)];
    
    cur_eye_data = find(eye_tx > cur_stim_taxis(1) & eye_tx < cur_stim_taxis(end));
    interp_eye_data = interp1(eye_tx(cur_eye_data),eyepos(cur_eye_data,:),cur_stim_taxis);
    good_eye_data = (interp_eye_data(:,1) >= target_window(1,1) & interp_eye_data(:,1) <= target_window(1,2) & ...
        interp_eye_data(:,2) >= target_window(2,1) & interp_eye_data(:,2) <= target_window(2,2));
    
    cur_sacs = sac_times(sac_times > cur_stim_taxis(1) & sac_times < cur_stim_taxis(end));
    post_sac_times = zeros(size(cur_stim_taxis));
    for s = 1:length(cur_sacs)
        post_sac_times(cur_stim_taxis >= cur_sacs(s)-sac_win(1) & cur_stim_taxis < cur_sacs(s) + sac_win(2)) = 1;
    end
    good_times = find(post_sac_times == 0 & good_eye_data');
    
% %just use all data for now
%     good_times = 1:length(cur_stim_taxis);    
    
    n_used_samps(i) = length(good_times);
    n_tot_samps(i) = length(cur_stim_taxis);
    
    binned_spks = histc(spiketimes,cur_bin_edges);
    binned_spks(end) = [];
    
    un_spk_cnts = unique(binned_spks);
    spk_bins = [];
    for tt = 2:length(un_spk_cnts)
        spk_bins = [spk_bins repmat(find(binned_spks==un_spk_cnts(tt)),1,un_spk_cnts(tt))];
    end
    spk_bins = spk_bins + delta_bin;
    spk_bins(spk_bins <= 0 | spk_bins > n_recon_samps(i)) = [];
    
    %use only times outside of post-saccade interval
    spk_bins = spk_bins(ismember(spk_bins,good_times));
    
    cur_spkavg_stim = squeeze(nansum(STIMrec(spk_bins,:,:),1));
    ov_spktrg_avg = squeeze(ov_spktrg_avg) + cur_spkavg_stim;
   
%     spk_bins2 = all_spike_bins{i};
%     cur_spkavg_stim = squeeze(nansum(STIMrec(spk_bins2,:,:),1));
%     ov_spktrg_avg2 = squeeze(ov_spktrg_avg2) + cur_spkavg_stim;

    n_tot_spks = n_tot_spks + length(spk_bins);
%     n_tot_spks2 = n_tot_spks2 + length(spk_bins2);

    cur_avg_stim = squeeze(nanmean(STIMrec(good_times,:,:)));
    ov_avg_stim = ov_avg_stim + cur_avg_stim*n_used_samps(i);
end

ov_avg_stim = ov_avg_stim/sum(n_used_samps);
ov_spktrg_avg = bsxfun(@rdivide,ov_spktrg_avg,n_tot_spks);
% ov_spktrg_avg2 = bsxfun(@rdivide,ov_spktrg_avg2,n_tot_spks2);
ov_spktrg_avg_cor = ov_spktrg_avg - ov_avg_stim;
% ov_spktrg_avg_cor2 = ov_spktrg_avg2 - ov_avg_stim;

%%

