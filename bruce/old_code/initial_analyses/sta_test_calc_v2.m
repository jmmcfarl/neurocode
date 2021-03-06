clear all
close all

cd ~/Data/bruce/2_27_12/
load Blocks.mat

block_id = 1;
block_times = Blocks{block_id}.blocktimes;
stim_times = Blocks{block_id}.stimtime;
stimids = Blocks{block_id}.stimids;
is_stim_filtered = mod(stimids-1,4) >= 2;
n_stims = length(stim_times);
Pix2Deg = 1.1279 / 60;
dsfrac = 2;

cd ~/Data/bruce/2_27_12/saccades/
load(sprintf('lemM232.5%d.em.sac.mat',block_id))
sac_times = Expt.Trials.EyeMsac.sacT;

cd ~/Data/bruce/2_27_12/stimrecon/
load all_image_mat


%%
n_units = length(Blocks{block_id}.spktimes);

stim_NY = size(all_Image_array,1);
stim_NX = size(all_Image_array,2);

stimres = 0.01; %resolution of time bins for STA calc
sac_win = [-0.1 0.2]; %time window around saccades to reject data
sac_win_bins = round(sac_win/stimres);

target_window = [-5 5; -5 5]; %acceptance window (in degrees) for gaze direction.
target_window_pix = target_window/Pix2Deg/dsfrac;

delta_bin = round(-0.06/stimres); %number of time bins to shift stim relative to spike data

%%
n_zpads = 350;
padded_image_array = cat(1,zeros(n_zpads,size(all_Image_array,2),size(all_Image_array,3)),all_Image_array);
padded_image_array = cat(1,padded_image_array,zeros(n_zpads,size(padded_image_array,2),size(all_Image_array,3)));
padded_image_array = cat(2,zeros(size(padded_image_array,1),n_zpads,size(all_Image_array,3)),padded_image_array);
padded_image_array = cat(2,padded_image_array,zeros(size(padded_image_array,1),n_zpads,size(all_Image_array,3)));

x_offset = round(stim_NX/2 + n_zpads);
y_offset = round(stim_NY/2 + n_zpads);

% padded_image_array = flipdim(padded_image_array,1); %flip y axis
% padded_image_array = flipdim(padded_image_array,2); %flip x axis

%%
% eyepos = [Expt.Trials.Eyevals.rh Expt.Trials.Eyevals.rv];% right eye
eyepos = [Expt.Trials.Eyevals.lh Expt.Trials.Eyevals.lv];%left eye
Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
eye_tx = Expt.Trials.Start/1e4:Eyedt:Expt.Trials.End/1e4;
eye_tx(end) = []; %throw out last bin so that aligns to eye data vector

time_axis = block_times(1,1):stimres:block_times(2,end);

eye_pix = interp1(eye_tx,eyepos,time_axis);
eye_pix = eye_pix/Pix2Deg/dsfrac;

stim_id_vec = nan(size(time_axis));
for i = 1:n_stims-1
    cur_range = find(time_axis >= stim_times(i) & time_axis < stim_times(i+1));
    stim_id_vec(cur_range) = i;
end
last_stim_onset = find(stim_times(end) <= time_axis,1,'first');
stim_id_vec(last_stim_onset:end) = n_stims;

%find eye samples within the target window
good_eye_inds = (eye_pix(:,1) >= target_window_pix(1,1) & eye_pix(:,1) <= target_window_pix(1,2) & ...
    eye_pix(:,2) >= target_window_pix(2,1) & eye_pix(:,2) <= target_window_pix(2,2));

%locate time samples within some window of a saccade
sac_inds = hist(sac_times,time_axis);
sac_inds = find(sac_inds > 0);
post_sac_inds = zeros(size(time_axis));
edge_sacs = find(sac_inds <= -sac_win_bins(1) | sac_inds >= length(time_axis) - sac_win_bins(2));
for i = 1:length(edge_sacs)
    el = sac_inds(edge_sacs(i))+sac_win_bins;
    el(el < 1 | el > length(time_axis)) = [];
    post_sac_inds(el) = 1;
end
sac_inds(edge_sacs) = [];
for i = 1:length(sac_win_bins)
    post_sac_inds(sac_inds+sac_win_bins(i)) = 1;
end
post_sac_inds = logical(post_sac_inds');

used_samples = find(~post_sac_inds & good_eye_inds & time_axis' >= stim_times(1));

% %use only data during filtered stims
% used_samples = used_samples(is_stim_filtered(stim_id_vec(used_samples)));

%use only data during unfiltered stims
used_samples = used_samples(~is_stim_filtered(stim_id_vec(used_samples)));

%%
% cell_id = 4; %cell to work on
for cell_id = 1:10
    
    fprintf('Analyzing cell %d of %d\n',cell_id,10);
    
    %find set of used spikes
    spike_times = Blocks{block_id}.spktimes{cell_id};
    spikes_binned = hist(spike_times,time_axis);
    
    un_spk_cnts = unique(spikes_binned);
    spk_bins = [];
    for tt = 2:length(un_spk_cnts)
        spk_bins = [spk_bins repmat(find(spikes_binned==un_spk_cnts(tt)),1,un_spk_cnts(tt))];
    end
    spk_bins = sort(spk_bins);
    spk_bins(spk_bins <= 0 | spk_bins > length(time_axis)) = [];
    
    shifted_spk_bins = spk_bins + delta_bin;
    shifted_spk_bins = shifted_spk_bins(ismember(shifted_spk_bins,used_samples));
    
    spike_eye_pos = round(eye_pix(shifted_spk_bins,:));
    spike_im_ids = stim_id_vec(shifted_spk_bins);
    
    %% compute STA
    window_x = x_offset+[-300:300];
    window_y = y_offset+[-300:300];
    
    spktrg_avg(cell_id,:,:) = zeros(length(window_y),length(window_x));
    n_spikes(cell_id) = length(shifted_spk_bins);
    for i = 1:length(shifted_spk_bins)
        %   fprintf('spike %d of %d\n',i,n_spikes(cell_id));
        cur_image = squeeze(padded_image_array(:,:,spike_im_ids(i)));
        spktrg_avg(cell_id,:,:) = squeeze(spktrg_avg(cell_id,:,:)) + cur_image(spike_eye_pos(i,2)+window_y,spike_eye_pos(i,1)+window_x);
    end
    spktrg_avg(cell_id,:,:) = spktrg_avg(cell_id,:,:)/n_spikes(cell_id);
    
end

%%
% n_spikes = 2000;
n_rand_spikes = 5000;
rand_bins = used_samples(ceil(rand(n_rand_spikes,1)*length(used_samples)));
rspike_eye_pos = round(eye_pix(rand_bins,:));
rspike_im_ids = stim_id_vec(rand_bins);

rand_spktrg_avg = zeros(length(window_y),length(window_x));
for i = 1:length(rand_bins)
    %     fprintf('spike %d of %d\n',i,n_rand_spikes);
    cur_image = squeeze(padded_image_array(:,:,rspike_im_ids(i)));
    rand_spktrg_avg = rand_spktrg_avg + cur_image(rspike_eye_pos(i,2)+window_y,rspike_eye_pos(i,1)+window_x);
end
rand_spktrg_avg = rand_spktrg_avg/n_rand_spikes;

%%
i = 3;
clf
subplot(2,1,1)
spktrg_avg_cor = squeeze(spktrg_avg(i,:,:));
imagesc(spktrg_avg_cor)
subplot(2,1,2)
spktrg_avg_cor = squeeze(spktrg_avg(i,:,:)) - rand_spktrg_avg;
imagesc(spktrg_avg_cor)


