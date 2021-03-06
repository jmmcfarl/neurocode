clear all
% close all

cd ~/Data/bruce/2_27_12/
load Blocks.mat

block_id = 1;
block_times = Blocks{block_id}.blocktimes;
stim_times = Blocks{block_id}.stimtime;
stimids = Blocks{block_id}.stimids;
is_stim_filtered = mod(stimids-1,4) >= 2;
n_stims = length(stim_times);
Pix2Deg = 1.1279 / 60;
dsfrac = 4;
Fsd = 1/Pix2Deg/dsfrac;


cd ~/Data/bruce/2_27_12/saccades/
load(sprintf('lemM232.5%d.em.sac.mat',block_id))
sac_times = Expt.Trials.EyeMsac.sacT;

cd ~/Data/bruce/2_27_12/stimrecon/
load all_image_mat_ds4
% load all_image_Energyarrays
% all_Image_array = all_Energy_array;

load surrogate_spikes4

load used_eyesig_p5
%%
stim_NY = size(all_Image_array,1);
stim_NX = size(all_Image_array,2);

stimres = 0.05; %resolution of time bins for STA calc
% sac_win = [-0.1 0.2]; %time window around saccades to reject data
sac_win = [0 0]; %time window around saccades to reject data
sac_win_bins = round(sac_win/stimres);

target_window = [-6 6; -6 6]; %acceptance window (in degrees) for gaze direction.
target_window_pix = target_window/Pix2Deg/dsfrac;

shift_time = 0;
% delta_bin = round(-0.06/stimres); %number of time bins to shift stim relative to spike data

%%
n_zpads = 200;
padded_image_array = cat(1,zeros(n_zpads,size(all_Image_array,2),size(all_Image_array,3)),all_Image_array);
padded_image_array = cat(1,padded_image_array,zeros(n_zpads,size(padded_image_array,2),size(all_Image_array,3)));
padded_image_array = cat(2,zeros(size(padded_image_array,1),n_zpads,size(all_Image_array,3)),padded_image_array);
padded_image_array = cat(2,padded_image_array,zeros(size(padded_image_array,1),n_zpads,size(all_Image_array,3)));

x_offset = round(stim_NX/2 + n_zpads);
y_offset = round(stim_NY/2 + n_zpads);

% padded_image_array = flipdim(padded_image_array,1); %flip y axis
% padded_image_array = flipdim(padded_image_array,2); %flip x axis

%%
eyepos = [Expt.Trials.Eyevals.rh Expt.Trials.Eyevals.rv];% right eye
% eyepos = [Expt.Trials.Eyevals.lh Expt.Trials.Eyevals.lv];%left eye
Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
eye_tx = Expt.Trials.Start/1e4:Eyedt:Expt.Trials.End/1e4;
eye_tx(end) = []; %throw out last bin so that aligns to eye data vector

eye_pix = interp1(eye_tx,eyepos,recon_t);
% eye_pix = ov_eyepos;
eye_pix = eye_pix/Pix2Deg/dsfrac;


%%
stim_id_vec = nan(size(recon_t));
for i = 1:n_stims-1
    cur_range = find(recon_t >= stim_times(i) & recon_t < stim_times(i+1));
    stim_id_vec(cur_range) = i;
end
last_stim_onset = find(stim_times(end) <= recon_t,1,'first');
stim_id_vec(last_stim_onset:end) = n_stims;

%find eye samples within the target window
good_eye_inds = (eye_pix(:,1) >= target_window_pix(1,1) & eye_pix(:,1) <= target_window_pix(1,2) & ...
    eye_pix(:,2) >= target_window_pix(2,1) & eye_pix(:,2) <= target_window_pix(2,2));

%locate time samples within some window of a saccade
sac_inds = hist(sac_times,recon_t);
sac_inds = find(sac_inds > 0);
post_sac_inds = zeros(size(recon_t));
edge_sacs = find(sac_inds <= -sac_win_bins(1) | sac_inds >= length(recon_t) - sac_win_bins(2));
for i = 1:length(edge_sacs)
    el = sac_inds(edge_sacs(i))+sac_win_bins;
    el(el < 1 | el > length(recon_t)) = [];
    post_sac_inds(el) = 1;
end
sac_inds(edge_sacs) = [];
for i = 1:length(sac_win_bins)
    post_sac_inds(sac_inds+sac_win_bins(i)) = 1;
end
post_sac_inds = logical(post_sac_inds');

used_samples = find(~post_sac_inds & good_eye_inds & recon_t' >= stim_times(1));

% %use only data during filtered stims
% used_samples = used_samples(is_stim_filtered(stim_id_vec(used_samples)));
%
% %use only data during unfiltered stims
% used_samples = used_samples(~is_stim_filtered(stim_id_vec(used_samples)));

%%
%find set of used spikes
spike_times_shifted = spiketimes+shift_time;

spikes_binned = hist(spike_times_shifted,recon_t);


%% compute STAs
window_x = x_offset+[-160:160];
window_y = y_offset+[-128:128];
spktrg_avg = zeros(length(window_y),length(window_x));
ov_avg = zeros(length(window_y),length(window_x));

n_spikes = sum(spikes_binned(used_samples));
f_rates = n_spikes/(length(used_samples)*stimres);

%% add noise to eye signal
% noise_vals = [0 0.1 0.25 0.5 1 1.5];
noise_vals = [0];
for n = 1:length(noise_vals)
        fprintf('%d of %d\n',n,length(noise_vals));
    
%     eye_pix = interp1(eye_tx,eyepos,recon_t);
%     eye_pix = eye_pix/Pix2Deg/dsfrac;
%     eye_pix = eye_pix + randn(size(eye_pix))*noise_vals(n)*Fsd;
%     
    eye_pix_round = round(eye_pix);
    for i = 1:length(used_samples)
%         fprintf('%d of %d\n',i,length(used_samples));
        cur_image = squeeze(padded_image_array(-eye_pix_round(used_samples(i),2)+window_y,eye_pix_round(used_samples(i),1)+window_x,stim_id_vec(used_samples(i))));
        ov_avg = ov_avg + cur_image;
        
        spktrg_avg = spktrg_avg + bsxfun(@times,spikes_binned(used_samples(i)),cur_image);
    end
    ov_avg = ov_avg/length(used_samples);
    spktrg_avg = bsxfun(@rdivide,spktrg_avg,n_spikes);
    
    spktrg_avg_cor(n,:,:) = spktrg_avg - ov_avg;
    
end

