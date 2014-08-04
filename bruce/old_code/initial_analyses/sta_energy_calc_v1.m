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
dsfrac = 4;

cd ~/Data/bruce/2_27_12/saccades/
load(sprintf('lemM232.5%d.em.sac.mat',block_id))
sac_times = Expt.Trials.EyeMsac.sacT;

cd ~/Data/bruce/2_27_12/stimrecon/
load all_image_Energyarrays_p5_3

%%
n_units = length(Blocks{block_id}.spktimes);

stim_NY = size(all_Energy_array,1);
stim_NX = size(all_Energy_array,2);

stimres = 0.05; %resolution of time bins for STA calc
sac_win = [-0.05 0.1]; %time window around saccades to reject data
sac_win_bins = round(sac_win/stimres);

target_window = [-5 5; -5 5]; %acceptance window (in degrees) for gaze direction.
target_window_pix = target_window/Pix2Deg/dsfrac;

shift_time = -0.05;
% delta_bin = round(-0.06/stimres); %number of time bins to shift stim relative to spike data

%%
n_zpads = 100;
padded_energy_array = cat(1,zeros(n_zpads,size(all_Energy_array,2),size(all_Energy_array,3),size(all_Energy_array,4)),...
    all_Energy_array,zeros(n_zpads,size(all_Energy_array,2),size(all_Energy_array,3),size(all_Energy_array,4)));
padded_energy_array = cat(2,zeros(size(padded_energy_array,1),n_zpads,size(padded_energy_array,3),size(padded_energy_array,4)),...
    padded_energy_array,zeros(size(padded_energy_array,1),n_zpads,size(padded_energy_array,3),size(padded_energy_array,4)));

x_offset = round(stim_NX/2 + n_zpads);
y_offset = round(stim_NY/2 + n_zpads);

padded_image_array = flipdim(padded_energy_array,1); %flip y axis
% padded_image_array = flipdim(padded_energy_array,2); %flip x axis

%%
eyepos = [Expt.Trials.Eyevals.rh Expt.Trials.Eyevals.rv];% right eye
% eyepos = [Expt.Trials.Eyevals.lh Expt.Trials.Eyevals.lv];%left eye
Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
eye_tx = Expt.Trials.Start/1e4:Eyedt:Expt.Trials.End/1e4;
eye_tx(end) = []; %throw out last bin so that aligns to eye data vector

% time_axis = block_times(1,1):stimres:block_times(2,end);
time_axis = [];
for i = 1:length(stim_times)-1
    start_T = stim_times(i);
    end_T = stim_times(i+1)-0.2;
    time_axis = [time_axis start_T:stimres:end_T];
end

between_blocks = zeros(size(time_axis));
for i = 1:size(block_times,2)-1
    cur_set = find(time_axis >= block_times(2,i) & time_axis < block_times(1,i+1));
    between_blocks(cur_set) = 1;
end

eye_pix = interp1(eye_tx,eyepos,time_axis);
eye_pix = eye_pix/Pix2Deg/dsfrac;
eye_grid = eye_pix/grid_dx;

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

used_samples = find(~post_sac_inds & good_eye_inds & time_axis' >= stim_times(1) & between_blocks' == 0);

% %use only data during filtered stims
% used_samples = used_samples(is_stim_filtered(stim_id_vec(used_samples)));
% 
% %use only data during unfiltered stims
% used_samples = used_samples(~is_stim_filtered(stim_id_vec(used_samples)));

%%
for cell_id = 1:10
    
    %find set of used spikes
    spike_times = Blocks{block_id}.spktimes{cell_id};
    spike_times_shifted = spike_times + shift_time;
    
    spikes_binned(cell_id,:) = hist(spike_times_shifted,time_axis);
    
end

%% compute STAs
window_x = x_offset+round([-300:300]/grid_dx);
window_y = y_offset+round([-300:300]/grid_dx);
spktrg_avg = zeros(10,length(window_y),length(window_x),4);
ov_avg = zeros(length(window_y),length(window_x),4);

n_spikes = sum(spikes_binned(:,used_samples),2);
f_rates = n_spikes/(length(used_samples)*stimres);

eye_bin_round = round(eye_grid);
for i = 1:length(used_samples)
    fprintf('%d of %d\n',i,length(used_samples));
    cur_image = squeeze(padded_energy_array(-eye_bin_round(used_samples(i),2)+window_y,...
        eye_bin_round(used_samples(i),1)+window_x,:,stim_id_vec(used_samples(i))));
    ov_avg = ov_avg + cur_image;
    
    spktrg_avg = spktrg_avg + bsxfun(@times,spikes_binned(:,used_samples(i)),...
        repmat(reshape(cur_image,1,length(window_y),length(window_x),4),[10 1 1 1]));
    
end
ov_avg = ov_avg/length(used_samples);
spktrg_avg = bsxfun(@rdivide,spktrg_avg,n_spikes);

spktrg_avg_cor = spktrg_avg - repmat(reshape(ov_avg,[1 length(window_y),length(window_x), 4]),[10 1 1 1]);

%%
cd ~/Data/bruce/2_27_12/
save spktrg_avg_energy_reye_flipy spktrg_avg_cor spktrg_avg ov_avg stimres 

%%
close all
i = 2;
for j = 1:4
subplot(2,2,j)
imagesc(squeeze(spktrg_avg_cor(i,:,:,j)/6.5));
caxis([-0.5 0.5])
end
