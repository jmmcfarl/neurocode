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

X = [];
spk_cnts = [];
all_sac_amps = [];
all_sac_vels = [];
all_sac_verspeeds = [];
all_sac_durs = [];
all_fix_start_times = [];
all_fix_stop_times = [];
all_stim_filtered = [];
all_sac_locs_left = [];
all_sac_locs_right = [];
blockids = [];
for blockid = 1:4
    
    fprintf('Processing block %d of %d...\n',blockid,4);
    
    block_times = Blocks{blockid}.blocktimes;
    stim_times = Blocks{blockid}.stimtime;
    stimID = Blocks{blockid}.stimids;
    
    n_stims = length(stim_times);
    
    cd ~/Data/bruce/2_27_12/saccades/
    load(sprintf('lemM232.5%d.em.sac.mat',blockid))
    
    % identify saccade start and stop times
    EyeStartT = Expt.Trials.Start/10000; % time of first eye sample
    EyeEndT = Expt.Trials.End/10000; % time of last eye sample
    Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
    eyets = EyeStartT:Eyedt:EyeEndT; %eye tracking time axis (sec)
    sac_buffer_inds = round(sac_buffer/Eyedt);
    
    reye_pos = [Expt.Trials.Eyevals.rh Expt.Trials.Eyevals.rv];
    leye_pos = [Expt.Trials.Eyevals.lh Expt.Trials.Eyevals.lv];
    
    %if we haven't already computed saccade times do so now
    sname = sprintf('sac_times_corf_block%d',blockid);
    %     if ~exist([sname '.mat'],'file')
    fprintf('Computing saccade times\n');
    avg_eyepos = (reye_pos + leye_pos)/2;
    clear sm_avg_eyepos eye_vel
    sm_avg_eyepos(:,1) = smooth(avg_eyepos(:,1),3);
    sm_avg_eyepos(:,2) = smooth(avg_eyepos(:,2),3);
    eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/Eyedt;
    eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/Eyedt;
    eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
    
    
    vergence = reye_pos - leye_pos;
    vervel = [0 0;diff(vergence)]/Eyedt;
    %         verspeed = abs(vervel);
    verspeed = sqrt(vervel(:,1).^2+vervel(:,2).^2);
    %         verspeed = max(vervel,[],2); %take the larger of the hor and ver vergence speeds (Read and Cummin 2003).
    blink_on = find(verspeed(2:end) > blink_thresh & verspeed(1:end-1) <= blink_thresh);
    blink_off = find(verspeed(2:end) <= blink_thresh & verspeed(1:end-1) > blink_thresh);
    blink_dur = eyets(blink_off) - eyets(blink_on);
    is_blink = find(blink_dur > min_blink_dur);
    blink_times = eyets(round((blink_on(is_blink)+blink_off(is_blink))/2));
    
    %find saccade start and stop indices
    sac_inds = find(eye_speed(1:end-1) < sac_eyespeed & eye_speed(2:end) > sac_eyespeed);
    
    sac_start_inds = nan(size(sac_inds));
    sac_stop_inds = nan(size(sac_inds));
    for i = 1:length(sac_inds)
        temp = find(eye_speed(1:sac_inds(i)) < thresh_eyespeed,1,'last');
        if ~isempty(temp)
            sac_start_inds(i) = temp;
        end
        temp = find(eye_speed(sac_inds(i)+1:end) < thresh_eyespeed,1,'first');
        if ~isempty(temp)
            sac_stop_inds(i) = sac_inds(i)+temp-1;
        end
    end
    
    %identify start and stop times of unique saccades
    sac_vec = zeros(size(reye_pos,1),1);
    for i = 1:length(sac_start_inds)
        if ~isnan(sac_start_inds(i)) & ~isnan(sac_stop_inds(i))
            sac_vec(sac_start_inds(i):sac_stop_inds(i)+sac_buffer_inds) = 1;
        end
    end
    sac_vec(length(eyets)+1:end) = [];
    sac_vec([1 end]) = 0;
    sac_start_indsn = find(sac_vec(1:end-1) == 0 & sac_vec(2:end) == 1);
    sac_stop_indsn = find(sac_vec(1:end-1) == 1 & sac_vec(2:end) == 0);
    if length(sac_start_indsn) ~= length(sac_stop_indsn)
        error('saccade mis-alignment');
    end
    
    sac_start_times = eyets(sac_start_indsn);
    sac_stop_times = eyets(sac_stop_indsn);
        
    %compute saccade amplitudes, velocities, and durations
    sac_dx = avg_eyepos(sac_stop_indsn,1) - avg_eyepos(sac_start_indsn,1);
    sac_dy = avg_eyepos(sac_stop_indsn,2) - avg_eyepos(sac_start_indsn,2);
    sac_rdx = reye_pos(sac_stop_indsn,1) - reye_pos(sac_start_indsn,1);
    sac_rdy = reye_pos(sac_stop_indsn,2) - reye_pos(sac_start_indsn,2);
    sac_ldx = leye_pos(sac_stop_indsn,1) - leye_pos(sac_start_indsn,1);
    sac_ldy = leye_pos(sac_stop_indsn,2) - leye_pos(sac_start_indsn,2);
    sac_amps = sqrt(sac_dx.^2+sac_dy.^2);
    sac_peakvel = zeros(size(sac_start_indsn));
    sac_peakverspeed = zeros(size(sac_start_indsn));
    for i = 1:length(sac_start_indsn)
        sac_peakvel(i) = max(eye_speed(sac_start_indsn(i):sac_stop_indsn(i)));
        sac_peakverspeed(i) = max(verspeed(sac_start_indsn(i):sac_stop_indsn(i)));
    end
    sac_durs = sac_stop_times-sac_start_times;
    
    cd ~/Data/bruce/2_27_12/saccades/
    save(sname,'sac_start_times','sac_stop_times','sac_start_inds','sac_stop_inds','sac_vec',...
        'sac_amps','sac_peakvel','sac_durs','blink_times','sac_peakverspeed')
    %     else
    %         load(sname)
    %     end
    
    % load patch video data
    cd ~/Data/bruce/2_27_12/stimrecon/
    %     sname = sprintf('image_patch_block%d_reye_50_dsf4',blockid);
    %     sname = sprintf('image_patch_block%d_leye_50_dsf4',blockid);
    sname = sprintf('image_patch_block%d_reye_25_dsf4_patch2',blockid);
    load(sname)
    
    % interpolate eye signal onto stimulus time axis
    eye_interp_right = interp1(eyets(1:end-1),reye_pos,recon_t);
    eye_interp_left = interp1(eyets(1:end-1),leye_pos,recon_t);
    %         eye_interp = interp1(eyets(1:end-1),leye_pos+reye_pos,recon_t)/2;
    
%     eye_interp = eye_interp_left;
    eye_interp = eye_interp_right;
    
    %create an interpolated 0-1 vector of saccade times
    sac_vec_interp = zeros(size(recon_t));
    for i = 1:length(sac_start_times)
        cur_set = find(recon_t >= sac_start_times(i) & recon_t <= sac_stop_times(i));
        sac_vec_interp(cur_set) = 1;
    end
    sac_vec_interp([1 end]) = 0;
    
    %indices of the start and stops of fixations
    fix_start_inds = 1+[0 find(sac_vec_interp(1:end-1) == 1 & sac_vec_interp(2:end) == 0)];
    fix_stop_inds = [find(sac_vec_interp(1:end-1) == 0 & sac_vec_interp(2:end) == 1) length(recon_t)-1];
    if length(fix_start_inds) ~= length(fix_stop_inds)
        error('Fixation mis-alignment');
    end
    fix_durs = (fix_stop_inds - fix_start_inds)*stimres;
    too_short = find(fix_durs < min_fix_dur); %only keep fixations that are minimum duration
    
    fprintf('%d of %d fixations too short\n',length(too_short),length(fix_durs));
    fix_start_inds(too_short) = []; fix_stop_inds(too_short) = []; fix_durs(too_short) = [];
    fix_start_times = recon_t(fix_start_inds);
    
    fix_dx = eye_interp(fix_stop_inds,1) - eye_interp(fix_start_inds,1);
    fix_dy = eye_interp(fix_stop_inds,2) - eye_interp(fix_start_inds,2);
    drift_amps = sqrt(fix_dx.^2+fix_dy.^2);
    
    %determine correspondence between fixation starts and previous saccades
    fix_prev_sac_inds = zeros(size(fix_start_inds));
    for i = 1:length(fix_start_inds)
        [~,fix_prev_sac_inds(i)] = min(abs(sac_stop_times-fix_start_times(i)));
    end
    
    % create vector determining whether current stimulus was high-pass filtered
    filt_stims = mod(0:500,4) + 1 >2;
    is_stim_filtered = filt_stims(stimID(stim_num));
    
    %find times where we skip stimulus presentation data
    %also find blinks
    bad_times = [0 diff(recon_t)] > 0.05;
    fix_skips = zeros(size(fix_start_inds));
    blinks = zeros(size(fix_start_inds));
    for i = 1:length(fix_start_inds)
        if any(bad_times(fix_start_inds(i):fix_stop_inds(i)))
            fix_skips(i) = 1;
        end
        if min(abs(blink_times - fix_start_times(i))) < 0.1
            blinks(i) = 1;
        end
    end
    
    % find times where eye-position is within central window
    in_window_left = (eye_interp_left(:,1) >= accept_window(1,1) & eye_interp_left(:,1) <= accept_window(1,2) & ...
        eye_interp_left(:,2) >= accept_window(2,1) & eye_interp_left(:,2) <= accept_window(2,2));
    in_window_right = (eye_interp_right(:,1) >= accept_window(1,1) & eye_interp_right(:,1) <= accept_window(1,2) & ...
        eye_interp_right(:,2) >= accept_window(2,1) & eye_interp_right(:,2) <= accept_window(2,2));
    in_window = in_window_left==1 & in_window_right==1;
    off_screen = max(any(isnan(ov_im_patch),3),[],2);
    
    
    used_fix_start_inds = find(in_window(fix_start_inds) & off_screen(fix_start_inds) == 0 & fix_skips'==0);
    %     used_fix_start_inds = find(in_window(fix_start_inds)==1 & is_stim_filtered(fix_start_inds)' == 0);
    n_used_fixs = length(used_fix_start_inds);
    fprintf('Out of window: %d,  off-screen:  %d, gray-stim: %d,  blinks: %d\n',sum(in_window(fix_start_inds)==0),sum(off_screen(fix_start_inds)==1),sum(fix_skips==1),sum(blinks==1));
    fprintf('analyzing %d of %d fixations\n',n_used_fixs,length(fix_start_inds));
    
    
    %compute binned spike count vectors at stimulus resolution
    spike_bin_vecs = [];
    for t = 1:n_stims
        onsetT = stim_times(t); %stimulus onset time (s)
        if t < n_stims %for the last image, take the end time as the block end time
            endT = stim_times(t+1)-stimres; %set end time as one time bin before the following stimulus presentation
            endT = endT - 0.2; %last 200ms doesn't have any image presentation!
        else
            endT = Blocks{blockid}.blocktimes(2,end);
        end
        Tax = onsetT:stimres:endT;
        time_bin_edges = [(Tax(1)-stimres/2) (Tax+stimres/2)];
        cur_spikes_binned = zeros(10,length(Tax));
        for cellid = 1:10
            spikes_binned = histc(Blocks{blockid}.spktimes{cellid},time_bin_edges);
            spikes_binned(end) = [];
            cur_spikes_binned(cellid,:) = spikes_binned;
        end
        spike_bin_vecs = [spike_bin_vecs; cur_spikes_binned'];
    end
    spike_bin_vecs = spike_bin_vecs';
    
    cur_fix_image_set = zeros(n_used_fixs,size(ov_im_patch,2),size(ov_im_patch,3));
    cur_fix_spkcnt_set = zeros(n_used_fixs,10);
    cur_fix_loc_left = zeros(n_used_fixs,2);
    cur_fix_loc_right = zeros(n_used_fixs,2);
    for i = 1:n_used_fixs
        %        cur_set = fix_start_inds(used_fix_start_inds(i)):fix_stop_inds(used_fix_start_inds(i)); %use entire fixation
        cur_set = fix_start_inds(used_fix_start_inds(i)) : (fix_start_inds(used_fix_start_inds(i))+round(min_fix_dur/stimres)); %use only window at begenning of fixation
        cur_avg_image = nanmean(ov_im_patch(cur_set,:,:));
        cur_fix_loc_right(i,:) = mean(eye_interp_right(cur_set,:));
        cur_fix_loc_left(i,:) = mean(eye_interp_left(cur_set,:));
        
        cur_fix_image_set(i,:,:) = cur_avg_image - nanmean(cur_avg_image(:)); %within fixation mean subtraction
        
        cur_fix_spkcnt_set(i,:) = sum(spike_bin_vecs(:,cur_set),2); %total spike counts within fixation window
    end
    
    vert_disp = cur_fix_loc_right(:,2) - cur_fix_loc_left(:,2);
    hor_disp = cur_fix_loc_right(:,1) - cur_fix_loc_left(:,1);
    vert_disp = vert_disp - median(vert_disp);
    hor_disp = hor_disp - median(hor_disp);
    good_fixs = find(abs(vert_disp) < 1 & abs(hor_disp) < 1);
    
    
    X = [X; cur_fix_image_set(good_fixs,:,:)];
    spk_cnts = [spk_cnts; cur_fix_spkcnt_set(good_fixs,:)];
    all_sac_locs_left = [all_sac_locs_left; cur_fix_loc_left(good_fixs,:)];
    all_sac_locs_right = [all_sac_locs_right; cur_fix_loc_right(good_fixs,:)];
    all_sac_amps = [all_sac_amps; sac_amps(fix_prev_sac_inds(used_fix_start_inds(good_fixs)))];
    all_sac_durs = [all_sac_durs; sac_durs(fix_prev_sac_inds(used_fix_start_inds(good_fixs)))'];
    all_sac_vels = [all_sac_vels; sac_peakvel(fix_prev_sac_inds(used_fix_start_inds(good_fixs)))];
    all_sac_verspeeds = [all_sac_verspeeds; sac_peakverspeed(fix_prev_sac_inds(used_fix_start_inds(good_fixs)))];
    all_fix_start_times = [all_fix_start_times; sac_stop_times(fix_prev_sac_inds(used_fix_start_inds(good_fixs)))'];
    all_fix_stop_times = [all_fix_stop_times; sac_start_times(fix_prev_sac_inds(used_fix_start_inds(good_fixs(1:end-1)))+1)'];
    all_fix_stop_times = [all_fix_stop_times; recon_t(end)];
    blockids = [blockids; ones(length(good_fixs),1)*blockid];
    all_stim_filtered = [all_stim_filtered; is_stim_filtered(fix_start_inds(used_fix_start_inds(good_fixs)))'];
end

%% convert spike data to bin indices
% spikebins = cell(10,1);
% for cellid = 1:10
%     un_spk_cnts = unique(spk_cnts(:,cellid));
%     cur_spikebins = [];
%     for i = 1:length(un_spk_cnts)
%         cur_set = find(spk_cnts(:,cellid) == un_spk_cnts(i));
%         cur_spikebins = [cur_spikebins; repmat(cur_set(:),un_spk_cnts(i),1)];
%     end
%     cur_spikebins = sort(cur_spikebins);
%     spikebins{cellid} = [spikebins{cellid}; cur_spikebins ];
% end

%%
cd /Users/James/Data/bruce/2_27_12/stimrecon
save fixation_stim_righteye_patch2 all_* block* spk_* X 