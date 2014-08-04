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
dsfrac = 4;
Fsd = 1/Pix2Deg/dsfrac;

% desired temporal resolution
stimres = 0.025; %in s


spk_cnts = [];
mua_cnts = [];
all_sac_amps = [];
all_sac_vels = [];
all_sac_verspeeds = [];
all_sac_durs = [];
all_fix_start_times = [];
all_fix_stop_times = [];
all_stim_filtered = [];
all_sac_locs_left = [];
all_sac_locs_right = [];
all_stim_nums = [];
blockids = [];
%%
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
%         if ~exist([sname '.mat'],'file')
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
%             verspeed = max(abs(vervel),[],2); %take the larger of the hor and ver vergence speeds (Read and Cummin 2003).
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
    
    n_sacs = length(sac_start_times);
    
    %compute saccade amplitudes, velocities, and durations
    sac_stats(blockid).sac_dx = avg_eyepos(sac_stop_indsn,1) - avg_eyepos(sac_start_indsn,1);
    sac_stats(blockid).sac_dy = avg_eyepos(sac_stop_indsn,2) - avg_eyepos(sac_start_indsn,2);
    sac_stats(blockid).sac_amps = sqrt(sac_stats(blockid).sac_dx.^2+sac_stats(blockid).sac_dy.^2);
    sac_stats(blockid).sac_peakvel = zeros(n_sacs,1);
    sac_stats(blockid).sac_avgvel = zeros(n_sacs,1);
    sac_stats(blockid).sac_peakverspeed = zeros(n_sacs,1);
    sac_stats(blockid).sac_rstddev = zeros(n_sacs,1);
    sac_stats(blockid).sac_lstddev = zeros(n_sacs,1);
    sac_stats(blockid).sac_rpos = zeros(n_sacs,2);
    sac_stats(blockid).sac_lpos = zeros(n_sacs,2);
    sac_stats(blockid).sac_blinkdur = zeros(n_sacs,1);
    for i = 1:n_sacs
        sac_stats(blockid).sac_peakvel(i) = max(eye_speed(sac_start_indsn(i):sac_stop_indsn(i)));
        sac_stats(blockid).sac_peakverspeed(i) = max(verspeed(sac_start_indsn(i):sac_stop_indsn(i)));
        sac_stats(blockid).sac_avgvel(i) = mean(eye_speed(sac_start_indsn(i):sac_stop_indsn(i)));
        sac_stats(blockid).sac_rpos(i,:) = mean(reye_pos(sac_start_indsn(i):sac_stop_indsn(i),:));
        sac_stats(blockid).sac_lpos(i,:) = mean(leye_pos(sac_start_indsn(i):sac_stop_indsn(i),:));
        sac_stats(blockid).sac_lstddev(i) = sqrt(sum(var(leye_pos(sac_start_indsn(i):sac_stop_indsn(i),:))));
        sac_stats(blockid).sac_rstddev(i) = sqrt(sum(var(reye_pos(sac_start_indsn(i):sac_stop_indsn(i),:))));
        if sac_stats(blockid).sac_peakverspeed(i) > blink_thresh
           [~,ploc] = max(eye_speed(sac_start_indsn(i):sac_stop_indsn(i)));
           ploc = ploc + sac_start_indsn(i);
           prev_blink = find(verspeed(1:ploc) < blink_thresh,1,'last');
           next_blink = find(verspeed(ploc:end) < blink_thresh,1,'first')+ploc;
           sac_stats(blockid).sac_blinkdur(i) = (next_blink-prev_blink)*Eyedt;
        end
    end
    sac_stats(blockid).sac_durs = (sac_stop_times-sac_start_times)*1000; %in ms
    sac_stats(blockid).sac_start_time = sac_start_times;
    sac_stats(blockid).sac_stop_time = sac_stop_times;  
    
end


%%
close all
figure
for i = 1:2
subplot(3,2,i)
hold on
plot(eyets(sac_start_indsn),verspeed(sac_start_indsn),'r.')
plot(eyets(sac_start_indsn(temp)),verspeed(sac_start_indsn(temp)),'k*')
plot(eyets(1:end-1),abs(vervel(:,1)),'r')
plot(eyets(1:end-1),abs(vervel(:,2)),'k')
plot(eyets(1:end-1),sqrt(smooth(vervel(:,1).^2,20)),'b','linewidth',1)
plot(eyets(1:end-1),sqrt(smooth(vervel(:,2).^2,20)),'g','linewidth',1)
ylim([0 10])
end
subplot(3,2,3)
plot(eyets(1:end-1),reye_pos(:,1),'k')
hold on
plot(eyets(1:end-1),leye_pos(:,1),'r')
subplot(3,2,4)
plot(eyets(1:end-1),reye_pos(:,1)-leye_pos(:,1))
subplot(3,2,5)
plot(eyets(1:end-1),reye_pos(:,2),'k')
hold on
plot(eyets(1:end-1),leye_pos(:,2),'r')
subplot(3,2,6)
plot(eyets(1:end-1),reye_pos(:,2)-leye_pos(:,2))
%%
xl = [434 437];
% subplot(3,2,[1 2])
% xlim(xl)
for i = 1:6
    subplot(3,2,i)
    xlim(xl)
end
shg


%%
cd /Users/James/Data/bruce/2_27_12/stimrecon
save saccade_data sac_stats