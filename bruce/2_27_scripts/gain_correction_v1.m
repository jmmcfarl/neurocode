clear all
% close all

cd ~/Data/bruce/2_27_12/
load Blocks.mat
accept_window = [-7 7;-7 7];
sac_eyespeed = 10;
thresh_eyespeed = 2.5;
sac_buffer = 0;
min_fix_dur = 0.2;
blink_thresh = 5;
min_blink_dur = 0.05;
% nlags = 4;

Pix2Deg = 0.018837;
% down-sampling fraction for image
dsfrac = 4;
Fsd = 1/Pix2Deg/dsfrac;

%%
ov_reye_pos = [];
ov_leye_pos = [];
for blockid = 1:4
    
    fprintf('Processing block %d of %d...\n',blockid,4);    
    cd ~/Data/bruce/2_27_12/saccades/
    load(sprintf('lemM232.5%d.em.sac.mat',blockid))
         EyeStartT = Expt.Trials.Start/10000; % time of first eye sample
    EyeEndT = Expt.Trials.End/10000; % time of last eye sample
    Eyedt = Expt.Header.CRrates(1); % resolution of eye signal (sec)
    eyets = EyeStartT:Eyedt:EyeEndT; %eye tracking time axis (sec)
    sac_buffer_inds = round(sac_buffer/Eyedt);
   
    reye_pos = [Expt.Trials.Eyevals.rh Expt.Trials.Eyevals.rv];
    leye_pos = [Expt.Trials.Eyevals.lh Expt.Trials.Eyevals.lv];
  
    fprintf('Computing saccade times\n');
    avg_eyepos = (reye_pos + leye_pos)/2;
    clear sm_avg_eyepos eye_vel
    sm_avg_eyepos(:,1) = smooth(avg_eyepos(:,1),3);
    sm_avg_eyepos(:,2) = smooth(avg_eyepos(:,2),3);
    eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/Eyedt;
    eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/Eyedt;
    eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
        
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

    avg_pos = 0.5*reye_pos + 0.5*leye_pos;
    in_window = (avg_pos(:,1) >= accept_window(1,1) & avg_pos(:,1) <= accept_window(1,2) & ...
        avg_pos(:,2) >= accept_window(2,1) & avg_pos(:,2) <= accept_window(2,2));

    used_times = find(sac_vec == 0 & in_window == 1);
    
    ov_reye_pos = [ov_reye_pos; reye_pos(used_times,:)];
    ov_leye_pos = [ov_leye_pos; leye_pos(used_times,:)];
end

%%
disparity = ov_reye_pos - ov_leye_pos;
gaze_loc = 0.5*ov_reye_pos+0.5*ov_leye_pos;

med_disparity = median(disparity);
outliers = find(max(abs(bsxfun(@minus,disparity,med_disparity)),[],2) > 1);
used_pts = setdiff(1:size(disparity,1),outliers);

%%
leye_pos = ov_leye_pos(used_pts,:);
reye_pos = ov_reye_pos(used_pts,:);
disparity = reye_pos - leye_pos;

%get rid of overall disparity
reye_pos = bsxfun(@minus,reye_pos,mean(disparity));
disparity = reye_pos - leye_pos;

%%
K0 = [1 0 0 1 1 0 0 1];
options.Display = 'iter';
scale_pen = 0.5;
[Kbest] = fminunc(@(K) disparity_error(K,reye_pos,leye_pos,scale_pen), K0,options);

%%
cur_r(:,1) = Kbest(1)*reye_pos(:,1) + Kbest(2)*reye_pos(:,2);
cur_r(:,2) = Kbest(3)*reye_pos(:,1) + Kbest(4)*reye_pos(:,2);

cur_l(:,1) = Kbest(5)*leye_pos(:,1) + Kbest(6)*leye_pos(:,2);
cur_l(:,2) = Kbest(7)*leye_pos(:,1) + Kbest(8)*leye_pos(:,2);

cur_disparity = cur_r - cur_l;

%%

[rbefore,b] = corr(disparity, reye_pos);
[lbefore,b] = corr(disparity, leye_pos);

[rafter,b] = corr(cur_disparity, cur_r);
[lafter,b] = corr(cur_disparity, cur_l);








































