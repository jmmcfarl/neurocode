function [sac_data,in_sac,eye_speed] = get_saccades(rEyeXY,lEyeXY,eye_ts,in_blink)

%%
sac_eyespeed = 10;
thresh_eyespeed = 4;
% sac_eyespeed = 5;
% thresh_eyespeed = 2;

eye_dt = median(diff(eye_ts));
eye_fs = 1/eye_dt;
min_isi = 0.05;
% min_isi = 0.01;

%%
avg_eyepos = 0.5*lEyeXY + 0.5*lEyeXY;
sm_avg_eyepos(:,1) = smooth(avg_eyepos(:,1),3);
sm_avg_eyepos(:,2) = smooth(avg_eyepos(:,2),3);
eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/eye_dt;
eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/eye_dt;


eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
eye_speed(in_blink==1) = 0;

%%
sac_data = struct('start_times',[],'stop_times',[],'amplitude',[],'dx',[],'dy',[],'dx_left',[],'dy_left',[],'dx_right',[],'dy_right',[]);
in_sac = zeros(size(eye_ts));
if max(eye_speed) > sac_eyespeed

    proposed_sac_starts = find(eye_speed(1:end-1) < thresh_eyespeed & eye_speed(2:end) >= thresh_eyespeed);
    proposed_sac_stops = find(eye_speed(1:end-1) >= thresh_eyespeed & eye_speed(2:end) < thresh_eyespeed);
    if eye_speed(1) >= thresh_eyespeed
        proposed_sac_starts = [1; proposed_sac_starts];
    end
    if eye_speed(end) >= thresh_eyespeed
        proposed_sac_stops = [proposed_sac_stops; length(eye_speed)];
    end
    fail_sacs = [];
    for i = 1:length(proposed_sac_starts)
        if max(eye_speed(proposed_sac_starts(i):proposed_sac_stops(i))) < sac_eyespeed
            fail_sacs = [fail_sacs i];
        end
    end
    proposed_sac_starts(fail_sacs) = [];
    proposed_sac_stops(fail_sacs) = [];
    
    n_sacs = length(proposed_sac_starts);
    for i = 1:n_sacs
        in_sac(proposed_sac_starts(i):proposed_sac_stops(i)) = 1;
    end
    
    sac_start_times = eye_ts(proposed_sac_starts);
    sac_stop_times = eye_ts(proposed_sac_stops);
    
    isi = diff(sac_start_times);
    too_short = find(isi < min_isi);
    proposed_sac_starts(too_short) = [];
    proposed_sac_stops(too_short) = [];
    
    n_sacs = length(proposed_sac_starts);
    sac_start_times = eye_ts(proposed_sac_starts);
    sac_stop_times = eye_ts(proposed_sac_stops);
    sac_amps_dx= avg_eyepos(proposed_sac_stops,1) - avg_eyepos(proposed_sac_starts,1);
    sac_amps_dy = avg_eyepos(proposed_sac_stops,2) - avg_eyepos(proposed_sac_starts,2);
    sac_amps_ldx= lEyeXY(proposed_sac_stops,1) - lEyeXY(proposed_sac_starts,1);
    sac_amps_ldy = lEyeXY(proposed_sac_stops,2) - lEyeXY(proposed_sac_starts,2);
    sac_amps_rdx= rEyeXY(proposed_sac_stops,1) - rEyeXY(proposed_sac_starts,1);
    sac_amps_rdy = rEyeXY(proposed_sac_stops,2) - rEyeXY(proposed_sac_starts,2);

    sac_amps = sqrt(sac_amps_dx.^2 + sac_amps_dy.^2);
    if n_sacs > 0
        sac_data = struct('start_inds',proposed_sac_starts,'stop_inds',proposed_sac_stops,...
            'start_time',mat2cell(sac_start_times(:),ones(n_sacs,1)),...
            'stop_time',mat2cell(sac_stop_times(:),ones(n_sacs,1)),'amplitude',mat2cell(sac_amps(:),ones(n_sacs,1)),...
            'dx',mat2cell(sac_amps_dx(:),ones(n_sacs,1)),'dy',mat2cell(sac_amps_dy(:),ones(n_sacs,1)),...
            'dx_left',mat2cell(sac_amps_ldx(:),ones(n_sacs,1)), 'dy_left',mat2cell(sac_amps_ldy(:),ones(n_sacs,1)),...
            'dx_right',mat2cell(sac_amps_rdx(:),ones(n_sacs,1)), 'dy_right',mat2cell(sac_amps_rdy(:),ones(n_sacs,1)));
    end
else
    n_sacs = 0;
end

% fprintf('%d sacs located\n',n_sacs);
    

%% compute saccade amplitudes, velocities, and durations
%     sac_stats(blockid).sac_dx = avg_eyepos(sac_stop_indsn,1) - avg_eyepos(sac_start_indsn,1);
%     sac_stats(blockid).sac_dy = avg_eyepos(sac_stop_indsn,2) - avg_eyepos(sac_start_indsn,2);
%     sac_stats(blockid).sac_amps = sqrt(sac_stats(blockid).sac_dx.^2+sac_stats(blockid).sac_dy.^2);
%     sac_stats(blockid).sac_peakvel = zeros(n_sacs,1);
%     sac_stats(blockid).sac_avgvel = zeros(n_sacs,1);
%     sac_stats(blockid).sac_peakverspeed = zeros(n_sacs,1);
%     sac_stats(blockid).sac_rstddev = zeros(n_sacs,1);
%     sac_stats(blockid).sac_lstddev = zeros(n_sacs,1);
%     sac_stats(blockid).sac_rpos = zeros(n_sacs,2);
%     sac_stats(blockid).sac_lpos = zeros(n_sacs,2);
%     sac_stats(blockid).sac_sacdur = zeros(n_sacs,1);
%     for i = 1:n_sacs
%         sac_stats(blockid).sac_peakvel(i) = max(eye_speed(sac_start_indsn(i):sac_stop_indsn(i)));
%         sac_stats(blockid).sac_peakverspeed(i) = max(verspeed(sac_start_indsn(i):sac_stop_indsn(i)));
%         sac_stats(blockid).sac_avgvel(i) = mean(eye_speed(sac_start_indsn(i):sac_stop_indsn(i)));
%         sac_stats(blockid).sac_rpos(i,:) = mean(reye_pos(sac_start_indsn(i):sac_stop_indsn(i),:));
%         sac_stats(blockid).sac_lpos(i,:) = mean(leye_pos(sac_start_indsn(i):sac_stop_indsn(i),:));
%         sac_stats(blockid).sac_lstddev(i) = sqrt(sum(var(leye_pos(sac_start_indsn(i):sac_stop_indsn(i),:))));
%         sac_stats(blockid).sac_rstddev(i) = sqrt(sum(var(reye_pos(sac_start_indsn(i):sac_stop_indsn(i),:))));
%         if sac_stats(blockid).sac_peakverspeed(i) > sac_thresh
%            [~,ploc] = max(eye_speed(sac_start_indsn(i):sac_stop_indsn(i)));
%            ploc = ploc + sac_start_indsn(i);
%            prev_sac = find(verspeed(1:ploc) < sac_thresh,1,'last');
%            next_sac = find(verspeed(ploc:end) < sac_thresh,1,'first')+ploc;
%            sac_stats(blockid).sac_sacdur(i) = (next_sac-prev_sac)*Eyedt;
%         end
%     end
