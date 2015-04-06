% 
% close all
% 
% blocknum = 7;
% emtimes = Expts{blocknum}.Header.emtimes;
% nTrials = length(Expts{blocknum}.Trials);
% 
% for ii = 1:nTrials
%     subplot(2,1,1); hold on
%     plot(emtimes/1e4,Expts{blocknum}.Trials(ii).EyeData(:,1));
%     plot(emtimes/1e4,Expts{blocknum}.Trials(ii).EyeData(:,2),'r');
%     plot(emtimes/1e4,Expts{blocknum}.Trials(ii).EyeData(:,5),'k');
%     xlim([0 4]);
%     %     ylim([-5 5]);
%     ylim([-1 1]);
%     
%     subplot(2,1,2); hold on
%     plot(emtimes/1e4,Expts{blocknum}.Trials(ii).EyeData(:,3));
%     plot(emtimes/1e4,Expts{blocknum}.Trials(ii).EyeData(:,4),'r');
%     plot(emtimes/1e4,Expts{blocknum}.Trials(ii).EyeData(:,6),'k');
%     xlim([0 4]);
%     % ylim([-1 1]);
%     ylim([-5 5]);
%     pause
%     clf
% end

clear all
close all
cd ~/Data/bruce/misc/
load ThreeEyeSaccades
%%
for blocknum = 4:11;
emtimes = Expts{blocknum}.Header.emtimes;
nTrials = length(Expts{blocknum}.Trials);

Fa = Expts{blocknum}.Stimvals.Fa
all_Fas(blocknum) = Fa;

all_Edata = cat(3,Expts{blocknum}.Trials(:).EyeData);
all_Edata = permute(all_Edata,[1 3 2]);
all_Edata = reshape(all_Edata,[],6);
bad_pts = find(isnan(all_Edata(:,1)));
all_Edata(bad_pts,:) = 0;

%subtract off median 
corrected_eye_vals = bsxfun(@minus,all_Edata,nanmedian(all_Edata));

%rotate data into bar-oriented coordinate system
Rmat = [cosd(Fa-90) sind(Fa-90); -sind(Fa-90) cosd(Fa-90)];
corrected_eye_vals(:,[1 3]) = corrected_eye_vals(:,[1 3])*Rmat';
corrected_eye_vals(:,[2 4]) = corrected_eye_vals(:,[2 4])*Rmat';
corrected_eye_vals(:,[5 6]) = corrected_eye_vals(:,[5 6])*Rmat';

vev = std(corrected_eye_vals(:,1:4));
if vev(3) > vev(1)
    nv = corrected_eye_vals;
    nv(:,[3 4 6]) = corrected_eye_vals(:,[1 2 5]);
    nv(:,[1 2 5]) = corrected_eye_vals(:,[3 4 6]);
    corrected_eye_vals = nv;
end

%positions
b = robustfit(corrected_eye_vals(:,1),corrected_eye_vals(:,3));
pred_eyevals = corrected_eye_vals(:,1)*b(2) + b(1);
corrected_eye_vals(:,3) = corrected_eye_vals(:,3) - pred_eyevals;

b = robustfit(corrected_eye_vals(:,2),corrected_eye_vals(:,4));
pred_eyevals = corrected_eye_vals(:,2)*b(2) + b(1);
corrected_eye_vals(:,4) = corrected_eye_vals(:,4) - pred_eyevals;

b = robustfit(corrected_eye_vals(:,5),corrected_eye_vals(:,6));
pred_eyevals = corrected_eye_vals(:,5)*b(2) + b(1);
corrected_eye_vals(:,6) = corrected_eye_vals(:,6) - pred_eyevals;

%%
vt_diff = abs([0; diff(corrected_eye_vals(:,6))]);
bad_pts = find(vt_diff > 1 | abs(corrected_eye_vals(:,6)) > 1);
good_pts = setdiff(1:length(vt_diff),bad_pts);
new_vt_orth = interp1(good_pts,corrected_eye_vals(good_pts,6),1:length(vt_diff));
corrected_eye_vals(:,6) = new_vt_orth;
corrected_eye_vals(isnan(corrected_eye_vals)) = 0;
%%
eye_dt = Expts{blocknum}.Header.CRrates(1);
eye_fs = 1/eye_dt;
lEyeXY = all_Edata(:,[1 3]);
rEyeXY = all_Edata(:,[2 4]);

eye_vel = lEyeXY; %initialization
%uses smoothed left eye coil signal to define instantaneous speed
eye_smooth = 3;
sm_avg_eyepos = lEyeXY;
% sm_avg_eyepos(:,1) = smooth(lEyeXY(:,1),eye_smooth);
% sm_avg_eyepos(:,2) = smooth(lEyeXY(:,2),eye_smooth);
eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/eye_dt;
eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/eye_dt;

eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);

eye_ts = (1:length(eye_speed))/eye_fs;
et_params.eye_fs = eye_fs;
et_params.eye_smooth = eye_smooth;
et_params.use_coils = [1 1];
[saccades,et_params] = detect_saccades_v2(corrected_eye_vals(:,[1 3 2 4]),all_Edata(:,[1 3 2 4]),eye_speed,eye_ts(:),et_params);
is_blink = detect_blinks(eye_ts(:),corrected_eye_vals(:,[1 3 2 4]),saccades,et_params);
[saccades,is_blink] = merge_blinks(saccades,is_blink);
%start and stop positions of each saccade
sac_prepos = reshape([saccades(:).pre_pos],[],length(saccades));
sac_postpos = reshape([saccades(:).post_pos],[],length(saccades));
%saccade amplitude along parallel axis
sac_deltaX = sac_postpos(1,:) - sac_prepos(1,:);

sac_start_inds = round(interp1(eye_ts,1:length(eye_ts),[saccades(:).start_time]));
sac_stop_inds = round(interp1(eye_ts,1:length(eye_ts),[saccades(:).stop_time]));

%saccade amplitude along orthogonal axis
sac_deltaY = sac_postpos(2,:) - sac_prepos(2,:);

gsacs = find(abs(sac_deltaX) > 2 & abs(sac_deltaX) < 5 & ~is_blink');

Yright_gsacs = gsacs(find(sac_deltaY(gsacs) > 0));
Yleft_gsacs = gsacs(find(sac_deltaY(gsacs) < 0));
Xright_gsacs = gsacs(find(sac_deltaX(gsacs) > 0));
Xleft_gsacs = gsacs(find(sac_deltaX(gsacs) < 0));

pre_inds = sac_start_inds;
pre_inds(pre_inds < 1) = 1;
post_inds = sac_stop_inds;
post_inds(post_inds > length(eye_ts)) = length(eye_ts);
sac_full_prepos = corrected_eye_vals(pre_inds,:);
sac_full_postpos = corrected_eye_vals(post_inds,:);

clear sac_full_deltaY
sac_full_deltaY(:,1) = sac_full_postpos(:,3) - sac_full_prepos(:,3);
sac_full_deltaY(:,2) = sac_full_postpos(:,4) - sac_full_prepos(:,4);
sac_full_deltaY(:,3) = sac_full_postpos(:,6) - sac_full_prepos(:,6);

%%
sm_eyepos = corrected_eye_vals;
for ii = 1:6
sm_eyepos(:,ii) = smooth(corrected_eye_vals(:,ii),3);
end
eyevel = ([zeros(1,6); diff(sm_eyepos)])/eye_dt;

sac_backlag = round(0.05*eye_fs);
sac_forlag = round(0.1*eye_fs);
% [trig_avg_eyepos_right,lags] = get_event_trig_avg_v3(corrected_eye_vals,sac_start_inds(right_gsacs),sac_backlag,sac_forlag);
% [trig_avg_eyepos_left,lags] = get_event_trig_avg_v3(corrected_eye_vals,sac_start_inds(left_gsacs),sac_backlag,sac_forlag);

[trig_avg_eyevel_right,lags] = get_event_trig_avg_v3(eyevel,sac_start_inds(Yright_gsacs),sac_backlag,sac_forlag);
[trig_avg_eyevel_left,lags] = get_event_trig_avg_v3(eyevel,sac_start_inds(Yleft_gsacs),sac_backlag,sac_forlag);

avg_trig_orthvel_Yright(blocknum,:,:) =(trig_avg_eyevel_right);
avg_trig_orthvel_Yleft(blocknum,:,:) = (trig_avg_eyevel_left);

[trig_avg_eyevel_right,lags] = get_event_trig_avg_v3(eyevel,sac_start_inds(Xright_gsacs),sac_backlag,sac_forlag);
[trig_avg_eyevel_left,lags] = get_event_trig_avg_v3(eyevel,sac_start_inds(Xleft_gsacs),sac_backlag,sac_forlag);

avg_trig_orthvel_Xright(blocknum,:,:) = trig_avg_eyevel_right;
avg_trig_orthvel_Xleft(blocknum,:,:) = trig_avg_eyevel_left;

end
%%

for blocknum = 4:11
  Fa = Expts{blocknum}.Stimvals.Fa
  close all
    subplot(2,2,1); hold on
    plot(lags,squeeze(avg_trig_orthvel_Yright(blocknum,:,1)))
    plot(lags,squeeze(avg_trig_orthvel_Yright(blocknum,:,2)),'r')
        plot(lags,squeeze(avg_trig_orthvel_Yright(blocknum,:,5)),'k')
    plot(lags,squeeze(avg_trig_orthvel_Yleft(blocknum,:,1)),'--')
    plot(lags,squeeze(avg_trig_orthvel_Yleft(blocknum,:,2)),'r--')
        plot(lags,squeeze(avg_trig_orthvel_Yleft(blocknum,:,5)),'k--')
        
        
    subplot(2,2,2); hold on
    plot(lags,squeeze(avg_trig_orthvel_Yright(blocknum,:,3)))
    plot(lags,squeeze(avg_trig_orthvel_Yright(blocknum,:,4)),'r')
        plot(lags,squeeze(avg_trig_orthvel_Yright(blocknum,:,6)),'k')
    plot(lags,squeeze(avg_trig_orthvel_Yleft(blocknum,:,3)),'--')
    plot(lags,squeeze(avg_trig_orthvel_Yleft(blocknum,:,4)),'r--')
        plot(lags,squeeze(avg_trig_orthvel_Yleft(blocknum,:,6)),'k--')
    subplot(2,2,3); hold on
    plot(lags,squeeze(avg_trig_orthvel_Xright(blocknum,:,1)))
    plot(lags,squeeze(avg_trig_orthvel_Xright(blocknum,:,2)),'r')
        plot(lags,squeeze(avg_trig_orthvel_Xright(blocknum,:,5)),'k')
    plot(lags,squeeze(avg_trig_orthvel_Xleft(blocknum,:,1)),'--')
    plot(lags,squeeze(avg_trig_orthvel_Xleft(blocknum,:,2)),'r--')
        plot(lags,squeeze(avg_trig_orthvel_Xleft(blocknum,:,5)),'k--')
    subplot(2,2,4); hold on
    plot(lags,squeeze(avg_trig_orthvel_Xright(blocknum,:,3)))
    plot(lags,squeeze(avg_trig_orthvel_Xright(blocknum,:,4)),'r')
        plot(lags,squeeze(avg_trig_orthvel_Xright(blocknum,:,6)),'k')
    plot(lags,squeeze(avg_trig_orthvel_Xleft(blocknum,:,3)),'--')
    plot(lags,squeeze(avg_trig_orthvel_Xleft(blocknum,:,4)),'r--')
        plot(lags,squeeze(avg_trig_orthvel_Xleft(blocknum,:,6)),'k--')
    pause
    clf
end

% emtimes = Expts{blocknum}.Header.emtimes;
% nTrials = length(Expts{blocknum}.Trials);
% 
% close all
% for ii = 1:nTrials
%     cur_inds = (1:length(emtimes)) + (ii-1)*length(emtimes);
%     subplot(2,1,1); hold on
%     plot(emtimes/1e4,corrected_eye_vals(cur_inds,1));
%     plot(emtimes/1e4,corrected_eye_vals(cur_inds,2),'r');
%     plot(emtimes/1e4,corrected_eye_vals(cur_inds,5),'k');
%     xlim([0 4]);
%         ylim([-5 5]);
% %     ylim([-1 1]);
%     
%     subplot(2,1,2); hold on
%     plot(emtimes/1e4,corrected_eye_vals(cur_inds,3));
%     plot(emtimes/1e4,corrected_eye_vals(cur_inds,4),'r');
%     plot(emtimes/1e4,corrected_eye_vals(cur_inds,6),'k');
%     xlim([0 4]);
%     ylim([-1 1]);
% %     ylim([-5 5]);
%     pause
%     clf
% end
