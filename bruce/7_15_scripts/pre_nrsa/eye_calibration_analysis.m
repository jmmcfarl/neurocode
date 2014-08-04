clear all
close all
cd ~/Data/bruce/7_15_12/G029
load ./jbeE005.eyecal.mat
%%
is_good = [Expt.Trials(:).good];
fx = [Expt.Trials(:).fx];
fy = [Expt.Trials(:).fy];
Trial_starts = [Expt.Trials(:).Start]/1e4;
Trial_ends = [Expt.Trials(:).End]/1e4;
good_trials = find(is_good==1);
eye_dt = Expt.Header.CRrates(1);
blink_thresh = 0.5;
blink_dur = 0.5;
sac_eyespeed = 10;
thresh_eyespeed = 2.5;
max_error = 1.5;
eye_fs = 1/eye_dt;
[b,a] = butter(2,[.5 10]/(eye_fs/2));

n_trials = length(fx);
lEye_avg_pos = nan(n_trials,2);
rEye_avg_pos = nan(n_trials,2);

all_tnum = [];
all_inblink = [];
all_insac = [];
all_infixation = [];
all_lEye = [];
all_rEye = [];
left_sac_dx = [];
left_sac_dy = [];
right_sac_dx = [];
right_sac_dy = [];
for i = 1:length(good_trials)
    
    lEyeXY = [Expt.Trials(good_trials(i)).Eyevals.lh Expt.Trials(good_trials(i)).Eyevals.lv];
    rEyeXY = [Expt.Trials(good_trials(i)).Eyevals.rh Expt.Trials(good_trials(i)).Eyevals.rv];
        
    trial_n = size(lEyeXY,1);
    eye_ts = Trial_starts(good_trials(i)):eye_dt:(Trial_starts(good_trials(i)) + eye_dt*(trial_n - 1));
    
    [blink_data,in_blink,tot_disp_f] = get_blinks(rEyeXY,lEyeXY,eye_ts);
    [sac_data,in_sac,eye_speed] = get_saccades(rEyeXY,lEyeXY,eye_ts,in_blink);
   
    all_insac = [all_insac in_sac];
    all_inblink = [all_inblink in_blink];
    
    [fixation_data,in_fixation] = parse_fixations(in_blink,in_sac,eye_ts);
    cur_lmean = mean(lEyeXY(in_fixation,:));
    cur_rmean = mean(rEyeXY(in_fixation,:));
    for j = 1:length(fixation_data)
        cur_inds = fixation_data(j).start_inds:fixation_data(j).stop_inds;
        cur_fix_lerror = cur_lmean - mean(lEyeXY(cur_inds,:));
        cur_fix_lerror = norm(cur_fix_lerror);
         cur_fix_rerror = cur_rmean - mean(rEyeXY(cur_inds,:));
        cur_fix_rerror = norm(cur_fix_rerror);
        if cur_fix_lerror > max_error | cur_fix_rerror > max_error
            in_fixation(cur_inds) = 0;
        end
    end
        
    lEye_avg_pos(good_trials(i),:) = mean(lEyeXY(in_fixation,:));
    rEye_avg_pos(good_trials(i),:) = mean(rEyeXY(in_fixation,:));
    
    all_lEye = [all_lEye; lEyeXY];
    all_rEye = [all_rEye; rEyeXY];
    
    left_sac_dx = [left_sac_dx sac_data(:).dx_left];
    left_sac_dy = [left_sac_dx sac_data(:).dy_left];
    right_sac_dx = [left_sac_dx sac_data(:).dx_right];
    right_sac_dy = [left_sac_dx sac_data(:).dy_right];
    
    all_infixation = [all_infixation in_fixation];
    all_tnum = [all_tnum good_trials(i)*ones(size(in_fixation))];
    
% %     figure
%     subplot(2,1,1)
%     plot(lEyeXY(:,1))
%     hold on
%     plot(rEyeXY(:,1),'r')
%     plot(find(in_fixation==1),lEyeXY(in_fixation==1,1),'c.')
%     xl = xlim();
%     line(xl,-[fx(good_trials(i)) fx(good_trials(i))],'color','k')
%     ylim([-10 10])
%     subplot(2,1,2)
%     plot(lEyeXY(:,2))
%     hold on
%     plot(rEyeXY(:,2),'r')
%     plot(find(in_fixation==1),lEyeXY(in_fixation==1,2),'c.')
%     xl = xlim();
%     line(xl,[fy(good_trials(i)) fy(good_trials(i))],'color','k')
%     ylim([-10 10])
%      
%     pause
%     clf
    
end

%%

%% FOR LEFT EYE
Y = [fx' fy'];
X = lEye_avg_pos(good_trials,:);
Xq = [lEye_avg_pos(good_trials,:) lEye_avg_pos(good_trials,:).^2 lEye_avg_pos(good_trials,1).*lEye_avg_pos(good_trials,2)];

bx_left = robustfit(X,Y(good_trials,1));
by_left = robustfit(X,Y(good_trials,2));
left_offset = [bx_left(1) by_left(1)];
left_gain = [bx_left(2:end) by_left(2:end)];
bx_left = robustfit(Xq,Y(good_trials,1));
by_left = robustfit(Xq,Y(good_trials,2));
corrected_lpos = bsxfun(@plus,lEye_avg_pos*left_gain,left_offset);
corrected_all_leye = bsxfun(@plus,all_lEye*left_gain,left_offset);
qcorrected_lpos(:,1) = bsxfun(@plus,Xq*bx_left(2:end),bx_left(1));
qcorrected_lpos(:,2) = bsxfun(@plus,Xq*by_left(2:end),by_left(1));
X = [all_lEye all_lEye.^2 all_lEye(:,1).*all_lEye(:,2)];
qcorrected_all_leye(:,1) = bsxfun(@plus,X*bx_left(2:end),bx_left(1));
qcorrected_all_leye(:,2) = bsxfun(@plus,X*by_left(2:end),by_left(1));

figure
hold on
plot(lEye_avg_pos(:,1),lEye_avg_pos(:,2),'k.')
plot(fx,fy,'r*','linewidth',2)
xlim([-12 12]); ylim([-12 12])

figure
hold on
plot(corrected_lpos(:,1),corrected_lpos(:,2),'k.')
plot(qcorrected_lpos(:,1),qcorrected_lpos(:,2),'b.')
plot(fx,fy,'r*','linewidth',2)
xlim([-12 12]); ylim([-12 12])

%% FOR RIGHT EYE
Y = [fx' fy'];
X = rEye_avg_pos(good_trials,:);
bx_right = robustfit(X,Y(good_trials,1));
by_right = robustfit(X,Y(good_trials,2));
right_offset = [bx_right(1) by_right(1)];
right_gain = [bx_right(2:end) by_right(2:end)];
Xq = [rEye_avg_pos(good_trials,:) rEye_avg_pos(good_trials,:).^2 rEye_avg_pos(good_trials,1).*rEye_avg_pos(good_trials,2)];
bx_right = robustfit(Xq,Y(good_trials,1));
by_right = robustfit(Xq,Y(good_trials,2));
corrected_rpos = bsxfun(@plus,rEye_avg_pos*right_gain,right_offset);
corrected_all_reye = bsxfun(@plus,all_rEye*right_gain,right_offset);
qcorrected_rpos(:,1) = bsxfun(@plus,Xq*bx_right(2:end),bx_right(1));
qcorrected_rpos(:,2) = bsxfun(@plus,Xq*by_right(2:end),by_right(1));
X = [all_rEye all_rEye.^2 all_rEye(:,1).*all_rEye(:,2)];
qcorrected_all_reye(:,1) = bsxfun(@plus,X*bx_right(2:end),bx_right(1));
qcorrected_all_reye(:,2) = bsxfun(@plus,X*by_right(2:end),by_right(1));

figure
hold on
plot(rEye_avg_pos(:,1),rEye_avg_pos(:,2),'k.')
hold on
plot(fx,fy,'r*','linewidth',2)
xlim([-12 12]); ylim([-12 12])

figure
hold on
plot(corrected_rpos(:,1),corrected_rpos(:,2),'k.')
plot(qcorrected_rpos(:,1),qcorrected_rpos(:,2),'b.')
plot(fx,fy,'r*','linewidth',2)
xlim([-12 12]); ylim([-12 12])

%%
figure
hold on
plot(rEye_avg_pos(:,1),rEye_avg_pos(:,2),'k.')
hold on
plot(lEye_avg_pos(:,1),lEye_avg_pos(:,2),'r.')
plot(fx,fy,'r*','linewidth',2)
xlim([-12 12]); ylim([-12 12])

figure
hold on
plot(corrected_rpos(:,1),corrected_rpos(:,2),'k.','markersize',15)
plot(corrected_lpos(:,1),corrected_lpos(:,2),'r.','markersize',15)
plot(fx,fy,'r*','linewidth',2)
xlim([-12 12]); ylim([-12 12])
for i = 1:length(good_trials)
    line([corrected_lpos(i,1) corrected_rpos(i,1)],[corrected_lpos(i,2) corrected_rpos(i,2)])
end

%%
clear diff_loc correct_loc
n_groups = 25;
poss_x = -8:4:8;
poss_y = -8:4:8;
% avg_cor_diff = corrected_rpos - corrected_lpos;
avg_qcor_diff = qcorrected_all_reye - qcorrected_all_leye;
avg_cor_diff = corrected_all_reye - corrected_all_leye;
avg_ucor_diff = all_rEye - all_lEye;
gcnt = 1;
for xx = 1:length(poss_x)
    for yy = 1:length(poss_y)
        
        all_cur = find(fx==poss_x(xx) & fy == poss_y(yy));
        all_cur = find(ismember(all_tnum,all_cur) & all_infixation==1);
        diff_loc(gcnt,:) = nanmean(avg_cor_diff(all_cur,:));
        diff_qloc(gcnt,:) = nanmean(avg_qcor_diff(all_cur,:));
        diff_uloc(gcnt,:) = nanmean(avg_ucor_diff(all_cur,:));
        diff_mloc(gcnt,:) = nanmedian(avg_cor_diff(all_cur,:));
        correct_loc(gcnt,:) = [-poss_x(xx) poss_y(yy)];
        
        gcnt = gcnt + 1;
    end
end

plot(correct_loc(:,1),correct_loc(:,2),'o')
xlim([-10 10])
ylim([-10 10])
hold on
plot(correct_loc(:,1)+diff_loc(:,1),correct_loc(:,2)+diff_loc(:,2),'r*')
plot(correct_loc(:,1)+diff_uloc(:,1),correct_loc(:,2)+diff_uloc(:,2),'k*')
plot(correct_loc(:,1)+diff_qloc(:,1),correct_loc(:,2)+diff_qloc(:,2),'g*')


%%
cd ~/Data/bruce/7_15_12/G029
save eye_calibration_data right_offset right_gain left_offset left_gain bx* by*
