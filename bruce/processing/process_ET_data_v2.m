function [all_eye_vals,all_eye_ts,all_eye_speed,et_params] = process_ET_data_v2(all_t_axis,all_blockvec,cur_block_set,Expt_name,trial_toffset,use_coils)

if Expt_name(1) == 'G'
    expt_type = 'UA';
elseif Expt_name(1) == 'M'
    expt_type = 'LP';
else
    error('Unrecognized experiment name');
end

if nargin < 5 || isempty(trial_toffset)
    trial_toffset = zeros(size(cur_block_set));
end
if nargin < 6 || isempty(use_coils)
    use_coils = [1 0];
end

eye_smooth = 3;

if strcmp(expt_type,'UA')
    emfile = ['jbe' Expt_name '.em.mat']; %will eventually need to specify the animal name here...
    load(emfile);
end

all_eye_vals = [];
all_eye_speed = [];
all_eye_ts = [];
all_eye_blockvec = [];
for ee = 1:length(cur_block_set);
    fprintf('Loading ET data for expt %s, block %d of %d\n',Expt_name,ee,length(cur_block_set));
    
    if strcmp(expt_type,'LP')
        emfile = sprintf('lem%s.%d.em.mat',Expt_name,cur_block_set(ee));
        load(emfile);
    end
    
    cur_set = find(all_blockvec==ee);
    if ee > 1
        cur_toffset = trial_toffset(ee-1);
    else
        cur_toffset = 0;
    end
    if ~isempty(cur_set)
        [eye_vals_interp,eye_ts_interp,eye_insample] = get_eye_tracking_data_new(Expt,all_t_axis(cur_set([1 end])) - cur_toffset);
        
        eye_dt = Expt.Header.CRrates(1);
        eye_fs = 1/eye_dt;
        lEyeXY = eye_vals_interp(:,1:2);
        rEyeXY = eye_vals_interp(:,3:4);
        
        eye_vel = lEyeXY; %initialization
%         %slight smoothing before computing speed
%         if use_coils(1) && ~use_coils(2)
%             sm_avg_eyepos = lEyeXY; 
%         elseif use_coils(2) && ~use_coils(1)
%             sm_avg_eyepos = rEyeXY; 
%         elseif use_coils(1) && use_coils(2)
%             sm_avg_eyepos = 0.5*lEyeXY + 0.5*rEyeXY;
%         else
%            error('Invalid setting for use_coils!'); 
%         end
%uses smoothed left eye coil signal to define instantaneous speed
sm_avg_eyepos = lEyeXY;
        sm_avg_eyepos(:,1) = smooth(lEyeXY(:,1),eye_smooth);
        sm_avg_eyepos(:,2) = smooth(lEyeXY(:,2),eye_smooth);
        eye_vel(:,1) = [0; diff(sm_avg_eyepos(:,1))]/eye_dt;
        eye_vel(:,2) = [0; diff(sm_avg_eyepos(:,2))]/eye_dt;
        
        eye_speed = sqrt(eye_vel(:,1).^2+eye_vel(:,2).^2);
        
        all_eye_vals = [all_eye_vals; lEyeXY rEyeXY];
        all_eye_speed = [all_eye_speed; eye_speed];
        all_eye_ts = [all_eye_ts; eye_ts_interp' + cur_toffset];
        all_eye_blockvec = [all_eye_blockvec; ee*ones(size(eye_speed))];
    end
end

back_pts = 1 + find(diff(all_eye_ts) <= 0);
double_samples = [];
for i = 1:length(back_pts)
    next_forward = find(all_eye_ts > all_eye_ts(back_pts(i)-1),1,'first');
    double_samples = [double_samples back_pts(i):next_forward];
end
all_eye_ts(double_samples) = [];
all_eye_speed(double_samples) = [];
all_eye_vals(double_samples,:) = [];

et_params.eye_fs = eye_fs;
et_params.eye_smooth = eye_smooth;
et_params.use_coils = use_coils;